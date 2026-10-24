"""P1 backward-Euler Burgers evolution with nonconservative rezoning."""

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from matplotlib.tri import LinearTriInterpolator, Triangulation
from scipy.sparse.linalg import spsolve
from fealpy.backend import backend_manager as bm
from fealpy.decorator import barycentric, cartesian
from fealpy.fem import (BilinearForm, LinearForm, ScalarMassIntegrator,
                        ScalarDiffusionIntegrator, ScalarConvectionIntegrator, ScalarSourceIntegrator)
from fealpy.functionspace import LagrangeFESpace
from fealpy.mesh import TriangleMesh

from src.geometry import signed_triangle_areas
from src.fealpy_mesh_evolution import solve as move_mesh


@dataclass(frozen=True)
class BurgersData:
    """u_t + u*(u_x+u_y) - nu*Delta(u) = f on [0,1]^2."""

    nu: float = .05
    constant: bool = False

    def fields(self, p, t):
        if self.constant:
            return np.ones(p.shape[:-1]), np.zeros_like(p), np.zeros(p.shape[:-1])
        phase = np.pi*(p[..., 0]-.4*t)
        y = np.pi*p[..., 1]
        scale = np.exp(-t)
        u = scale*np.sin(phase)*np.sin(y)
        grad = scale*np.pi*np.stack((np.cos(phase)*np.sin(y), np.sin(phase)*np.cos(y)), axis=-1)
        ut = -u-.4*grad[..., 0]
        source = ut+u*grad.sum(axis=-1)+2*self.nu*np.pi**2*u
        return u, grad, source


def transfer(old_nodes, cells, values, new_nodes):
    """P1 interpolation on the old physical mesh; exact on constants/affines."""
    if signed_triangle_areas(old_nodes, cells).min() <= 0 or signed_triangle_areas(new_nodes, cells).min() <= 0:
        raise ValueError("invalid_mesh: nonpositive rezoning triangle")
    result = LinearTriInterpolator(Triangulation(old_nodes[:, 0], old_nodes[:, 1], cells), values)(
        new_nodes[:, 0], new_nodes[:, 1])
    if np.ma.getmaskarray(result).any():
        raise ValueError("Rezoning nodes lie outside the old mesh")
    return np.asarray(result)


def mass(nodes, cells, values):
    return float(np.sum(signed_triangle_areas(nodes, cells)*values[cells].mean(axis=1)))


def numerical_metric(nodes, cells, values, strength=2.):
    triangles = nodes[cells]
    edges = triangles[:, 1:]-triangles[:, :1]
    difference = values[cells[:, 1:]]-values[cells[:, :1]]
    grad = np.linalg.solve(edges, difference[..., None])[..., 0]
    norm2 = np.sum(grad**2, axis=1)
    return np.eye(2)[None]+strength*(grad[:, :, None]*grad[:, None, :])/(1+norm2[:, None, None])


def solve_step(nodes, cells, old_values, time, dt, data=BurgersData(), *, tolerance=1e-9, max_picard=30):
    """Nonlinear Galerkin step on the current, fixed physical mesh."""
    if signed_triangle_areas(nodes, cells).min() <= 0:
        raise ValueError("invalid_mesh: Burgers triangles must be positively oriented")
    mesh = TriangleMesh(bm.asarray(nodes), bm.asarray(cells))
    space = LagrangeFESpace(mesh, p=1)

    def matrix(integrator):
        form = BilinearForm(space)
        form.add_integrator(integrator)
        return form.assembly().to_scipy()

    M = matrix(ScalarMassIntegrator(q=4))
    K = matrix(ScalarDiffusionIntegrator(q=4, method=None))
    boundary = np.asarray(mesh.boundary_node_flag())
    free = ~boundary
    exact_boundary = data.fields(nodes[boundary], time)[0]

    @cartesian
    def forcing(p):
        return data.fields(p, time)[2]

    load = LinearForm(space)
    load.add_integrator(ScalarSourceIntegrator(forcing, q=6))
    rhs = M@old_values/dt+np.asarray(load.assembly())

    def operator(values):
        uh = space.function(bm.asarray(values))

        @barycentric
        def convection(bcs, index):
            value = uh(bcs, index)
            return np.stack((value, value), axis=-1)

        return M/dt+data.nu*K+matrix(ScalarConvectionIntegrator(convection, q=4))

    values = old_values.copy()
    values[boundary] = exact_boundary
    residual_history = []
    status = "nonlinear_iteration_cap"
    for iteration in range(1, max_picard+1):
        A = operator(values)
        values[free] = spsolve(A[free][:, free], rhs[free]-A[free][:, boundary]@exact_boundary)
        residual = operator(values)@values-rhs
        relative = float(np.linalg.norm(residual[free])/max(np.linalg.norm(rhs[free]), 1.))
        residual_history.append(relative)
        if relative <= tolerance:
            status = "converged"
            break
    uh = space.function(bm.asarray(values))

    @cartesian
    def exact(p):
        return data.fields(p, time)[0]

    return values, {"status": status, "nonlinear_iterations": iteration,
                    "relative_nonlinear_residual": residual_history[-1], "residual_history": residual_history,
                    "boundary_error": float(np.max(np.abs(values[boundary]-exact_boundary))),
                    "l2_error": float(mesh.error(exact, uh, q=6)), "mass": mass(nodes, cells, values)}


def run(*, n=8, steps=3, dt=.02, inner_cap=3., max_outer=6, moving=True,
        static=False, perturbation=0., data=BurgersData(), save_inner_steps=None,
        inner_frames=8, on_step=None, stop_on_outer_cap=False):
    """Alternate PDE/metric/mesh at t_(n+1); transfer from the fixed t_n state.

    Each middle iteration solves the same physical step, rather than advancing
    time repeatedly. A final PDE solve on the retained mesh supplies its field.
    No ALE transport term is used: transfer and PDE evolution are split.
    """
    if n < 2 or steps < 1 or dt <= 0 or inner_cap <= 0 or max_outer < 1:
        raise ValueError("n >= 2, steps/max_outer >= 1 and positive time steps/caps are required")
    bm.set_backend("numpy")
    started = perf_counter()
    mesh = TriangleMesh.from_box([0., 1., 0., 1.], nx=n, ny=n)
    reference, cells = np.array(mesh.node), np.array(mesh.cell)
    boundary = np.asarray(mesh.boundary_node_flag())
    nodes = reference.copy()
    nodes[~boundary] += np.random.default_rng(42).uniform(-perturbation/n, perturbation/n, (sum(~boundary), 2))
    if signed_triangle_areas(nodes, cells).min() <= 0:
        raise ValueError("invalid_mesh: perturbed initial physical mesh")
    values = data.fields(nodes, 0.)[0]
    states = [{"physical_time": 0., "physical_step_index": 0, "nodes": nodes.tolist(), "values": values.tolist(),
               "metric": numerical_metric(nodes, cells, values).tolist(), "outer_iterations": []}]
    states[0].update(t_start=0., t_end=0., physical_dt=0., outer_status="not_applicable",
                     min_signed_area=float(signed_triangle_areas(nodes,cells).min()),
                     input_physical=nodes.tolist(), k_actual=0,k_max=max_outer,outer_residual=None,outer_tol=1.)
    if on_step:
        on_step(states[0])
    metric_version = 0
    for step in range(1 if static else steps):
        time = 0. if static else (step+1)*dt
        base_nodes, base_values = nodes.copy(), values.copy()
        iterations = []
        previous_solution = None
        for k in range(max_outer if moving and not static else 1):
            old_on_current = transfer(base_nodes, cells, base_values, nodes)
            if static:
                provisional, pde = old_on_current, {"status": "not_advanced_static_input", "l2_error": None}
            else:
                provisional, pde = solve_step(nodes, cells, old_on_current, time, dt, data)
                if pde["status"] != "converged":
                    raise RuntimeError("Burgers nonlinear solve did not converge")
            frozen_metric = numerical_metric(nodes, cells, provisional)
            metric_version += 1
            old_mesh = nodes.copy()
            if moving:
                inner = move_mesh(reference, nodes, cells, boundary, frozen_metric,
                                  duration=inner_cap, frames=inner_frames,
                                  save_inner_meshes=save_inner_steps is None or step+1 in save_inner_steps,
                                  epoch_id=f"physical_step_{0 if static else step+1}_outer_{k+1}_metric_{metric_version}",
                                  metric_provenance="I + 2*grad(U_h)*grad(U_h)^T/(1+|grad(U_h)|^2), numerical P1 cell gradient")
                if inner["stop_reason"] in ("invalid_mesh", "solver_failure", "max_accepted_steps"):
                    raise RuntimeError(f"Mesh step failed: {inner['stop_reason']}: {inner['solver_message']}")
                nodes = np.array(inner["frames"][-1]["physical"])
                inner["endpoint_states"] = [inner["frames"][0], inner["frames"][-1]]
            else:
                inner = None
            mesh_change = float(np.max(np.linalg.norm(nodes-old_mesh, axis=1))*n)
            # Compare the same physical-time solution on the current mesh.
            field_change = None if previous_solution is None else float(
                np.linalg.norm(provisional-transfer(previous_solution[0], cells, previous_solution[1], old_mesh))
                /max(np.linalg.norm(provisional), 1.))
            consistent = bool(mesh_change <= 1e-3 and field_change is not None and field_change <= 1e-4
                              and inner is not None and inner["stop_reason"] == "converged")
            outer_residual = None if field_change is None else max(mesh_change/1e-3, field_change/1e-4)
            outer_status = ("not_applicable" if static or not moving else "not_checked" if field_change is None else
                            "converged" if consistent else "cap_reached" if k+1==max_outer else "not_converged")
            unmet = (["首次耦合轮尚无前一轮同时间解，场一致性未检查"] if field_change is None else [])
            if mesh_change>1e-3:
                unmet.append(f"网格变化 {mesh_change:.6g} > 0.001")
            if field_change is not None and field_change>1e-4:
                unmet.append(f"同时间场变化 {field_change:.6g} > 0.0001")
            if inner is not None and not inner["criterion_met"]:
                unmet.append("内层网格判据未满足")
            iterations.append({"outer_iteration": k+1, "physical_time": time, "metric_version": metric_version,
                               "frozen_physical": old_mesh.tolist(), "frozen_metric": frozen_metric.tolist(),
                               "provisional_values": provisional.tolist(), "pde": pde,
                               "mesh_inner": inner, "relative_mesh_change": mesh_change,
                               "relative_field_change": field_change, "outer_converged": consistent,
                               "outer_status": outer_status, "outer_residual": outer_residual, "outer_tol": 1.,
                               "outer_unmet_reasons": unmet})
            previous_solution = old_mesh, provisional.copy()
            if consistent or not moving or static:
                break
        transferred = transfer(base_nodes, cells, base_values, nodes)
        transfer_mass_change = mass(nodes, cells, transferred)-mass(base_nodes, cells, base_values)
        if static:
            values, final_pde = transferred, {"status": "not_advanced_static_input", "l2_error": None}
        else:
            values, final_pde = solve_step(nodes, cells, transferred, time, dt, data)
            if final_pde["status"] != "converged":
                raise RuntimeError("Final Burgers nonlinear solve did not converge")
        states.append({"physical_time": time, "physical_step_index": 0 if static else step+1,
                       "nodes": nodes.tolist(), "values": values.tolist(),
                       "metric": numerical_metric(nodes, cells, values).tolist(), "outer_iterations": iterations,
                       "pde": final_pde, "outer_converged": iterations[-1]["outer_converged"],
                       "outer_stop_reason": ("static_single_optimization" if static else
                                             "fixed_grid" if not moving else
                                             "converged" if iterations[-1]["outer_converged"] else "outer_iteration_cap"),
                       "transfer_mass_change": transfer_mass_change,
                       "min_signed_area": float(signed_triangle_areas(nodes, cells).min()),
                       "physical_dt": 0. if static else dt})
        states[-1].update(t_start=0. if static else step*dt, t_end=time,
                          input_physical=base_nodes.tolist(), k_actual=len(iterations), k_max=max_outer,
                          outer_residual=iterations[-1]["outer_residual"], outer_tol=1.,
                          outer_status=iterations[-1]["outer_status"],
                          overall_success=bool(static or not moving or iterations[-1]["outer_converged"]))
        if on_step:
            on_step(states[-1])
        if stop_on_outer_cap and moving and not static and not iterations[-1]["outer_converged"]:
            break
    return {"config": {"n": n, "steps": steps, "physical_dt": dt, "inner_pseudo_time_cap": inner_cap,
                        "max_outer": max_outer, "moving": moving, "static": static,
                        "perturbation": perturbation, "nu": data.nu, "constant_solution": data.constant},
            "reference": reference.tolist(), "cells": cells.tolist(),
            "edges": np.asarray(mesh.entity("edge")).tolist(), "states": states,
            "wall_seconds": perf_counter()-started,
            "method": "P1 nonlinear Galerkin / backward Euler / Picard / nonconservative P1 rezoning",
            "field_origin": ("Initial FE data transferred by P1; PDE time is not advanced" if static else
                             "Numerical PDE solve; exact solution is used only for initial/boundary/source and error"),
            "transfer": "Interpolate the fixed t_n FE field onto each trial mesh; not mass conservative"}
