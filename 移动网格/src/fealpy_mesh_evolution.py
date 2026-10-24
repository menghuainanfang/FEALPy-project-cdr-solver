"""Frozen-metric evolution using FEALPy's MetricTensorAdaptive numerical core."""

from __future__ import annotations

import numpy as np
from dataclasses import asdict
from time import perf_counter, process_time
from scipy.integrate import BDF, solve_ivp
from scipy.sparse import coo_matrix

from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh
from fealpy.mmesh import MetricTensorAdaptive
from fealpy.mmesh.metric_shared import GeometricDiscreteCore

from src.geometry import signed_triangle_areas
from src.mesh_convergence import AcceptedMeshState, Convergence, MeshConvergenceMonitor
from src.mesh_trace import epoch_metadata
from src.variational_mesh import reconstruct_physical, interpolation_rmse


class FealpyMeshEvolution(MetricTensorAdaptive):
    """Two-dimensional unit-square driver on FEALPy's schema-based mesh.

    Reuses the library's energy, local derivative, assembly, balancing and
    tangential boundary velocity. Geometry is frozen for one coordinate solve.
    """

    def __init__(self, physical, cell, metric, *, boundary_mode="tangent", gamma=2., tau=.1):
        self.mesh = TriangleMesh(bm.asarray(physical), bm.asarray(cell))
        self.cell = self.mesh.cell
        self.node = self.mesh.node
        self.NN = self.mesh.number_of_nodes()
        self.NC = self.mesh.number_of_cells()
        self.GD = 2
        self.kwargs0 = bm.context(self.node)
        self.gamma = gamma
        self.tau = tau
        self.boundary_mode = boundary_mode
        self.cm = self.mesh.entity_measure("cell")
        self.geo_core = GeometricDiscreteCore(self.mesh)
        self.R = self.geo_core.R_matrix()
        self.I_p = bm.broadcast_to(bm.eye(2), (self.NC, 2, 2)).copy()
        self.boundary = self.mesh.boundary_node_flag()
        on_x = (bm.abs(self.node[:, 0]) < 1e-12) | (bm.abs(self.node[:, 0]-1) < 1e-12)
        on_y = (bm.abs(self.node[:, 1]) < 1e-12) | (bm.abs(self.node[:, 1]-1) < 1e-12)
        self.Vertices_idx = bm.where(on_x & on_y)[0]
        self.Bdinnernode_idx = bm.where(self.boundary & ~(on_x & on_y))[0]
        normal = bm.zeros((self.NN, 2), **self.kwargs0)
        normal = bm.set_at(normal, bm.where(on_x)[0], bm.array([1., 0.]))
        normal = bm.set_at(normal, bm.where(on_y & ~on_x)[0], bm.array([0., 1.]))
        self.Bi_Lnode_normal = normal[self.Bdinnernode_idx]
        self.BD_projector()
        self._build_jac_pattern()
        self.metric = bm.asarray(metric)
        self.M_node = bm.zeros((self.NN, 2, 2), **self.kwargs0)
        weights = bm.zeros(self.NN, **self.kwargs0)
        self.M_node = bm.index_add(self.M_node, self.cell,
                                  self.cm[:, None, None, None] * self.metric[:, None])
        weights = bm.index_add(weights, self.cell, bm.broadcast_to(self.cm[:, None], self.cell.shape))
        self.M_node /= weights[:, None, None]
        self.metric_inverse = bm.linalg.inv(self.metric)
        self.metric_theta = self.theta(self.metric)
        self._prepare_ivp_cache(self.node, self.metric, self.metric_theta)

    def vector_construction(self, A, g, trA, E_hat, return_local=False):
        result = super().vector_construction(A, g, trA, E_hat, return_local=return_local)
        velocity = result[0] if return_local else result
        if self.boundary_mode == "fixed":
            velocity = bm.set_at(velocity, self.boundary, 0.)
        return (velocity, result[1]) if return_local else velocity

    def quantities(self, computational):
        E_hat = self.edge_matrix(computational)
        A = self.A(self._ivp_cache["E_K"], E_hat, self.metric_inverse)
        return A, bm.linalg.det(A), bm.trace(A, axis1=1, axis2=2), E_hat

    def energy_velocity(self, computational):
        A, g, trA, E_hat = self.quantities(bm.asarray(computational))
        value = self.I_func(self.metric_theta, self.lam(self.metric_theta),
                            trA, g, self.rho(self.metric))
        velocity = self.vector_construction(A, g, trA, E_hat)
        return float(value), np.asarray(bm.to_numpy(velocity))


def boundary_diagnostics(nodes, original):
    """Displacements relative to each node's original unit-square side."""
    left = np.isclose(original[:, 0], 0.)
    right = np.isclose(original[:, 0], 1.)
    bottom = np.isclose(original[:, 1], 0.)
    top = np.isclose(original[:, 1], 1.)
    corners = (left | right) & (bottom | top)
    boundary = left | right | bottom | top
    delta = nodes - original
    normal_error = max(np.max(np.abs(delta[left | right, 0]), initial=0.),
                       np.max(np.abs(delta[bottom | top, 1]), initial=0.))
    return {"boundary_normal_displacement": float(normal_error),
            "boundary_tangent_displacement": float(np.max(np.linalg.norm(delta[boundary & ~corners], axis=1), initial=0.)),
            "corner_displacement": float(np.max(np.linalg.norm(delta[corners], axis=1), initial=0.))}


def quality_metrics(nodes, reference, cells, metric):
    """Metric quality relative to the given reference triangles; ideal is one."""
    physical_mesh = TriangleMesh(bm.asarray(nodes), bm.asarray(cells))
    core = GeometricDiscreteCore(physical_mesh)
    E = core.edge_matrix(bm.asarray(nodes))
    E_ref = core.edge_matrix(bm.asarray(reference))
    F = E @ bm.linalg.inv(E_ref)
    gram = bm.swapaxes(F, -1, -2) @ bm.asarray(metric) @ F
    alignment = bm.trace(gram, axis1=1, axis2=2) / (2 * bm.sqrt(bm.linalg.det(gram)))
    volume = physical_mesh.entity_measure("cell") * bm.sqrt(bm.linalg.det(bm.asarray(metric)))
    relative_volume = volume / bm.mean(volume)
    return {"alignment_rms": float(bm.sqrt(bm.mean(alignment ** 2))),
            "equidistribution_rms": float(bm.sqrt(bm.mean(relative_volume ** 2)))}


def solve_ivp_baseline(reference, physical, cells, boundary, metric, *, duration=.5, frames=61,
          boundary_mode="tangent", metric_function=None):
    """Evolve one frozen-coordinate problem with sparse finite-difference BDF."""
    model = FealpyMeshEvolution(physical, cells, metric, boundary_mode=boundary_mode)
    size = model.NN
    rows = np.asarray(bm.to_numpy(model.I))
    cols = np.asarray(bm.to_numpy(model.J))
    sparsity = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(2*size, 2*size)).tocsr()

    def rhs(t, y):
        return model.energy_velocity(y.reshape(2, size).T)[1].ravel(order="F")

    output_times = np.linspace(0., duration, frames)
    solution = solve_ivp(rhs, (0., duration), np.asarray(reference).ravel(order="F"),
                         method="BDF", t_eval=output_times, jac_sparsity=sparsity,
                         rtol=1e-6, atol=1e-8, max_step=duration/10)
    if not solution.success:
        raise RuntimeError(solution.message)
    samples_1d = (np.arange(100) + .5) / 100
    xx, yy = np.meshgrid(samples_1d, samples_1d)
    samples = np.column_stack((xx.ravel(), yy.ravel()))
    history, saved = [], []
    for i, (time, y) in enumerate(zip(solution.t, solution.y.T)):
        computational = y.reshape(2, size).T.copy()
        on_x = np.isclose(reference[:, 0], 0.) | np.isclose(reference[:, 0], 1.)
        on_y = np.isclose(reference[:, 1], 0.) | np.isclose(reference[:, 1], 1.)
        computational[on_x, 0] = reference[on_x, 0]
        computational[on_y, 1] = reference[on_y, 1]
        if boundary_mode == "fixed":
            computational[boundary] = reference[boundary]
        moved = reconstruct_physical(computational, physical, cells, reference)
        computational_area = signed_triangle_areas(computational, cells)
        physical_area = signed_triangle_areas(moved, cells)
        if min(computational_area.min(), physical_area.min()) <= 0:
            raise ValueError("FEALPy evolution produced a nonpositive displayed triangle")
        energy, velocity = model.energy_velocity(computational)
        actual_metric = metric_function(moved[cells].mean(axis=1)) if metric_function else metric
        row = {"iteration": i, "pseudo_time": float(time), "energy": energy,
               "computational_min_area": float(computational_area.min()),
               "physical_min_area": float(physical_area.min()),
               "gradient_norm": float(np.linalg.norm(velocity)),
               "step": float(time - solution.t[i-1]) if i else 0., "backtracks": 0,
               **boundary_diagnostics(moved, physical),
               **quality_metrics(moved, reference, cells, actual_metric)}
        history.append(row)
        saved.append({**row, "computational": computational.tolist(), "physical": moved.tolist(),
                      "interpolation_rmse": interpolation_rmse(moved, cells, samples)})
    return {"history": history, "frames": saved, "stop_reason": "pseudo_time_end",
            "boundary_mode": boundary_mode, "integrator": "SciPy BDF",
            "algorithm": "FEALPy MetricTensorAdaptive",
            "solver_stats": {"nfev": solution.nfev, "njev": solution.njev, "nlu": solution.nlu},
            "rmse_sampling": "100x100 midpoint interpolation samples, not PDE error"}


def solve(reference, physical, cells, boundary, metric, *, duration=.5, frames=61,
          boundary_mode="tangent", metric_function=None, rtol=1e-6, atol=1e-8,
          tau=.1, convergence=Convergence(), max_internal_steps=2000,
          save_inner_meshes=False, epoch_id="frozen_mesh_0", stop_on_convergence=True,
          metric_provenance="caller-supplied SPD cell tensors, frozen during this solve"):
    """Step BDF with geometry checks and independent replay sampling.

    R = tau*v/h is the dimensionless balanced, projected residual, using the
    minimum incident reference edge length h at each node. The raw projected
    energy derivative is -tau*v/P, where P is FEALPy's frozen balance diagonal.
    The independent monitor checks fixed-window displacement and frozen energy
    change on accepted states. Replay sampling does not enter its state machine.
    """
    started, cpu_started = perf_counter(), process_time()
    reference, physical, cells = np.asarray(reference), np.asarray(physical), np.asarray(cells)
    boundary = np.asarray(boundary, dtype=bool)
    if duration <= 0 or frames < 2 or tau <= 0 or max_internal_steps < 1:
        raise ValueError("Positive duration/tau/step limit and at least two output frames are required")
    if boundary_mode not in ("fixed", "tangent"):
        raise ValueError("boundary_mode must be fixed or tangent")
    if min(signed_triangle_areas(reference, cells).min(), signed_triangle_areas(physical, cells).min()) <= 0:
        raise ValueError("invalid_mesh: initial triangles must have positive signed area")
    if not np.allclose(metric, np.swapaxes(metric, -1, -2)) or np.linalg.eigvalsh(metric).min() <= 0:
        raise ValueError("The frozen metric must be symmetric positive definite")
    model = FealpyMeshEvolution(physical, cells, metric, boundary_mode=boundary_mode, tau=tau)
    size = model.NN
    sparsity = coo_matrix((np.ones(len(model.I)), (model.I, model.J)),
                          shape=(2*size, 2*size)).tocsr()
    local_h = np.full(size, np.inf)
    for a, b in ((0, 1), (1, 2), (2, 0)):
        length = np.linalg.norm(reference[cells[:, a]]-reference[cells[:, b]], axis=1)
        np.minimum.at(local_h, cells[:, a], length)
        np.minimum.at(local_h, cells[:, b], length)
    on_x = np.isclose(reference[:, 0], 0.) | np.isclose(reference[:, 0], 1.)
    on_y = np.isclose(reference[:, 1], 0.) | np.isclose(reference[:, 1], 1.)
    free = np.ones(reference.shape, dtype=bool)
    free[on_x, 0] = False
    free[on_y, 1] = False
    if boundary_mode == "fixed":
        free[boundary] = False
    balance = np.asarray(model._ivp_cache["P_diag"])
    rhs_calls = 0

    def rhs(t, y):
        nonlocal rhs_calls
        rhs_calls += 1
        return model.energy_velocity(y.reshape(2, size).T)[1].ravel(order="F")

    def state(t, q):
        q = q.copy()
        bc = boundary_diagnostics(q, reference)
        if bc["boundary_normal_displacement"] > 1e-10 or bc["corner_displacement"] > 1e-10:
            raise ValueError("invalid_mesh: computational boundary constraint violated")
        if boundary_mode == "fixed" and np.max(np.abs(q[boundary]-reference[boundary])) > 1e-10:
            raise ValueError("invalid_mesh: fixed computational boundary moved")
        if signed_triangle_areas(q, cells).min() <= 0:
            raise ValueError("invalid_mesh: computational triangle lost positive orientation")
        q[on_x, 0], q[on_y, 1] = reference[on_x, 0], reference[on_y, 1]
        if boundary_mode == "fixed":
            q[boundary] = reference[boundary]
        moved = reconstruct_physical(q, physical, cells, reference)
        physical_bc = boundary_diagnostics(moved, physical)
        if max(physical_bc["boundary_normal_displacement"], physical_bc["corner_displacement"]) > 1e-10:
            raise ValueError("invalid_mesh: physical boundary constraint violated")
        if boundary_mode == "fixed" and physical_bc["boundary_tangent_displacement"] > 1e-10:
            raise ValueError("invalid_mesh: fixed physical boundary moved")
        moved[on_x, 0], moved[on_y, 1] = physical[on_x, 0], physical[on_y, 1]
        if boundary_mode == "fixed":
            moved[boundary] = physical[boundary]
        area = float(signed_triangle_areas(moved, cells).min())
        if area <= 0:
            raise ValueError("invalid_mesh: reconstructed physical triangle lost positive orientation")
        energy, velocity = model.energy_velocity(q)
        residual = (tau*velocity/local_h[:, None])[free]
        projected_gradient = (-tau*velocity/balance[:, None])[free]
        row = {"pseudo_time": float(t), "energy": energy,
               "velocity_norm": float(np.linalg.norm(velocity[free])),
               "scaled_residual_rms": float(np.sqrt(np.mean(residual**2))),
               "projected_gradient_norm": float(np.linalg.norm(projected_gradient)),
               "computational_min_area": float(signed_triangle_areas(q, cells).min()),
               "physical_min_area": area, **physical_bc}
        if not np.all(np.isfinite(list(row.values()))):
            raise ValueError("invalid_mesh: nonfinite geometry or energy")
        return row, q, moved

    first, q0, moved0 = state(0., reference)
    monitor = MeshConvergenceMonitor(local_h, tau=tau, pseudo_time_cap=duration,
                                     config=convergence, max_accepted_steps=max_internal_steps)

    def diagnose(row, q, j, *, geometry_valid=True, geometry_message="", solver_success=True):
        return monitor.check(AcceptedMeshState(row["pseudo_time"], q, row["energy"],
                                              row["scaled_residual_rms"], geometry_valid,
                                              geometry_message, epoch_id), j, solver_success=solver_success)

    convergence_trace = [{**diagnose(first, q0, 0), "pseudo_step_size": 0.,
                          "computational_min_area": first["computational_min_area"],
                          "physical_min_area": first["physical_min_area"]}]
    threshold = convergence_trace[0]["residual_threshold"]
    traces = [{**first, "accepted_step_index": 0, "accepted_dt": 0.,
               "relative_displacement": 0., "relative_energy_rate": 0.}]
    samples = [(0., q0)]
    if save_inner_meshes:
        traces[0].update(computational=q0.tolist(), physical=moved0.tolist())
    sampling_times = np.linspace(0., duration, frames)
    cursor, accepted_steps, reason, message = 1, 0, "pseudo_time_cap", ""
    invalid_accepted_states, invalid_output_samples = 0, 0
    solver = BDF(rhs, 0., reference.ravel(order="F"), duration, max_step=duration/10,
                 jac_sparsity=sparsity, rtol=rtol, atol=atol)
    convergence_trace[0].update(nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu)
    previous_q, previous_row = q0, first
    last_valid_t = 0.
    for _ in range(max_internal_steps):
        try:
            message = solver.step() or ""
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            reason, message = "solver_failure", str(error)
            convergence_trace.append(diagnose(previous_row, previous_q, accepted_steps, solver_success=False))
            convergence_trace[-1]["solver_status"] = "failed"
            break
        if solver.status == "failed":
            reason = "solver_failure"
            convergence_trace.append(diagnose(previous_row, previous_q, accepted_steps, solver_success=False))
            convergence_trace[-1]["solver_status"] = "failed"
            break
        accepted_steps += 1
        try:
            row, q, moved = state(solver.t, solver.y.reshape(2, size).T)
        except ValueError as error:
            reason, message = "invalid_mesh", str(error)
            invalid_accepted_states += 1
            failed_row = {**previous_row, "pseudo_time": float(solver.t)}
            convergence_trace.append(diagnose(failed_row, previous_q, accepted_steps,
                                              geometry_valid=False, geometry_message=message))
            convergence_trace[-1].update(mesh_energy=None, mesh_residual_rms=None,
                                         candidate_values_available=False)
            break
        dt = float(solver.step_size)
        displacement = float(np.max(np.linalg.norm((q-previous_q)/local_h[:, None], axis=1)))
        energy_rate = tau*abs(row["energy"]-previous_row["energy"])/(dt*max(abs(first["energy"]), 1.))
        row.update(accepted_step_index=len(traces), accepted_dt=dt, t_old=float(solver.t_old),
                   pseudo_step_size=dt,
                   relative_displacement=displacement, relative_energy_rate=energy_rate,
                   nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu)
        if save_inner_meshes:
            row.update(computational=q.tolist(), physical=moved.tolist())
        traces.append(row)
        diagnostic = diagnose(row, q, accepted_steps)
        diagnostic.update(pseudo_step_size=dt, nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
                          computational_min_area=row["computational_min_area"], physical_min_area=row["physical_min_area"])
        convergence_trace.append(diagnostic)
        last_valid_t = float(solver.t)
        dense = solver.dense_output()
        while cursor < len(sampling_times) and sampling_times[cursor] <= solver.t:
            t = float(sampling_times[cursor])
            samples.append((t, dense(t).reshape(2, size).T))
            cursor += 1
        if diagnostic["stop_reason"] is not None:
            if diagnostic["stop_reason"] != "converged" or stop_on_convergence:
                reason = diagnostic["stop_reason"]
                break
        if solver.status == "finished":
            break
        previous_q, previous_row = q, row
    else:
        reason = "max_accepted_steps"
    if samples[-1][0] < last_valid_t:
        samples.append((last_valid_t, q))
    samples_1d = (np.arange(100)+.5)/100
    xx, yy = np.meshgrid(samples_1d, samples_1d)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    history, saved = [], []
    for i, (t, sampled_q) in enumerate(samples):
        try:
            row, sampled_q, moved = state(t, sampled_q)
        except ValueError as error:
            reason, message = "invalid_mesh", str(error)
            invalid_output_samples += 1
            break
        actual_metric = metric_function(moved[cells].mean(axis=1)) if metric_function else metric
        row.update(output_frame_index=i, sample_dt=t-samples[i-1][0] if i else 0.,
                   display_frame_index=i, output_sample_dt=t-samples[i-1][0] if i else 0.,
                   **quality_metrics(moved, reference, cells, actual_metric))
        history.append(row)
        saved.append({**row, "computational": sampled_q.tolist(), "physical": moved.tolist(),
                      "interpolation_rmse": interpolation_rmse(moved, cells, points)})
    terminal = convergence_trace[-1]
    if not stop_on_convergence and reason not in ("invalid_mesh", "solver_failure", "max_accepted_steps"):
        reason = "pseudo_time_cap" if terminal["cap_reached"] else reason
    if invalid_output_samples:
        terminal.update(mesh_status="invalid_mesh", geometry_valid=False, criterion_met=False,
                        stop_reason="invalid_mesh", unmet_reasons=[message])
    terminal["stop_reason"] = reason
    wall_seconds = perf_counter()-started
    return {**epoch_metadata(monitor, convergence_trace, free_scalar_dofs=int(free.sum()), tau=tau,
                            local_h=local_h, frames=saved, wall_seconds=wall_seconds),
            "history": history, "frames": saved, "accepted_step_trace": traces,
            "stop_reason": reason, "solver_message": message,
            "terminal_pseudo_time": float(solver.t), "last_valid_pseudo_time": last_valid_t,
            "boundary_mode": boundary_mode, "integrator": "SciPy BDF.step",
            "pseudo_time_cap": duration, "metric_provenance": metric_provenance,
            "solver_config": {"rtol": rtol, "atol": atol, "max_step": duration/10,
                              "max_accepted_steps": max_internal_steps, "jacobian": "sparse finite difference"},
            "algorithm": "FEALPy MetricTensorAdaptive", "tau": tau,
            "solver_stats": {"accepted_bdf_steps": accepted_steps, "valid_accepted_states": len(traces)-1,
                             "invalid_accepted_states": invalid_accepted_states,
                             "invalid_output_samples": invalid_output_samples, "nfev": solver.nfev,
                             "njev": solver.njev, "nlu": solver.nlu,
                             "rhs_calls_including_jacobian_fd": rhs_calls,
                             "rejected_bdf_attempts": None,
                             "rejected_error_test_steps": None},
            "convergence": None if convergence is None else asdict(convergence),
            "residual_threshold": threshold, "free_scalar_dofs": int(free.sum()),
            "velocity_definition": "Projected FEALPy v = -P*projected(dI/dxi)/tau; Euclidean norm on free scalar dofs",
            "residual_definition": "R_rms = rms(tau*v_i/h_i) on free scalar dofs; h_i=min incident reference edge length",
            "wall_seconds": wall_seconds, "cpu_seconds": process_time()-cpu_started,
            "stop_on_convergence": stop_on_convergence,
            "rmse_sampling": "100x100 midpoint interpolation samples, not PDE error"}
