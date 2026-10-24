"""Tangential mesh optimization with nested frozen-metric inner solves."""

from dataclasses import asdict, dataclass
from time import perf_counter

from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh
import numpy as np

from src.fealpy_mesh_evolution import boundary_diagnostics, solve
from src.geometry import signed_triangle_areas
from src.poisson_mesh_comparison import PoissonCase


@dataclass(frozen=True)
class OuterConvergence:
    mesh_tolerance: float = 1e-3
    metric_tolerance: float = 1e-3
    consecutive_rounds: int = 2


def run(case="band", *, n=8, perturbation=.2, inner_cap=40., max_outer=60,
        outer=OuterConvergence(), metric_relaxation=.5, initial_physical=None,
        prescribed_metric=None, epoch_prefix=None, save_inner_meshes=True, verbose=True):
    """Rebuild the physical mesh and resample its metric until a fixed point.

    Each round starts computational coordinates at the same reference mesh.
    Its physical mesh and cell metric remain frozen during the inner solve.
    Only the accepted final map becomes the next round's physical mesh.
    """
    if case not in ("smooth", "band", "peak"):
        raise ValueError("case must be smooth, band or peak")
    if n<2 or max_outer<1 or inner_cap<=0 or not 0<metric_relaxation<=1:
        raise ValueError("n>=2, positive caps and 0<metric_relaxation<=1 are required")
    bm.set_backend("numpy")
    started = perf_counter()
    mesh = TriangleMesh.from_box(nx=n, ny=n)
    reference = np.asarray(mesh.node).copy()
    cells = np.asarray(mesh.cell)
    boundary = np.asarray(mesh.boundary_node_flag())
    nodes = reference.copy() if initial_physical is None else np.array(initial_physical).copy()
    if initial_physical is None:
        nodes[~boundary] += np.random.default_rng(42).uniform(-perturbation/n, perturbation/n, (sum(~boundary), 2))
    initial = nodes.copy()
    metric_function = prescribed_metric or PoissonCase(case).metric
    metric = metric_function(nodes[cells].mean(axis=1))
    rounds, consecutive = [], 0
    stop_reason = "outer_iteration_cap"
    for k in range(1, max_outer+1):
        physical = nodes.copy()
        inner = solve(reference, physical, cells, boundary, metric, duration=inner_cap,
                      frames=2, save_inner_meshes=save_inner_meshes, boundary_mode="tangent",
                      epoch_id=f"{epoch_prefix or case}_outer_{k}", metric_function=metric_function,
                      metric_provenance=f"prescribed analytic {case} metric, resampled at physical cell centers each outer round")
        nodes = np.asarray(inner["frames"][-1]["physical"])
        sampled_metric = metric_function(nodes[cells].mean(axis=1))
        next_metric = (1-metric_relaxation)*metric+metric_relaxation*sampled_metric
        mesh_change = float(n*np.linalg.norm(nodes-physical, axis=1).max())
        metric_change = float(np.linalg.norm(sampled_metric-metric)/max(np.linalg.norm(metric), 1.))
        geometry_valid = bool(signed_triangle_areas(nodes, cells).min()>0 and
                              boundary_diagnostics(nodes, initial)["boundary_normal_displacement"]<=1e-10 and
                              boundary_diagnostics(nodes, initial)["corner_displacement"]<=1e-10)
        inner_converged = bool(inner["criterion_met"] and inner["stop_reason"] == "converged")
        conditions_met = bool(inner_converged and geometry_valid and mesh_change<=outer.mesh_tolerance
                              and metric_change<=outer.metric_tolerance)
        consecutive = consecutive+1 if conditions_met else 0
        converged = consecutive>=outer.consecutive_rounds
        reasons = []
        if not inner_converged:
            reasons.append("本轮内层尚未完整收敛")
        if not geometry_valid:
            reasons.append("物理网格几何检查未通过")
        if mesh_change>outer.mesh_tolerance:
            reasons.append(f"网格变化 {mesh_change:.6g} > {outer.mesh_tolerance:g}")
        if metric_change>outer.metric_tolerance:
            reasons.append(f"度量变化 {metric_change:.6g} > {outer.metric_tolerance:g}")
        if not converged:
            reasons.append(f"外层连续达标 {consecutive}/{outer.consecutive_rounds}")
        rounds.append({"outer_iteration": k, "input_physical": physical.tolist(),
                       "frozen_metric": metric.tolist(), "output_physical": nodes.tolist(),
                       "sampled_metric": sampled_metric.tolist(), "updated_metric": next_metric.tolist(),
                       "mesh_change": mesh_change, "metric_change": metric_change,
                       "inner_converged": inner_converged, "geometry_valid": geometry_valid,
                       "conditions_met": conditions_met, "consecutive_success": consecutive,
                       "outer_converged": converged, "unmet_reasons": reasons,
                       "inner": {**{key: value for key, value in inner.items() if key not in ("frames", "history", "field_dictionary")},
                                 "endpoint_states": [inner["frames"][0], inner["frames"][-1]]}})
        if verbose:
            print(f"{case} k={k} inner={inner['stop_reason']} j={inner['solver_stats']['accepted_bdf_steps']} dx={mesh_change:.3g} dM={metric_change:.3g} consecutive={consecutive}", flush=True)
        if not inner_converged or not geometry_valid:
            stop_reason = "inner_not_converged" if not inner_converged else "invalid_mesh"
            break
        if converged:
            stop_reason = "converged"
            break
        metric = next_metric
    return {"case": case, "config": {"n": n, "seed": 42, "perturbation": perturbation,
                                      "inner_cap": inner_cap, "max_outer": max_outer,
                                      "boundary_mode": "tangent", "metric_relaxation": metric_relaxation},
            "outer_criterion_version": "mesh-metric-fixed-point-v1",
            "outer_convergence": {**asdict(outer),
                                  "mesh_definition": "max_i |x_i^(k+1)-x_i^k| / h_ref, h_ref=1/n",
                                  "metric_definition": "Frobenius(M_sample(x^(k+1))-M_frozen^k)/max(Frobenius(M_frozen^k),1)"},
            "reference": reference.tolist(), "initial_physical": initial.tolist(),
            "cells": cells.tolist(), "edges": np.asarray(mesh.entity("edge")).tolist(),
            "boundary": boundary.tolist(), "rounds": rounds, "stop_reason": stop_reason,
            "outer_converged": rounds[-1]["outer_converged"], "final_physical": nodes.tolist(),
            "wall_seconds": perf_counter()-started}
