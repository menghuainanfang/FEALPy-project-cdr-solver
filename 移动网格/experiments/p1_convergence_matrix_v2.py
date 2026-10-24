"""Bounded P1-v2 checks of caps, solver accuracy, sampling and continuation."""

import csv
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import sys

import fealpy
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh
import numpy as np
import scipy

from experiments.e_variational_mesh_visual import git_state
from src.fealpy_mesh_evolution import Convergence, solve
from src.variational_mesh import band_metric

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/mentor_2026-10-23/p1_tangent/matrix"


def problem(n, name, amplitude=.2):
    mesh = TriangleMesh.from_box(nx=n, ny=n)
    x, c, b = np.array(mesh.node), np.array(mesh.cell), np.array(mesh.boundary_node_flag())
    physical = x.copy()
    physical[~b] += np.random.default_rng(42).uniform(-amplitude/n, amplitude/n, (sum(~b), 2))
    metric = np.broadcast_to(np.eye(2), (len(c), 2, 2)).copy() if name == "smoothing" else band_metric(physical[c].mean(axis=1))
    return x, physical, c, b, metric, "tangent"


def main():
    bm.set_backend("numpy")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("p1_legacy_matrix", ROOT / "results/mentor_2026-10-23/p1_v2/legacy_evolution.py")
    legacy = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = legacy
    spec.loader.exec_module(legacy)
    rows, sampling, equivalence = [], {}, []

    def run_case(n, scene, cap, rtol, atol, *, amplitude=.2, frames=61, group="core", continuation=False):
        x, physical, c, b, m, mode = problem(n, scene, amplitude)
        result = solve(x, physical, c, b, m, boundary_mode=mode, duration=cap, frames=frames,
                       rtol=rtol, atol=atol, stop_on_convergence=not continuation,
                       epoch_id=f"{scene}_n{n}_a{amplitude}_frozen")
        trace = result["convergence_trace"]
        steps = [r["accepted_dt"] for r in result["accepted_step_trace"][1:]]
        first = next((r for r in trace if r["criterion_met"]), None)
        row = {"group": group, "scenario": scene, "nx": n, "amplitude": amplitude,
               "threshold": result["residual_threshold"], "rtol": rtol, "atol": atol,
               "pseudo_time_cap": cap, "first_converged_s": result["first_converged_s"],
               "final_s": result["terminal_pseudo_time"], "stop_reason": result["stop_reason"],
               "mesh_status": result["mesh_status"], "solver_status": result["solver_status"],
               "criterion_met": result["criterion_met"], "cap_reached": result["cap_reached"],
               "requested_frames": frames, "display_frame_count": result["saved_frames"],
               "bdf_accepted_steps": result["solver_stats"]["accepted_bdf_steps"],
               "nfev": result["solver_stats"]["nfev"], "njev": result["solver_stats"]["njev"], "nlu": result["solver_stats"]["nlu"],
               "rejected_bdf_attempts": None,
               "min_signed_area": min(min(r["physical_min_area"], r["computational_min_area"]) for r in result["accepted_step_trace"]),
               "residual_initial": trace[0]["mesh_residual_rms"], "residual_final": trace[-1]["mesh_residual_rms"],
               "energy_initial": trace[0]["mesh_energy"], "energy_final": trace[-1]["mesh_energy"],
               "first_converged_energy": None if first is None else first["mesh_energy"],
               "criterion_losses_after_first": 0 if first is None else sum(not r["criterion_met"] for r in trace if r["mesh_pseudo_time"] > first["mesh_pseudo_time"]),
               "accepted_dt_min": min(steps), "accepted_dt_median": float(np.median(steps)), "accepted_dt_max": max(steps),
               "free_scalar_dofs": result["free_scalar_dofs"], "wall_seconds": result["wall_seconds"]}
        key = f"{group}_{scene}_n{n}_a{amplitude}_cap{cap}_r{rtol}_f{frames}"
        endpoint = result["frames"][-1]
        record = {"config": {"scenario": scene, "n": n, "amplitude": amplitude, "seed": 42, "frames": frames,
                              "pseudo_time_cap": cap, "rtol": rtol, "atol": atol, "max_step": cap/10,
                              "stop_on_convergence": not continuation, "convergence": asdict(Convergence())},
                  "reference": x.tolist(), "fixed_physical": physical.tolist(), "cells": c.tolist(),
                  "frozen_metric": m.tolist(), "endpoint": endpoint,
                  "result": {k: v for k, v in result.items() if k not in ("frames", "history")},
                  "output_samples": [{k: v for k, v in f.items() if k not in ("physical", "computational")} for f in result["frames"]]}
        (OUTPUT / f"{key}.json").write_text(json.dumps(record, allow_nan=False), encoding="utf-8")
        with (OUTPUT / f"{key}_trace.csv").open("w", newline="", encoding="utf-8-sig") as f:
            flat = [{k: json.dumps(v, ensure_ascii=False) if isinstance(v, list) else v for k, v in r.items()} for r in trace]
            writer = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in flat for k in r)))
            writer.writeheader()
            writer.writerows(flat)
        rows.append(row)
        if cap == 3. and rtol == 1e-6 and amplitude == .2 and not continuation:
            sampling.setdefault((n, scene), {})[frames] = {"row": row, "times": [r["pseudo_time"] for r in result["accepted_step_trace"]],
                                                        "endpoint": endpoint["computational"]}
        print(key, result["stop_reason"], result["first_converged_s"], flush=True)
        return result

    for n in (8, 16):
        for scene in ("smoothing", "adaptive"):
            for cap in (.5, 3., 5.):
                for rtol, atol in ((1e-4, 1e-6), (1e-6, 1e-8)):
                    run_case(n, scene, cap, rtol, atol)
            for amplitude in (.05, .3):
                run_case(n, scene, 3., 1e-6, 1e-8, amplitude=amplitude, group="perturbation")
            for frames in (21, 121):
                run_case(n, scene, 3., 1e-6, 1e-8, frames=frames, group="sampling")
            continued = run_case(n, scene, 5., 1e-6, 1e-8, group="continuation", continuation=True)
            x, physical, c, b, m, mode = problem(n, scene)
            old = legacy.solve(x, physical, c, b, m, boundary_mode=mode, duration=5., frames=2, convergence=None)
            equivalence.append({"n": n, "scenario": scene, "config": "rtol=1e-6 atol=1e-8 cap=5 max_step=.5, no convergence stop",
                                "endpoint_max_difference": float(np.max(np.abs(np.array(old["frames"][-1]["computational"])-np.array(continued["frames"][-1]["computational"])))),
                                "nfev_old": old["solver_stats"]["nfev"], "nfev_new": continued["solver_stats"]["nfev"],
                                "accepted_steps_old": old["solver_stats"]["accepted_bdf_steps"], "accepted_steps_new": continued["solver_stats"]["accepted_bdf_steps"]})
    sampling_report = []
    for (n, scene), records in sampling.items():
        base = records[61]
        for frames in (21, 121):
            row = records[frames]
            sampling_report.append({"nx": n, "scenario": scene, "frames": frames, "reference_frames": 61,
                                    "accepted_times_equal": row["times"] == base["times"],
                                    "endpoint_max_difference": float(np.max(np.abs(np.array(row["endpoint"])-np.array(base["endpoint"])))),
                                    "stop_reason_equal": row["row"]["stop_reason"] == base["row"]["stop_reason"],
                                    "first_converged_s_equal": row["row"]["first_converged_s"] == base["row"]["first_converged_s"]})
    report = {"project": git_state(ROOT), "python": sys.executable, "numpy": np.__version__, "scipy": scipy.__version__,
              "fealpy_import": fealpy.__file__, "fealpy": git_state(Path(fealpy.__file__).resolve().parents[1]),
              "rows": rows, "sampling_comparisons": sampling_report, "integrator_equivalence": equivalence,
              "criterion_version": "projected-balanced-residual-window-v2"}
    (OUTPUT / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    with (OUTPUT / "summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
