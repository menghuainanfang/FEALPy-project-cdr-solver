"""Frozen-metric comparisons across mesh size, perturbation and BDF accuracy."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh

from experiments.e_variational_mesh_visual import git_state
from src.fealpy_mesh_evolution import solve
from src.variational_mesh import band_metric

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/mentor_2026-10-23/comparison"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cap", type=float, default=5.)
    args = parser.parse_args()
    bm.set_backend("numpy")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for n in (8, 16):
        mesh = TriangleMesh.from_box([0., 1., 0., 1.], nx=n, ny=n)
        reference, cells, boundary = np.asarray(mesh.node), np.asarray(mesh.cell), np.asarray(mesh.boundary_node_flag())
        for amplitude in (.05, .30):
            physical = reference.copy()
            physical[~boundary] += np.random.default_rng(42).uniform(-amplitude/n, amplitude/n, (sum(~boundary), 2))
            for name, mode in (("smoothing", "tangent"), ("adaptive_fixed", "fixed"), ("adaptive", "tangent")):
                metric = (np.broadcast_to(np.eye(2), (len(cells), 2, 2)).copy() if name == "smoothing"
                          else band_metric(physical[cells].mean(axis=1)))
                for rtol, atol in ((1e-4, 1e-6), (1e-6, 1e-8)):
                    for policy in ("fixed_cap", "convergence"):
                        config = {"n": n, "amplitude": amplitude, "scene": name, "rtol": rtol,
                                  "atol": atol, "policy": policy, "seed": 42}
                        options = {"convergence": None} if policy == "fixed_cap" else {}
                        result = solve(reference, physical, cells, boundary, metric,
                                       duration=.5 if policy == "fixed_cap" else args.cap,
                                       frames=9, rtol=rtol, atol=atol, boundary_mode=mode,
                                       metric_function=band_metric if name != "smoothing" else None, **options)
                        trace = result["accepted_step_trace"]
                        final = result["frames"][-1]
                        key = f"n{n}_{name}_a{amplitude}_rtol{rtol}_{policy}"
                        record = {"config": config, "project": git_state(ROOT), "reference": reference.tolist(),
                                  "physical": physical.tolist(), "cells": cells.tolist(), "metric": metric.tolist(),
                                  "result": result}
                        (OUTPUT / f"{key}.json").write_text(json.dumps(record, allow_nan=False), encoding="utf-8")
                        with (OUTPUT / f"{key}_accepted.csv").open("w", newline="", encoding="utf-8") as f:
                            fields = list(dict.fromkeys(k for row in trace for k in row))
                            writer = csv.DictWriter(f, fieldnames=fields)
                            writer.writeheader()
                            writer.writerows(trace)
                        row = {**config, "stop_reason": result["stop_reason"],
                               "pseudo_time": result["terminal_pseudo_time"], **result["solver_stats"],
                               "wall_seconds": result["wall_seconds"], "cpu_seconds": result["cpu_seconds"],
                               "energy_initial": trace[0]["energy"], "energy_final": final["energy"],
                               "velocity_norm": final["velocity_norm"],
                               "scaled_residual_rms": final["scaled_residual_rms"],
                               "residual_threshold": result["residual_threshold"],
                               "energy_rate": trace[-1]["relative_energy_rate"],
                               "relative_displacement": trace[-1]["relative_displacement"],
                               "min_area": min(min(t["physical_min_area"], t["computational_min_area"]) for t in trace),
                               "invalid_accepted_states": sum(min(t["physical_min_area"], t["computational_min_area"]) <= 0 for t in trace),
                               "accepted_energy_increases_over_1e-10": int(sum(np.diff([t["energy"] for t in trace]) > 1e-10)),
                               "interpolation_rmse": final["interpolation_rmse"]}
                        rows.append(row)
                        print(key, row["stop_reason"], row["accepted_bdf_steps"], flush=True)
                        if amplitude == .30 and rtol == 1e-6 and policy == "convergence":
                            col = ("smoothing", "adaptive_fixed", "adaptive").index(name)
                            times = [t["pseudo_time"] for t in trace]
                            axes[0, col].semilogy(times, [max(t["scaled_residual_rms"], 1e-16) for t in trace], label=f"n={n}")
                            axes[1, col].plot(times, [t["energy"] for t in trace], label=f"n={n}")
                            axes[0, col].axhline(result["residual_threshold"], ls=":", alpha=.5)
    for col, name in enumerate(("smoothing", "adaptive_fixed", "adaptive")):
        axes[0, col].set(title=name, xlabel="Mesh pseudo-time s", ylabel="Scaled residual RMS")
        axes[1, col].set(xlabel="Mesh pseudo-time s", ylabel="Frozen-metric energy")
        for ax in axes[:, col]:
            ax.grid(alpha=.25)
            ax.legend()
    fig.savefig(OUTPUT / "convergence.png", dpi=150)
    plt.close(fig)
    (OUTPUT / "summary.json").write_text(json.dumps(rows, indent=2, allow_nan=False), encoding="utf-8")
    with (OUTPUT / "summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
