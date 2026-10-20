"""F: equal-DOF Poisson comparisons with FEALPy boundary-constrained meshes."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from importlib.metadata import version
import json
from pathlib import Path
import sys
from time import perf_counter

import fealpy
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

from experiments.e_variational_mesh_visual import git_state
from src.poisson_mesh_comparison import PoissonCase, solve_poisson
from src.fealpy_mesh_evolution import solve as move_mesh, quality_metrics, boundary_diagnostics

ROOT = Path(__file__).resolve().parents[1]
LABELS = {"smooth": "光滑正弦解", "band": "弯曲陡变带", "peak": "局部高斯尖峰"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--levels", type=int, nargs="+", default=[8, 16, 32])
    parser.add_argument("--duration", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if min(args.levels) < 2 or args.duration <= 0:
        parser.error("Mesh levels must be >= 2 and duration > 0")
    bm.set_backend("numpy")
    rows, datasets, audit = [], {}, []
    for n in args.levels:
        started = perf_counter()
        mesh = TriangleMesh.from_box([0., 1., 0., 1.], nx=n, ny=n)
        uniform = np.asarray(mesh.entity("node")).copy()
        cells = np.asarray(mesh.entity("cell"))
        boundary = np.asarray(mesh.boundary_node_flag())
        generation_seconds = perf_counter() - started
        perturbed = uniform.copy()
        rng = np.random.default_rng(args.seed)
        perturbed[~boundary] += rng.uniform(-0.25 / n, 0.25 / n,
                                           (np.count_nonzero(~boundary), 2))
        started = perf_counter()
        smoothed_run = move_mesh(uniform, perturbed, cells, boundary,
                                 np.broadcast_to(np.eye(2), (len(cells), 2, 2)),
                                 duration=args.duration, frames=2)
        smoothing_seconds = perf_counter() - started
        smoothed = np.asarray(smoothed_run["frames"][-1]["physical"])
        for name in LABELS:
            problem = PoissonCase(name)
            started = perf_counter()
            metric = problem.metric(uniform[cells].mean(axis=1))
            adapted_run = move_mesh(uniform, uniform, cells, boundary, metric,
                                    duration=args.duration, frames=2,
                                    metric_function=problem.metric)
            adaptation_seconds = perf_counter() - started
            adapted = np.asarray(adapted_run["frames"][-1]["physical"])
            started = perf_counter()
            fixed_run = move_mesh(uniform, uniform, cells, boundary, metric,
                                  duration=args.duration, frames=2, boundary_mode="fixed",
                                  metric_function=problem.metric)
            fixed_seconds = perf_counter() - started
            fixed = np.asarray(fixed_run["frames"][-1]["physical"])
            variants = {"uniform": (uniform, 0., None),
                        "perturbed": (perturbed, 0., None),
                        "smoothed": (smoothed, smoothing_seconds, smoothed_run),
                        "adaptive_fixed": (fixed, fixed_seconds, fixed_run),
                        "adaptive": (adapted, adaptation_seconds, adapted_run)}
            key = f"{name}-{n}"
            datasets[key] = {"case": name, "n": n, "cells": cells.tolist(), "variants": {}}
            for method, (nodes, moving_seconds, movement) in variants.items():
                result = solve_poisson(nodes, cells, problem)
                row = {"case": name, "n": n, "method": method,
                       **result["metrics"], "generation_seconds": generation_seconds,
                       "moving_seconds": moving_seconds,
                       "total_seconds": generation_seconds + moving_seconds + result["metrics"]["pde_seconds"],
                       "mesh_stop": movement["stop_reason"] if movement else "not_moved",
                       "mesh_rhs_evaluations": movement["solver_stats"]["nfev"] if movement else 0}
                if movement:
                    history = movement["history"]
                    row["mesh_initial_energy"] = history[0]["energy"]
                    row["mesh_final_energy"] = history[-1]["energy"]
                    row["min_computational_area"] = min(h["computational_min_area"] for h in history)
                    row["max_energy_increase"] = float(max(np.diff([h["energy"] for h in history]), default=0.))
                    last = history[-1]
                    for metric_key in ("boundary_normal_displacement", "boundary_tangent_displacement",
                                "corner_displacement", "alignment_rms", "equidistribution_rms"):
                        row[metric_key] = last[metric_key]
                else:
                    row.update(mesh_initial_energy=None, mesh_final_energy=None,
                               min_computational_area=None, max_energy_increase=None)
                    for metric_key in ("boundary_normal_displacement", "boundary_tangent_displacement",
                                "corner_displacement", "alignment_rms", "equidistribution_rms"):
                        row[metric_key] = None
                    row.update(boundary_diagnostics(nodes, nodes))
                row.update(quality_metrics(nodes, uniform, cells, problem.metric(nodes[cells].mean(axis=1))))
                assert row["boundary_error"] < 1e-11
                assert row["relative_residual"] < 1e-10
                assert abs(row["total_area"] - 1) < 1e-12
                rows.append(row)
                result["metrics"] = row
                datasets[key]["variants"][method] = result
            baseline = datasets[key]["variants"]["uniform"]["metrics"]
            for method in variants:
                row = datasets[key]["variants"][method]["metrics"]
                row["l2_ratio_to_uniform"] = row["l2_error"] / baseline["l2_error"]
                row["h1_ratio_to_uniform"] = row["h1_seminorm_error"] / baseline["h1_seminorm_error"]
            if n == min(args.levels):
                for method in ("uniform", "adaptive"):
                    base = datasets[key]["variants"][method]
                    high = solve_poisson(np.asarray(base["nodes"]), cells, problem,
                                         quadrature=16, error_quadrature=20)
                    relative_l2_change = abs(high["metrics"]["l2_error"] / base["metrics"]["l2_error"] - 1)
                    audit.append({"case": name, "method": method, "n": n,
                                  "load_q": [12, 16], "error_q": [16, 20],
                                  "relative_l2_change": relative_l2_change})
                    assert relative_l2_change < 0.01, audit[-1]
            final = datasets[key]["variants"]["adaptive"]["metrics"]
            print(json.dumps({"case": name, "n": n, "uniform_l2": baseline["l2_error"],
                              "adaptive_l2": final["l2_error"], "ratio": final["l2_ratio_to_uniform"],
                              "adaptive_total_seconds": final["total_seconds"]}), flush=True)
    for name in LABELS:
        for method in ("uniform", "perturbed", "smoothed", "adaptive_fixed", "adaptive"):
            subset = sorted((r for r in rows if r["case"] == name and r["method"] == method), key=lambda r: r["n"])
            for i, row in enumerate(subset):
                previous = subset[i-1] if i else None
                row["l2_order"] = (float(np.log(previous["l2_error"] / row["l2_error"]) /
                                         np.log(row["n"] / previous["n"])) if previous else None)
                row["h1_order"] = (float(np.log(previous["h1_seminorm_error"] / row["h1_seminorm_error"]) /
                                         np.log(row["n"] / previous["n"])) if previous else None)
    smooth_rows = sorted((r for r in rows if r["case"] == "smooth" and r["method"] == "uniform"), key=lambda r:r["n"])
    if len(smooth_rows) > 1:
        assert min(r["l2_order"] for r in smooth_rows[1:]) > 1.8
        assert min(r["h1_order"] for r in smooth_rows[1:]) > 0.9
    result = {"config": vars(args), "labels": LABELS, "datasets": datasets,
              "rows": rows, "quadrature_audit": audit,
              "provenance": {"created_at": datetime.now().astimezone().isoformat(),
                             "python": sys.executable, "fealpy": version("fealpy"),
                             "fealpy_source": str(Path(fealpy.__file__).resolve()),
                             "fealpy_git": git_state(Path(fealpy.__file__).resolve().parents[1]),
                             "numpy": np.__version__, "scipy": scipy.__version__,
                             "project": git_state(ROOT), "backend": "numpy", "dtype": "float64"}}
    output = ROOT / "results"
    output.mkdir(exist_ok=True)
    encoded = json.dumps(result, ensure_ascii=False, allow_nan=False)
    (output / "f_poisson_comparison.json").write_text(encoded, encoding="utf-8")
    with (output / "f_poisson_comparison.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    template = Path(__file__).with_name("poisson_viewer.html").read_text(encoding="utf-8")
    (output / "f_poisson_comparison.html").write_text(template.replace("__PDE_DATA__", encoded), encoding="utf-8")
    fig, axes = plt.subplots(3, 3, figsize=(13, 12), constrained_layout=True)
    finest = max(args.levels)
    for i, name in enumerate(LABELS):
        ds = datasets[f"{name}-{finest}"]
        common_max = max(max(ds["variants"][m]["cell_l2_error"]) for m in ("uniform", "adaptive"))
        for j, method in enumerate(("uniform", "adaptive")):
            variant = ds["variants"][method]
            nodes = np.asarray(variant["nodes"])
            collection = axes[i,j].tripcolor(nodes[:,0], nodes[:,1], np.asarray(ds["cells"]),
                                             facecolors=variant["cell_l2_error"], cmap="magma", shading="flat",
                                             vmin=0, vmax=common_max)
            axes[i,j].triplot(nodes[:,0],nodes[:,1],np.asarray(ds["cells"]),color="white",lw=.15,alpha=.5)
            axes[i,j].set_aspect("equal")
            axes[i,j].set_title(f"{name} / {method}: cell L2 error")
            fig.colorbar(collection, ax=axes[i,j], shrink=.7)
        for method in ("uniform", "perturbed", "smoothed", "adaptive_fixed", "adaptive"):
            subset = sorted((r for r in rows if r["case"] == name and r["method"] == method), key=lambda r:r["n"])
            axes[i,2].loglog([r["dofs"] for r in subset], [r["l2_error"] for r in subset], "o-", label=method)
        axes[i,2].set(xlabel="P1 degrees of freedom", ylabel="L2 PDE solution error")
        axes[i,2].grid(alpha=.25)
        axes[i,2].legend()
    fig.savefig(output / "f_poisson_comparison.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
