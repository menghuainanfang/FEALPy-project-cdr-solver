"""E: visualize FEALPy coordinate evolution with fixed and tangential boundaries."""

from __future__ import annotations

import argparse
from datetime import datetime
from importlib.metadata import version
import json
import os
from pathlib import Path
import subprocess
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

from src.geometry import signed_triangle_areas
from src.variational_mesh import band_metric
from src.fealpy_mesh_evolution import solve
from src.mesh_trace import FIELD_DICTIONARY, SCHEMA_VERSION, render_diagnostic_assets

ROOT = Path(__file__).resolve().parents[1]


def write_viewer(result, output):
    result = {**result, "scenes": {name: scene for name, scene in result["scenes"].items()
                                  if scene["boundary_mode"] == "tangent"}}
    encoded = json.dumps(result, ensure_ascii=False, allow_nan=False)
    template = Path(__file__).with_name("mesh_viewer.html").read_text(encoding="utf-8")
    (output / "e_variational_mesh.html").write_text(
        render_diagnostic_assets(template).replace("__EXPERIMENT_DATA__", encoded)
        .replace("../demos/index.html", Path(os.path.relpath(ROOT / "demos/index.html", output)).as_posix())
        .replace("../index.html", Path(os.path.relpath(ROOT / "index.html", output)).as_posix()), encoding="utf-8")


def git_state(folder):
    try:
        revision = subprocess.run(["git", "-C", str(folder), "rev-parse", "HEAD"],
                                  capture_output=True, text=True)
        if revision.returncode:
            return {"sha": None, "dirty": None}
        status = subprocess.run(["git", "-C", str(folder), "status", "--porcelain"],
                                capture_output=True, text=True, check=True)
        return {"sha": revision.stdout.strip(), "dirty": bool(status.stdout.strip())}
    except FileNotFoundError:
        return {"sha": None, "dirty": None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=16)
    parser.add_argument("--duration", type=float, default=0.5)
    parser.add_argument("--frames", type=int, default=61)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--rtol", type=float, default=1e-6)
    parser.add_argument("--atol", type=float, default=1e-8)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    if args.nx < 2 or args.duration <= 0 or args.frames < 2:
        parser.error("nx >= 2, duration > 0 and frames >= 2 are required")
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box([0., 1., 0., 1.], nx=args.nx, ny=args.nx)
    reference = np.asarray(bm.to_numpy(mesh.entity("node"))).copy()
    cells = np.asarray(bm.to_numpy(mesh.entity("cell")), dtype=np.int64)
    boundary = np.asarray(bm.to_numpy(mesh.boundary_node_flag()), dtype=bool)
    if np.any(signed_triangle_areas(reference, cells) <= 0):
        raise ValueError("FEALPy mesh orientation is not positive")
    rng = np.random.default_rng(args.seed)
    perturbed = reference.copy()
    perturbed[~boundary] += rng.uniform(-0.20 / args.nx, 0.20 / args.nx,
                                       size=(np.count_nonzero(~boundary), 2))
    identity = np.broadcast_to(np.eye(2), (len(cells), 2, 2)).copy()
    centers = reference[cells].mean(axis=1)
    started = perf_counter()
    scenes = {}
    for name, physical, metric, mode in (
            ("smoothing", perturbed, identity, "tangent"),
            ("adaptive", reference, band_metric(centers), "tangent")):
        metric_function = band_metric if name.startswith("adaptive") else None
        scenes[name] = solve(reference, physical, cells, boundary, metric,
                             duration=args.duration, frames=args.frames,
                             rtol=args.rtol, atol=args.atol,
                             epoch_id=f"E_{name}_frozen_metric",
                             metric_provenance="identity cell metric" if name == "smoothing" else "src.variational_mesh.band_metric(strength=30), prescribed analytic band",
                             boundary_mode=mode, metric_function=metric_function)
        scenes[name]["fixed_physical"] = physical.tolist()
        scenes[name]["frozen_metric"] = metric.tolist()
    elapsed = perf_counter() - started
    provenance = {"created_at": datetime.now().astimezone().isoformat(),
                  "python": sys.executable, "fealpy_version": version("fealpy"),
                  "fealpy_import": fealpy.__file__, "numpy": np.__version__,
                  "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
                  "project": git_state(ROOT),
                  "fealpy": git_state(Path(fealpy.__file__).resolve().parents[1]),
                  "elapsed_seconds": elapsed, "backend": "numpy", "dtype": "float64"}
    config = {**vars(args), "output_dir": str(args.output_dir)}
    result = {"schema_version": SCHEMA_VERSION, "field_dictionary": FIELD_DICTIONARY,
              "config": config, "provenance": provenance,
              "reference": reference.tolist(), "cells": cells.tolist(),
              "edges": np.asarray(bm.to_numpy(mesh.entity("edge"))).tolist(),
              "boundary": boundary.tolist(), "scenes": scenes}
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(result, ensure_ascii=False, allow_nan=False)
    (output / "e_variational_mesh.json").write_text(encoded, encoding="utf-8")
    write_viewer(result, output)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    for row, (name, scene) in enumerate(scenes.items()):
        first, last = scene["frames"][0], scene["frames"][-1]
        for col, frame in enumerate((first, last)):
            nodes = np.asarray(frame["physical"])
            axes[row, col].triplot(nodes[:, 0], nodes[:, 1], cells, color="#087f82", lw=0.6)
            axes[row, col].set_aspect("equal")
            axes[row, col].set_title(f"{name}: {'initial' if col == 0 else 'final'}")
        history = scene["history"]
        axes[row, 2].plot([h["pseudo_time"] for h in history], [h["energy"] for h in history])
        axes[row, 2].set(xlabel="Pseudo-time", ylabel="Frozen-metric FEALPy energy")
        axes[row, 2].grid(alpha=0.25)
        summary = {"scene": name, "energy": [first["energy"], last["energy"]],
                   "interpolation_rmse": [first["interpolation_rmse"], last["interpolation_rmse"]],
                   "minimum_computational_area": min(h["computational_min_area"] for h in history),
                   "minimum_physical_area": min(h["physical_min_area"] for h in history),
                   "boundary_normal_displacement": last["boundary_normal_displacement"],
                   "boundary_tangent_displacement": last["boundary_tangent_displacement"],
                   "corner_displacement": last["corner_displacement"],
                   "alignment_rms": last["alignment_rms"],
                   "equidistribution_rms": last["equidistribution_rms"],
                   "saved_frames": len(scene["frames"]), "stop_reason": scene["stop_reason"]}
        print(json.dumps(summary, ensure_ascii=False))
    fig.savefig(output / "e_variational_mesh.png", dpi=160)
    plt.close(fig)
    print(f"Viewer: {output / 'e_variational_mesh.html'}")


if __name__ == "__main__":
    main()
