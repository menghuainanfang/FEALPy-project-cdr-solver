"""H: physical-time Burgers, PDE/metric alternation and accepted mesh steps."""

import argparse
import csv
import json
from pathlib import Path
import sys
from datetime import datetime

import fealpy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

from experiments.e_variational_mesh_visual import git_state
from src.burgers_rezone import run, transfer

ROOT = Path(__file__).resolve().parents[1]


def baseline_refinement(output):
    spatial, temporal, values = [], [], []
    for n in (4, 8, 16):
        result = run(n=n, steps=3, moving=False)
        pde = result["states"][-1]["pde"]
        spatial.append({"n": n, "physical_dt": .02, "time": .06, "l2_error": pde["l2_error"],
                        "nonlinear_residual": pde["relative_nonlinear_residual"], "wall_seconds": result["wall_seconds"]})
    for steps in (3, 6, 12):
        result = run(n=16, steps=steps, dt=.06/steps, moving=False)
        values.append(np.array(result["states"][-1]["values"]))
        temporal.append({"n": 16, "steps": steps, "physical_dt": .06/steps,
                         "l2_error": result["states"][-1]["pde"]["l2_error"], "wall_seconds": result["wall_seconds"]})
    differences = [float(np.linalg.norm(a-b)) for a, b in zip(values, values[1:])]
    record = {"spatial": spatial, "temporal": temporal, "temporal_nodal_self_differences": differences,
              "self_difference_ratio": differences[0]/differences[1], "project": git_state(ROOT)}
    output.mkdir(parents=True, exist_ok=True)
    (output / "pde_baseline_refinement.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=8)
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument("--dt", type=float, default=.02)
    parser.add_argument("--inner-cap", type=float, default=3.)
    parser.add_argument("--max-outer", type=int, default=6)
    parser.add_argument("--perturbation", type=float, default=.2)
    parser.add_argument("--static", action="store_true")
    parser.add_argument("--baseline-refinement", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/mentor_2026-10-23/burgers")
    args = parser.parse_args()
    if args.baseline_refinement:
        baseline_refinement(args.output_dir)
        return
    result = run(n=args.nx, steps=args.steps, dt=args.dt, inner_cap=args.inner_cap,
                 max_outer=args.max_outer, perturbation=args.perturbation, static=args.static)
    result["fixed_grid_baseline"] = run(n=args.nx, steps=args.steps, dt=args.dt, moving=False,
                                         static=args.static, perturbation=args.perturbation)
    result["provenance"] = {"created_at": datetime.now().astimezone().isoformat(), "python": sys.executable,
                            "python_version": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
                            "project": git_state(ROOT), "fealpy": git_state(Path(fealpy.__file__).resolve().parents[1]),
                            "fealpy_import": fealpy.__file__, "command": [sys.executable, "-m", "experiments.h_burgers_physical_time", *sys.argv[1:]]}
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(result, allow_nan=False)
    (output / "h_burgers.json").write_text(encoded, encoding="utf-8")
    template = Path(__file__).with_name("physical_time_viewer.html").read_text(encoding="utf-8")
    import os
    (output / "h_burgers.html").write_text(template.replace("__EXPERIMENT_DATA__", encoded).replace(
        "__DEMO_LINK__", Path(os.path.relpath(ROOT / "demos/index.html", output)).as_posix()), encoding="utf-8")
    rows = []
    cells = np.array(result["cells"])
    fig, axes = plt.subplots(len(result["states"]), 3, figsize=(12, 3.5*len(result["states"])), squeeze=False, constrained_layout=True)
    for i, state in enumerate(result["states"]):
        nodes, values = np.array(state["nodes"]), np.array(state["values"])
        axes[i, 0].triplot(nodes[:, 0], nodes[:, 1], cells, lw=.6)
        field = axes[i, 1].tripcolor(nodes[:, 0], nodes[:, 1], cells, values, shading="gouraud")
        metric = np.trace(np.array(state["metric"]), axis1=1, axis2=2)
        axes[i, 2].tripcolor(nodes[:, 0], nodes[:, 1], cells, facecolors=metric)
        axes[i, 0].set_title(f"Physical mesh, t={state['physical_time']:.2f}")
        axes[i, 1].set_title("Initial data interpolation" if i == 0 else "Numerical Burgers u")
        axes[i, 2].set_title("trace(M(u)) on final mesh")
        fig.colorbar(field, ax=axes[i, 1])
        for ax in axes[i]:
            ax.set_aspect("equal")
        if i:
            previous = result["states"][i-1]
            rows.append({"physical_time": state["physical_time"], "outer_iterations": len(state["outer_iterations"]),
                         "mesh_inner_accepted_steps": sum(k["mesh_inner"]["solver_stats"]["accepted_bdf_steps"] for k in state["outer_iterations"]),
                         "pde_status": state["pde"]["status"], "l2_error": state["pde"]["l2_error"],
                         "fixed_grid_l2_error": result["fixed_grid_baseline"]["states"][i]["pde"]["l2_error"],
                         "relative_nonlinear_residual": state["pde"].get("relative_nonlinear_residual"),
                         "outer_stop_reason": state["outer_stop_reason"], "min_signed_area": state["min_signed_area"],
                         "transfer_mass_change": state["transfer_mass_change"],
                         "mesh_change_max": float(np.max(np.linalg.norm(nodes-np.array(previous["nodes"]), axis=1))),
                         "field_change_on_current_nodes": float(np.linalg.norm(values-transfer(np.array(previous["nodes"]), cells, np.array(previous["values"]), nodes))),
                         "metric_cellwise_change_norm": float(np.linalg.norm(np.array(state["metric"])-np.array(previous["metric"])))})
    fig.savefig(output / "physical_time.png", dpi=140)
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    for i in range(len(result["states"])):
        box = Bbox.union([ax.get_tightbbox(fig.canvas.get_renderer()) for ax in axes[i]])
        box = box.transformed(fig.dpi_scale_trans.inverted()).expanded(1.02, 1.05)
        fig.savefig(output / f"physical_time_step_{i}.png", dpi=140, bbox_inches=box)
    plt.close(fig)
    with (output / "summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"wall_seconds": result["wall_seconds"], "summary": rows}, indent=2))


if __name__ == "__main__":
    main()
