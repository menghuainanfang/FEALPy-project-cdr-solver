"""Reproduce the original solve_ivp runs and expose BDF accepted steps."""

import csv
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import BDF
from scipy.sparse import coo_matrix

from src.fealpy_mesh_evolution import FealpyMeshEvolution
from fealpy.backend import backend_manager as bm

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/mentor_2026-10-23"


def trace(reference, physical, cells, metric, mode):
    model = FealpyMeshEvolution(physical, cells, metric, boundary_mode=mode)
    size = model.NN
    pattern = coo_matrix((np.ones(len(model.I)), (model.I, model.J)),
                         shape=(2*size, 2*size)).tocsr()

    def rhs(t, y):
        return model.energy_velocity(y.reshape(2, size).T)[1].ravel(order="F")

    started = perf_counter()
    solver = BDF(rhs, 0., reference.ravel(order="F"), .5,
                 rtol=1e-6, atol=1e-8, max_step=.05, jac_sparsity=pattern)
    rows = []
    while solver.status == "running":
        message = solver.step()
        if solver.status == "failed":
            raise RuntimeError(message)
        rows.append({"accepted_step_index": len(rows)+1, "pseudo_time": solver.t,
                     "t_old": solver.t_old, "accepted_dt": solver.step_size,
                     "nfev": solver.nfev, "njev": solver.njev, "nlu": solver.nlu})
    return rows, solver.y.reshape(2, size).T, perf_counter()-started


def main():
    bm.set_backend("numpy")
    summary = []
    for n in (8, 16):
        old = json.loads((OUTPUT / f"baseline_n{n}/e_variational_mesh.json").read_text(encoding="utf-8"))
        reference, cells = np.array(old["reference"]), np.array(old["cells"])
        for name, scene in old["scenes"].items():
            rows, last, elapsed = trace(reference, np.array(scene["fixed_physical"]),
                                        cells, np.array(scene["frozen_metric"]), scene["boundary_mode"])
            first, final = scene["frames"][0], scene["frames"][-1]
            row = {"n": n, "scene": name, "output_frames": len(scene["frames"]),
                   "t_eval_count": len(scene["history"]), "pseudo_time_end": final["pseudo_time"],
                   "accepted_bdf_steps": len(rows), **scene["solver_stats"],
                   "traced_wall_seconds": elapsed,
                   "original_three_scene_wall_seconds": old["provenance"]["elapsed_seconds"],
                   "velocity_initial": first["gradient_norm"], "velocity_final": final["gradient_norm"],
                   "energy_initial": first["energy"], "energy_final": final["energy"],
                   "min_computational_area": min(h["computational_min_area"] for h in scene["history"]),
                   "min_physical_area": min(h["physical_min_area"] for h in scene["history"]),
                   "stop_reason": scene["stop_reason"], "solver_success": True,
                   "trace_endpoint_max_difference": float(np.max(np.abs(last-np.array(final["computational"]))))}
            for counter in ("nfev", "njev", "nlu"):
                if rows[-1][counter] != row[counter]:
                    raise AssertionError(f"BDF trace differs from solve_ivp: {counter}")
            summary.append(row)
            (OUTPUT / f"baseline_n{n}/{name}_accepted_steps.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    (OUTPUT / "baseline_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (OUTPUT / "baseline_summary.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    for row in summary:
        print(row)


if __name__ == "__main__":
    main()
