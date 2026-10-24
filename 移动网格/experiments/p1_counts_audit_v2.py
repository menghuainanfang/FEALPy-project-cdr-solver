"""Trace displayed counts to archived fields and rerun the pre-v2 driver."""

import importlib.util
import json
from pathlib import Path
import sys

import fealpy
import numpy as np
import scipy
from fealpy.backend import backend_manager as bm

from experiments.e_variational_mesh_visual import git_state

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/mentor_2026-10-23"
OUTPUT = BASE / "p1_v2"


def main():
    bm.set_backend("numpy")
    spec = importlib.util.spec_from_file_location("p1_legacy_evolution", OUTPUT / "legacy_evolution.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    e = json.loads((BASE / "e_convergence/e_variational_mesh.json").read_text(encoding="utf-8"))
    h = json.loads((BASE / "burgers/h_burgers.json").read_text(encoding="utf-8"))
    records = []
    cases = []
    for name, scene in e["scenes"].items():
        cases.append((f"E_{name}", e, scene, scene["fixed_physical"], scene["frozen_metric"],
                      5., 61, scene["boundary_mode"], f"scenes.{name}", "e_convergence/e_variational_mesh.json"))
    for state_index in (1, 2):
        k = h["states"][state_index]["outer_iterations"][0]
        cases.append((f"H_t{state_index}_k1", h, k["mesh_inner"], k["frozen_physical"], k["frozen_metric"],
                      3., 8, "tangent", f"states[{state_index}].outer_iterations[0].mesh_inner", "burgers/h_burgers.json"))
    for name, data, archived, physical, metric, cap, frames, mode, field_path, filename in cases:
        reference, cells = np.array(data["reference"]), np.array(data["cells"])
        boundary = ((reference[:, 0] == 0) | (reference[:, 0] == 1) |
                    (reference[:, 1] == 0) | (reference[:, 1] == 1))
        result = module.solve(reference, np.array(physical), cells, boundary, np.array(metric),
                              duration=cap, frames=frames, boundary_mode=mode, save_inner_meshes=True)
        (OUTPUT / f"legacy_rerun_{name}.json").write_text(json.dumps(result, allow_nan=False), encoding="utf-8")
        config = archived["convergence"]
        first_residual_met = None
        first_all_met = None
        first_converged = None
        consecutive = 0
        for row in result["accepted_step_trace"][1:]:
            residual_met = row["scaled_residual_rms"] <= result["residual_threshold"]
            if residual_met and first_residual_met is None:
                first_residual_met = row["pseudo_time"]
            met = (residual_met and row["relative_displacement"] <= config["displacement_atol"] and
                   row["relative_energy_rate"] <= config["energy_rate_atol"] and
                   row["pseudo_time"] >= config["min_pseudo_time"])
            if met and first_all_met is None:
                first_all_met = row["pseudo_time"]
            consecutive = consecutive+1 if met else 0
            if consecutive >= config["consecutive_steps"] and first_converged is None:
                first_converged = row["pseudo_time"]
        steps = [row["accepted_dt"] for row in result["accepted_step_trace"][1:]]
        record = {"case": name, "source_file": filename, "json_path": field_path,
                  "archived_git": e["provenance"]["project"] if name.startswith("E") else h["provenance"]["project"],
                  "display_frame_count": len(result["frames"]), "bdf_accepted_steps": result["solver_stats"]["accepted_bdf_steps"],
                  "nfev": result["solver_stats"]["nfev"], "pseudo_time_cap": cap,
                  "final_s": result["terminal_pseudo_time"], "stop_reason": result["stop_reason"],
                  "first_residual_met_s": first_residual_met, "first_all_conditions_met_s": first_all_met,
                  "first_converged_s": first_converged, "consecutive_success_at_stop": consecutive,
                  "cap_reached": result["terminal_pseudo_time"] >= cap,
                  "residual_final": result["history"][-1]["scaled_residual_rms"],
                  "residual_threshold": result["residual_threshold"],
                  "accepted_dt_min": min(steps), "accepted_dt_max": max(steps), "accepted_dt_median": float(np.median(steps)),
                  "endpoint_max_difference": float(np.max(np.abs(np.array(result["frames"][-1]["computational"])-
                                                                np.array(archived["frames"][-1]["computational"])))),
                  "archived_counts_match": result["solver_stats"]["accepted_bdf_steps"] == archived["solver_stats"]["accepted_bdf_steps"]}
        records.append(record)
        print(json.dumps(record, ensure_ascii=False), flush=True)
    report = {"source_snapshot_commit": "ba0d1aa", "python": sys.executable, "python_version": sys.version,
              "numpy": np.__version__, "scipy": scipy.__version__, "project": git_state(ROOT),
              "fealpy_import": fealpy.__file__, "fealpy": git_state(Path(fealpy.__file__).resolve().parents[1]),
              "records": records}
    (OUTPUT / "counts_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
