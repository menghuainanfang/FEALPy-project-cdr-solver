"""JSON metadata and shared field definitions for mesh diagnostic viewers."""

import json
from pathlib import Path

FIELD_DICTIONARY = json.loads(Path(__file__).with_name("mesh_field_dictionary.json").read_text(encoding="utf-8"))
SCHEMA_VERSION = "mesh-evolution-v2.0"


def render_diagnostic_assets(template):
    script = Path(__file__).resolve().parents[1].joinpath("experiments/mesh_diagnostics_view.js").read_text(encoding="utf-8")
    return template.replace("__MESH_DIAGNOSTICS_JS__", "globalThis.MESH_FIELD_DICTIONARY="+
                            json.dumps(FIELD_DICTIONARY, ensure_ascii=False)+";\n"+script)


def epoch_metadata(monitor, convergence_trace, *, free_scalar_dofs, tau, local_h, frames, wall_seconds):
    terminal = convergence_trace[-1]
    return {"schema_version": SCHEMA_VERSION, "field_dictionary": FIELD_DICTIONARY,
            "epoch_id": terminal["epoch_id"], "criterion_version": terminal["criterion_version"],
            "convergence_trace": convergence_trace, "convergence_definition": monitor.parameters(),
            "criterion_met": terminal["criterion_met"], "cap_reached": terminal["cap_reached"],
            "solver_status": terminal["solver_status"], "mesh_status": terminal["mesh_status"],
            "first_converged_s": terminal["first_converged_s"],
            "last_unmet_reasons": terminal["unmet_reasons"], "saved_frames": len(frames),
            "wall_clock_seconds": wall_seconds,
            "norm_definition": {"free_scalar_dofs": free_scalar_dofs, "residual_reference": monitor.r_ref,
                                "tau": tau, "h_definition": "minimum incident reference edge length per node",
                                "h_min": float(min(local_h)), "h_max": float(max(local_h)),
                                "r_ref_definition": "initial R_rms of this frozen metric epoch",
                                "residual_kind": "projected balanced velocity residual, not raw gradient"}}
