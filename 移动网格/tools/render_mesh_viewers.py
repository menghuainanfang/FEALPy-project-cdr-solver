"""Render E/H numerical JSON with the common field-explanation panels."""

import json
import os
from pathlib import Path

from experiments.e_variational_mesh_visual import write_viewer
from experiments.f_poisson_mesh_comparison import write_viewer as write_poisson_viewer
from src.mesh_trace import render_diagnostic_assets

ROOT = Path(__file__).resolve().parents[1]


def main():
    write_poisson_viewer(json.loads((ROOT / "results/f_poisson_comparison.json").read_text(encoding="utf-8")), ROOT / "results")
    for relative in ("results", "results/mentor_2026-10-23/e_convergence", "results/mentor_2026-10-23/p1_v2/e"):
        folder = ROOT / relative
        data = json.loads((folder / "e_variational_mesh.json").read_text(encoding="utf-8"))
        write_viewer(data, folder)
    template = (ROOT / "experiments/physical_time_viewer.html").read_text(encoding="utf-8")
    for relative in ("results/mentor_2026-10-23/burgers", "results/mentor_2026-10-23/static",
                     "results/mentor_2026-10-23/p1_v2/h_preview", "results/mentor_2026-10-23/p1_v2/static"):
        folder = ROOT / relative
        encoded = (folder / "h_burgers.json").read_text(encoding="utf-8")
        (folder / "h_burgers.html").write_text(render_diagnostic_assets(template).replace(
            "__EXPERIMENT_DATA__", encoded).replace("__DEMO_LINK__", Path(os.path.relpath(
                ROOT / "demos/index.html", folder)).as_posix()), encoding="utf-8")
    print("Rendered E/H explanations; numerical JSON files retained.")
    from experiments.mesh_optimization_flow import OUTPUT, write_viewer as write_flow_viewer
    write_flow_viewer(json.loads((OUTPUT / "optimization_flow.json").read_text(encoding="utf-8")), OUTPUT)


if __name__ == "__main__":
    main()
