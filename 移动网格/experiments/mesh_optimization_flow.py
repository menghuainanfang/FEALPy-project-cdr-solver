"""Export complete converged mesh optimization trajectories and a compact viewer."""

import argparse
from datetime import datetime
from importlib.metadata import version
import json
from pathlib import Path
import numpy as np

from experiments.e_variational_mesh_visual import git_state
from src.mesh_optimization import run
from src.fealpy_mesh_evolution import quality_metrics
from src.poisson_mesh_comparison import PoissonCase

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/mentor_2026-10-23/optimization_flow"
LABELS = {"smooth": "扰动网格平滑", "band": "弯曲陡变带自适应", "peak": "局部尖峰自适应"}


def add_display_metrics(result):
    for key, case in result["cases"].items():
        cells, reference = np.array(case["cells"]), np.array(case["reference"])
        for row in case["rounds"]:
            for state in row["inner"]["accepted_step_trace"]:
                nodes = np.array(state["physical"])
                state.update(quality_metrics(nodes, reference, cells, PoissonCase(key).metric(nodes[cells].mean(axis=1))))
                state["energy_relative_change"]=state["energy"]/row["inner"]["accepted_step_trace"][0]["energy"]-1.


def write_viewer(result, output):
    template = Path(__file__).with_name("mesh_optimization_viewer.html").read_text(encoding="utf-8")
    script = Path(__file__).with_name("mesh_optimization_view.js").read_text(encoding="utf-8")
    encoded = json.dumps(result, ensure_ascii=False, allow_nan=False)
    (output / "index.html").write_text(template.replace("__FLOW_DATA__", encoded).replace(
        "__FLOW_SCRIPT__", script), encoding="utf-8")
    if output.resolve() == OUTPUT.resolve():
        (output.parent / "p1_complete/e_variational_mesh.html").write_text(
            '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta http-equiv="refresh" content="0;url=../optimization_flow/index.html">'
            '<title>切向移动 · 网格优化全流程</title><a href="../optimization_flow/index.html">打开网格优化全流程</a></html>', encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=8)
    parser.add_argument("--inner-cap", type=float, default=40.)
    parser.add_argument("--max-outer", type=int, default=60)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.render_only:
        result = json.loads((args.output_dir / "optimization_flow.json").read_text(encoding="utf-8"))
    else:
        cases = {}
        for key, label in LABELS.items():
            case = run(key, n=args.nx, inner_cap=args.inner_cap, max_outer=args.max_outer)
            case["label"] = label
            (args.output_dir / f"{key}.json").write_text(json.dumps(case, ensure_ascii=False, allow_nan=False), encoding="utf-8")
            if not case["outer_converged"]:
                raise RuntimeError(f"{key}: {case['stop_reason']}; full-flow demonstration requires outer convergence")
            cases[key] = case
        result = {"schema_version": "mesh-optimization-flow-v1", "cases": cases,
                  "provenance": {"created_at": datetime.now().astimezone().isoformat(),
                                 "project": git_state(ROOT), "fealpy": version("fealpy")}}
        (args.output_dir / "optimization_flow.json").write_text(
            json.dumps(result, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    add_display_metrics(result)
    (args.output_dir / "optimization_flow.json").write_text(json.dumps(result, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    for key, case in result["cases"].items():
        (args.output_dir / f"{key}.json").write_text(json.dumps(case, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    write_viewer(result, args.output_dir)


if __name__ == "__main__":
    main()
