"""Build a recorded case library with complete accepted-state replay."""

import argparse
from datetime import datetime
from importlib.metadata import version
import json
from pathlib import Path
from time import perf_counter

import fealpy
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh
import numpy as np

from experiments.e_variational_mesh_visual import ROOT, git_state, write_viewer
from src.fealpy_mesh_evolution import solve, quality_metrics
from src.mesh_trace import FIELD_DICTIONARY, SCHEMA_VERSION
from src.variational_mesh import band_metric, interpolation_rmse

OUTPUT = ROOT / "results/mentor_2026-10-23/p1_complete"


def replay_frames(scene, reference, cells, metric_function, metric):
    """Evaluate display quantities on the recorded accepted meshes themselves."""
    points_1d = (np.arange(100) + .5) / 100
    xx, yy = np.meshgrid(points_1d, points_1d)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    frames = []
    for j, row in enumerate(scene["accepted_step_trace"]):
        nodes = np.asarray(row["physical"])
        actual_metric = metric_function(nodes[cells].mean(axis=1)) if metric_function else metric
        frames.append({**row, "output_frame_index": j, "display_frame_index": j,
                       "sample_dt": row["accepted_dt"], "output_sample_dt": row["accepted_dt"],
                       **quality_metrics(nodes, reference, cells, actual_metric),
                       "interpolation_rmse": interpolation_rmse(nodes, cells, points)})
    scene["sampled_output_history"] = scene["history"]
    scene["frames"] = frames
    scene["history"] = [{k: v for k, v in f.items() if k not in ("physical", "computational")}
                        for f in frames]
    scene["saved_frames"] = len(frames)
    scene["replay_mode"] = "accepted_states"
    scene["replay_definition"] = "One displayed mesh per accepted BDF state, including j=0; sampled_output_history retains requested output diagnostics."
    compact_replay_storage(scene)


def compact_replay_storage(scene):
    """Store each accepted mesh once, in the display frame with the same j."""
    scene["accepted_step_trace"] = [{k: v for k, v in row.items() if k not in ("physical", "computational")}
                                    for row in scene["accepted_step_trace"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--long-cap", type=float, default=20.)
    args = parser.parse_args()
    if args.long_cap <= .5:
        parser.error("long-cap must exceed the short cap 0.5")
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box([0., 1., 0., 1.], nx=16, ny=16)
    reference = np.asarray(bm.to_numpy(mesh.entity("node"))).copy()
    cells = np.asarray(bm.to_numpy(mesh.entity("cell")), dtype=np.int64)
    boundary = np.asarray(bm.to_numpy(mesh.boundary_node_flag()), dtype=bool)
    perturbed = reference.copy()
    perturbed[~boundary] += np.random.default_rng(42).uniform(-.2/16, .2/16, (sum(~boundary), 2))
    strong = reference.copy()
    strong[~boundary] += np.random.default_rng(42).uniform(-.3/16, .3/16, (sum(~boundary), 2))
    identity = np.broadcast_to(np.eye(2), (len(cells), 2, 2)).copy()
    band = band_metric(reference[cells].mean(axis=1))
    cases = []
    for key, label, physical, metric, mode in (
            ("adaptive", "陡变带 · 切向边界", reference, band, "tangent"),
            ("adaptive_fixed", "陡变带 · 固定边界", reference, band, "fixed"),
            ("smoothing", "扰动平滑", perturbed, identity, "tangent")):
        cases.append((key, label + " · 长时间", physical, metric, mode, args.long_cap, True, 2000,
                      "延长伪时间上限，观察残量、固定窗口和连续检查何时共同达标。"))
        cases.append((key + "_short", label + " · 短时截断", physical, metric, mode, .5, True, 2000,
                      "与对应长时间案例使用相同初值和度量，展示上限到达但判据尚未满足。"))
    cases.extend([
        ("equilibrium", "初始平衡 · 仍需窗口检查", reference, identity, "tangent", 2., True, 2000,
         "初始残量接近零，也要等待固定窗口就绪并连续达标；j=0不直接判收敛。"),
        ("strong_smoothing", "较强扰动 · 平滑恢复", strong, identity, "tangent", args.long_cap, True, 2000,
         "内部节点扰动幅度由0.2/n增加至0.3/n，边界与角点保持约束。"),
        ("smoothing_continue", "平滑达标后 · 继续积分", perturbed, identity, "tangent", 8., False, 2000,
         "明确关闭达标早停，比较首次达标与继续积分后的能量；最终按时间上限结束。"),
        ("adaptive_budget", "计算预算 · 12个接受步", reference, band, "tangent", args.long_cap, True, 12,
         "实际在第12个接受步结束，展示计算预算耗尽与时间上限、网格收敛的区别。"),
    ])
    scenes, summary = {}, []
    started = perf_counter()
    for key, label, physical, metric, mode, cap, early_stop, budget, description in cases:
        metric_function = band_metric if key.startswith("adaptive") else None
        scene = solve(reference, physical, cells, boundary, metric, duration=cap, frames=121,
                      boundary_mode=mode, metric_function=metric_function,
                      save_inner_meshes=True, stop_on_convergence=early_stop,
                      max_internal_steps=budget, epoch_id=f"complete_{key}",
                      metric_provenance="prescribed band(strength=30)" if metric_function else "identity")
        scene.update(label=label, description=description, fixed_physical=physical.tolist(),
                     frozen_metric=metric.tolist(), case_key=key)
        replay_frames(scene, reference, cells, metric_function, metric)
        scenes[key] = scene
        summary.append({"case": key, "label": label, "cap": cap, "stop_reason": scene["stop_reason"],
                        "terminal_s": scene["terminal_pseudo_time"], "first_converged_s": scene["first_converged_s"],
                        "criterion_met": scene["criterion_met"], "cap_reached": scene["cap_reached"],
                        "accepted_bdf_steps": scene["solver_stats"]["accepted_bdf_steps"],
                        "replay_frames": len(scene["frames"]),
                        "minimum_physical_area": min(r["physical_min_area"] for r in scene["convergence_trace"])})
        print(json.dumps(summary[-1], ensure_ascii=False), flush=True)
    result = {"schema_version": SCHEMA_VERSION, "field_dictionary": FIELD_DICTIONARY,
              "config": {"seed": 42, "nx": 16, "long_cap": args.long_cap},
              "provenance": {"created_at": datetime.now().astimezone().isoformat(),
                             "fealpy_version": version("fealpy"),
                             "fealpy_import": fealpy.__file__, "project": git_state(ROOT),
                             "elapsed_seconds": perf_counter()-started},
              "reference": reference.tolist(), "cells": cells.tolist(),
              "edges": np.asarray(bm.to_numpy(mesh.entity("edge"))).tolist(),
              "boundary": boundary.tolist(), "scenes": scenes,
              "demo_library": True, "summary": summary}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "e_variational_mesh.json").write_text(
        json.dumps(result, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    (args.output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_viewer(result, args.output_dir)


if __name__ == "__main__":
    main()
