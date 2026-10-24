"""Record whether the installed upstream Burgers moving-mesh API initializes."""

import importlib.util
import json
from pathlib import Path
import traceback

import fealpy
from fealpy.backend import backend_manager as bm
from experiments.e_variational_mesh_visual import git_state


def main():
    bm.set_backend("numpy")
    root = Path(fealpy.__file__).resolve().parents[1]
    path = root / "example/meshopt/mmesh/scalar_burgers_mm_solver.py"
    record = {"source": str(path), "fealpy": git_state(root), "phases": []}
    try:
        spec = importlib.util.spec_from_file_location("upstream_burgers", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        record["phases"].append("module import succeeded")
        pde = module.ScalarBurgersData("exp(-t)*sin(pi*x)*sin(pi*y)", ["x", "y", "t"], T=[0., .02], Re=20.)
        pde.set_mesh(nx=4, ny=4)
        record["phases"].append("ScalarBurgersData.set_mesh succeeded")
        solver = module.Burgers_MMsolver(pde, nt=1, method="MetricTensorAdaptive", plot=False, mm_steps=1)
        solver.linear_system()
        record["phases"].append("linear_system succeeded")
        solver.moving_mesher()
        record["phases"].append("moving_mesher succeeded")
        record["status"] = "initialized_not_time_stepped"
    except Exception as error:
        record.update(status="initialization_failed", exception=repr(error), traceback=traceback.format_exc())
    output = Path(__file__).resolve().parents[1] / "results/mentor_2026-10-23/upstream_probe.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
