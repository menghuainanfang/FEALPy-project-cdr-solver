import json
from pathlib import Path

import numpy as np
import pytest

from src.geometry import signed_triangle_areas
from src.mesh_optimization import run
from src.poisson_mesh_comparison import PoissonCase
from src.variational_mesh import reconstruct_physical

OUTPUT = Path(__file__).resolve().parents[1] / "results/mentor_2026-10-23/optimization_flow"


def test_live_driver_stops_only_after_outer_convergence_and_reports_budget():
    result = run("smooth", n=4)
    assert result["outer_converged"]
    assert result["stop_reason"] == "converged"
    assert result["rounds"][-1]["consecutive_success"] == 2
    limited = run("band", n=4, max_outer=1)
    assert not limited["outer_converged"]
    assert limited["stop_reason"] == "outer_iteration_cap"


@pytest.mark.parametrize("key", ["smooth", "band", "peak"])
def test_published_flow_preserves_all_inner_states_and_every_outer_map(key):
    result = json.loads((OUTPUT / f"{key}.json").read_text(encoding="utf-8"))
    reference, cells = np.array(result["reference"]), np.array(result["cells"])
    previous = np.array(result["initial_physical"])
    previous_metric = None
    consecutive = 0
    for k, row in enumerate(result["rounds"], 1):
        assert row["outer_iteration"] == k
        np.testing.assert_array_equal(row["input_physical"], previous)
        frozen = np.array(row["frozen_metric"])
        if previous_metric is not None:
            np.testing.assert_array_equal(frozen, previous_metric)
        inner = row["inner"]
        accepted, diagnostics = inner["accepted_step_trace"], inner["convergence_trace"]
        assert len(accepted) == len(diagnostics) == inner["solver_stats"]["accepted_bdf_steps"]+1
        assert inner["criterion_met"] and inner["stop_reason"] == "converged"
        np.testing.assert_allclose(accepted[0]["computational"], reference, atol=1e-14)
        for j, (state, diagnostic) in enumerate(zip(accepted, diagnostics)):
            assert state["accepted_step_index"] == diagnostic["accepted_step_index"] == j
            assert state["pseudo_time"] == diagnostic["mesh_pseudo_time"]
            assert state["scaled_residual_rms"] == diagnostic["mesh_residual_rms"]
            assert diagnostic["geometry_valid"]
            assert signed_triangle_areas(np.array(state["physical"]), cells).min()>0
            assert signed_triangle_areas(np.array(state["computational"]), cells).min()>0
        end = accepted[-1]
        mapped = reconstruct_physical(np.array(end["computational"]), previous, cells, reference)
        np.testing.assert_allclose(mapped, row["output_physical"], atol=1e-13)
        np.testing.assert_array_equal(end["physical"], row["output_physical"])
        physical = np.array(row["output_physical"])
        raw_metric = PoissonCase(key).metric(physical[cells].mean(axis=1))
        np.testing.assert_array_equal(raw_metric, row["sampled_metric"])
        delta_x = result["config"]["n"]*np.linalg.norm(physical-previous, axis=1).max()
        delta_m = np.linalg.norm(raw_metric-frozen)/max(np.linalg.norm(frozen), 1.)
        assert delta_x == row["mesh_change"]
        assert delta_m == row["metric_change"]
        cfg = result["outer_convergence"]
        met = delta_x<=cfg["mesh_tolerance"] and delta_m<=cfg["metric_tolerance"]
        consecutive = consecutive+1 if met else 0
        assert consecutive == row["consecutive_success"]
        assert row["outer_converged"] == (consecutive>=cfg["consecutive_rounds"])
        if k<len(result["rounds"]):
            assert not row["outer_converged"]
        previous = physical
        alpha = result["config"]["metric_relaxation"]
        previous_metric = (1-alpha)*frozen+alpha*raw_metric
        np.testing.assert_array_equal(row["updated_metric"], previous_metric)
    assert result["outer_converged"] and result["stop_reason"] == "converged"
