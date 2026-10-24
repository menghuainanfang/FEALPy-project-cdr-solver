import gzip
import json
from pathlib import Path

import numpy as np
import pytest
import src.burgers_rezone as rezone

from src.burgers_rezone import run as burgers, numerical_metric
from src.time_mesh_demo import run as analytic, target_field
from src.geometry import signed_triangle_areas

ROOT=Path(__file__).resolve().parents[1]/"results/mentor_2026-10-23/p2_v2"


def decoded(inner):
    if "convergence_trace" in inner:return inner["convergence_trace"]
    return [{**inner["trace_constants"],**dict(zip(inner["trace_columns"],row))} for row in inner["trace_rows"]]


def test_saved_coordinates_and_output_sampling_do_not_change_physical_solution():
    a=burgers(n=4,steps=2,inner_cap=10,max_outer=8,inner_frames=2,save_inner_steps=None)
    b=burgers(n=4,steps=2,inner_cap=10,max_outer=8,inner_frames=13,save_inner_steps=set())
    for x,y in zip(a["states"],b["states"]):
        np.testing.assert_array_equal(x["nodes"],y["nodes"])
        np.testing.assert_array_equal(x["values"],y["values"])
        assert x["outer_status"]==y["outer_status"]
        for k,z in zip(x["outer_iterations"],y["outer_iterations"]):
            assert k["mesh_inner"]["solver_stats"]==z["mesh_inner"]["solver_stats"]
            assert k["mesh_inner"]["first_converged_s"]==z["mesh_inner"]["first_converged_s"]
    assert a["states"][1]["outer_iterations"][0]["outer_status"]=="not_checked"
    assert a["states"][1]["outer_iterations"][0]["outer_residual"] is None
    limited=burgers(n=4,steps=10,inner_cap=10,max_outer=1,stop_on_outer_cap=True)
    assert len(limited["states"])==2
    assert limited["states"][-1]["outer_status"]=="not_checked"
    assert not limited["states"][-1]["overall_success"]


def test_static_limit_has_no_physical_advance():
    r=analytic(n=4,steps=10,static=True)
    assert len(r["states"])==2
    assert all(s["physical_step_index"]==0 and s["physical_time"]==0 for s in r["states"])
    for s in r["states"]:
        np.testing.assert_array_equal(s["values"],target_field(np.array(s["nodes"]),0.))
        assert s["pde"]["status"]=="not_applicable"


def test_pde_failure_never_records_a_successful_physical_step(monkeypatch):
    saved=[]
    def failed_step(nodes,cells,old_values,*args,**kwargs):
        return old_values,{"status":"nonlinear_iteration_cap","relative_nonlinear_residual":1.}
    monkeypatch.setattr(rezone,"solve_step",failed_step)
    with pytest.raises(RuntimeError,match="nonlinear solve did not converge"):
        rezone.run(n=4,steps=10,on_step=saved.append)
    assert len(saved)==1 and saved[0]["physical_step_index"]==0


@pytest.mark.parametrize("scenario,steps",[("analytic",100),("burgers",60)])
def test_long_measured_trajectory_inherits_mesh_and_has_independent_epochs(scenario,steps):
    path=ROOT/f"{scenario}_n8_steps{steps}_dt0.02/dataset.json"
    d=json.loads(path.read_text(encoding="utf-8"));cells=np.array(d["cells"])
    assert len(d["states"])==steps+1
    epochs=set()
    for index,(previous,s) in enumerate(zip(d["states"],d["states"][1:]),1):
        assert s["physical_step_index"]==index
        assert s["physical_time"]==pytest.approx(index*.02)
        np.testing.assert_array_equal(s["input_physical"],previous["nodes"])
        np.testing.assert_array_equal(s["outer_iterations"][0]["frozen_physical"],previous["nodes"])
        assert s["outer_status"]=="converged"
        assert signed_triangle_areas(np.array(s["nodes"]),cells).min()>0
        assert not np.array_equal(s["values"],previous["values"])
        if scenario=="burgers":
            assert s["pde"]["status"]=="converged"
            assert s["pde"]["relative_nonlinear_residual"]<1e-9
            np.testing.assert_allclose(s["metric"],numerical_metric(np.array(s["nodes"]),cells,np.array(s["values"])))
        for k in s["outer_iterations"]:
            inner=k["mesh_inner"];trace=decoded(inner)
            assert trace[0]["mesh_pseudo_time"]==0.
            assert inner["epoch_id"] not in epochs;epochs.add(inner["epoch_id"])
            assert trace[0]["epoch_id"]==trace[-1]["epoch_id"]==inner["epoch_id"]
            assert trace[0]["consecutive_success"]==0
            assert len(trace)==inner["solver_stats"]["accepted_bdf_steps"]+1
            assert trace[-1]["criterion_met"]
            for coordinates in inner["mesh_states"]:
                assert signed_triangle_areas(np.array(coordinates["physical"]),cells).min()>0
                assert signed_triangle_areas(np.array(coordinates["computational"]),cells).min()>0
    assert len(epochs)>steps


def test_compact_trace_decodes_to_actual_raw_records():
    folder=ROOT/"burgers_n8_steps30_dt0.02"
    compact=json.loads((folder/"dataset.json").read_text(encoding="utf-8"))["states"][1]
    raw=json.loads(gzip.decompress((folder/"steps/step_0001.json.gz").read_bytes()))
    for a,b in zip(compact["outer_iterations"],raw["outer_iterations"]):
        for actual,packed in zip(b["mesh_inner"]["convergence_trace"],decoded(a["mesh_inner"])):
            assert all(packed[key]==value for key,value in actual.items())
