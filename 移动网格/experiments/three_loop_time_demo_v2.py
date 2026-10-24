"""Generate long physical-time mesh records with bounded detailed inner traces."""

import argparse
from copy import deepcopy
import csv
from datetime import datetime
import json
import gzip
import base64
from pathlib import Path
from time import perf_counter
import sys

import fealpy
from experiments.e_variational_mesh_visual import git_state
from src.burgers_rezone import run as burgers
from src.time_mesh_demo import KEY_STEPS, run as analytic
from src.mesh_trace import FIELD_DICTIONARY, render_diagnostic_assets

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results/mentor_2026-10-23/p2_v2"


def compact_trace(inner):
    trace=inner.pop("convergence_trace")
    keys=list(dict.fromkeys(k for row in trace for k in row))
    constants={key:trace[0].get(key) for key in keys if all(row.get(key)==trace[0].get(key) for row in trace)}
    columns=[key for key in keys if key not in constants]
    inner.update(trace_constants=constants,trace_columns=columns,trace_rows=[[row.get(key) for key in columns] for row in trace])


def pack_state(state):
    state=deepcopy(state)
    for iteration in state["outer_iterations"]:
        inner=iteration["mesh_inner"]
        inner["mesh_states"]=[{"accepted_step_index":r["accepted_step_index"],"pseudo_time":r["pseudo_time"],
                               "computational":r["computational"],"physical":r["physical"]}
                              for r in inner["accepted_step_trace"] if "physical" in r]
        for key in ("accepted_step_trace","frames","history","field_dictionary"):
            inner.pop(key,None)
        inner["coordinate_coverage"]="all_accepted" if len(inner["mesh_states"])==len(inner["convergence_trace"]) else "endpoints_only"
        compact_trace(inner)
    return state


def render_site():
    paths=list(OUTPUT.glob("*/dataset.json"))
    # For each scenario/resolution retain the longest measured trajectory in the page.
    models={}
    for path in paths:
        data=json.loads(path.read_text(encoding="utf-8"))
        key=f"{data['scenario']}_n{data['config']['n']}_dt{data['config']['physical_dt']}"
        if key not in models or len(data["states"])>len(models[key]["states"]):models[key]=data
    template=(ROOT/"experiments/time_evolution_viewer.html").read_text(encoding="utf-8")
    script=(ROOT/"experiments/time_evolution_view.js").read_text(encoding="utf-8")
    for model in models.values():
        blobs=[];story=[]
        for state in model["states"]:
            blobs.append(base64.b64encode(gzip.compress(json.dumps(state,ensure_ascii=False,allow_nan=False).encode("utf-8"),compresslevel=6)).decode("ascii"))
            light={**state,"outer_iterations":[]}
            for k in state["outer_iterations"]:
                light["outer_iterations"].append({key:({inner_key:value for inner_key,value in k["mesh_inner"].items()
                                                       if inner_key not in ("trace_rows","trace_columns","trace_constants","convergence_trace","mesh_states","endpoint_states")}
                                                       if key=="mesh_inner" else value)
                                                   for key,value in k.items() if key not in ("frozen_physical","frozen_metric","output_physical","provisional_values")})
            story.append(light)
        model["states"]=story;model["algorithm_blobs"]=blobs
    encoded=json.dumps({"models":models,"field_dictionary":FIELD_DICTIONARY},ensure_ascii=False,allow_nan=False)
    if len(encoded.encode("utf-8"))>80*1024**2:raise RuntimeError("Offline page dataset exceeds 80 MiB; reduce detailed inner traces")
    (OUTPUT/"index.html").write_text(render_diagnostic_assets(template).replace("__TIME_DATA__",encoded).replace("__TIME_SCRIPT__",script),encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario",choices=["analytic","burgers","static"],default="analytic")
    parser.add_argument("--physical-steps",type=int,choices=[10,30,60,100],default=60)
    parser.add_argument("--dt-physical",type=float,default=.02)
    parser.add_argument("--nx",type=int,default=8)
    parser.add_argument("--inner-cap",type=float)
    parser.add_argument("--max-outer",type=int)
    parser.add_argument("--save-inner-traces",choices=["key","all","none"],default="key")
    parser.add_argument("--render-only",action="store_true")
    args=parser.parse_args()
    if args.render_only:render_site();return
    if not 2<=args.nx<=32 or args.dt_physical<=0:parser.error("2<=nx<=32 and dt-physical>0 required")
    name=f"{args.scenario}_n{args.nx}_steps{0 if args.scenario=='static' else args.physical_steps}_dt{args.dt_physical}"
    folder=OUTPUT/name;(folder/"steps").mkdir(parents=True,exist_ok=True)
    rows=[];bytes_written=0;raw_bytes=0;start=perf_counter()
    def save_step(state):
        nonlocal bytes_written,raw_bytes
        index=len(rows);encoded=json.dumps(state,ensure_ascii=False,allow_nan=False)
        content=encoded.encode("utf-8");compressed=gzip.compress(content,compresslevel=6)
        (folder/"steps"/f"step_{index:04d}.json.gz").write_bytes(compressed)
        bytes_written+=len(compressed);raw_bytes+=len(content)
        if bytes_written>200*1024**2:raise RuntimeError("Run exceeded 200 MiB raw-state budget")
        iterations=state["outer_iterations"];inners=[k["mesh_inner"] for k in iterations]
        rows.append({"n":state["physical_step_index"],"time":state["physical_time"],"dt":state.get("physical_dt",0),
                     "k_actual":len(iterations),"outer_status":state.get("outer_status"),
                     "outer_residual":state.get("outer_residual"),"outer_tol":state.get("outer_tol"),
                     "pde_status":state.get("pde",{}).get("status","initial_interpolation"),
                     "pde_residual":state.get("pde",{}).get("relative_nonlinear_residual"),
                     "pde_l2_error":state.get("pde",{}).get("l2_error"),
                     "accepted_bdf_steps":sum(i["solver_stats"]["accepted_bdf_steps"] for i in inners),
                     "inner_all_converged":all(i["criterion_met"] for i in inners) if inners else None,
                     "minimum_area":state.get("min_signed_area"),"transfer_mass_change":state.get("transfer_mass_change"),
                     "target_rmse":state.get("target_interpolation_rmse"),"fixed_target_rmse":state.get("fixed_grid_target_rmse")})
        if index:print(f"{name} n={state['physical_step_index']} t={state['physical_time']:.3f} k={len(iterations)} outer={state.get('outer_status')}",flush=True)
    try:
        options=dict(n=args.nx,steps=args.physical_steps,dt=args.dt_physical,
                     inner_cap=args.inner_cap or (10. if args.scenario=="burgers" else 40.),
                     max_outer=args.max_outer or (8 if args.scenario=="burgers" else 30),on_step=save_step)
        if args.scenario=="burgers":
            result=burgers(**options,perturbation=.2,inner_frames=2,
                           save_inner_steps=None if args.save_inner_traces=="all" else KEY_STEPS if args.save_inner_traces=="key" else set())
            baseline=burgers(n=args.nx,steps=args.physical_steps,dt=args.dt_physical,moving=False,perturbation=.2)
            result["fixed_grid_baseline"]=[{"physical_time":s["physical_time"],"nodes":s["nodes"],"values":s["values"],"pde":s.get("pde")} for s in baseline["states"]]
            for row,s in zip(rows,baseline["states"]):row["fixed_pde_l2_error"]=s.get("pde",{}).get("l2_error")
        else:
            result=analytic(**options,static=args.scenario=="static",save_inner=args.save_inner_traces)
        result["states"]=[pack_state(s) for s in result["states"]]
        result.update(scenario=args.scenario,schema_version="three-loop-time-v2",field_dictionary=FIELD_DICTIONARY,
                      provenance={"created_at":datetime.now().astimezone().isoformat(),"project":git_state(ROOT),
                                  "fealpy":git_state(Path(fealpy.__file__).resolve().parents[1]),"python":sys.executable,
                                  "fealpy_import":fealpy.__file__,"command":sys.argv},
                      performance={"wall_seconds":perf_counter()-start,"compressed_state_bytes":bytes_written,"raw_bytes":raw_bytes,"rejected_bdf_attempts":None,
                                   "max_inner_steps":max((k["mesh_inner"]["solver_stats"]["accepted_bdf_steps"] for s in result["states"] for k in s["outer_iterations"]),default=0)})
        encoded=json.dumps(result,ensure_ascii=False,allow_nan=False)
        (folder/"dataset.json").write_text(encoded,encoding="utf-8")
        (folder/"performance.json").write_text(json.dumps({**result["performance"],"dataset_bytes":len(encoded.encode("utf-8")),"actual_physical_steps":max(s["physical_step_index"] for s in result["states"]),"all_outer_converged":all(s["outer_converged"] for s in result["states"][1:])},indent=2),encoding="utf-8")
    except Exception as error:
        (folder/"failure.json").write_text(json.dumps({"error":str(error),"completed_states":len(rows),"wall_seconds":perf_counter()-start},ensure_ascii=False,indent=2),encoding="utf-8")
        raise
    finally:
        if rows:
            with (folder/"states.csv").open("w",newline="",encoding="utf-8-sig") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for row in rows for k in row)));writer.writeheader();writer.writerows(rows)
    print(json.dumps(result["performance"],indent=2),flush=True)
    # Rendering is independent so long runs can finish before UI preparation.


if __name__=="__main__":main()
