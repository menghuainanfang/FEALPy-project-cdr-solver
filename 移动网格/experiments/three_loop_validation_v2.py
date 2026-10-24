"""Summarize measured long trajectories, comparison errors and storage costs."""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.burgers_rezone import run as burgers

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results/mentor_2026-10-23/p2_v2"


def main():
    summary=[];comparisons=[]
    for p in sorted(OUTPUT.glob("*/dataset.json")):
        d=json.loads(p.read_text(encoding="utf-8"));perf=json.loads(p.with_name("performance.json").read_text(encoding="utf-8"))
        rows=d["states"];continuity=max((float(np.max(np.abs(np.array(s["input_physical"])-np.array(old["nodes"])))) for old,s in zip(rows,rows[1:])),default=0.)
        summary.append({"run":p.parent.name,**perf,"continuity_max_difference":continuity,
                        "validation_role":"controlled_budget_exit" if p.parent.name=="burgers_n4_steps10_dt0.02_outer1" else "measured_run",
                        "minimum_physical_area":min(s["min_signed_area"] for s in rows),
                        "max_outer_iterations":max(len(s["outer_iterations"]) for s in rows),
                        "all_pde_converged":all(s.get("pde",{}).get("status")=="converged" for s in rows[1:]) if d["scenario"]=="burgers" else None})
        if d["scenario"]=="burgers" and perf["actual_physical_steps"]>=10:
            index=min(range(len(rows)),key=lambda i:abs(rows[i]["physical_time"]-.2));s=rows[index]
            comparisons.append({"n":d["config"]["n"],"dt":d["config"]["physical_dt"],"time":s["physical_time"],
                                "moving_l2":s["pde"]["l2_error"],"fixed_l2":d["fixed_grid_baseline"][index]["pde"]["l2_error"],
                                "pde_residual":s["pde"]["relative_nonlinear_residual"]})
    # Temporal consistency on the same fixed mesh, same initial field and same final time.
    temporal=[];values=[]
    for steps in (10,20,40):
        r=burgers(n=8,steps=steps,dt=.2/steps,moving=False,perturbation=.2)
        s=r["states"][-1];values.append(np.array(s["values"]))
        temporal.append({"n":8,"steps":steps,"dt":.2/steps,"time":.2,"l2_error":s["pde"]["l2_error"],"pde_residual":s["pde"]["relative_nonlinear_residual"]})
    differences=[float(np.linalg.norm(x-y)) for x,y in zip(values,values[1:])]
    comparisons=list({(r["n"],r["dt"],r["time"]):r for r in comparisons}.values())
    record={"runs":summary,"coupled_refinement_at_t02":comparisons,"fixed_mesh_temporal":temporal,
            "temporal_self_differences":differences,"temporal_self_difference_ratio":differences[0]/differences[1],
            "failure_records":[str(p.relative_to(OUTPUT)) for p in OUTPUT.glob("*/failure.json")]}
    (OUTPUT/"validation.json").write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding="utf-8")
    for name,rows in (("performance",summary),("refinement",comparisons)):
        with (OUTPUT/f"{name}.csv").open("w",newline="",encoding="utf-8-sig") as f:
            writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));writer.writeheader();writer.writerows(rows)
    print(json.dumps(record,ensure_ascii=False,indent=2))
    figures()


def figures():
    a=json.loads((OUTPUT/"analytic_n8_steps100_dt0.02/dataset.json").read_text(encoding="utf-8"))
    b=json.loads((OUTPUT/"burgers_n8_steps60_dt0.02/dataset.json").read_text(encoding="utf-8"))
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    t=[s["physical_time"] for s in a["states"]]
    axes[0,0].plot(t,[s["target_position"] for s in a["states"]],label="analytic interface center")
    axes[0,0].plot(t,[s["nodes"][36][0] for s in a["states"]],label="physical node 36 x")
    axes[0,0].set(ylabel="x position",title="Prescribed moving interface and inherited mesh")
    axes[0,1].plot(t,[s["target_interpolation_rmse"] for s in a["states"]],label="moving mesh")
    axes[0,1].plot(t,[s["fixed_grid_target_rmse"] for s in a["states"]],label="uniform reference mesh")
    axes[0,1].set(ylabel="sampled target RMSE",title="Analytic target interpolation; no PDE solve")
    t=[s["physical_time"] for s in b["states"][1:]]
    axes[1,0].plot(t,[s["pde"]["l2_error"] for s in b["states"][1:]],label="moving mesh Burgers")
    axes[1,0].plot(t,[s["pde"]["l2_error"] for s in b["fixed_grid_baseline"][1:]],label="same-initial fixed mesh")
    axes[1,0].set(ylabel="PDE L2 error",title="Manufactured Burgers PDE error")
    axes[1,1].plot(t,[s["transfer_mass_change"] for s in b["states"][1:]],label="P1 transfer mass difference")
    axes[1,1].axhline(0,color="gray",ls=":");axes[1,1].set(ylabel="mass change from transfer",title="Nonconservative rezoning diagnostic")
    for ax in axes.flat:ax.set_xlabel("physical time t");ax.grid(alpha=.25);ax.legend()
    fig.savefig(OUTPUT/"long_time_diagnostics.png",dpi=150);plt.close(fig)


if __name__=="__main__":main()
