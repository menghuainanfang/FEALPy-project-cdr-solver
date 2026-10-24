"""Plot accepted-state residual, frozen energy and real BDF step lengths."""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/mentor_2026-10-23/p1_v2/e"


def main():
    result = json.loads((OUTPUT / "e_variational_mesh.json").read_text(encoding="utf-8"))
    scenes = {k: s for k, s in result["scenes"].items() if s["boundary_mode"] == "tangent"}
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), constrained_layout=True)
    for i, (name, scene) in enumerate(scenes.items()):
        trace = scene["convergence_trace"]
        with (OUTPUT / f"{name}_convergence_trace.csv").open("w", newline="", encoding="utf-8-sig") as f:
            flat = [{k: json.dumps(v, ensure_ascii=False) if isinstance(v, list) else v for k, v in row.items()} for row in trace]
            writer = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in flat for k in r)))
            writer.writeheader()
            writer.writerows(flat)
        times = [r["mesh_pseudo_time"] for r in trace]
        axes[i, 0].semilogy(times, [max(r["mesh_residual_rms"], 1e-16) for r in trace])
        axes[i, 0].axhline(scene["residual_threshold"], color="orange", ls="--", label="residual threshold")
        axes[i, 1].plot(times, [r["mesh_energy"] for r in trace])
        axes[i, 2].semilogy([r["accepted_step_index"] for r in trace[1:]], [r["pseudo_step_size"] for r in trace[1:]])
        if scene["first_converged_s"] is not None:
            for ax in axes[i, :2]:
                ax.axvline(scene["first_converged_s"], color="orange", ls=":")
        axes[i, 0].set(title=f"{name}: {scene['stop_reason']}", xlabel="Mesh pseudo-time s", ylabel="R_rms")
        axes[i, 1].set(title=f"Frozen epoch: {scene['epoch_id']}", xlabel="Mesh pseudo-time s", ylabel="I_h")
        axes[i, 2].set(title=f"criterion={scene['criterion_met']}, cap={scene['cap_reached']}", xlabel="Accepted BDF step j", ylabel="Accepted step size delta_s")
        for ax in axes[i]:
            ax.grid(alpha=.25)
    fig.savefig(OUTPUT / "convergence_diagnostics.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
