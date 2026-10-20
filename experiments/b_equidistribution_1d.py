"""B: one-dimensional equidistribution teaching experiment (E0)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad


ROOT = Path(__file__).resolve().parents[1]
ALPHA = 9.0
BETA = 100.0
N = 40


def monitor(x: np.ndarray | float) -> np.ndarray | float:
    return 1.0 + ALPHA * np.exp(-BETA * (np.asarray(x) - 0.5) ** 2)


def target(x: np.ndarray) -> np.ndarray:
    return np.exp(-BETA * (x - 0.5) ** 2)


def interpolation_errors(nodes: np.ndarray) -> tuple[float, float]:
    dense = np.linspace(0.0, 1.0, 20001)
    exact = target(dense)
    interpolated = np.interp(dense, nodes, target(nodes))
    error = exact - interpolated
    return float(np.sqrt(np.trapezoid(error * error, dense))), float(np.max(np.abs(error)))


def main() -> None:
    dense = np.linspace(0.0, 1.0, 200001)
    cumulative = np.concatenate(([0.0], cumulative_trapezoid(monitor(dense), dense)))
    cumulative /= cumulative[-1]
    adapted = np.interp(np.linspace(0.0, 1.0, N + 1), cumulative, dense)
    uniform = np.linspace(0.0, 1.0, N + 1)

    masses = np.array([quad(monitor, adapted[i], adapted[i + 1], epsabs=1e-13)[0] for i in range(N)])
    mean_mass = float(masses.mean())
    relative_deviation = np.abs(masses / mean_mass - 1.0)
    uniform_l2, uniform_linf = interpolation_errors(uniform)
    adapted_l2, adapted_linf = interpolation_errors(adapted)

    result = {
        "alpha": ALPHA,
        "beta": BETA,
        "segments": N,
        "strictly_increasing": bool(np.all(np.diff(adapted) > 0.0)),
        "left_endpoint": float(adapted[0]),
        "right_endpoint": float(adapted[-1]),
        "max_monitor_mass_relative_deviation": float(relative_deviation.max()),
        "uniform_l2_interpolation_error": uniform_l2,
        "adapted_l2_interpolation_error": adapted_l2,
        "l2_error_ratio_adapted_over_uniform": adapted_l2 / uniform_l2,
        "uniform_linf_interpolation_error": uniform_linf,
        "adapted_linf_interpolation_error": adapted_linf,
    }
    assert result["strictly_increasing"]
    assert adapted[0] == 0.0 and adapted[-1] == 1.0
    assert result["max_monitor_mass_relative_deviation"] < 1e-6
    assert adapted_l2 < uniform_l2

    with (ROOT / "results" / "b_equidistribution_segments.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["segment", "x_left", "x_right", "monitor_mass", "relative_deviation"])
        for i, mass in enumerate(masses):
            writer.writerow([i, adapted[i], adapted[i + 1], mass, relative_deviation[i]])

    (ROOT / "results" / "b_equidistribution.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.2), constrained_layout=True)
    xplot = np.linspace(0.0, 1.0, 1000)
    axes[0].plot(xplot, monitor(xplot), color="black", label=r"$\rho(x)$")
    axes[0].plot(adapted, np.zeros_like(adapted), "|", ms=10, label="equidistributed nodes")
    axes[0].set_ylabel("monitor")
    axes[0].legend()
    axes[1].plot(np.arange(N), relative_deviation, "o-", ms=3)
    axes[1].set_xlabel("segment")
    axes[1].set_ylabel("relative mass deviation")
    fig.savefig(ROOT / "results" / "b_equidistribution.png", dpi=180)
    plt.close(fig)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

