"""Illustrate the P1 Galerkin increment recurrence for three Peclet numbers."""

from pathlib import Path
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run():
    output = Path(__file__).resolve().parent / "results_final"
    output.mkdir(exist_ok=True)
    n = 16
    h = 1 / n
    nodes = np.linspace(0, 1, n + 1)
    dense = np.linspace(0, 1, 4001)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    records = []
    for ax, eps in zip(axes, (0.1, 0.03125, 0.005)):
        lower, diagonal, upper = -eps / h - 0.5, 2 * eps / h, -eps / h + 0.5
        matrix = (
            np.diag(np.full(n - 1, diagonal))
            + np.diag(np.full(n - 2, lower), -1)
            + np.diag(np.full(n - 2, upper), 1)
        )
        rhs = np.zeros(n - 1)
        rhs[-1] = -upper
        values = np.r_[0.0, np.linalg.solve(matrix, rhs), 1.0]
        increments = np.diff(values)
        pe = h / (2 * eps)
        # The derived recurrence must hold even at Pe=1, where division fails.
        recurrence = (eps - h / 2) * increments[1:] - (eps + h / 2) * increments[:-1]
        assert np.max(np.abs(recurrence)) < 1e-13
        if pe <= 1:
            assert values.min() >= -1e-13 and values.max() <= 1 + 1e-13
        else:
            assert values.min() < -0.1
        exact = (np.exp((dense - 1) / eps) - np.exp(-1 / eps)) / (-np.expm1(-1 / eps))
        ax.plot(dense, exact, color="#292929", lw=2, label="Exact continuum solution")
        ax.plot(nodes, values, "o-", ms=4, color="#167d9a", label="P1 Galerkin")
        ax.set(title=f"Pe = {pe:g}   /   epsilon = {eps:g}", xlabel="x", ylabel="u")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8, loc="upper left")
        records.append(
            dict(
                epsilon=eps,
                n=n,
                Pe=pe,
                minimum=float(values.min()),
                maximum=float(values.max()),
                x=nodes.tolist(),
                u=values.tolist(),
            )
        )
    fig.suptitle("A unique monotone PDE solution can have oscillating discrete approximations")
    fig.savefig(output / "peclet_mechanism.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    (output / "peclet_mechanism.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    # Symmetric reaction eigenvalues used in the equation analysis.
    eigenvalues = np.linalg.eigvalsh([[2, -3 / 8], [-3 / 8, 3]])
    np.testing.assert_allclose(eigenvalues, [15 / 8, 25 / 8])
    print("PASS: increment recurrence, Peclet regimes and reaction eigenvalues")


if __name__ == "__main__":
    run()
