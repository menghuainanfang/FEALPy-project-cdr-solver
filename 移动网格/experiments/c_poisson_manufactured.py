"""C: P1 FEALPy solution of a manufactured Poisson problem."""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import spsolve

from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from fealpy.fem import BilinearForm, DirichletBC, LinearForm
from fealpy.fem import ScalarDiffusionIntegrator, ScalarSourceIntegrator
from fealpy.functionspace import LagrangeFESpace
from fealpy.mesh import TriangleMesh


ROOT = Path(__file__).resolve().parents[1]


@cartesian
def solution(p):
    x, y = p[..., 0], p[..., 1]
    return bm.sin(bm.pi * x) * bm.sin(bm.pi * y)


@cartesian
def source(p):
    return 2.0 * bm.pi**2 * solution(p)


@cartesian
def gradient(p):
    x, y = p[..., 0], p[..., 1]
    return bm.stack(
        [bm.pi * bm.cos(bm.pi * x) * bm.sin(bm.pi * y),
         bm.pi * bm.sin(bm.pi * x) * bm.cos(bm.pi * y)],
        axis=-1,
    )


def solve(n: int) -> dict[str, float | int]:
    started = time.perf_counter()
    mesh = TriangleMesh.from_box([0.0, 1.0, 0.0, 1.0], nx=n, ny=n)
    space = LagrangeFESpace(mesh, p=1)

    bform = BilinearForm(space)
    bform.add_integrator(ScalarDiffusionIntegrator(q=4, method=None))
    matrix = bform.assembly()

    lform = LinearForm(space)
    lform.add_integrator(ScalarSourceIntegrator(source, q=5))
    rhs = lform.assembly()
    matrix, rhs = DirichletBC(space, gd=solution).apply(matrix, rhs)

    scipy_matrix = matrix.to_scipy()
    rhs_numpy = bm.to_numpy(rhs)
    values = spsolve(scipy_matrix, rhs_numpy)
    uh = space.function()
    uh[:] = bm.from_numpy(np.asarray(values))

    residual = scipy_matrix @ values - rhs_numpy
    boundary = bm.to_numpy(space.is_boundary_dof())
    exact_boundary = bm.to_numpy(solution(mesh.entity("node")[boundary]))
    boundary_error = float(np.max(np.abs(values[boundary] - exact_boundary)))
    l2_error = float(mesh.error(solution, uh, q=5))
    h1_error = float(mesh.error(gradient, uh.grad_value, q=5))

    return {
        "n": n,
        "cells": int(mesh.number_of_cells()),
        "dofs": int(space.number_of_global_dofs()),
        "l2_error": l2_error,
        "h1_seminorm_error": h1_error,
        "boundary_max_error": boundary_error,
        "relative_linear_residual": float(np.linalg.norm(residual) / np.linalg.norm(rhs_numpy)),
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    bm.set_backend("numpy")
    rows = [solve(n) for n in (8, 16, 32, 64)]
    for i, row in enumerate(rows):
        if i == 0:
            row["l2_order"] = None
            row["h1_order"] = None
        else:
            row["l2_order"] = float(np.log2(rows[i - 1]["l2_error"] / row["l2_error"]))
            row["h1_order"] = float(np.log2(rows[i - 1]["h1_seminorm_error"] / row["h1_seminorm_error"]))

    assert max(row["boundary_max_error"] for row in rows) < 1e-12
    assert max(row["relative_linear_residual"] for row in rows) < 1e-10
    assert min(row["l2_order"] for row in rows[1:]) > 1.8
    assert min(row["h1_order"] for row in rows[1:]) > 0.9

    with (ROOT / "results" / "c_poisson_convergence.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    (ROOT / "results" / "c_poisson_convergence.json").write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

