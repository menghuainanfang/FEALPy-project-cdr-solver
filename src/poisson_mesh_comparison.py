"""Manufactured Poisson problems and a shared FEALPy P1 solve."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy.sparse.linalg import spsolve
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from fealpy.fem import (BilinearForm, LinearForm, DirichletBC,
                        ScalarDiffusionIntegrator, ScalarSourceIntegrator)
from fealpy.functionspace import LagrangeFESpace
from fealpy.mesh import TriangleMesh

from src.geometry import signed_triangle_areas
from src.variational_mesh import band_metric


@dataclass(frozen=True)
class PoissonCase:
    name: str

    def fields(self, p):
        x, y = p[..., 0], p[..., 1]
        if self.name == "smooth":
            u = np.sin(np.pi * x) * np.sin(np.pi * y)
            grad = np.stack((np.pi * np.cos(np.pi * x) * np.sin(np.pi * y),
                             np.pi * np.sin(np.pi * x) * np.cos(np.pi * y)), axis=-1)
            return u, grad, 2 * np.pi ** 2 * u
        if self.name == "band":
            s = x - 0.5 - 0.15 * np.sin(2 * np.pi * y)
            sy = -0.3 * np.pi * np.cos(2 * np.pi * y)
            syy = 0.6 * np.pi ** 2 * np.sin(2 * np.pi * y)
            u = np.tanh(25 * s)
            first = 25 * (1 - u ** 2)
            second = -1250 * u * (1 - u ** 2)
            grad = np.stack((first, first * sy), axis=-1)
            return u, grad, -(second * (1 + sy ** 2) + first * syy)
        if self.name == "peak":
            delta = p - np.array([0.5, 0.5])
            radius_squared = np.sum(delta ** 2, axis=-1)
            u = np.exp(-150 * radius_squared)
            grad = -300 * delta * u[..., None]
            return u, grad, (600 - 90000 * radius_squared) * u
        if self.name == "affine":
            u = 1 + x - 2 * y
            return u, np.broadcast_to([1., -2.], p.shape), np.zeros_like(x)
        raise ValueError(f"Unknown problem: {self.name}")

    def metric(self, centers):
        if self.name == "band":
            return band_metric(centers)
        if self.name == "peak":
            delta = centers - np.array([0.5, 0.5])
            u = self.fields(centers)[0]
            hessian = (90000 * delta[:, :, None] * delta[:, None, :]
                       - 300 * np.eye(2)[None]) * u[:, None, None]
            eigenvalues, eigenvectors = np.linalg.eigh(hessian)
            absolute = (eigenvectors * np.abs(eigenvalues)[:, None, :]) @ eigenvectors.transpose(0, 2, 1)
            return np.eye(2)[None] + 30 * absolute / 300.
        return np.broadcast_to(np.eye(2), (len(centers), 2, 2)).copy()


def solve_poisson(nodes, cells, case, *, quadrature=12, error_quadrature=16):
    """Solve -Delta u=f, u=g on the entire unit-square boundary."""
    started = perf_counter()
    area = signed_triangle_areas(nodes, cells)
    if area.min() <= 0:
        raise ValueError("Poisson mesh has a nonpositive signed triangle area")
    mesh = TriangleMesh(bm.from_numpy(np.asarray(nodes)), bm.from_numpy(np.asarray(cells)))
    space = LagrangeFESpace(mesh, p=1)

    @cartesian
    def exact(p):
        return case.fields(p)[0]

    @cartesian
    def exact_gradient(p):
        return case.fields(p)[1]

    @cartesian
    def source(p):
        return case.fields(p)[2]

    @cartesian
    def zero(p):
        return np.zeros(p.shape[:-1])

    stiffness = BilinearForm(space)
    stiffness.add_integrator(ScalarDiffusionIntegrator(q=3, method=None))
    matrix = stiffness.assembly()
    load = LinearForm(space)
    load.add_integrator(ScalarSourceIntegrator(source, q=quadrature))
    rhs = load.assembly()
    matrix, rhs = DirichletBC(space, gd=exact).apply(matrix, rhs)
    assembled = perf_counter()
    scipy_matrix = matrix.to_scipy()
    values = np.asarray(spsolve(scipy_matrix, bm.to_numpy(rhs)))
    solved = perf_counter()
    uh = space.function()
    uh[:] = bm.from_numpy(values)
    boundary = np.asarray(bm.to_numpy(space.is_boundary_dof()), dtype=bool)
    boundary_error = float(np.max(np.abs(values[boundary] - exact(nodes[boundary]))))
    residual = scipy_matrix @ values - bm.to_numpy(rhs)
    relative_residual = float(np.linalg.norm(residual) / max(np.linalg.norm(bm.to_numpy(rhs)), 1.))
    l2 = float(mesh.error(exact, uh, q=error_quadrature))
    h1 = float(mesh.error(exact_gradient, uh.grad_value, q=error_quadrature))
    exact_l2 = float(mesh.error(exact, zero, q=error_quadrature))
    cell_l2 = np.asarray(mesh.error(exact, uh, q=error_quadrature, cell_axis=True))
    centers = np.asarray(nodes)[cells].mean(axis=1)
    center_values = values[cells].mean(axis=1)
    metrics = {"dofs": int(space.number_of_global_dofs()),
               "free_dofs": int(np.count_nonzero(~boundary)), "cells": len(cells),
               "l2_error": l2, "h1_seminorm_error": h1,
               "relative_l2_error": l2 / exact_l2,
               "boundary_error": boundary_error, "relative_residual": relative_residual,
               "assembly_seconds": assembled - started,
               "linear_solve_seconds": solved - assembled,
               "pde_seconds": solved - started,
               "diagnostic_seconds": perf_counter() - solved,
               "min_signed_area": float(area.min()), "total_area": float(area.sum())}
    return {"metrics": metrics, "nodes": np.asarray(nodes).tolist(), "values": values.tolist(),
            "center_values": center_values.tolist(),
            "center_absolute_error": np.abs(center_values - exact(centers)).tolist(),
            "cell_l2_error": cell_l2.tolist()}
