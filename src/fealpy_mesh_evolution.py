"""Frozen-metric evolution using FEALPy's MetricTensorAdaptive numerical core."""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix

from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh
from fealpy.mmesh import MetricTensorAdaptive
from fealpy.mmesh.metric_shared import GeometricDiscreteCore

from src.geometry import signed_triangle_areas
from src.variational_mesh import reconstruct_physical, interpolation_rmse


class FealpyMeshEvolution(MetricTensorAdaptive):
    """Two-dimensional unit-square driver on FEALPy's schema-based mesh.

    Reuses the library's energy, local derivative, assembly, balancing and
    tangential boundary velocity. Geometry is frozen for one coordinate solve.
    """

    def __init__(self, physical, cell, metric, *, boundary_mode="tangent", gamma=2., tau=.1):
        self.mesh = TriangleMesh(bm.asarray(physical), bm.asarray(cell))
        self.cell = self.mesh.cell
        self.node = self.mesh.node
        self.NN = self.mesh.number_of_nodes()
        self.NC = self.mesh.number_of_cells()
        self.GD = 2
        self.kwargs0 = bm.context(self.node)
        self.gamma = gamma
        self.tau = tau
        self.boundary_mode = boundary_mode
        self.cm = self.mesh.entity_measure("cell")
        self.geo_core = GeometricDiscreteCore(self.mesh)
        self.R = self.geo_core.R_matrix()
        self.I_p = bm.broadcast_to(bm.eye(2), (self.NC, 2, 2)).copy()
        self.boundary = self.mesh.boundary_node_flag()
        on_x = (bm.abs(self.node[:, 0]) < 1e-12) | (bm.abs(self.node[:, 0]-1) < 1e-12)
        on_y = (bm.abs(self.node[:, 1]) < 1e-12) | (bm.abs(self.node[:, 1]-1) < 1e-12)
        self.Vertices_idx = bm.where(on_x & on_y)[0]
        self.Bdinnernode_idx = bm.where(self.boundary & ~(on_x & on_y))[0]
        normal = bm.zeros((self.NN, 2), **self.kwargs0)
        normal = bm.set_at(normal, bm.where(on_x)[0], bm.array([1., 0.]))
        normal = bm.set_at(normal, bm.where(on_y & ~on_x)[0], bm.array([0., 1.]))
        self.Bi_Lnode_normal = normal[self.Bdinnernode_idx]
        self.BD_projector()
        self._build_jac_pattern()
        self.metric = bm.asarray(metric)
        self.M_node = bm.zeros((self.NN, 2, 2), **self.kwargs0)
        weights = bm.zeros(self.NN, **self.kwargs0)
        self.M_node = bm.index_add(self.M_node, self.cell,
                                  self.cm[:, None, None, None] * self.metric[:, None])
        weights = bm.index_add(weights, self.cell, bm.broadcast_to(self.cm[:, None], self.cell.shape))
        self.M_node /= weights[:, None, None]
        self.metric_inverse = bm.linalg.inv(self.metric)
        self.metric_theta = self.theta(self.metric)
        self._prepare_ivp_cache(self.node, self.metric, self.metric_theta)

    def vector_construction(self, A, g, trA, E_hat, return_local=False):
        result = super().vector_construction(A, g, trA, E_hat, return_local=return_local)
        velocity = result[0] if return_local else result
        if self.boundary_mode == "fixed":
            velocity = bm.set_at(velocity, self.boundary, 0.)
        return (velocity, result[1]) if return_local else velocity

    def quantities(self, computational):
        E_hat = self.edge_matrix(computational)
        A = self.A(self._ivp_cache["E_K"], E_hat, self.metric_inverse)
        return A, bm.linalg.det(A), bm.trace(A, axis1=1, axis2=2), E_hat

    def energy_velocity(self, computational):
        A, g, trA, E_hat = self.quantities(bm.asarray(computational))
        value = self.I_func(self.metric_theta, self.lam(self.metric_theta),
                            trA, g, self.rho(self.metric))
        velocity = self.vector_construction(A, g, trA, E_hat)
        return float(value), np.asarray(bm.to_numpy(velocity))


def boundary_diagnostics(nodes, original):
    """Displacements relative to each node's original unit-square side."""
    left = np.isclose(original[:, 0], 0.)
    right = np.isclose(original[:, 0], 1.)
    bottom = np.isclose(original[:, 1], 0.)
    top = np.isclose(original[:, 1], 1.)
    corners = (left | right) & (bottom | top)
    boundary = left | right | bottom | top
    delta = nodes - original
    normal_error = max(np.max(np.abs(delta[left | right, 0]), initial=0.),
                       np.max(np.abs(delta[bottom | top, 1]), initial=0.))
    return {"boundary_normal_displacement": float(normal_error),
            "boundary_tangent_displacement": float(np.max(np.linalg.norm(delta[boundary & ~corners], axis=1), initial=0.)),
            "corner_displacement": float(np.max(np.linalg.norm(delta[corners], axis=1), initial=0.))}


def quality_metrics(nodes, reference, cells, metric):
    """Metric quality relative to the given reference triangles; ideal is one."""
    physical_mesh = TriangleMesh(bm.asarray(nodes), bm.asarray(cells))
    core = GeometricDiscreteCore(physical_mesh)
    E = core.edge_matrix(bm.asarray(nodes))
    E_ref = core.edge_matrix(bm.asarray(reference))
    F = E @ bm.linalg.inv(E_ref)
    gram = bm.swapaxes(F, -1, -2) @ bm.asarray(metric) @ F
    alignment = bm.trace(gram, axis1=1, axis2=2) / (2 * bm.sqrt(bm.linalg.det(gram)))
    volume = physical_mesh.entity_measure("cell") * bm.sqrt(bm.linalg.det(bm.asarray(metric)))
    relative_volume = volume / bm.mean(volume)
    return {"alignment_rms": float(bm.sqrt(bm.mean(alignment ** 2))),
            "equidistribution_rms": float(bm.sqrt(bm.mean(relative_volume ** 2)))}


def solve(reference, physical, cells, boundary, metric, *, duration=.5, frames=61,
          boundary_mode="tangent", metric_function=None):
    """Evolve one frozen-coordinate problem with sparse finite-difference BDF."""
    model = FealpyMeshEvolution(physical, cells, metric, boundary_mode=boundary_mode)
    size = model.NN
    rows = np.asarray(bm.to_numpy(model.I))
    cols = np.asarray(bm.to_numpy(model.J))
    sparsity = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(2*size, 2*size)).tocsr()

    def rhs(t, y):
        return model.energy_velocity(y.reshape(2, size).T)[1].ravel(order="F")

    output_times = np.linspace(0., duration, frames)
    solution = solve_ivp(rhs, (0., duration), np.asarray(reference).ravel(order="F"),
                         method="BDF", t_eval=output_times, jac_sparsity=sparsity,
                         rtol=1e-6, atol=1e-8, max_step=duration/10)
    if not solution.success:
        raise RuntimeError(solution.message)
    samples_1d = (np.arange(100) + .5) / 100
    xx, yy = np.meshgrid(samples_1d, samples_1d)
    samples = np.column_stack((xx.ravel(), yy.ravel()))
    history, saved = [], []
    for i, (time, y) in enumerate(zip(solution.t, solution.y.T)):
        computational = y.reshape(2, size).T.copy()
        on_x = np.isclose(reference[:, 0], 0.) | np.isclose(reference[:, 0], 1.)
        on_y = np.isclose(reference[:, 1], 0.) | np.isclose(reference[:, 1], 1.)
        computational[on_x, 0] = reference[on_x, 0]
        computational[on_y, 1] = reference[on_y, 1]
        if boundary_mode == "fixed":
            computational[boundary] = reference[boundary]
        moved = reconstruct_physical(computational, physical, cells, reference)
        computational_area = signed_triangle_areas(computational, cells)
        physical_area = signed_triangle_areas(moved, cells)
        if min(computational_area.min(), physical_area.min()) <= 0:
            raise ValueError("FEALPy evolution produced a nonpositive displayed triangle")
        energy, velocity = model.energy_velocity(computational)
        actual_metric = metric_function(moved[cells].mean(axis=1)) if metric_function else metric
        row = {"iteration": i, "pseudo_time": float(time), "energy": energy,
               "computational_min_area": float(computational_area.min()),
               "physical_min_area": float(physical_area.min()),
               "gradient_norm": float(np.linalg.norm(velocity)),
               "step": float(time - solution.t[i-1]) if i else 0., "backtracks": 0,
               **boundary_diagnostics(moved, physical),
               **quality_metrics(moved, reference, cells, actual_metric)}
        history.append(row)
        saved.append({**row, "computational": computational.tolist(), "physical": moved.tolist(),
                      "interpolation_rmse": interpolation_rmse(moved, cells, samples)})
    return {"history": history, "frames": saved, "stop_reason": "pseudo_time_end",
            "boundary_mode": boundary_mode, "integrator": "SciPy BDF",
            "algorithm": "FEALPy MetricTensorAdaptive",
            "solver_stats": {"nfev": solution.nfev, "njev": solution.njev, "nlu": solution.nlu},
            "rmse_sampling": "100x100 midpoint interpolation samples, not PDE error"}
