"""Fixed-physical-mesh Winslow discretization in computational coordinates."""

from __future__ import annotations

import numpy as np
from matplotlib.tri import LinearTriInterpolator, Triangulation

from src.geometry import signed_triangle_areas


def edge_matrices(nodes, cells):
    triangles = np.asarray(nodes)[cells]
    return np.stack((triangles[:, 1] - triangles[:, 0],
                     triangles[:, 2] - triangles[:, 0]), axis=-1)


class WinslowEnergy:
    """I_h = sum |K| tr(J M_K^-1 J^T), J = hat(E) E^-1.

    Physical vertices, physical areas, and cell metrics stay fixed throughout
    each inner solve. Gradients are with respect to computational vertices.
    """

    def __init__(self, physical, cells, metric):
        self.cells = np.asarray(cells, dtype=np.int64)
        self.area = signed_triangle_areas(physical, self.cells)
        if np.any(self.area <= 0):
            raise ValueError("The fixed physical mesh must be positively oriented.")
        inverse = np.linalg.inv(edge_matrices(physical, self.cells))
        self.coefficient = inverse @ np.linalg.inv(metric) @ inverse.transpose(0, 2, 1)

    def value_gradient(self, computational):
        edges = edge_matrices(computational, self.cells)
        product = edges @ self.coefficient
        value = float(np.sum(self.area * np.sum(edges * product, axis=(1, 2))))
        derivative = 2 * self.area[:, None, None] * product
        gradient = np.zeros_like(computational)
        np.add.at(gradient, self.cells[:, 1], derivative[:, :, 0])
        np.add.at(gradient, self.cells[:, 2], derivative[:, :, 1])
        np.add.at(gradient, self.cells[:, 0], -derivative.sum(axis=2))
        return value, gradient


def reconstruct_physical(computational, fixed_physical, cells, reference):
    """Evaluate Phi_h: current computational mesh -> fixed physical mesh."""
    triangulation = Triangulation(computational[:, 0], computational[:, 1], cells)
    coordinates = []
    for axis in range(2):
        values = LinearTriInterpolator(triangulation, fixed_physical[:, axis])(
            reference[:, 0], reference[:, 1])
        if np.any(np.ma.getmaskarray(values)):
            raise ValueError("The reference grid falls outside the computational mesh.")
        coordinates.append(np.asarray(values))
    return np.column_stack(coordinates)


def target_function(points):
    x, y = points[..., 0], points[..., 1]
    return np.tanh(25 * (x - 0.5 - 0.15 * np.sin(2 * np.pi * y)))


def band_metric(centers, strength=30.0):
    """Prescribed smooth normal-direction monitor for the teaching experiment."""
    x, y = centers[:, 0], centers[:, 1]
    distance = x - 0.5 - 0.15 * np.sin(2 * np.pi * y)
    normal = np.column_stack((np.ones_like(x), -0.3 * np.pi * np.cos(2 * np.pi * y)))
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    weight = strength * np.exp(-(distance / 0.09) ** 2)
    return np.eye(2)[None] + weight[:, None, None] * normal[:, :, None] * normal[:, None, :]


def interpolation_rmse(nodes, cells, samples):
    triangulation = Triangulation(nodes[:, 0], nodes[:, 1], cells)
    values = LinearTriInterpolator(triangulation, target_function(nodes))(
        samples[:, 0], samples[:, 1])
    if np.any(np.ma.getmaskarray(values)):
        raise ValueError("Interpolation samples fall outside the physical mesh.")
    return float(np.sqrt(np.mean((np.asarray(values) - target_function(samples)) ** 2)))


def solve(reference, physical, cells, boundary, metric, *, iterations=500,
          step=0.1, frame_stride=5, gradient_tolerance=1e-7):
    """Backtracked explicit gradient descent with signed-area acceptance.

    Both computational and reconstructed physical triangles must remain
    positive. This is a discrete endpoint check, not a nondegeneracy theorem.
    """
    energy = WinslowEnergy(physical, cells, metric)
    q = reference.copy()
    value, gradient = energy.value_gradient(q)
    samples_1d = (np.arange(100) + 0.5) / 100
    xx, yy = np.meshgrid(samples_1d, samples_1d)
    samples = np.column_stack((xx.ravel(), yy.ravel()))
    history, frames = [], []
    moved = physical.copy()
    accepted_step, backtracks = 0.0, 0
    stop_reason = "iteration_limit"

    def record(iteration, save_frame):
        projected = gradient.copy()
        projected[boundary] = 0.0
        row = {"iteration": iteration, "energy": value,
               "computational_min_area": float(signed_triangle_areas(q, cells).min()),
               "physical_min_area": float(signed_triangle_areas(moved, cells).min()),
               "gradient_norm": float(np.linalg.norm(projected)),
               "step": accepted_step, "backtracks": backtracks}
        history.append(row)
        if save_frame:
            frames.append({**row, "computational": q.tolist(), "physical": moved.tolist(),
                           "interpolation_rmse": interpolation_rmse(moved, cells, samples)})

    record(0, True)
    for iteration in range(1, iterations + 1):
        gradient[boundary] = 0.0
        norm_squared = float(np.sum(gradient ** 2))
        if np.sqrt(norm_squared) <= gradient_tolerance:
            stop_reason = "gradient_tolerance"
            break
        accepted_step = step
        for backtracks in range(40):
            candidate = q - accepted_step * gradient
            if np.all(signed_triangle_areas(candidate, cells) > 0):
                candidate_value, candidate_gradient = energy.value_gradient(candidate)
                if candidate_value <= value - 1e-4 * accepted_step * norm_squared:
                    candidate_physical = reconstruct_physical(candidate, physical, cells, reference)
                    if np.all(signed_triangle_areas(candidate_physical, cells) > 0):
                        break
            accepted_step *= 0.5
        else:
            stop_reason = "line_search_stalled"
            break
        q, moved = candidate, candidate_physical
        value, gradient = candidate_value, candidate_gradient
        record(iteration, iteration % frame_stride == 0)
    if frames[-1]["iteration"] != history[-1]["iteration"]:
        frames.append({**history[-1], "computational": q.tolist(), "physical": moved.tolist(),
                       "interpolation_rmse": interpolation_rmse(moved, cells, samples)})
    return {"history": history, "frames": frames, "stop_reason": stop_reason,
            "rmse_sampling": "100x100 fixed midpoint grid; sampled interpolation RMSE, not PDE error"}
