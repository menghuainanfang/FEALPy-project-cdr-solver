"""Orientation-sensitive geometry helpers for simplex meshes."""

from __future__ import annotations

import numpy as np


def edge_matrix(vertices: np.ndarray) -> np.ndarray:
    """Return [v1-v0, ..., vd-v0] with edge vectors stored as columns."""
    vertices = np.asarray(vertices, dtype=np.float64)
    if vertices.ndim != 2 or vertices.shape[0] != vertices.shape[1] + 1:
        raise ValueError("vertices must have shape (d + 1, d)")
    return (vertices[1:] - vertices[0]).T


def signed_simplex_measure(vertices: np.ndarray) -> float:
    """Return signed length/area/volume for a simplex in dimensions 1--3."""
    matrix = edge_matrix(vertices)
    dimension = matrix.shape[0]
    factorial = (1, 1, 2, 6)[dimension]
    return float(np.linalg.det(matrix) / factorial)


def signed_triangle_areas(nodes: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Vectorized signed areas; negative values identify reversed cells."""
    nodes = np.asarray(nodes, dtype=np.float64)
    cells = np.asarray(cells, dtype=np.int64)
    tri = nodes[cells]
    e1 = tri[:, 1] - tri[:, 0]
    e2 = tri[:, 2] - tri[:, 0]
    return 0.5 * (e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])

