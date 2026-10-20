import numpy as np

from src.geometry import edge_matrix, signed_simplex_measure, signed_triangle_areas


def test_edge_matrix_and_signed_area():
    vertices = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]])
    np.testing.assert_allclose(edge_matrix(vertices), [[2.0, 0.0], [0.0, 1.0]])
    assert signed_simplex_measure(vertices) == 1.0
    assert signed_simplex_measure(vertices[[0, 2, 1]]) == -1.0


def test_degenerate_triangle_is_not_hidden():
    vertices = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    assert signed_simplex_measure(vertices) == 0.0


def test_vectorized_signed_areas_preserve_orientation():
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    cells = np.array([[0, 1, 2], [0, 2, 1]])
    np.testing.assert_allclose(signed_triangle_areas(nodes, cells), [0.5, -0.5])

