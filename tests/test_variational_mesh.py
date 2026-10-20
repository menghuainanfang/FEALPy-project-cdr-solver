import numpy as np
import pytest
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh

from src.geometry import signed_triangle_areas
from src.variational_mesh import WinslowEnergy, band_metric, reconstruct_physical, solve


def mesh_data(n=4):
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box([0., 1., 0., 1.], nx=n, ny=n)
    return (np.asarray(mesh.entity("node")).copy(), np.asarray(mesh.entity("cell")),
            np.asarray(mesh.boundary_node_flag()))


def test_single_triangle_matches_hand_energy_and_coordinate_derivative():
    physical = np.array([[0., 0.], [2., 0.], [0., 1.]])
    cells = np.array([[0, 1, 2]])
    energy = WinslowEnergy(physical, cells, np.eye(2)[None])
    for a in (0.7, 1., 1.3):
        computational = np.array([[0., 0.], [a, 0.], [0., 1.]])
        value, gradient = energy.value_gradient(computational)
        assert value == pytest.approx(1 + a * a / 4)
        assert gradient[1, 0] == pytest.approx(a / 2)
        np.testing.assert_allclose(gradient.sum(axis=0), 0, atol=1e-14)


@pytest.mark.parametrize("adaptive", [False, True])
def test_assembled_gradient_matches_independent_directional_difference(adaptive):
    physical, cells, boundary = mesh_data()
    metric = band_metric(physical[cells].mean(axis=1)) if adaptive else np.broadcast_to(np.eye(2), (len(cells), 2, 2))
    energy = WinslowEnergy(physical, cells, metric)
    rng = np.random.default_rng(3)
    computational = physical + rng.normal(0, 0.01, physical.shape)
    direction = rng.normal(size=physical.shape)
    direction[boundary] = 0
    _, gradient = energy.value_gradient(computational)
    h = 1e-5
    difference = (energy.value_gradient(computational + h * direction)[0]
                  - energy.value_gradient(computational - h * direction)[0]) / (2 * h)
    assert difference == pytest.approx(np.sum(gradient * direction), rel=1e-8, abs=1e-8)


def test_reconstruction_is_identity_for_unchanged_correspondence():
    reference, cells, _ = mesh_data()
    np.testing.assert_allclose(reconstruct_physical(reference, reference, cells, reference), reference, atol=1e-14)


def test_backtracking_preserves_energy_orientation_and_fixed_boundaries():
    reference, cells, boundary = mesh_data()
    physical = reference.copy()
    physical[~boundary] += np.random.default_rng(4).uniform(-0.025, 0.025, (np.count_nonzero(~boundary), 2))
    metric = np.broadcast_to(np.eye(2), (len(cells), 2, 2))
    result = solve(reference, physical, cells, boundary, metric, iterations=50, step=10.)
    history = result["history"]
    assert history[-1]["energy"] < history[0]["energy"]
    assert np.all(np.diff([h["energy"] for h in history]) <= 1e-12)
    assert any(h["backtracks"] > 0 for h in history)
    for frame in result["frames"]:
        for key in ("physical", "computational"):
            nodes = np.asarray(frame[key])
            assert np.min(signed_triangle_areas(nodes, cells)) > 0
            assert signed_triangle_areas(nodes, cells).sum() == pytest.approx(1.)
            np.testing.assert_allclose(nodes[boundary], reference[boundary], atol=1e-14)
