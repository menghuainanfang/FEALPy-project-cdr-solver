import numpy as np
import pytest
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh

from src.fealpy_mesh_evolution import FealpyMeshEvolution, solve, quality_metrics
from src.variational_mesh import band_metric


def mesh_data():
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box([0., 1., 0., 1.], nx=4, ny=4)
    return np.asarray(mesh.node), np.asarray(mesh.cell), np.asarray(mesh.boundary_node_flag())


@pytest.mark.parametrize("anisotropic", [False, True])
def test_fealpy_coordinate_gradient_matches_energy_difference(anisotropic):
    nodes, cells, _ = mesh_data()
    metric = band_metric(nodes[cells].mean(axis=1)) if anisotropic else np.broadcast_to(np.eye(2), (len(cells),2,2))
    model = FealpyMeshEvolution(nodes, cells, metric)
    direction = np.random.default_rng(42).normal(size=nodes.shape)
    A, g, trA, E_hat = model.quantities(nodes)
    local = model.Idxi_from_Ehat(A, g, trA, E_hat, model.rho(model.metric), model.metric_theta)
    gradient = np.asarray(model.geo_core.assemble_vector(local))
    h = 1e-7
    difference = (model.energy_velocity(nodes+h*direction)[0]
                  - model.energy_velocity(nodes-h*direction)[0]) / (2*h)
    np.testing.assert_allclose(difference, np.sum(gradient*direction), rtol=2e-6, atol=1e-7)


def test_fealpy_tangent_velocity_has_no_normal_component_and_fixed_corners():
    nodes, cells, boundary = mesh_data()
    model = FealpyMeshEvolution(nodes, cells, band_metric(nodes[cells].mean(axis=1)))
    velocity = model.energy_velocity(nodes)[1]
    np.testing.assert_allclose(np.sum(velocity[model.Bdinnernode_idx]*model.Bi_Lnode_normal, axis=1), 0., atol=1e-13)
    np.testing.assert_allclose(velocity[model.Vertices_idx], 0., atol=1e-13)
    assert np.linalg.norm(velocity[boundary]) > .01


@pytest.mark.parametrize("mode", ["fixed", "tangent"])
def test_fealpy_replay_preserves_geometry_and_boundary_constraints(mode):
    nodes, cells, boundary = mesh_data()
    result = solve(nodes, nodes, cells, boundary, band_metric(nodes[cells].mean(axis=1)),
                   duration=.1, frames=5, boundary_mode=mode, metric_function=band_metric)
    rows = result["history"]
    assert np.max(np.diff([r["energy"] for r in rows])) < 1e-7
    assert min(r["physical_min_area"] for r in rows) > 0
    assert min(r["computational_min_area"] for r in rows) > 0
    assert max(r["boundary_normal_displacement"] for r in rows) < 1e-12
    assert max(r["corner_displacement"] for r in rows) < 1e-12
    tangent = rows[-1]["boundary_tangent_displacement"]
    assert tangent < 1e-12 if mode == "fixed" else tangent > .01


def test_reference_grid_has_unit_metric_quality():
    nodes, cells, _ = mesh_data()
    result = quality_metrics(nodes, nodes, cells, np.broadcast_to(np.eye(2), (len(cells),2,2)))
    np.testing.assert_allclose(list(result.values()), [1.,1.], atol=1e-12)
