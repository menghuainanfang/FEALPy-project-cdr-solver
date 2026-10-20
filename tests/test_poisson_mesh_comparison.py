import numpy as np
import pytest
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh

from src.poisson_mesh_comparison import PoissonCase, solve_poisson


def test_peak_metric_is_independent_of_sampling_batch():
    case = PoissonCase("peak")
    points = np.array([[.43,.5],[.5,.6],[.7,.7]])
    with_center = np.vstack((points, [[.5,.5]]))
    np.testing.assert_allclose(case.metric(points), case.metric(with_center)[:len(points)], atol=1e-14)


@pytest.mark.parametrize("name", ["smooth", "band", "peak"])
def test_manufactured_source_matches_independent_gradient_divergence(name):
    case = PoissonCase(name)
    points = np.array([[.51,.5],[.54,.46],[.2,.3],[.7,.6]])
    h = 1e-6
    divergence = np.zeros(len(points))
    for axis in range(2):
        direction = np.zeros(2)
        direction[axis] = h
        divergence += (case.fields(points + direction)[1][:,axis]
                       - case.fields(points - direction)[1][:,axis]) / (2*h)
    np.testing.assert_allclose(case.fields(points)[2], -divergence, rtol=2e-7, atol=2e-5)


def test_p1_reproduces_affine_solution_on_a_perturbed_mesh():
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box([0.,1.,0.,1.],nx=4,ny=4)
    nodes = np.asarray(mesh.entity("node")).copy()
    boundary = np.asarray(mesh.boundary_node_flag())
    nodes[~boundary] += np.random.default_rng(2).uniform(-.03,.03,(np.count_nonzero(~boundary),2))
    result = solve_poisson(nodes,np.asarray(mesh.entity("cell")),PoissonCase("affine"))
    assert result["metrics"]["l2_error"] < 1e-12
    assert result["metrics"]["h1_seminorm_error"] < 1e-11
    assert result["metrics"]["boundary_error"] < 1e-13
    assert result["metrics"]["relative_residual"] < 1e-12
    assert len(result["cell_l2_error"]) == mesh.number_of_cells()
    np.testing.assert_allclose(np.linalg.norm(result["cell_l2_error"]),
                               result["metrics"]["l2_error"], atol=1e-14)
