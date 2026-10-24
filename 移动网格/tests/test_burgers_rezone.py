import numpy as np
import pytest
from fealpy.mesh import TriangleMesh
from fealpy.backend import backend_manager as bm

from src.burgers_rezone import BurgersData, numerical_metric, run, solve_step, transfer


def test_manufactured_source_matches_independent_differences():
    p = np.random.default_rng(42).uniform(.1, .9, (15, 2))
    data, t, h = BurgersData(), .07, 1e-4
    u, grad, f = data.fields(p, t)
    time_derivative = (data.fields(p, t+h)[0]-data.fields(p, t-h)[0])/(2*h)
    laplacian = np.zeros(len(p))
    for axis in range(2):
        delta = np.eye(2)[axis]*h
        laplacian += (data.fields(p+delta, t)[0]-2*u+data.fields(p-delta, t)[0])/h**2
    np.testing.assert_allclose(time_derivative+u*grad.sum(axis=1)-data.nu*laplacian, f, atol=3e-7)


def test_constant_and_affine_transfer_and_constant_pde_are_preserved():
    bm.set_backend("numpy")
    mesh = TriangleMesh.from_box(nx=4, ny=4)
    x, c = np.array(mesh.node), np.array(mesh.cell)
    y = x.copy()
    b = np.array(mesh.boundary_node_flag())
    y[~b] += np.random.default_rng(42).uniform(-.03, .03, (sum(~b), 2))
    np.testing.assert_allclose(transfer(x, c, np.ones(len(x)), y), 1., atol=1e-14)
    np.testing.assert_allclose(transfer(x, c, x[:, 0]+2*x[:, 1], y), y[:, 0]+2*y[:, 1], atol=1e-14)
    values, pde = solve_step(y, c, np.ones(len(x)), .02, .02, BurgersData(constant=True))
    np.testing.assert_allclose(values, 1., atol=1e-12)
    assert pde["status"] == "converged"
    assert pde["relative_nonlinear_residual"] < 1e-10


def test_fixed_grid_pde_error_decreases_with_refinement():
    coarse = run(n=4, steps=2, moving=False)
    fine = run(n=8, steps=2, moving=False)
    assert fine["states"][-1]["pde"]["l2_error"] < .4*coarse["states"][-1]["pde"]["l2_error"]


def test_time_mesh_and_metric_continue_from_previous_state():
    result = run(n=4, steps=2, dt=.02, inner_cap=.5, max_outer=2, perturbation=.3)
    cells = np.array(result["cells"])
    versions = []
    for previous, current in zip(result["states"], result["states"][1:]):
        assert current["physical_time"]-previous["physical_time"] == pytest.approx(.02)
        np.testing.assert_array_equal(current["outer_iterations"][0]["frozen_physical"], previous["nodes"])
        assert not np.allclose(current["values"], previous["values"])
        assert not np.allclose(current["metric"], previous["metric"])
        np.testing.assert_allclose(current["metric"], numerical_metric(np.array(current["nodes"]), cells, np.array(current["values"])))
        assert current["pde"]["status"] == "converged"
        for k in current["outer_iterations"]:
            versions.append(k["metric_version"])
            assert k["physical_time"] == current["physical_time"]
            assert k["mesh_inner"]["accepted_step_trace"][0]["pseudo_time"] == 0.
            assert min(r["physical_min_area"] for r in k["mesh_inner"]["accepted_step_trace"]) > 0
            expected = numerical_metric(np.array(k["frozen_physical"]), cells, np.array(k["provisional_values"]))
            np.testing.assert_allclose(k["frozen_metric"], expected)
    assert versions == list(range(1, len(versions)+1))


def test_static_input_optimizes_once_and_bad_mesh_is_rejected():
    static = run(n=4, static=True, inner_cap=.5, perturbation=.2)
    assert len(static["states"]) == 2
    final = static["states"][-1]
    assert final["physical_time"] == 0.
    assert final["physical_step_index"] == 0
    assert len(final["outer_iterations"]) == 1
    assert final["pde"]["status"] == "not_advanced_static_input"
    assert final["outer_stop_reason"] == "static_single_optimization"
    with pytest.raises(ValueError, match="invalid_mesh"):
        run(n=4, steps=1, perturbation=2.)
