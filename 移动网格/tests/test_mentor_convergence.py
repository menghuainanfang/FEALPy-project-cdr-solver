import numpy as np
import pytest

from src.fealpy_mesh_evolution import Convergence, FealpyMeshEvolution, solve, solve_ivp_baseline
from src.variational_mesh import band_metric
from tests.test_fealpy_mesh_evolution import mesh_data


def problem():
    nodes, cells, boundary = mesh_data()
    return nodes, cells, boundary, band_metric(nodes[cells].mean(axis=1))


def test_output_sampling_does_not_change_bdf_steps_or_endpoint():
    x, c, b, m = problem()
    a = solve(x, x, c, b, m, frames=3, convergence=None)
    z = solve(x, x, c, b, m, frames=17, convergence=None)
    assert a["solver_stats"] == z["solver_stats"]
    np.testing.assert_array_equal([r["accepted_dt"] for r in a["accepted_step_trace"]],
                                  [r["accepted_dt"] for r in z["accepted_step_trace"]])
    np.testing.assert_array_equal(a["frames"][-1]["physical"], z["frames"][-1]["physical"])
    baseline = solve_ivp_baseline(x, x, c, b, m, frames=3)
    np.testing.assert_allclose(a["frames"][-1]["computational"], baseline["frames"][-1]["computational"], atol=1e-14)
    assert a["solver_stats"]["nfev"] == baseline["solver_stats"]["nfev"]
    assert "gradient_norm" not in a["history"][0]
    assert "backtracks" not in a["history"][0]


def test_equilibrium_converges_before_cap_and_strict_problem_reaches_cap():
    x, c, b, m = problem()
    identity = np.broadcast_to(np.eye(2), m.shape).copy()
    equilibrium = solve(x, x, c, b, identity, frames=4)
    assert equilibrium["stop_reason"] == "converged"
    assert .05 <= equilibrium["terminal_pseudo_time"] < .5
    assert equilibrium["accepted_step_trace"][-1]["scaled_residual_rms"] <= equilibrium["residual_threshold"]
    strict = solve(x, x, c, b, m, duration=.01, frames=3,
                   convergence=Convergence(residual_atol=1e-14, residual_rtol=0.))
    assert strict["stop_reason"] == "pseudo_time_cap"
    assert strict["terminal_pseudo_time"] == .01


def test_nonpositive_initial_mesh_fails_and_step_limit_is_explicit():
    x, c, b, m = problem()
    inverted = c[:, [0, 2, 1]]
    with pytest.raises(ValueError, match="invalid_mesh"):
        solve(x, x, inverted, b, m)
    damaged = x.copy()
    damaged[6] = [2., 2.]
    with pytest.raises(ValueError, match="invalid_mesh"):
        solve(x, damaged, c, b, m)
    limited = solve(x, x, c, b, m, frames=3, max_internal_steps=1)
    assert limited["stop_reason"] == "max_internal_steps"
    assert limited["solver_stats"]["accepted_bdf_steps"] == 1


def test_scaled_residual_does_not_depend_on_tau():
    x, c, b, m = problem()
    a = solve(x, x, c, b, m, duration=.001, frames=2, tau=.1)
    z = solve(x, x, c, b, m, duration=.001, frames=2, tau=1.)
    np.testing.assert_allclose(a["history"][0]["scaled_residual_rms"], z["history"][0]["scaled_residual_rms"], rtol=1e-12)
    np.testing.assert_allclose(a["history"][0]["projected_gradient_norm"], z["history"][0]["projected_gradient_norm"], rtol=1e-12)
    np.testing.assert_allclose(a["history"][0]["velocity_norm"], 10*z["history"][0]["velocity_norm"], rtol=1e-12)


@pytest.mark.parametrize("mode", ["fixed", "tangent"])
def test_reported_gradient_matches_projected_energy_derivative(mode):
    x, c, b, m = problem()
    model = FealpyMeshEvolution(x, c, m, boundary_mode=mode)
    A, g, trA, edges = model.quantities(x)
    local = model.Idxi_from_Ehat(A, g, trA, edges, model.rho(model.metric), model.metric_theta)
    gradient = np.asarray(model.geo_core.assemble_vector(local))
    on_x = np.isclose(x[:, 0], 0.) | np.isclose(x[:, 0], 1.)
    on_y = np.isclose(x[:, 1], 0.) | np.isclose(x[:, 1], 1.)
    gradient[on_x, 0], gradient[on_y, 1] = 0., 0.
    if mode == "fixed":
        gradient[b] = 0.
    result = solve(x, x, c, b, m, duration=.001, frames=2, boundary_mode=mode)
    np.testing.assert_allclose(result["history"][0]["projected_gradient_norm"], np.linalg.norm(gradient), rtol=1e-12)
