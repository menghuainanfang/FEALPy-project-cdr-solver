import numpy as np
import pytest

from src.mesh_convergence import AcceptedMeshState, Convergence, MeshConvergenceMonitor
from src.fealpy_mesh_evolution import solve
from tests.test_fealpy_mesh_evolution import mesh_data
from src.variational_mesh import band_metric


def state(s, *, x=0., energy=1., residual=0., epoch="one", valid=True):
    return AcceptedMeshState(s, np.array([[x, 0.]]), energy, residual, valid, "翻转" if not valid else "", epoch)


def monitor(*, cap=1., **config):
    return MeshConvergenceMonitor(np.ones(1), tau=1., pseudo_time_cap=cap,
                                   config=Convergence(window_scaled_time=.1, **config))


def test_cap_reached_without_residual_is_not_converged():
    m = monitor(consecutive_steps=1)
    m.check(state(0., residual=1.), 0)
    row = m.check(state(1., residual=1.), 1)
    assert row["cap_reached"] and not row["criterion_met"]
    assert row["mesh_status"] == "not_converged"
    assert row["solver_status"] == "success"
    assert row["stop_reason"] == "pseudo_time_cap"
    assert any("残量" in reason for reason in row["unmet_reasons"])


def test_criterion_and_cap_at_same_endpoint_preserve_both_facts():
    m = monitor(consecutive_steps=1)
    m.check(state(0.), 0)
    row = m.check(state(1.), 1)
    assert row["criterion_met"] and row["cap_reached"]
    assert row["stop_reason"] == "converged"
    assert row["first_converged_s"] == 1.


def test_zero_initial_motion_still_needs_window_and_consecutive_checks():
    m = monitor()
    assert not m.check(state(0.), 0)["criterion_met"]
    assert not m.check(state(.01), 1)["window_ready"]
    assert not m.check(state(.1), 2)["criterion_met"]
    assert not m.check(state(.11), 3)["criterion_met"]
    assert m.check(state(.12), 4)["criterion_met"]


def test_tiny_step_cannot_hide_motion_over_the_window():
    m = monitor(consecutive_steps=1)
    m.check(state(0.), 0)
    m.check(state(.099999999, x=.01), 1)
    row = m.check(state(.1, x=.01000000001), 2)
    assert row["displacement"] > .009
    assert not row["displacement_met"] and not row["criterion_met"]


def test_geometry_or_solver_failure_cannot_be_converged():
    for valid, solver_success, reason in ((False, True, "invalid_mesh"), (True, False, "solver_failure")):
        m = monitor(consecutive_steps=1)
        m.check(state(0.), 0)
        row = m.check(state(1., valid=valid), 1, solver_success=solver_success)
        assert not row["criterion_met"] and row["stop_reason"] == reason


def test_energy_scaling_and_epoch_reset_do_not_mix_energy_curves():
    rates = []
    for scale in (1., 10.):
        m = monitor(consecutive_steps=1, energy_scale=scale)
        m.check(state(0., energy=2*scale), 0)
        row = m.check(state(.1, energy=(2+1e-8)*scale), 1)
        rates.append(row["energy_change_rate"])
        reset = m.check(state(0., energy=1000., epoch="two"), 0)
        assert reset["epoch_reset"] and reset["energy_change"] is None
        assert not reset["criterion_met"] and reset["consecutive_success"] == 0
    assert rates[0] == pytest.approx(rates[1], rel=1e-6)


@pytest.mark.parametrize("mode,expected", [("fixed", 18), ("tangent", 30)])
def test_free_scalar_dof_definition_matches_boundary_constraints(mode, expected):
    x, c, b = mesh_data()
    result = solve(x, x, c, b, band_metric(x[c].mean(axis=1)), boundary_mode=mode, duration=.01, frames=2)
    assert result["free_scalar_dofs"] == result["norm_definition"]["free_scalar_dofs"] == expected
    assert result["solver_status"] == "success"
    assert result["schema_version"] == "mesh-evolution-v2.0"
