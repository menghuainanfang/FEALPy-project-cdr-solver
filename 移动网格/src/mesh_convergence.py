"""Accepted-state mesh convergence diagnostics for one frozen metric epoch."""

from collections import deque
from dataclasses import asdict, dataclass

import numpy as np

CRITERION_VERSION = "projected-balanced-residual-window-v2"


@dataclass(frozen=True)
class Convergence:
    residual_atol: float = 1e-6
    residual_rtol: float = 1e-4
    displacement_atol: float = 1e-5
    energy_rate_atol: float = 1e-7
    min_pseudo_time: float = .05
    consecutive_steps: int = 3
    window_scaled_time: float = 1.
    energy_scale: float = 1.


@dataclass(frozen=True)
class AcceptedMeshState:
    pseudo_time: float
    coordinates: np.ndarray
    energy: float
    residual_rms: float
    geometry_valid: bool = True
    geometry_message: str = ""
    epoch_id: str = "frozen_mesh_0"


class MeshConvergenceMonitor:
    """Monitor R=rms(tau*v/h) and changes over fixed scaled time s/tau.

    The window anchor is linearly interpolated between accepted states at
    s-window_scaled_time*tau. Dense replay samples never enter this monitor.
    A new epoch resets energy, residual reference and consecutive successes.
    """

    def __init__(self, local_h, *, tau, pseudo_time_cap, config=Convergence(), max_accepted_steps=2000):
        self.h = np.asarray(local_h)
        self.tau = tau
        self.cap = pseudo_time_cap
        self.config = config
        self.max_steps = max_accepted_steps
        if config is not None and (config.window_scaled_time <= 0 or config.energy_scale <= 0 or config.consecutive_steps < 1):
            raise ValueError("Positive diagnostic window/energy scale and at least one consecutive success are required")
        self.epoch = None
        self.history = deque()
        self.consecutive = 0
        self.first_residual = self.first_conditions = self.first_converged = None
        self.last = None

    def check(self, state, accepted_step_index, *, solver_success=True):
        cfg = self.config
        reset = state.epoch_id != self.epoch
        if reset:
            self.epoch = state.epoch_id
            self.history.clear()
            self.initial_s, self.r_ref = state.pseudo_time, state.residual_rms
            self.consecutive = 0
            self.first_residual = self.first_conditions = self.first_converged = None
        finite = bool(np.isfinite(state.energy) and np.isfinite(state.residual_rms) and
                      np.isfinite(state.coordinates).all())
        geometry_valid = bool(state.geometry_valid and finite)
        residual_limit = None if cfg is None else cfg.residual_atol+cfg.residual_rtol*self.r_ref
        window_s = None if cfg is None else cfg.window_scaled_time*self.tau
        ready = bool(cfg is not None and self.history and
                     state.pseudo_time-self.initial_s >= window_s)
        displacement = energy_change = energy_rate = energy_reference = anchor_s = None
        if ready:
            anchor_s = state.pseudo_time-window_s
            while len(self.history) >= 2 and self.history[1][0] <= anchor_s:
                self.history.popleft()
            left = self.history[0]
            right = self.history[1] if len(self.history) >= 2 else (state.pseudo_time, state.coordinates, state.energy)
            alpha = (anchor_s-left[0])/(right[0]-left[0])
            anchor_q = left[1]+alpha*(right[1]-left[1])
            anchor_energy = left[2]+alpha*(right[2]-left[2])
            displacement = float(np.max(np.linalg.norm((state.coordinates-anchor_q)/self.h[:, None], axis=1)))
            energy_reference = max(abs(anchor_energy), cfg.energy_scale)
            energy_change = abs(state.energy-anchor_energy)/energy_reference
            energy_rate = energy_change/cfg.window_scaled_time
        residual_met = bool(cfg is not None and state.residual_rms <= residual_limit)
        displacement_met = bool(ready and displacement <= cfg.displacement_atol)
        energy_met = bool(ready and energy_rate <= cfg.energy_rate_atol)
        time_met = bool(cfg is not None and state.pseudo_time-self.initial_s >= cfg.min_pseudo_time)
        conditions_met = bool(residual_met and displacement_met and energy_met and time_met and
                              geometry_valid and solver_success and accepted_step_index > 0 and not reset)
        self.consecutive = self.consecutive+1 if conditions_met else 0
        criterion_met = bool(cfg is not None and self.consecutive >= cfg.consecutive_steps)
        if residual_met and self.first_residual is None:
            self.first_residual = state.pseudo_time
        if conditions_met and self.first_conditions is None:
            self.first_conditions = state.pseudo_time
        if criterion_met and self.first_converged is None:
            self.first_converged = state.pseudo_time
        cap_reached = bool(state.pseudo_time >= self.cap)
        step_limit_reached = bool(accepted_step_index >= self.max_steps)
        reason = ("invalid_mesh" if not geometry_valid else "solver_failure" if not solver_success else
                  "converged" if criterion_met else "pseudo_time_cap" if cap_reached else
                  "max_accepted_steps" if step_limit_reached else None)
        unmet = []
        if not geometry_valid:
            unmet.append(state.geometry_message or "几何无效或状态含非有限数")
        if not solver_success:
            unmet.append("数值积分器未成功推进")
        if cfg is None:
            unmet.append("未启用网格收敛停止，仅记录轨迹")
        else:
            if not residual_met:
                unmet.append(f"残量 {state.residual_rms:.6g} > 阈值 {residual_limit:.6g}")
            if not ready:
                unmet.append(f"诊断窗口不足 {window_s:.6g} 网格伪时间")
            elif not displacement_met:
                unmet.append(f"窗口相对位移 {displacement:.6g} > 阈值 {cfg.displacement_atol:.6g}")
            if ready and not energy_met:
                unmet.append(f"窗口能量变化率 {energy_rate:.6g} > 阈值 {cfg.energy_rate_atol:.6g}")
            if not time_met:
                unmet.append(f"未达到最小伪时间 {cfg.min_pseudo_time:.6g}")
            if not criterion_met:
                unmet.append(f"连续达标 {self.consecutive}/{cfg.consecutive_steps}")
        row = {"criterion_version": CRITERION_VERSION, "epoch_id": self.epoch,
               "mesh_pseudo_time": float(state.pseudo_time), "accepted_step_index": accepted_step_index,
               "mesh_energy": float(state.energy), "mesh_residual_rms": float(state.residual_rms),
               "residual_reference": float(self.r_ref), "residual_threshold": residual_limit,
               "window_pseudo_time": window_s, "window_scaled_time": None if cfg is None else cfg.window_scaled_time,
               "window_anchor_s": anchor_s, "window_ready": ready,
               "displacement": displacement, "displacement_threshold": None if cfg is None else cfg.displacement_atol,
               "energy_change": energy_change, "energy_reference": energy_reference,
               "energy_change_rate": energy_rate, "energy_rate_threshold": None if cfg is None else cfg.energy_rate_atol,
               "residual_met": residual_met, "displacement_met": displacement_met, "energy_met": energy_met,
               "minimum_time_met": time_met, "geometry_valid": geometry_valid,
               "conditions_met": conditions_met, "consecutive_success": self.consecutive,
               "required_successes": None if cfg is None else cfg.consecutive_steps,
               "criterion_met": criterion_met, "cap_reached": cap_reached,
               "step_limit_reached": step_limit_reached, "stop_reason": reason,
               "solver_status": "failed" if not solver_success else "success" if accepted_step_index > 0 else "not_started",
               "mesh_status": "converged" if criterion_met else "invalid_mesh" if not geometry_valid else "not_converged",
               "first_residual_met_s": self.first_residual, "first_conditions_met_s": self.first_conditions,
               "first_converged_s": self.first_converged, "unmet_reasons": unmet,
               "epoch_reset": reset}
        if geometry_valid and solver_success:
            self.history.append((state.pseudo_time, state.coordinates.copy(), state.energy))
        self.last = row
        return row

    def parameters(self):
        return {"criterion_version": CRITERION_VERSION, "config": None if self.config is None else asdict(self.config),
                "pseudo_time_cap": self.cap, "tau": self.tau, "max_accepted_steps": self.max_steps,
                "precedence": ["invalid_mesh", "solver_failure", "converged", "pseudo_time_cap", "max_accepted_steps"],
                "cap_policy": "criterion_met wins at the same accepted endpoint; cap_reached remains true",
                "residual_formula": "sqrt(sum_free((tau*v_i/h_i)^2)/N_free_scalar)",
                "residual_kind": "dimensionless projected balanced velocity residual; not raw gradient",
                "displacement_formula": "max_i norm((xi(s)-xi(s-window))/h_i)",
                "energy_formula": "abs(I(s)-I(s-window))/max(abs(I(s-window)),I_scale)/(window/tau)",
                "window_anchor": "linear interpolation between accepted states; no output frames used"}
