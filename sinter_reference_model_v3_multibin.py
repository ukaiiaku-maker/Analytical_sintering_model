#!/usr/bin/env python3
"""
sinter_reference_model_v3_multibin.py

Compact multi-bin pore-population reference model for pressureless sintering.

This is a reduced-order exploratory model for Codex-driven mechanism search. It couples:
- two-barrier renewal densification: nucleation waiting + sink completion;
- surface-diffusion and grain-boundary grain-growth channels;
- pore-size-dependent coverage, removability, drag, and coalescence;
- explicit pore volume bins so small pores close first and large pores remain.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, Optional

import numpy as np

kB = 1.380649e-23
Rgas = 8.31446261815324


def sigmoid(x: float) -> float:
    if x > 60:
        return 1.0
    if x < -60:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def arrhenius(A: float, Q: float, T: float) -> float:
    return A * math.exp(-Q / (Rgas * T))


@dataclass
class Params:
    gamma_s: float = 1.0
    gamma_gb: float = 0.03
    b: float = 3.5e-10
    Omega: float = 4.3174e-29
    rho0: float = 0.74
    G0: float = 200e-9
    Q_nuc: float = 500e3
    Vact: float = 1.4e-28
    D0_gb: float = 7.85113e-4
    Q_gb: float = 446.219e3
    C1: float = 1.702e-6
    sink_radius_factor: float = 0.5
    sigma_factor: float = 2.0
    sigma_cap: float = 3.0e8
    anchor_T_C: float = 1375.0
    anchor_edot: float = 1.5e-3
    phi_max: float = 1.0
    eta_floor: float = 0.05
    Q_surf: float = 539e3
    Q_SD_growth: float = 539e3
    kSD0_m3_s: float = 2.5e-8
    Q_GG: float = 425e3
    kGG0_m3_s: float = 6.1e-7
    T_gbm_mid_C: float = 1180.0
    T_gbm_width_C: float = 55.0
    n_bins: int = 9
    rp0: float = 25e-9
    rmax_factor: float = 20.0
    pore_ln_sigma: float = 0.55
    phi_tail0_frac: float = 0.01
    tail0_center_frac: float = 0.75
    tail0_width_bins: float = 1.2
    coverage_chi: float = 2.2
    zener_k: float = 6.0
    mobile_T_mid_C: float = 1210.0
    mobile_T_width_C: float = 70.0
    highT_damage_rho_mid: float = 0.835
    highT_damage_rho_width: float = 0.020
    loss_SD_coeff: float = 0.08
    loss_drag_coeff: float = 4.0
    loss_clean_coeff: float = 2.0
    loss_GB_pore_coeff: float = 0.05
    removability_q: float = 2.5
    eta_pore_alpha: float = 3.0
    coverage_eta_exp: float = 0.5
    eta_pore_floor: float = 0.06
    PR_coeff: float = 2.0e-3
    PR_activity_exp: float = 1.2
    PR_rho_mid: float = 0.82
    PR_rho_width: float = 0.045
    PR_small_pore_exp: float = 1.2
    pore_drag_coeff: float = 0.6
    drag_damage_coeff: float = 0.12
    clean_damage_coeff: float = 0.05
    coalesce_coeff: float = 0.03
    tail_drag_radius_exp: float = 1.0
    tail_mobile_boost_C: float = 100.0
    drag_transfer_radius_exp: float = 0.8
    rho_cap: float = 0.985
    t_max_s: float = 1.0e6
    dt_min_s: float = 0.5
    dt_max_s: float = 240.0
    dT_max_C: float = 2.0
    drho_max: float = 0.0025
    dG_frac_max: float = 0.015
    dphi_frac_max: float = 0.03


@dataclass
class State:
    t: float
    G: float
    r_bins: np.ndarray
    phi_bins: np.ndarray

    @property
    def rho(self) -> float:
        return 1.0 - float(np.sum(self.phi_bins))


def pore_volume(r: np.ndarray) -> np.ndarray:
    return (4.0 / 3.0) * math.pi * r**3


def initial_bins(p: Params) -> tuple[np.ndarray, np.ndarray]:
    n = int(max(3, p.n_bins))
    r = np.geomspace(p.rp0, p.rp0 * p.rmax_factor, n)
    x = np.log(r / p.rp0)
    core = np.exp(-0.5 * (x / max(p.pore_ln_sigma, 1e-9)) ** 2)
    core /= max(np.sum(core), 1e-300)
    center = p.tail0_center_frac * (n - 1)
    idx = np.arange(n)
    tail = np.exp(-0.5 * ((idx - center) / max(p.tail0_width_bins, 1e-9)) ** 2)
    tail /= max(np.sum(tail), 1e-300)
    fv0 = 1.0 - p.rho0
    phi_tail = min(max(p.phi_tail0_frac, 0.0), 0.8) * fv0
    phi_core = fv0 - phi_tail
    return r, phi_core * core + phi_tail * tail


def initial_state(p: Params) -> State:
    r, phi = initial_bins(p)
    return State(t=0.0, G=p.G0, r_bins=r, phi_bins=phi)


def number_bins(phi: np.ndarray, r: np.ndarray) -> np.ndarray:
    return np.maximum(phi, 0.0) / np.maximum(pore_volume(r), 1e-300)


def active_radius(p: Params, phi: np.ndarray, r: np.ndarray) -> float:
    w = np.maximum(phi, 0.0) * (p.rp0 / np.maximum(r, 1e-30)) ** max(p.removability_q, 0.0)
    if np.sum(w) <= 0:
        return float(p.rp0)
    return float(np.sum(w * r) / np.sum(w))


def volume_radius(phi: np.ndarray, r: np.ndarray) -> float:
    if np.sum(phi) <= 0:
        return float(r[0])
    return float(np.sum(phi * r) / np.sum(phi))


def large_fraction(p: Params, phi: np.ndarray, r: np.ndarray) -> float:
    if np.sum(phi) <= 0:
        return 0.0
    return float(np.sum(phi[r >= 4.0 * p.rp0]) / np.sum(phi))


def coverage_fraction(p: Params, s: State) -> float:
    N = number_bins(s.phi_bins, s.r_bins)
    A_pore_V = float(np.sum(N * math.pi * s.r_bins**2))
    A_gb_V = max(2.0 / max(s.G, 1e-30), 1e-30)
    return 1.0 - math.exp(-max(p.coverage_chi * A_pore_V / A_gb_V, 0.0))


def calibrate_nu(p: Params) -> float:
    s = initial_state(p)
    ract = active_radius(p, s.phi_bins, s.r_bins)
    T = p.anchor_T_C + 273.15
    sigma = min(p.sigma_factor * 2.0 * p.gamma_s / max(ract, 1e-30), p.sigma_cap)
    Dgb = arrhenius(p.D0_gb, p.Q_gb, T)
    Rs = p.sink_radius_factor * s.G
    tau = p.C1 * (kB * T / max(sigma * p.Omega, 1e-300)) * (Rs**2 / max(Dgb, 1e-300))
    ceiling = p.phi_max * (p.b / s.G) / max(tau, 1e-300)
    target = min(p.anchor_edot, 0.75 * ceiling)
    denom = max(p.phi_max * (p.b / s.G) - target * tau, 1e-300)
    r_req = target / denom
    exponent = -p.Q_nuc / (Rgas * T) + p.Vact * sigma / (kB * T)
    return max(r_req / max(math.exp(np.clip(exponent, -700, 700)), 1e-300), 1e-300)


def compute_rates(p: Params, s: State, T_C: float, nu_eff: float) -> Dict[str, float]:
    T = T_C + 273.15
    r = s.r_bins
    phi = np.maximum(s.phi_bins, 0.0)
    fv = max(np.sum(phi), 1e-12)
    rho = s.rho
    ract = active_radius(p, phi, r)
    rvol = volume_radius(phi, r)
    sigma = min(p.sigma_factor * 2.0 * p.gamma_s / max(ract, 1e-30), p.sigma_cap)
    Dgb = arrhenius(p.D0_gb, p.Q_gb, T)
    Rs = p.sink_radius_factor * s.G
    tau = p.C1 * (kB * T / max(sigma * p.Omega, 1e-300)) * (Rs**2 / max(Dgb, 1e-300))
    rnuc = nu_eff * math.exp(np.clip(-p.Q_nuc / (Rgas * T) + p.Vact * sigma / (kB * T), -700, 700))
    Lam = rnuc * tau
    activity = Lam / (1.0 + Lam)
    edot0 = p.phi_max * (p.b / max(s.G, 1e-30)) * activity / max(tau, 1e-300)
    P_dens0 = fv * sigma * edot0
    k_sd = arrhenius(p.kSD0_m3_s, p.Q_SD_growth, T)
    k_gg = arrhenius(p.kGG0_m3_s, p.Q_GG, T)
    dG_sd = k_sd / max(3.0 * s.G**2, 1e-300)
    dG_gb_clean = k_gg / max(3.0 * s.G**2, 1e-300)
    w_gbm = sigmoid((T_C - p.T_gbm_mid_C) / max(p.T_gbm_width_C, 1e-9))
    f_pore = coverage_fraction(p, s)
    f_clean = max(0.0, 1.0 - f_pore)
    Rz = p.zener_k * ract / max(fv, 1e-12)
    S_zener = max(0.0, 1.0 - s.G / max(2.0 * Rz, 1e-300))
    g_mobile = sigmoid((T_C - p.mobile_T_mid_C) / max(p.mobile_T_width_C, 1e-9))
    gamma_pore = max(0.0, min(1.0, g_mobile * (0.25 + 0.75 * S_zener)))
    dG_gb_pore = dG_gb_clean * gamma_pore
    dGdt = (1.0 - w_gbm) * dG_sd + w_gbm * (f_pore * dG_gb_pore + f_clean * dG_gb_clean)
    pref = 4.0 / max(s.G**2, 1e-300)
    P_SD = p.gamma_s * pref * max((1.0 - w_gbm) * dG_sd, 0.0)
    P_GB_pore = p.gamma_gb * pref * max(w_gbm * f_pore * dG_gb_pore, 0.0)
    P_GB_clean = p.gamma_gb * pref * max(w_gbm * f_clean * dG_gb_clean, 0.0)
    rel = np.maximum(r / max(p.rp0, 1e-30), 1.0)
    g_mobile_i = 1.0 / (1.0 + np.exp(-np.clip((T_C - (p.mobile_T_mid_C - p.tail_mobile_boost_C * np.log2(rel))) / max(p.mobile_T_width_C, 1e-9), -60, 60)))
    drag_weight = rel ** p.tail_drag_radius_exp
    P_drag_bins = p.pore_drag_coeff * p.gamma_gb * phi * max(w_gbm * dG_gb_clean, 0.0) / np.maximum(r**2, 1e-300) * drag_weight * g_mobile_i
    P_drag = float(np.sum(P_drag_bins))
    P_loss = (p.loss_SD_coeff * (1.0 - activity) * P_SD + p.loss_drag_coeff * P_drag + p.loss_clean_coeff * P_GB_clean + p.loss_GB_pore_coeff * P_GB_pore)
    eta_Ons = max(p.eta_floor, P_dens0 / max(P_dens0 + P_loss, 1e-300))
    remov = np.sum(phi * (p.rp0 / np.maximum(r, 1e-30)) ** p.removability_q) / max(np.sum(phi), 1e-300)
    eta_pore = max(p.eta_pore_floor, min(1.0, math.exp(-p.eta_pore_alpha * (1.0 - remov)) * max(f_pore, 1e-12) ** p.coverage_eta_exp))
    rho_dot = fv * edot0 * eta_Ons * eta_pore
    f_1gb = 1.0 / (1.0 + math.exp(np.clip((rho - p.PR_rho_mid) / max(p.PR_rho_width, 1e-9), -60, 60)))
    pr_rate = p.PR_coeff * P_SD * (1.0 - activity) ** p.PR_activity_exp * f_1gb
    highT_stage = sigmoid((rho - p.highT_damage_rho_mid) / max(p.highT_damage_rho_width, 1e-9))
    drag_transfer_rate = highT_stage * (p.drag_damage_coeff * P_drag + p.clean_damage_coeff * P_GB_clean)
    coalesce_rate = highT_stage * p.coalesce_coeff * P_drag
    return {"T_C": T_C, "T_K": T, "rho": rho, "rr": ract, "r_active": ract, "r_vol": rvol, "sigma": sigma, "Dgb": Dgb, "tau_sink": tau, "rnuc": rnuc, "Lambda": Lam, "activity": activity, "edot0": edot0, "P_dens0": P_dens0, "P_SD": P_SD, "P_GB_pore": P_GB_pore, "P_GB_clean": P_GB_clean, "P_drag": P_drag, "P_loss": P_loss, "eta_Ons": eta_Ons, "eta_pore": eta_pore, "rho_dot": rho_dot, "dGdt": dGdt, "dG_sd": dG_sd, "dG_gb_clean": dG_gb_clean, "dG_gb_pore": dG_gb_pore, "w_gbm": w_gbm, "f_pore": f_pore, "f_clean": f_clean, "S_zener": S_zener, "g_mobile": g_mobile, "f_large": large_fraction(p, phi, r), "pr_rate": pr_rate, "drag_transfer_rate": drag_transfer_rate, "coalesce_rate": coalesce_rate}


class RampHoldCool:
    def __init__(self, heating_rate_C_min: float, soak_T_C: float = 1500.0, soak_time_s: float = 3600.0, T_start_C: float = 25.0, T_end_C: float = 25.0, cooling_rate_C_min: Optional[float] = None):
        self.hr = heating_rate_C_min / 60.0
        self.cr = (cooling_rate_C_min if cooling_rate_C_min is not None else heating_rate_C_min) / 60.0
        self.soak_T_C = soak_T_C
        self.soak_time_s = soak_time_s
        self.T_start_C = T_start_C
        self.T_end_C = T_end_C
        self.t_ramp = max(soak_T_C - T_start_C, 0.0) / max(self.hr, 1e-300)
        self.t_cool = max(soak_T_C - T_end_C, 0.0) / max(self.cr, 1e-300)
        self.t_end = self.t_ramp + soak_time_s + self.t_cool

    def T(self, t: float, rho: float) -> float:
        if t <= self.t_ramp:
            return self.T_start_C + self.hr * t
        if t <= self.t_ramp + self.soak_time_s:
            return self.soak_T_C
        if t <= self.t_end:
            return self.soak_T_C - self.cr * (t - self.t_ramp - self.soak_time_s)
        return self.T_end_C


class Iso:
    def __init__(self, T_C: float, t_max_s: float):
        self.T_C = T_C
        self.t_end = t_max_s

    def T(self, t: float, rho: float) -> float:
        return self.T_C


class TwoStep:
    def __init__(self, T1_C: float = 1350.0, T2_C: float = 1300.0, rho_switch: float = 0.83, t_max_s: float = 96*3600):
        self.T1_C = T1_C
        self.T2_C = T2_C
        self.rho_switch = rho_switch
        self.t_end = t_max_s

    def T(self, t: float, rho: float) -> float:
        return self.T2_C if rho >= self.rho_switch else self.T1_C


def choose_dt(p: Params, proto, s: State, rates: Dict[str, float]) -> float:
    dt = min(p.dt_max_s, max(getattr(proto, "t_end", p.t_max_s) - s.t, 0.0))
    dTdt = abs(proto.T(s.t + 1.0, s.rho) - proto.T(s.t, s.rho))
    if dTdt > 1e-12:
        dt = min(dt, p.dT_max_C / dTdt)
    if rates["rho_dot"] > 0:
        dt = min(dt, p.drho_max / max(rates["rho_dot"], 1e-300))
    if rates["dGdt"] > 0:
        dt = min(dt, p.dG_frac_max * s.G / max(rates["dGdt"], 1e-300))
    transfer = rates["pr_rate"] + rates["drag_transfer_rate"] + rates["rho_dot"]
    phi_tot = max(float(np.sum(s.phi_bins)), 1e-12)
    if transfer > 0:
        dt = min(dt, p.dphi_frac_max * phi_tot / transfer)
    return max(p.dt_min_s, dt)


def move_phi_up(phi: np.ndarray, amount: float, weights: np.ndarray) -> np.ndarray:
    phi = phi.copy()
    if amount <= 0 or np.sum(phi) <= 0:
        return phi
    weights = np.maximum(weights, 0.0) * (phi > 0)
    if np.sum(weights) <= 0:
        weights = phi.copy()
    weights = weights / max(np.sum(weights), 1e-300)
    moved = np.minimum(phi, amount * weights)
    rem = amount - float(np.sum(moved))
    if rem > 1e-300:
        w2 = np.maximum(phi - moved, 0.0)
        if np.sum(w2) > 0:
            moved += np.minimum(phi - moved, rem * w2 / np.sum(w2))
    phi -= moved
    up = np.zeros_like(phi)
    up[1:] = moved[:-1]
    up[-1] += moved[-1]
    return np.maximum(phi + up, 0.0)


def remove_for_densification(p: Params, phi: np.ndarray, r: np.ndarray, amount: float) -> np.ndarray:
    phi = phi.copy()
    if amount <= 0 or np.sum(phi) <= 0:
        return phi
    w = phi * (p.rp0 / np.maximum(r, 1e-30)) ** p.removability_q
    if np.sum(w) <= 0:
        w = phi.copy()
    rem = np.minimum(phi, amount * w / max(np.sum(w), 1e-300))
    leftover = amount - float(np.sum(rem))
    if leftover > 1e-300 and np.sum(phi - rem) > 0:
        rem += np.minimum(phi - rem, leftover * (phi - rem) / np.sum(phi - rem))
    return np.maximum(phi - rem, 0.0)


def run(p: Params, protocol, stop_at_rho: Optional[float] = None) -> Dict[str, np.ndarray]:
    nu_eff = calibrate_nu(p)
    s = initial_state(p)
    keys = ["t", "T_C", "rho", "G", "sigma", "Lambda", "activity", "edot0", "rho_dot", "P_dens0", "P_SD", "P_GB_pore", "P_GB_clean", "P_drag", "P_loss", "eta_Ons", "eta_pore", "f_pore", "f_clean", "S_zener", "g_mobile", "w_gbm", "f_large", "r_active", "r_vol", "dGdt", "pr_rate", "drag_transfer_rate", "coalesce_rate"]
    hist = {k: [] for k in keys}
    t_end = min(getattr(protocol, "t_end", p.t_max_s), p.t_max_s)
    while s.t < t_end - 1e-12 and s.rho < p.rho_cap:
        T_C = protocol.T(s.t, s.rho)
        rates = compute_rates(p, s, T_C, nu_eff)
        for k in keys:
            if k == "t": hist[k].append(s.t)
            elif k == "G": hist[k].append(s.G)
            elif k == "rho": hist[k].append(s.rho)
            else: hist[k].append(rates.get(k, np.nan))
        if stop_at_rho is not None and s.rho >= stop_at_rho:
            break
        dt = choose_dt(p, protocol, s, rates)
        if dt <= 0:
            break
        dens_amt = min(float(np.sum(s.phi_bins)), rates["rho_dot"] * dt)
        phi_new = remove_for_densification(p, s.phi_bins, s.r_bins, dens_amt)
        w_pr = phi_new * (p.rp0 / np.maximum(s.r_bins, 1e-30)) ** p.PR_small_pore_exp
        phi_new = move_phi_up(phi_new, rates["pr_rate"] * dt, w_pr)
        rel = np.maximum(s.r_bins / max(p.rp0, 1e-30), 1.0)
        w_drag = phi_new * rel ** p.drag_transfer_radius_exp
        phi_new = move_phi_up(phi_new, rates["drag_transfer_rate"] * dt, w_drag)
        phi_new = move_phi_up(phi_new, rates["coalesce_rate"] * dt, w_drag)
        s.phi_bins = phi_new
        s.G = max(1e-9, s.G + rates["dGdt"] * dt)
        s.t += dt
    T_C = protocol.T(s.t, s.rho)
    rates = compute_rates(p, s, T_C, nu_eff)
    for k in keys:
        if k == "t": hist[k].append(s.t)
        elif k == "G": hist[k].append(s.G)
        elif k == "rho": hist[k].append(s.rho)
        else: hist[k].append(rates.get(k, np.nan))
    out = {k: np.asarray(v, dtype=float) for k, v in hist.items()}
    out["nu_eff"] = np.asarray([nu_eff])
    return out


if __name__ == "__main__":
    for name, proto in {"slow": RampHoldCool(0.2, 1500, 3600), "fast": RampHoldCool(20.0, 1500, 3600), "high": Iso(1350, 96*3600), "two": TwoStep(1350, 1300, 0.83, 96*3600)}.items():
        pp = Params(G0=200e-9 if name in ["slow", "fast"] else 50e-9)
        r = run(pp, proto, stop_at_rho=0.92)
        print(name, "rho", r["rho"][-1], "G_nm", r["G"][-1]*1e9, "t_h", r["t"][-1]/3600)
