#!/usr/bin/env python3
"""Local, conservative early-stage PR/de-sintering competition closure."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import math
import numpy as np

import agentic_mechanism_model as base
import pore_location_agentic_model as action
import pore_location_topology_model as location
import topology_constrained_sintering as aggregate

R = aggregate.R
EARLY_MEMORY_MODES = ("disabled", "PR_desintering_competition", "PR_plus_connected_fine_attrition")


@dataclass
class PRMemoryParams:
    base: base.DiscoveryParams
    early_memory_mode: str = "disabled"
    k_PR_ref_s: float = 2.0e-5
    Q_PR_J_mol: float = 180e3
    T_PR_ref_C: float = 1050.0
    renewal_gate_mid: float = 0.35
    renewal_gate_width: float = 0.10
    renewal_power: float = 2.0
    fine_radius_mid_ratio: float = 1.0
    fine_radius_width_ratio: float = 0.25
    q_PR: float = 1.0
    topology_power: float = 1.0
    smoothing_share: float = 0.65
    GB_to_TJ_share: float = 0.25
    TJ_to_iso_share: float = 0.10
    attrition_factor: float = 1.8


@dataclass
class PRMemoryState:
    base: base.DiscoveryState
    cumulative_PR_desintering_work: float = 0.0
    cumulative_densifying_work: float = 0.0
    cumulative_non_densifying_work: float = 0.0


def validate(p: PRMemoryParams):
    if p.early_memory_mode not in EARLY_MEMORY_MODES:
        raise ValueError("invalid early_memory_mode")
    if min(p.k_PR_ref_s, p.Q_PR_J_mol, p.renewal_gate_width, p.fine_radius_width_ratio) <= 0:
        raise ValueError("PR rates and widths must be positive")
    shares = p.smoothing_share + p.GB_to_TJ_share + p.TJ_to_iso_share
    if abs(shares - 1.0) > 1e-12 or min(p.smoothing_share, p.GB_to_TJ_share, p.TJ_to_iso_share) < 0:
        raise ValueError("PR flux shares must be nonnegative and sum to one")
    base.validate(p.base)


def initial_state(p: PRMemoryParams):
    validate(p)
    return PRMemoryState(base.initial_state(p.base))


def clone_state(s: PRMemoryState):
    return PRMemoryState(base.clone_state(s.base), s.cumulative_PR_desintering_work,
                         s.cumulative_densifying_work, s.cumulative_non_densifying_work)


def _sig(x):
    return 1.0 / (1.0 + math.exp(float(np.clip(-x, -60, 60))))


def local_competition(s: PRMemoryState, T_C: float, p: PRMemoryParams):
    """Compute local fluxes from instantaneous state only."""
    validate(p)
    d = base.local_mechanism(s.base, T_C, p.base)
    pore = s.base.pore
    radii = pore.pore_radii
    r_ref = p.base.action.location.base.pore_radius0
    fine_weight = 1.0 / (1.0 + np.exp(np.clip((radii / r_ref - p.fine_radius_mid_ratio) /
                                               p.fine_radius_width_ratio, -60, 60)))
    activity = float(np.clip(d["activity"], 0, 1))
    low_activity = _sig((p.renewal_gate_mid - activity) / p.renewal_gate_width) * (1 - activity) ** p.renewal_power
    T = T_C + 273.15
    Tref = p.T_PR_ref_C + 273.15
    thermal = math.exp(float(np.clip(-p.Q_PR_J_mol / R * (1 / T - 1 / Tref), -40, 40)))
    connected = pore.phi_GBseg + pore.phi_TJ
    connected_total = max(float(np.sum(connected)), 1e-300)
    fine_connected = float(np.sum(connected * fine_weight))
    topology_gate = (fine_connected / connected_total) ** p.topology_power
    attrition = p.attrition_factor if p.early_memory_mode == "PR_plus_connected_fine_attrition" else 1.0
    propensity_bins = connected * fine_weight * (r_ref / radii) ** p.q_PR
    J_total = p.k_PR_ref_s * thermal * low_activity * topology_gate * attrition
    source = J_total * propensity_bins

    pr_smooth = np.zeros_like(radii)
    gb_source = source * pore.phi_GBseg / np.maximum(connected, 1e-300)
    move = p.smoothing_share * gb_source[:-1]
    pr_smooth[:-1] -= move
    pr_smooth[1:] += move
    pr_gb_to_tj = p.GB_to_TJ_share * source * pore.phi_GBseg / np.maximum(connected, 1e-300)
    pr_tj_to_iso = p.TJ_to_iso_share * source * pore.phi_TJ / np.maximum(connected, 1e-300)
    pr_flux = float(np.sum(move) + np.sum(pr_gb_to_tj) + np.sum(pr_tj_to_iso))

    if p.early_memory_mode == "disabled":
        pr_smooth[:] = 0
        pr_gb_to_tj[:] = 0
        pr_tj_to_iso[:] = 0
        pr_flux = 0.0

    H_dens = max(float(d["rho_dot"]), 0.0)
    H_PR = pr_flux
    w_dens = H_dens / max(H_dens + H_PR, 1e-300) if p.early_memory_mode != "disabled" else 1.0
    old_gb = d["GBseg_remove"].copy()
    old_tj = d["TJ_remove"].copy()
    d["GBseg_remove"] = old_gb * w_dens
    d["TJ_remove"] = old_tj * w_dens
    d["rho_dot"] = float(np.sum(d["GBseg_remove"] + d["TJ_remove"]))
    d["GB_smooth"] = d["GB_smooth"] + pr_smooth
    d["GB_to_TJ"] = d["GB_to_TJ"] + pr_gb_to_tj
    d["TJ_to_iso"] = d["TJ_to_iso"] + pr_tj_to_iso

    gamma = p.base.action.location.base.gamma_s
    capillary_work = gamma * pr_flux / max(r_ref, 1e-30)
    dens_work = d["sigma_act_total"] * d["rho_dot"]
    return {**d, "H_dens": H_dens, "H_PR": H_PR, "w_dens_competition": w_dens,
            "w_PR_competition": 1 - w_dens, "PR_desintering_flux": pr_flux,
            "PR_smoothing_flux": float(np.sum(move)) if p.early_memory_mode != "disabled" else 0.0,
            "PR_GB_to_TJ_flux": float(np.sum(pr_gb_to_tj)),
            "PR_TJ_to_iso_flux": float(np.sum(pr_tj_to_iso)),
            "connected_fine_pore_fraction": fine_connected / connected_total,
            "pore_mean_radius": float(np.sum(connected * radii) / connected_total),
            "large_pore_fraction": float(np.sum(connected[radii > r_ref]) / connected_total),
            "stress_generation_from_PR": capillary_work,
            "stress_release_by_densification": dens_work,
            "stress_release_by_desintering": capillary_work,
            "PR_work_dot": capillary_work, "densifying_work_dot": dens_work,
            "non_densifying_work_dot": capillary_work}


def run(p: PRMemoryParams, protocol, stop_at_rho: Optional[float] = None,
        initial: Optional[PRMemoryState] = None):
    validate(p)
    if p.early_memory_mode == "disabled":
        return base.run(p.base, protocol, stop_at_rho=stop_at_rho,
                        initial=None if initial is None else initial.base)
    s = initial_state(p) if initial is None else clone_state(initial)
    scalar = "t T_C rho G C_GBseg C_TJ f_clean_GB f_iso activity rho_dot G_dot E_G sigma_base sigma_GBseg_pore sigma_TJ_pore sigma_act_total X_J C_TJ_total C_TJ_pore C_TJ_structural C_TJ_constraint C_TJ_relaxed C_TJ_pinned R_TJ_pore_drag Lambda_TJ Lambda_TJ_structural K_TJ K_TJ_structural Lambda_over_K_TJ P_comp_TJ P_comp_TJ_structural H_dens H_PR w_dens_competition w_PR_competition PR_desintering_flux PR_smoothing_flux PR_GB_to_TJ_flux PR_TJ_to_iso_flux connected_fine_pore_fraction pore_mean_radius large_pore_fraction stress_generation_from_PR stress_release_by_densification stress_release_by_desintering cumulative_PR_desintering_work cumulative_densifying_work cumulative_non_densifying_work".split()
    powers = "P_GBseg_dens P_TJ_dens P_clean_GB P_persistent_junction_drag P_TJ_multihit P_TJ_pore_drag P_TJ_assisted_densification".split()
    h = {k: [] for k in scalar + powers}
    h.update(phi_GBseg=[], phi_TJ=[], phi_iso=[], N_GBseg=[], N_TJ=[], N_iso=[])
    lp = action.effective_location_params(p.base.action)
    while s.base.pore.t < min(protocol.t_end, lp.base.t_max_s) and s.base.pore.rho < lp.base.rho_cap:
        pore = s.base.pore
        T_C = protocol.T(pore.t, pore.rho)
        d = local_competition(s, T_C, p)
        vals = {"t": pore.t, "T_C": T_C, "rho": pore.rho, "G": pore.G,
                "E_G": d["rho_dot"] / (d["G_dot"] / max(pore.G, 1e-30) + 1e-30),
                "cumulative_PR_desintering_work": s.cumulative_PR_desintering_work,
                "cumulative_densifying_work": s.cumulative_densifying_work,
                "cumulative_non_densifying_work": s.cumulative_non_densifying_work, **d}
        for k in scalar + powers:
            h[k].append(vals[k])
        for k in ("phi_GBseg", "phi_TJ", "phi_iso", "N_GBseg", "N_TJ", "N_iso"):
            h[k].append(getattr(pore, k).copy())
        if stop_at_rho is not None and pore.rho >= stop_at_rho:
            break
        dt = min(lp.base.dt_max_s, protocol.t_end - pore.t)
        dT = abs(protocol.T(pore.t + 1, pore.rho) - T_C)
        if dT:
            dt = min(dt, lp.base.dT_max_C / dT)
        if d["rho_dot"] > 0:
            dt = min(dt, lp.base.drho_max / d["rho_dot"])
        if d["G_dot"] > 0:
            dt = min(dt, lp.base.dG_fraction_max * pore.G / d["G_dot"])
        outgb = -d["GBseg_remove"] + d["GB_smooth"] - d["GB_to_TJ"] + d["TJ_to_GBseg_capture"]
        outtj = -d["TJ_remove"] + d["GB_to_TJ"] - d["TJ_to_GBseg_capture"] - d["TJ_to_iso"]
        outiso = d["TJ_to_iso"]
        loss = max(float(np.max(np.maximum(-outgb, 0) / np.maximum(pore.phi_GBseg, 1e-300))),
                   float(np.max(np.maximum(-outtj, 0) / np.maximum(pore.phi_TJ, 1e-300))),
                   max(-d["X_J_dot"] / max(s.base.X_J, 1e-300), 0))
        if loss > 0:
            dt = min(dt, .2 / loss)
        dt = max(lp.base.dt_min_s, dt)
        pore.phi_GBseg = np.maximum(pore.phi_GBseg + outgb * dt, 0)
        pore.phi_TJ = np.maximum(pore.phi_TJ + outtj * dt, 0)
        pore.phi_iso = np.maximum(pore.phi_iso + outiso * dt, 0)
        pore.rho = 1 - float(np.sum(pore.phi_total))
        pore.N_GBseg = location._number(pore.phi_GBseg, pore.pore_radii)
        pore.N_TJ = location._number(pore.phi_TJ, pore.pore_radii)
        pore.N_iso = location._number(pore.phi_iso, pore.pore_radii)
        pore.G = max(pore.G + d["G_dot"] * dt, 1e-9)
        s.base.X_J = float(np.clip(s.base.X_J + d["X_J_dot"] * dt, 0, p.base.XJ_capacity))
        s.cumulative_PR_desintering_work += d["PR_work_dot"] * dt
        s.cumulative_densifying_work += d["densifying_work_dot"] * dt
        s.cumulative_non_densifying_work += d["non_densifying_work_dot"] * dt
        pore.t += dt
    out = {k: np.asarray(v, float) for k, v in h.items()}
    out["pore_radii"] = s.base.pore.pore_radii.copy()
    return out


def final_state(h, p: PRMemoryParams, index=-1):
    b = base.final_state(h, p.base, index)
    def val(name):
        return float(h[name][index]) if name in h else 0.0
    return PRMemoryState(b, val("cumulative_PR_desintering_work"),
                         val("cumulative_densifying_work"), val("cumulative_non_densifying_work"))


LOCAL_FUNCTIONS = (local_competition,)
