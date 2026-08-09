#!/usr/bin/env python3
"""
sinter_reference_model_v6_grainstress_multibin.py

Wrapper around sinter_reference_model_v3_multibin.py adding two mechanism-specific controls:

1. Lambda_ref: decouples renewal activity scale from Q_nuc.
2. grain-size baseline stress: sigma_base = baseline_stress_scale * 2*gamma_s/(grain_radius_factor*G), so nanograins have high capillary stress and coarse grains can coarsen out of the high-driving-force activity window.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

import sinter_reference_model_v3_multibin as base

kB = 1.380649e-23


@dataclass
class Params(base.Params):
    Lambda_ref: float = 10.0
    grain_drive_coeff: float = 0.0
    grain_drive_exp: float = 1.0
    grain_drive_ref: float = 100e-9
    lambda_floor: float = 1e-12
    baseline_stress_scale: float = 1.0
    concentration_stress_scale: float = 0.0
    baseline_radius_mode: str = "grain"
    grain_radius_factor: float = 0.5
    baseline_stress_cap: float = 1.0e12


RampHoldCool = base.RampHoldCool
Iso = base.Iso
TwoStep = base.TwoStep


def _activity_from_lambda(Lam):
    Lam = np.maximum(Lam, 1e-300)
    return Lam / (1.0 + Lam)


def _lambda_from_activity(a):
    a = np.clip(a, 1e-12, 1.0 - 1e-12)
    return a / (1.0 - a)


def grain_baseline_radius(p: Params, s) -> float:
    G = max(float(getattr(s, "G", getattr(p, "G0", 1.0))), 1e-30)
    return max(float(getattr(p, "grain_radius_factor", 0.5)) * G, 1e-30)


def grain_sigma_base(p: Params, G: float) -> float:
    gamma_s = float(getattr(p, "gamma_s", getattr(p, "gamma", 1.0)))
    Rg = max(float(getattr(p, "grain_radius_factor", 0.5)) * max(G, 1e-30), 1e-30)
    sig = float(getattr(p, "baseline_stress_scale", 1.0)) * 2.0 * gamma_s / Rg
    return min(max(sig, 0.0), float(getattr(p, "baseline_stress_cap", 1e12)))


def grain_lambda_factor(p: Params, G: float) -> float:
    c = max(float(getattr(p, "grain_drive_coeff", 0.0)), 0.0)
    if c <= 0:
        return 1.0
    m = max(float(getattr(p, "grain_drive_exp", 1.0)), 1e-12)
    Gref = max(float(getattr(p, "grain_drive_ref", 100e-9)), 1e-30)
    G = max(float(G), 1e-30)
    return max((1.0 + c * (Gref / G) ** m) / (1.0 + c), 1e-12)


def _pore_radius_from_state_or_dict(p, s, d):
    for key in ["r_active", "r_eff", "r_remove", "r_removable", "r_mean", "r_vol"]:
        if key in d and np.isfinite(d[key]) and d[key] > 0:
            return float(d[key])
    if hasattr(s, "r_bins"):
        rb = np.asarray(s.r_bins, dtype=float)
        ph = np.asarray(getattr(s, "phi_bins", np.ones_like(rb)), dtype=float)
        if rb.size and np.sum(ph) > 0:
            return float(np.sum(ph * rb) / max(np.sum(ph), 1e-300))
    return max(float(getattr(p, "rp0", 25e-9)), 1e-30)


def _baseline_radius(p, s, d):
    mode = str(getattr(p, "baseline_radius_mode", "grain")).lower()
    if mode in ["grain", "g", "grainsize", "grain_size"]:
        return grain_baseline_radius(p, s)
    if mode in ["pore", "pore_active", "active"]:
        return _pore_radius_from_state_or_dict(p, s, d)
    return max(float(getattr(p, "rp0", 25e-9)), 1e-30)


def _wrap_compute_rates(original_compute_rates):
    def wrapped_compute_rates(p, s, T_C, nu_eff):
        d = original_compute_rates(p, s, T_C, nu_eff)
        sigma_old = max(float(d.get("sigma", 0.0)), 1e-12)
        T_K = float(d.get("T_K", T_C + 273.15))
        Vact = float(getattr(p, "Vact", 0.0))
        sigma_cap = float(getattr(p, "sigma_cap", 1e99))
        r_base = _baseline_radius(p, s, d)
        gamma_s = float(getattr(p, "gamma_s", getattr(p, "gamma", 1.0)))
        sigma_base = float(getattr(p, "baseline_stress_scale", 1.0)) * 2.0 * gamma_s / max(r_base, 1e-30)
        sigma_base = min(max(sigma_base, 0.0), float(getattr(p, "baseline_stress_cap", 1e12)))
        sigma_conc = max(float(getattr(p, "concentration_stress_scale", 0.0)), 0.0) * sigma_old
        sigma_eff_uncapped = sigma_base + sigma_conc
        sigma_eff = min(max(sigma_eff_uncapped, 1e-12), sigma_cap)

        if "Lambda" in d:
            Lam_old = float(d["Lambda"])
        elif "activity" in d:
            Lam_old = float(_lambda_from_activity(d["activity"]))
        else:
            return d
        a_old = float(d.get("activity", _activity_from_lambda(Lam_old)))

        dSigma = sigma_eff - sigma_old
        if Vact != 0:
            expo = np.clip(Vact * dSigma / (kB * T_K), -80.0, 80.0)
            nuc_stress_factor = float(np.exp(expo))
        else:
            nuc_stress_factor = 1.0
        tau_factor = sigma_old / max(sigma_eff, 1e-12)
        lambda_scale = max(float(getattr(p, "Lambda_ref", 10.0)), 1e-12) / 10.0
        G_now = max(float(getattr(s, "G", getattr(p, "G0", 1.0))), 1e-30)
        gfac = grain_lambda_factor(p, G_now)
        Lam_new = max(Lam_old * nuc_stress_factor * tau_factor * lambda_scale * gfac, float(getattr(p, "lambda_floor", 1e-12)))
        a_new = float(_activity_from_lambda(Lam_new))

        edot_old = float(d.get("edot0", 0.0))
        edot_scale = (sigma_eff / max(sigma_old, 1e-12)) * (a_new / max(a_old, 1e-300))
        edot_new = edot_old * edot_scale
        Pdens_old = float(d.get("P_dens0", 0.0))
        Pdens_new = Pdens_old * (sigma_eff / max(sigma_old, 1e-12)) * edot_scale
        Ploss = float(d.get("P_loss", 0.0))
        eta_floor = float(getattr(p, "eta_floor", 0.0))
        eta_new = max(eta_floor, Pdens_new / max(Pdens_new + Ploss, 1e-300))
        eta_old = float(d.get("eta_Ons", eta_new))
        rho_dot_old = float(d.get("rho_dot", 0.0))
        rho_dot_new = rho_dot_old * edot_scale * eta_new / max(eta_old, 1e-300)

        d.update({
            "sigma_old_v3": sigma_old,
            "sigma_base": sigma_base,
            "sigma_concentration": sigma_conc,
            "sigma_eff_uncapped": sigma_eff_uncapped,
            "sigma": sigma_eff,
            "baseline_radius": r_base,
            "grain_lambda_factor": gfac,
            "Lambda_raw_v3": Lam_old,
            "activity_raw_v3": a_old,
            "Lambda": Lam_new,
            "activity": a_new,
            "edot0": edot_new,
            "P_dens0": Pdens_new,
            "eta_Ons": eta_new,
            "rho_dot": rho_dot_new,
        })
        return d
    return wrapped_compute_rates


def run(params: Params, protocol, stop_at_rho=None):
    original_compute_rates = base.compute_rates
    base.compute_rates = _wrap_compute_rates(original_compute_rates)
    try:
        if stop_at_rho is None:
            return base.run(params, protocol)
        return base.run(params, protocol, stop_at_rho=stop_at_rho)
    finally:
        base.compute_rates = original_compute_rates
