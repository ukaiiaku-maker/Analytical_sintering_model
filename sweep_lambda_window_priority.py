#!/usr/bin/env python3
"""
sweep_lambda_window_priority.py

Automated random search for a single reduced-order sintering formulation that can reproduce:
1. fast-heating-rate improvement at fixed density; and
2. two-step improvement at fixed density.

The score is percent based, not raw nm difference, and includes a renewal activity-window diagnostic.
"""
from __future__ import annotations

from dataclasses import fields, is_dataclass
import argparse
import csv
import importlib
import math
import random
from pathlib import Path
from typing import Any, Dict

import matplotlib.pyplot as plt
import numpy as np

R = 8.31446261815324

RANGES = {
    "Q_nuc_kJ": (400.0, 720.0, "linear"),
    "Q_SD_kJ": (450.0, 720.0, "linear"),
    "Q_GG_kJ": (360.0, 625.0, "linear"),
    "Lambda_ref": (0.3, 50.0, "log"),
    "grain_drive_coeff": (0.0, 80.0, "linear"),
    "grain_drive_exp": (0.5, 3.0, "linear"),
    "grain_drive_ref_nm": (40.0, 250.0, "linear"),
    "baseline_stress_scale": (0.0, 2.5, "linear"),
    "concentration_stress_scale": (0.0, 1.5, "linear"),
    "grain_radius_factor": (0.35, 0.75, "linear"),
    "anchor_edot_scale": (0.15, 12.0, "log"),
    "Vact_scale": (0.4, 3.0, "log"),
    "sigma_factor": (0.8, 4.5, "linear"),
    "sigma_cap_MPa": (120.0, 900.0, "linear"),
    "kSD_scale": (0.03, 30.0, "log"),
    "kGG_scale": (0.03, 30.0, "log"),
    "rp0_nm": (10.0, 70.0, "linear"),
    "rmax_factor": (8.0, 45.0, "linear"),
    "pore_ln_sigma": (0.20, 1.25, "linear"),
    "phi_tail0_frac": (2e-4, 0.15, "log"),
    "coverage_chi": (0.5, 7.0, "linear"),
    "zener_k": (1.5, 18.0, "linear"),
    "mobile_T_mid_C": (1100.0, 1320.0, "linear"),
    "mobile_T_width_C": (30.0, 140.0, "linear"),
    "loss_SD_coeff": (0.0, 0.45, "linear"),
    "loss_drag_coeff": (0.3, 50.0, "log"),
    "loss_clean_coeff": (0.1, 35.0, "log"),
    "PR_coeff": (1e-5, 2e-1, "log"),
    "PR_activity_exp": (0.35, 4.0, "linear"),
    "PR_rho_mid": (0.74, 0.89, "linear"),
    "pore_drag_coeff": (0.03, 4.0, "log"),
    "drag_damage_coeff": (0.005, 5.0, "log"),
    "clean_damage_coeff": (0.002, 2.5, "log"),
    "coalesce_coeff": (5e-4, 1.0, "log"),
    "tail_drag_radius_exp": (0.0, 3.5, "linear"),
    "tail_mobile_boost_C": (0.0, 260.0, "linear"),
}


def sample_value(lo: float, hi: float, mode: str, rng: random.Random) -> float:
    if mode == "log":
        return math.exp(rng.uniform(math.log(lo), math.log(hi)))
    return rng.uniform(lo, hi)


def allowed_fields(model) -> set[str]:
    if is_dataclass(model.Params):
        return {f.name for f in fields(model.Params)}
    return set(model.Params.__init__.__annotations__.keys())


def maybe(kwargs: Dict[str, Any], allowed: set[str], key: str, value: Any) -> None:
    if key in allowed:
        kwargs[key] = value


def anchored_prefactor(A_base: float, Q_base: float, Q_new: float, T_ref_K: float = 1375.0 + 273.15) -> float:
    return A_base * math.exp((Q_new - Q_base) / (R * T_ref_K))


def make_params(model, allowed: set[str], s: Dict[str, float], G0_nm: float, rho0: float):
    qsd = s["Q_SD_kJ"] * 1000.0
    qgg = s["Q_GG_kJ"] * 1000.0
    kwargs: Dict[str, Any] = {}
    maybe(kwargs, allowed, "rho0", rho0)
    maybe(kwargs, allowed, "G0", G0_nm * 1e-9)
    maybe(kwargs, allowed, "Q_nuc", s["Q_nuc_kJ"] * 1000.0)
    maybe(kwargs, allowed, "Q_SD_growth", qsd)
    maybe(kwargs, allowed, "Q_surf", qsd)
    maybe(kwargs, allowed, "Q_GG", qgg)
    maybe(kwargs, allowed, "kSD0_m3_s", anchored_prefactor(2.5e-8, 539e3, qsd) * s["kSD_scale"])
    maybe(kwargs, allowed, "kGG0_m3_s", anchored_prefactor(6.1e-7, 425e3, qgg) * s["kGG_scale"])
    maybe(kwargs, allowed, "anchor_edot", 1.5e-3 * s["anchor_edot_scale"])
    maybe(kwargs, allowed, "Vact", 1.4e-28 * s["Vact_scale"])
    maybe(kwargs, allowed, "sigma_factor", s["sigma_factor"])
    maybe(kwargs, allowed, "sigma_cap", s["sigma_cap_MPa"] * 1e6)
    maybe(kwargs, allowed, "rp0", s["rp0_nm"] * 1e-9)
    for key in ["rmax_factor", "pore_ln_sigma", "phi_tail0_frac", "coverage_chi", "zener_k", "mobile_T_mid_C", "mobile_T_width_C", "loss_SD_coeff", "loss_drag_coeff", "loss_clean_coeff", "PR_coeff", "PR_activity_exp", "PR_rho_mid", "pore_drag_coeff", "drag_damage_coeff", "clean_damage_coeff", "coalesce_coeff", "tail_drag_radius_exp", "tail_mobile_boost_C", "Lambda_ref", "grain_drive_coeff", "grain_drive_exp", "baseline_stress_scale", "concentration_stress_scale", "grain_radius_factor"]:
        maybe(kwargs, allowed, key, s[key])
    maybe(kwargs, allowed, "grain_drive_ref", s["grain_drive_ref_nm"] * 1e-9)
    if "baseline_radius_mode" in allowed:
        kwargs["baseline_radius_mode"] = "grain"
    return model.Params(**kwargs)


def g_at_rho(run: Dict[str, np.ndarray], target: float):
    idx = np.where(run["rho"] >= target)[0]
    if len(idx) == 0:
        return np.nan, np.nan
    i = int(idx[0])
    return float(run["G"][i] * 1e9), float(run["t"][i] / 3600.0)


def activity_array(run: Dict[str, np.ndarray]) -> np.ndarray:
    if "activity" in run:
        return np.asarray(run["activity"], dtype=float)
    L = np.asarray(run["Lambda"], dtype=float)
    return L / (1.0 + L)


def window_fraction(run: Dict[str, np.ndarray], lo: float, hi: float, activity_min: float) -> float:
    rho = np.asarray(run["rho"], dtype=float)
    a = activity_array(run)
    mask = (rho >= lo) & (rho <= hi)
    if not np.any(mask):
        return np.nan
    return float(np.mean(a[mask] >= activity_min))


def pct_gain(reference: float, improved: float) -> float:
    if not np.isfinite(reference) or not np.isfinite(improved) or reference <= 0:
        return np.nan
    return 100.0 * (reference - improved) / reference


def run_one(model, allowed: set[str], sample: Dict[str, float], args, sample_id: int) -> Dict[str, Any]:
    p_heat = make_params(model, allowed, sample, args.heat_G0_nm, args.rho0)
    slow = model.run(p_heat, model.RampHoldCool(0.2, soak_T_C=args.hr_soak_T_C, soak_time_s=args.hr_soak_h * 3600), stop_at_rho=args.rho_target)
    fast = model.run(p_heat, model.RampHoldCool(20.0, soak_T_C=args.hr_soak_T_C, soak_time_s=args.hr_soak_h * 3600), stop_at_rho=args.rho_target)
    G_slow, t_slow = g_at_rho(slow, args.rho_target)
    G_fast, t_fast = g_at_rho(fast, args.rho_target)
    HR_pct = pct_gain(G_slow, G_fast)

    p_two = make_params(model, allowed, sample, args.two_G0_nm, args.rho0)
    high = model.run(p_two, model.Iso(args.two_high_C, args.two_budget_h * 3600), stop_at_rho=args.rho_target)
    two = model.run(p_two, model.TwoStep(args.two_high_C, args.two_low_C, args.two_switch_rho, args.two_budget_h * 3600), stop_at_rho=args.rho_target)
    G_high, t_high = g_at_rho(high, args.rho_target)
    G_two, t_two = g_at_rho(two, args.rho_target)
    TS_pct = pct_gain(G_high, G_two)

    valid = np.all(np.isfinite([G_slow, G_fast, G_high, G_two])) and max(G_slow, G_fast, G_high, G_two) <= args.max_g_nm
    F_two = window_fraction(two, args.two_switch_rho, args.rho_target, args.activity_min)
    F_high = window_fraction(high, args.two_switch_rho, args.rho_target, args.activity_min)
    base_score = min(HR_pct, TS_pct) if valid and np.isfinite(HR_pct) and np.isfinite(TS_pct) else np.nan
    lambda_contrast = F_two - F_high if np.isfinite(F_two) and np.isfinite(F_high) else np.nan
    priority_score = base_score + args.lambda_weight * 100.0 * lambda_contrast if np.isfinite(base_score) and np.isfinite(lambda_contrast) else base_score

    row = {
        "sample_id": sample_id,
        "G_slow_nm": G_slow,
        "G_fast_nm": G_fast,
        "HR_gain_pct": HR_pct,
        "t_slow_h": t_slow,
        "t_fast_h": t_fast,
        "G_high_nm": G_high,
        "G_two_nm": G_two,
        "TS_gain_pct": TS_pct,
        "t_high_h": t_high,
        "t_two_h": t_two,
        "base_combined_score_pct": base_score,
        "F_lambda_two": F_two,
        "F_lambda_high": F_high,
        "lambda_window_contrast": lambda_contrast,
        "lambda_priority_score": priority_score,
        "passes_both": bool(valid and np.isfinite(HR_pct) and np.isfinite(TS_pct) and HR_pct > 0 and TS_pct > 0),
    }
    row.update(sample)
    return row


def save_csv(rows, path: Path):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def plot_summary(rows, outdir: Path):
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    x = np.array([r["HR_gain_pct"] for r in rows], dtype=float)
    y = np.array([r["TS_gain_pct"] for r in rows], dtype=float)
    c = np.array([r["lambda_window_contrast"] for r in rows], dtype=float)
    sc = ax[0].scatter(x, y, c=c, s=12, alpha=0.6)
    ax[0].axhline(0, color="k", ls="--", lw=1)
    ax[0].axvline(0, color="k", ls="--", lw=1)
    ax[0].set_xlabel("Heating-rate gain [%]")
    ax[0].set_ylabel("Two-step gain [%]")
    fig.colorbar(sc, ax=ax[0], label="F_two - F_high")
    ax[1].scatter([r["F_lambda_high"] for r in rows], [r["F_lambda_two"] for r in rows], s=12, alpha=0.6)
    ax[1].plot([0, 1], [0, 1], "k--", lw=1)
    ax[1].set_xlabel("High-T activity-window fraction")
    ax[1].set_ylabel("Two-step activity-window fraction")
    fig.tight_layout()
    fig.savefig(outdir / "lambda_search_summary.png", dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="sinter_reference_model_v6_grainstress_multibin")
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--outdir", default="sweep_lambda_window_priority")
    ap.add_argument("--rho0", type=float, default=0.74)
    ap.add_argument("--rho-target", type=float, default=0.92)
    ap.add_argument("--max-g-nm", type=float, default=4000.0)
    ap.add_argument("--heat-G0-nm", type=float, default=200.0)
    ap.add_argument("--two-G0-nm", type=float, default=50.0)
    ap.add_argument("--hr-soak-T-C", type=float, default=1500.0)
    ap.add_argument("--hr-soak-h", type=float, default=1.0)
    ap.add_argument("--two-high-C", type=float, default=1350.0)
    ap.add_argument("--two-low-C", type=float, default=1300.0)
    ap.add_argument("--two-switch-rho", type=float, default=0.83)
    ap.add_argument("--two-budget-h", type=float, default=96.0)
    ap.add_argument("--activity-min", type=float, default=0.80)
    ap.add_argument("--lambda-weight", type=float, default=0.50)
    args = ap.parse_args()

    model = importlib.import_module(args.model)
    allowed = allowed_fields(model)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    rows = []
    for i in range(args.n):
        sample = {k: sample_value(lo, hi, mode, rng) for k, (lo, hi, mode) in RANGES.items()}
        print(f"[{i + 1}/{args.n}]")
        try:
            rows.append(run_one(model, allowed, sample, args, i))
        except Exception as exc:
            row = {"sample_id": i, "error": repr(exc)}
            row.update(sample)
            rows.append(row)
            print("ERROR", repr(exc))
    full = [r for r in rows if "lambda_priority_score" in r]
    save_csv(rows, outdir / "summary_all.csv")
    save_csv(full, outdir / "summary_lambda_priority.csv")
    good = sorted([r for r in full if r["passes_both"]], key=lambda r: r["lambda_priority_score"], reverse=True)
    save_csv(good, outdir / "good_lambda_cases.csv")
    plot_summary(full, outdir)
    print("completed", len(full), "good", len(good), "outputs", outdir)


if __name__ == "__main__":
    main()
