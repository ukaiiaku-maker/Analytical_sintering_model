#!/usr/bin/env python3
"""Deterministic diagnostic for the grain-stress Lambda-window model."""
from __future__ import annotations

from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

import sinter_reference_model_v6_grainstress_multibin as sm

OUTDIR = Path("debug_lambda_v6_grainstress")
OUTDIR.mkdir(exist_ok=True)


def make_params(rho0=0.83, G0_nm=100, **kw):
    vals = dict(
        rho0=rho0,
        G0=G0_nm * 1e-9,
        Q_nuc=500e3,
        Q_SD_growth=539e3,
        Q_surf=539e3,
        Q_GG=425e3,
        Lambda_ref=40.0,
        grain_drive_coeff=12.0,
        grain_drive_exp=1.2,
        grain_drive_ref=100e-9,
        baseline_stress_scale=1.0,
        concentration_stress_scale=0.25,
        baseline_radius_mode="grain",
        grain_radius_factor=0.5,
        sigma_factor=2.0,
        sigma_cap=3e8,
        anchor_edot=4.0e-3,
        Vact=1.4e-28,
        rp0=25e-9,
        rmax_factor=20.0,
        pore_ln_sigma=0.55,
        phi_tail0_frac=0.01,
        tail0_center_frac=0.75,
        tail0_width_bins=1.2,
    )
    vals.update(kw)
    return sm.Params(**vals)


def idx_at(r, rho):
    idx = np.where(r["rho"] >= rho)[0]
    return int(idx[0]) if len(idx) else None


def med(r, key, lo, hi):
    if key not in r:
        return np.nan
    mask = (r["rho"] >= lo) & (r["rho"] <= hi)
    if not np.any(mask):
        return np.nan
    return float(np.nanmedian(r[key][mask]))


def sigma_base_from_G(p, G):
    return np.array([sm.grain_sigma_base(p, float(g)) for g in G])


def grain_factor_from_G(p, G):
    return np.array([sm.grain_lambda_factor(p, float(g)) for g in G])


def main():
    rows = []
    results = {}
    for G0_nm in [50, 100, 200, 500, 1000]:
        for T_C in [1150, 1250, 1300, 1350]:
            label = f"G0={G0_nm:g} nm, T={T_C:g} C"
            p = make_params(rho0=0.83, G0_nm=G0_nm)
            r = sm.run(p, sm.Iso(T_C, 48 * 3600), stop_at_rho=0.92)
            results[label] = (p, r)
            i = idx_at(r, 0.92)
            sig_base = sigma_base_from_G(p, r["G"])
            gfac = grain_factor_from_G(p, r["G"])
            rows.append({
                "G0_nm": G0_nm,
                "T_C": T_C,
                "reached_0p92": i is not None,
                "rho_final": float(r["rho"][-1]),
                "G_final_nm": float(r["G"][-1] * 1e9),
                "G_at_0p92_nm": float(r["G"][i] * 1e9) if i is not None else np.nan,
                "activity_med_083_092": med(r, "activity", 0.83, 0.92),
                "Lambda_med_083_092": med(r, "Lambda", 0.83, 0.92),
                "sigma_eff_med_MPa": med(r, "sigma", 0.83, 0.92) / 1e6,
                "sigma_base_from_G_med_MPa": float(np.nanmedian(sig_base)) / 1e6,
                "grain_factor_from_G_med": float(np.nanmedian(gfac)),
            })
    with open(OUTDIR / "second_step_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    selected = [k for k in results if k.startswith("G0=50") or k.startswith("G0=200") or k.startswith("G0=1000")]
    for label in selected:
        p, r = results[label]
        sig_base = sigma_base_from_G(p, r["G"])
        gfac = grain_factor_from_G(p, r["G"])
        ax[0,0].plot(r["t"]/3600, r["rho"], label=label)
        ax[0,1].plot(r["G"]*1e9, r["rho"], label=label)
        ax[0,2].plot(r["rho"], r["activity"], label=label)
        ax[1,0].plot(r["rho"], np.clip(r["Lambda"], 1e-12, 1e12), label=label)
        ax[1,1].plot(r["rho"], sig_base/1e6, label=label)
        ax[1,2].plot(r["rho"], gfac, label=label)
    titles = ["rho(t)", "trajectory", "activity", "Lambda", "sigma_base=2 gamma/R_G", "grain Lambda factor"]
    for a, title in zip(ax.flat, titles):
        a.set_title(title); a.grid(True, alpha=0.25)
    ax[1,0].set_yscale("log")
    ax[0,0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTDIR / "second_step_grainstress_window.png", dpi=170)
    print("Wrote", OUTDIR)


if __name__ == "__main__":
    main()
