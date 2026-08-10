#!/usr/bin/env python3
"""Bounded comparison of density and observable topology smoothing gates."""
from __future__ import annotations

import argparse
from dataclasses import fields, replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from density_window_processing_map import sample_at_density, write_csv

GATE_MODES = ("density", "fine_pore", "connectivity", "hybrid_topology")
RHO0_VALUES = (0.65, 0.70, 0.75, 0.80, 0.85)
TARGETS = (0.85, 0.88, 0.90, 0.92)
VARYING_FIELDS = {"rho0", "smoothing_gate_mode"}


def base_params() -> model.Params:
    return model.Params(memory_model="pore_bin_redistribution", rho0=.75, G0=150e-9,
                        pore_radius0=25e-9, pore_ln_sigma=.65,
                        smoothing_gate_mode="density")


def condition_params(base: model.Params, mode: str, rho0: float) -> model.Params:
    if mode not in GATE_MODES:
        raise ValueError(mode)
    return replace(base, smoothing_gate_mode=mode, rho0=rho0)


def assert_only_mode_and_rho_vary(parameter_sets: list[model.Params], base: model.Params) -> None:
    fixed = {f.name: getattr(base, f.name) for f in fields(base) if f.name not in VARYING_FIELDS}
    for params in parameter_sets:
        current = {f.name: getattr(params, f.name) for f in fields(params) if f.name not in VARYING_FIELDS}
        if current != fixed:
            raise AssertionError("mechanism parameters changed across gate modes")


def initial_observables(params: model.Params) -> dict:
    state = model.initial_state(params)
    pores = model.pore_distribution_diagnostics(state.pore_phi, state.pore_radii, params)
    return {
        "initial_removable_fine_pore_fraction": pores["removable_fine_pore_fraction"],
        "initial_f_pore": state.topology.f_pore,
        "initial_connectivity": state.topology.connectivity,
        "initial_connected_coverage": state.topology.f_pore * state.topology.connectivity,
        "initial_isolated_pore_fraction": state.topology.isolated_pore_fraction,
        "initial_large_pore_fraction": pores["large_pore_fraction"],
        "initial_gate_value": model.smoothing_topology_gate(state, params),
    }


def run_protocols(params: model.Params) -> dict[str, dict]:
    return {
        "slow": model.run(params, model.RampHold(.2)),
        "fast": model.run(params, model.RampHold(20.)),
        "high": model.run(params, model.Iso(1350.)),
        "two_step": model.run(params, model.TwoStep(1350., 1250., .85)),
    }


def rows_for(mode: str, rho0: float, params: model.Params, runs: dict[str, dict]) -> list[dict]:
    initial = initial_observables(params)
    rows = []
    for target in TARGETS:
        values = {name: sample_at_density(result, target, rho0) for name, result in runs.items()}
        hr_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("slow", "fast"))
        ts_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("high", "two_step"))
        row = {
            "condition_id": f"{mode}_rho0{rho0:.2f}", "smoothing_gate_mode": mode,
            "rho0": rho0, "target_density": target, "eligible_target": target > rho0 + 1e-12,
            "HR_scored": hr_scored, "TS_scored": ts_scored,
            "HR_pct": model.percent_gain(values["slow"]["G_at_target_nm"], values["fast"]["G_at_target_nm"]) if hr_scored else np.nan,
            "TS_pct": model.percent_gain(values["high"]["G_at_target_nm"], values["two_step"]["G_at_target_nm"]) if ts_scored else np.nan,
            **initial,
        }
        for name, sampled in values.items():
            row[f"{name}_reached"] = sampled["reached_target"]
            for metric, value in sampled.items():
                if metric not in ("eligible_target", "reached_target"):
                    row[f"{name}_{metric}"] = value
        rows.append(row)
    crossover = min((row["target_density"] for row in rows if row["HR_scored"] and row["HR_pct"] > 0), default=np.nan)
    for row in rows:
        row["HR_crossover_density"] = crossover
    return rows


def run_audit(base: model.Params) -> tuple[list[dict], dict]:
    settings = [(mode, rho0) for mode in GATE_MODES for rho0 in RHO0_VALUES]
    params_list = [condition_params(base, *setting) for setting in settings]
    assert_only_mode_and_rho_vary(params_list, base)
    rows, cache = [], {}
    for index, ((mode, rho0), params) in enumerate(zip(settings, params_list), 1):
        runs = run_protocols(params)
        cache[(mode, rho0)] = runs
        rows.extend(rows_for(mode, rho0, params, runs))
        print(f"[{index}/{len(settings)}] {mode} rho0={rho0:.2f}", flush=True)
    return rows, cache


def summary_rows(rows: list[dict]) -> list[dict]:
    output = []
    for mode in GATE_MODES:
        for target in TARGETS:
            group = [r for r in rows if r["smoothing_gate_mode"] == mode and r["target_density"] == target]
            hr = [r["HR_pct"] for r in group if r["HR_scored"]]
            ts = [r["TS_pct"] for r in group if r["TS_scored"]]
            output.append({
                "smoothing_gate_mode": mode, "target_density": target, "n_conditions": len(group),
                "fraction_HR_attainable": len(hr)/len(group), "fraction_TS_attainable": len(ts)/len(group),
                "fraction_HR_positive": float(np.mean(np.asarray(hr) > 0)) if hr else np.nan,
                "fraction_TS_positive": float(np.mean(np.asarray(ts) > 0)) if ts else np.nan,
                "mean_HR_pct": float(np.mean(hr)) if hr else np.nan,
                "mean_TS_pct": float(np.mean(ts)) if ts else np.nan,
            })
    return output


def make_plots(outdir: Path, rows: list[dict]) -> None:
    at90 = [r for r in rows if r["target_density"] == .90 and r["HR_scored"]]
    variables = (("initial_removable_fine_pore_fraction", "Initial removable fine-pore fraction"),
                 ("initial_connected_coverage", "Initial connected pore-boundary coverage"),
                 ("initial_gate_value", "Initial gate value"))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
    for axis, (key, label) in zip(axes, variables):
        for mode in GATE_MODES:
            group = [r for r in at90 if r["smoothing_gate_mode"] == mode]
            axis.plot([r[key] for r in group], [r["HR_pct"] for r in group], "o-", label=mode)
        axis.axhline(0, color="black", linestyle="--", linewidth=1)
        axis.set(xlabel=label, ylabel="HR_pct at rho=0.90")
        axis.grid(alpha=.25)
    axes[-1].legend(fontsize=8); fig.tight_layout()
    fig.savefig(outdir / "HR_pct_vs_topology_gate_variables.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for mode in GATE_MODES:
        group = [r for r in at90 if r["smoothing_gate_mode"] == mode]
        y = [r["slow_cumulative_redistributed_pore_volume"] for r in group]
        axes[0].plot([r["slow_removable_fine_pore_fraction"] for r in group], y, "o-", label=mode)
        axes[1].plot([r["slow_f_pore"] for r in group], y, "o-", label=mode)
    axes[0].set(xlabel="Slow-ramp removable fine-pore fraction at rho=0.90", ylabel="Cumulative redistribution")
    axes[1].set(xlabel="Slow-ramp f_pore at rho=0.90", ylabel="Cumulative redistribution")
    for axis in axes: axis.grid(alpha=.25)
    axes[1].legend(fontsize=8); fig.tight_layout()
    fig.savefig(outdir / "redistribution_vs_fine_pore_and_f_pore.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for mode in GATE_MODES:
        group = [r for r in rows if r["smoothing_gate_mode"] == mode]
        for rho0 in RHO0_VALUES:
            subset = [r for r in group if r["rho0"] == rho0]
            axes[0].plot([r["target_density"] for r in subset], [r["HR_pct"] for r in subset], marker="o", alpha=.65,
                         label=f"{mode}, rho0={rho0:.2f}" if rho0 == .75 else None)
            axes[1].plot([r["target_density"] for r in subset], [r["TS_pct"] for r in subset], marker="o", alpha=.65,
                         label=f"{mode}, rho0={rho0:.2f}" if rho0 == .75 else None)
    for axis, metric in zip(axes, ("HR_pct", "TS_pct")):
        axis.axhline(0, color="black", linestyle="--", linewidth=1); axis.set(xlabel="Target density", ylabel=metric); axis.grid(alpha=.25)
    axes[1].legend(fontsize=7); fig.tight_layout()
    fig.savefig(outdir / "HR_TS_density_window_by_gate_mode.png", dpi=150); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/topology_gate_identifiability")
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rows, _ = run_audit(base_params())
    unattainable = [r for r in rows if r["eligible_target"] and (not r["HR_scored"] or not r["TS_scored"])]
    write_csv(outdir / "gate_mode_results.csv", rows)
    write_csv(outdir / "gate_mode_summary.csv", summary_rows(rows))
    write_csv(outdir / "unattainable_cases.csv", unattainable)
    make_plots(outdir, rows)
    print(f"conditions={len(rows)//len(TARGETS)} rows={len(rows)} unattainable={len(unattainable)}")


if __name__ == "__main__":
    main()
