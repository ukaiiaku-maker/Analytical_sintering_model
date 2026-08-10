#!/usr/bin/env python3
"""Audit identifiability and robustness of the smoothing density gate."""
from __future__ import annotations

import argparse
import csv
from dataclasses import fields, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from density_window_processing_map import sample_at_density, write_csv

RHO_MIDS = (0.72, 0.76, 0.79, 0.82, 0.86)
WIDTHS = {"narrow": 0.0075, "baseline": 0.015, "broad": 0.03}
RHO0_VALUES = (0.65, 0.70, 0.75, 0.80, 0.85)
TARGETS = (0.85, 0.88, 0.90, 0.92)
GATE_FIELDS = {"rho0", "smoothing_rho_mid", "smoothing_rho_width", "smoothing_gate_form"}


def base_params() -> model.Params:
    return model.Params(memory_model="pore_bin_redistribution", rho0=0.75,
                        G0=150e-9, pore_radius0=25e-9, pore_ln_sigma=0.65)


def condition_params(base: model.Params, rho0: float, rho_mid: float, rho_width: float, gate_form: str="logistic") -> model.Params:
    return replace(base, rho0=rho0, smoothing_rho_mid=rho_mid,
                   smoothing_rho_width=rho_width, smoothing_gate_form=gate_form)


def assert_only_gate_and_rho_vary(parameter_sets: list[model.Params], base: model.Params) -> None:
    reference = {field.name: getattr(base, field.name) for field in fields(base) if field.name not in GATE_FIELDS}
    for params in parameter_sets:
        current = {field.name: getattr(params, field.name) for field in fields(params) if field.name not in GATE_FIELDS}
        if current != reference:
            raise AssertionError("non-gate mechanism parameter changed")


def run_protocols(params: model.Params) -> dict[str, dict]:
    return {
        "slow": model.run(params, model.RampHold(0.2)),
        "fast": model.run(params, model.RampHold(20.0)),
        "high": model.run(params, model.Iso(1350.0)),
        "two_step": model.run(params, model.TwoStep(1350.0, 1250.0, 0.85)),
    }


def condition_rows(condition_id: str, params: model.Params, width_label: str, runs: dict[str, dict]) -> list[dict]:
    rows = []
    for target in TARGETS:
        values = {name: sample_at_density(result, target, params.rho0) for name, result in runs.items()}
        hr_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("slow", "fast"))
        ts_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("high", "two_step"))
        row = {
            "condition_id": condition_id, "gate_form": params.smoothing_gate_form,
            "rho0": params.rho0, "smoothing_rho_mid": params.smoothing_rho_mid,
            "rho0_minus_rho_mid": params.rho0 - params.smoothing_rho_mid,
            "rho_width_label": width_label, "smoothing_rho_width": params.smoothing_rho_width,
            "target_density": target, "eligible_target": target > params.rho0 + 1e-12,
            "HR_scored": hr_scored, "TS_scored": ts_scored,
            "HR_pct": model.percent_gain(values["slow"]["G_at_target_nm"], values["fast"]["G_at_target_nm"]) if hr_scored else np.nan,
            "TS_pct": model.percent_gain(values["high"]["G_at_target_nm"], values["two_step"]["G_at_target_nm"]) if ts_scored else np.nan,
            "slow_time_budget_h": min(model.RampHold(0.2).t_end, params.t_max_s) / 3600,
            "fast_time_budget_h": min(model.RampHold(20.0).t_end, params.t_max_s) / 3600,
            "high_time_budget_h": min(model.Iso(1350.0).t_end, params.t_max_s) / 3600,
            "two_step_time_budget_h": min(model.TwoStep(1350.0, 1250.0, 0.85).t_end, params.t_max_s) / 3600,
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


def run_logistic_grid(base: model.Params) -> tuple[list[dict], dict]:
    rows, cache = [], {}
    settings = [(mid, label, width, rho0) for mid in RHO_MIDS for label, width in WIDTHS.items() for rho0 in RHO0_VALUES]
    parameter_sets = [condition_params(base, rho0, mid, width) for mid, _, width, rho0 in settings]
    assert_only_gate_and_rho_vary(parameter_sets, base)
    for index, ((mid, label, width, rho0), params) in enumerate(zip(settings, parameter_sets), 1):
        condition_id = f"logistic_mid{mid:.2f}_{label}_rho0{rho0:.2f}"
        runs = run_protocols(params)
        cache[(mid, label, rho0)] = runs
        rows.extend(condition_rows(condition_id, params, label, runs))
        print(f"[{index}/{len(settings)}] {condition_id}", flush=True)
    return rows, cache


def run_form_comparison(base: model.Params, logistic_rows: list[dict]) -> list[dict]:
    rows = [dict(row) for row in logistic_rows if row["smoothing_rho_mid"] == 0.79 and row["rho_width_label"] == "baseline"]
    parameter_sets = [condition_params(base, rho0, 0.79, WIDTHS["baseline"], "linear_clipped") for rho0 in RHO0_VALUES]
    assert_only_gate_and_rho_vary(parameter_sets, base)
    for params in parameter_sets:
        condition_id = f"linear_mid0.79_baseline_rho0{params.rho0:.2f}"
        rows.extend(condition_rows(condition_id, params, "baseline", run_protocols(params)))
    return rows


def sensitivity_summary(rows: list[dict]) -> list[dict]:
    summary = []
    for mid in RHO_MIDS:
        for label in WIDTHS:
            for target in TARGETS:
                group = [row for row in rows if row["smoothing_rho_mid"] == mid and row["rho_width_label"] == label and row["target_density"] == target]
                hr = [row["HR_pct"] for row in group if row["HR_scored"]]
                ts = [row["TS_pct"] for row in group if row["TS_scored"]]
                crossovers = [row["HR_crossover_density"] for row in group if np.isfinite(row["HR_crossover_density"])]
                redistributed = [row["slow_cumulative_redistributed_pore_volume"] for row in group if row["slow_reached"]]
                summary.append({
                    "smoothing_rho_mid": mid, "rho_width_label": label,
                    "smoothing_rho_width": WIDTHS[label], "target_density": target,
                    "n_rho0": len(group), "fraction_HR_attainable": len(hr) / len(group),
                    "fraction_TS_attainable": len(ts) / len(group),
                    "fraction_HR_positive": np.mean(np.asarray(hr) > 0) if hr else np.nan,
                    "fraction_TS_positive": np.mean(np.asarray(ts) > 0) if ts else np.nan,
                    "mean_HR_pct": np.mean(hr) if hr else np.nan,
                    "mean_TS_pct": np.mean(ts) if ts else np.nan,
                    "median_HR_crossover_density": np.median(crossovers) if crossovers else np.nan,
                    "mean_slow_cumulative_redistribution": np.mean(redistributed) if redistributed else np.nan,
                })
    return summary


def plot_calibration_views(outdir: Path, rows: list[dict], cache: dict) -> None:
    at_090 = [row for row in rows if row["target_density"] == 0.90 and row["eligible_target"]]
    styles = {"narrow": "^", "baseline": "o", "broad": "s"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for label in WIDTHS:
        group = [row for row in at_090 if row["rho_width_label"] == label]
        axes[0].scatter([row["rho0_minus_rho_mid"] for row in group], [row["slow_cumulative_redistributed_pore_volume"] for row in group], marker=styles[label], label=label, alpha=.75)
        axes[1].scatter([row["rho0_minus_rho_mid"] for row in group], [row["HR_pct"] for row in group], marker=styles[label], label=label, alpha=.75)
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[0].set(xlabel="rho0 - smoothing_rho_mid", ylabel="Slow cumulative redistribution", title="Redistribution relative to gate")
    axes[1].set(xlabel="rho0 - smoothing_rho_mid", ylabel="HR_pct at rho=0.90", title="Fast-heating response relative to gate")
    for axis in axes: axis.grid(alpha=.25); axis.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outdir / "gate_relative_redistribution_and_HR.png", dpi=150); plt.close(fig)

    fig, axis = plt.subplots(figsize=(8, 5.5))
    ramp = model.RampHold(0.2)
    for rho0 in RHO0_VALUES:
        result = cache[(0.79, "baseline", rho0)]["slow"]
        mask = result["t"] <= ramp.ramp_s + 1e-9
        axis.plot(result["T_C"][mask], result["pore_mean_radius"][mask] * 1e9, label=f"rho0={rho0:.2f}")
    axis.set(xlabel="Temperature [C]", ylabel="Predicted pore mean radius [nm]", title="Interrupted slow-ramp observable")
    axis.grid(alpha=.25); axis.legend(fontsize=8); fig.tight_layout()
    fig.savefig(outdir / "interrupted_ramp_mean_radius_vs_temperature.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for axis, label in zip(axes, WIDTHS):
        matrix = np.full((len(RHO_MIDS), len(RHO0_VALUES)), np.nan)
        for i, mid in enumerate(RHO_MIDS):
            for j, rho0 in enumerate(RHO0_VALUES):
                row = next(row for row in at_090 if row["smoothing_rho_mid"] == mid and row["rho_width_label"] == label and row["rho0"] == rho0)
                matrix[i, j] = row["HR_pct"]
        image = axis.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-20, vmax=20, interpolation="nearest")
        axis.set_xticks(range(len(RHO0_VALUES)), [f"{value:.2f}" for value in RHO0_VALUES])
        axis.set_yticks(range(len(RHO_MIDS)), [f"{value:.2f}" for value in RHO_MIDS])
        axis.set(xlabel="rho0", title=label)
    axes[0].set_ylabel("smoothing_rho_mid")
    fig.colorbar(image, ax=axes, label="HR_pct at rho=0.90")
    fig.subplots_adjust(wspace=.18, right=.90)
    fig.savefig(outdir / "HR_gate_center_width_map.png", dpi=150); plt.close(fig)

    fig, axis = plt.subplots(figsize=(8, 5.5))
    for label in WIDTHS:
        group = [row for row in at_090 if row["rho_width_label"] == label]
        axis.scatter([row["rho0_minus_rho_mid"] for row in group], [row["TS_pct"] for row in group], marker=styles[label], label=label, alpha=.75)
    axis.axhline(0, color="black", linestyle="--", linewidth=1)
    axis.set(xlabel="rho0 - smoothing_rho_mid", ylabel="TS_pct at rho=0.90", title="Two-step response relative to gate")
    axis.grid(alpha=.25); axis.legend(fontsize=8); fig.tight_layout()
    fig.savefig(outdir / "TS_vs_gate_relative_density.png", dpi=150); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/smoothing_gate_identifiability")
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    base = base_params()
    rows, cache = run_logistic_grid(base)
    form_rows = run_form_comparison(base, rows)
    summary = sensitivity_summary(rows)
    unattainable = [row for row in rows if row["eligible_target"] and (not row["HR_scored"] or not row["TS_scored"])]
    write_csv(outdir / "gate_sensitivity_results.csv", rows)
    write_csv(outdir / "gate_form_comparison.csv", form_rows)
    write_csv(outdir / "gate_sensitivity_summary.csv", summary)
    write_csv(outdir / "unattainable_cases.csv", unattainable)
    plot_calibration_views(outdir, rows, cache)
    print(f"gate conditions={len(rows)//len(TARGETS)} rows={len(rows)} unattainable={len(unattainable)}")


if __name__ == "__main__":
    main()
