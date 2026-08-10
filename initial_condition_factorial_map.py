#!/usr/bin/env python3
"""Bounded factorial/OAT audit of initial microstructure descriptors."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, fields, replace
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from density_window_processing_map import sample_at_density, write_csv

TARGETS = (0.80, 0.85, 0.88, 0.90, 0.92)
INITIAL_FIELDS = {"rho0", "G0", "pore_radius0", "pore_ln_sigma"}


@dataclass(frozen=True)
class DesignPoint:
    design_id: str
    design_type: str
    varied_descriptor: str
    rho0: float
    G0_nm: float
    pore_scale_nm: float
    log_width: float

    def params(self, base: model.Params) -> model.Params:
        return replace(base, rho0=self.rho0, G0=self.G0_nm * 1e-9,
                       pore_radius0=self.pore_scale_nm * 1e-9,
                       pore_ln_sigma=self.log_width)


BASELINE = DesignPoint("oat_baseline", "one_at_a_time", "baseline", 0.75, 150.0, 25.0, 0.65)


def oat_design() -> list[DesignPoint]:
    points = [BASELINE]
    for rho0 in (0.65, 0.70, 0.80, 0.85):
        points.append(DesignPoint(f"oat_rho0_{rho0:.2f}", "one_at_a_time", "rho0", rho0, 150, 25, 0.65))
    for grain in (75.0, 300.0):
        points.append(DesignPoint(f"oat_G0_{grain:g}", "one_at_a_time", "G0_nm", 0.75, grain, 25, 0.65))
    for pore in (15.0, 40.0):
        points.append(DesignPoint(f"oat_pore_{pore:g}", "one_at_a_time", "pore_scale_nm", 0.75, 150, pore, 0.65))
    for width in (0.45, 0.85):
        points.append(DesignPoint(f"oat_width_{width:.2f}", "one_at_a_time", "log_width", 0.75, 150, 25, width))
    return points


def factorial_design() -> list[DesignPoint]:
    points = []
    for index, (rho0, grain, pore, width) in enumerate(product((0.65, 0.85), (75.0, 300.0), (15.0, 40.0), (0.45, 0.85))):
        points.append(DesignPoint(f"factorial_{index:02d}", "corner_factorial", "factorial", rho0, grain, pore, width))
    return points


def all_design_points() -> list[DesignPoint]:
    return oat_design() + factorial_design()


def assert_shared_mechanism_parameters(points: list[DesignPoint], base: model.Params) -> None:
    reference = {field.name: getattr(base, field.name) for field in fields(base) if field.name not in INITIAL_FIELDS}
    for point in points:
        params = point.params(base)
        current = {field.name: getattr(params, field.name) for field in fields(params) if field.name not in INITIAL_FIELDS}
        if current != reference:
            raise AssertionError(f"mechanism parameter drift in {point.design_id}")


def run_point(point: DesignPoint, base: model.Params) -> dict[str, dict]:
    params = point.params(base)
    return {
        "slow": model.run(params, model.RampHold(0.2)),
        "fast": model.run(params, model.RampHold(20.0)),
        "high": model.run(params, model.Iso(1350.0)),
        "two_step": model.run(params, model.TwoStep(1350.0, 1250.0, 0.85)),
    }


def initial_diagnostics(point: DesignPoint, base: model.Params) -> dict:
    params = point.params(base)
    state = model.initial_state(params)
    pores = model.pore_distribution_diagnostics(state.pore_phi, state.pore_radii, params)
    return {
        "initial_f_pore": state.topology.f_pore,
        "initial_connectivity": state.topology.connectivity,
        "initial_isolated_fraction": state.topology.isolated_pore_fraction,
        "initial_pore_mean_radius_nm": pores["pore_mean_radius"] * 1e9,
        "initial_large_pore_fraction": pores["large_pore_fraction"],
        "initial_removable_fine_fraction": pores["removable_fine_pore_fraction"],
    }


def result_rows(point: DesignPoint, runs: dict[str, dict], base: model.Params) -> list[dict]:
    rows = []
    initial = initial_diagnostics(point, base)
    for target in TARGETS:
        values = {name: sample_at_density(result, target, point.rho0) for name, result in runs.items()}
        hr_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("slow", "fast"))
        ts_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("high", "two_step"))
        row = {
            "design_id": point.design_id, "design_type": point.design_type,
            "varied_descriptor": point.varied_descriptor,
            "rho0": point.rho0, "G0_nm": point.G0_nm,
            "pore_scale_nm": point.pore_scale_nm, "log_width": point.log_width,
            "target_density": target, "eligible_target": target > point.rho0 + 1e-12,
            "high_step_executed": point.rho0 < 0.85 - 1e-12,
            "slow_time_budget_h": min(model.RampHold(0.2).t_end, base.t_max_s) / 3600.0,
            "fast_time_budget_h": min(model.RampHold(20.0).t_end, base.t_max_s) / 3600.0,
            "high_time_budget_h": min(model.Iso(1350.0).t_end, base.t_max_s) / 3600.0,
            "two_step_time_budget_h": min(model.TwoStep(1350.0, 1250.0, 0.85).t_end, base.t_max_s) / 3600.0,
            "HR_pct": model.percent_gain(values["slow"]["G_at_target_nm"], values["fast"]["G_at_target_nm"]) if hr_scored else np.nan,
            "TS_pct": model.percent_gain(values["high"]["G_at_target_nm"], values["two_step"]["G_at_target_nm"]) if ts_scored else np.nan,
            "HR_scored": hr_scored, "TS_scored": ts_scored,
            **initial,
        }
        for name, sampled in values.items():
            row[f"{name}_reached"] = sampled["reached_target"]
            for metric, value in sampled.items():
                if metric not in ("eligible_target", "reached_target"):
                    row[f"{name}_{metric}"] = value
        rows.append(row)
    return rows


def run_design(points: list[DesignPoint], base: model.Params) -> list[dict]:
    rows = []
    for index, point in enumerate(points, 1):
        rows.extend(result_rows(point, run_point(point, base), base))
        print(f"[{index}/{len(points)}] {point.design_id}", flush=True)
    return rows


def crossover_density(rows: list[dict], design_id: str) -> float:
    candidates = [row["target_density"] for row in rows if row["design_id"] == design_id and row["HR_scored"] and row["HR_pct"] > 0]
    return min(candidates) if candidates else np.nan


def level_for(point: DesignPoint, descriptor: str) -> float:
    return getattr(point, descriptor)


def sensitivity_summary(rows: list[dict], points: list[DesignPoint]) -> list[dict]:
    summary = []
    descriptor_levels = {
        "rho0": sorted({point.rho0 for point in points}),
        "G0_nm": sorted({point.G0_nm for point in points}),
        "pore_scale_nm": sorted({point.pore_scale_nm for point in points}),
        "log_width": sorted({point.log_width for point in points}),
    }
    for design_type in ("one_at_a_time", "corner_factorial"):
        design_points = [point for point in points if point.design_type == design_type]
        for descriptor, levels in descriptor_levels.items():
            eligible_points = [point for point in design_points if design_type != "one_at_a_time" or point.varied_descriptor in ("baseline", descriptor)]
            for level in levels:
                selected = [point for point in eligible_points if level_for(point, descriptor) == level]
                if not selected:
                    continue
                ids = {point.design_id for point in selected}
                crossovers = [crossover_density(rows, design_id) for design_id in ids]
                finite_crossovers = [value for value in crossovers if np.isfinite(value)]
                for target in TARGETS:
                    group = [row for row in rows if row["design_id"] in ids and row["target_density"] == target]
                    hr = [row["HR_pct"] for row in group if row["HR_scored"]]
                    ts = [row["TS_pct"] for row in group if row["TS_scored"]]
                    summary.append({
                        "analysis_design": design_type, "descriptor": descriptor, "level": level,
                        "target_density": target, "n_conditions": len(selected),
                        "fraction_HR_attainable": len(hr) / len(group) if group else np.nan,
                        "fraction_TS_attainable": len(ts) / len(group) if group else np.nan,
                        "fraction_HR_positive": np.mean(np.asarray(hr) > 0) if hr else np.nan,
                        "fraction_TS_positive": np.mean(np.asarray(ts) > 0) if ts else np.nan,
                        "mean_HR_pct": np.mean(hr) if hr else np.nan,
                        "mean_TS_pct": np.mean(ts) if ts else np.nan,
                        "median_HR_crossover_density": np.median(finite_crossovers) if finite_crossovers else np.nan,
                        "fraction_with_HR_crossover": len(finite_crossovers) / len(crossovers),
                    })
    return summary


def oat_plot_rows(rows: list[dict], descriptor: str) -> list[dict]:
    allowed = {"baseline", descriptor}
    return [row for row in rows if row["design_type"] == "one_at_a_time" and row["varied_descriptor"] in allowed]


def plot_oat_metric(outdir: Path, rows: list[dict], descriptor: str, metric: str, filename: str) -> None:
    subset = oat_plot_rows(rows, descriptor)
    fig, axis = plt.subplots(figsize=(8, 5.5))
    for design_id in sorted({row["design_id"] for row in subset}):
        group = [row for row in subset if row["design_id"] == design_id]
        point = group[0]
        label = f"{descriptor}={point[descriptor]:g}"
        axis.plot([row["target_density"] for row in group], [row[metric] for row in group], marker="o", label=label)
    axis.axhline(0, color="black", linestyle="--", linewidth=1)
    axis.set(xlabel="Target density", ylabel=metric)
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=150)
    plt.close(fig)


def make_plots(outdir: Path, rows: list[dict], points: list[DesignPoint]) -> None:
    for descriptor, filename in (
        ("rho0", "HR_pct_by_rho0.png"), ("G0_nm", "HR_pct_by_G0.png"),
        ("pore_scale_nm", "HR_pct_by_pore_scale.png"), ("log_width", "HR_pct_by_log_width.png"),
    ):
        plot_oat_metric(outdir, rows, descriptor, "HR_pct", filename)
    plot_oat_metric(outdir, rows, "rho0", "TS_pct", "TS_pct_by_rho0.png")

    subset = oat_plot_rows(rows, "rho0")
    fig, axis = plt.subplots(figsize=(8, 5.5))
    for design_id in sorted({row["design_id"] for row in subset}):
        group = [row for row in subset if row["design_id"] == design_id]
        axis.plot([row["target_density"] for row in group], [row["slow_cumulative_redistributed_pore_volume"] for row in group], marker="o", label=f"rho0={group[0]['rho0']:g}")
    axis.set(xlabel="Target density", ylabel="Slow-ramp cumulative redistribution")
    axis.grid(alpha=0.25); axis.legend(fontsize=8); fig.tight_layout()
    fig.savefig(outdir / "cumulative_redistribution_vs_density.png", dpi=150); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for axis, descriptor in zip(axes.flat, ("rho0", "G0_nm", "pore_scale_nm", "log_width")):
        selected = [point for point in oat_design() if point.varied_descriptor in ("baseline", descriptor)]
        x = [level_for(point, descriptor) for point in selected]
        y = [crossover_density(rows, point.design_id) for point in selected]
        axis.scatter([a for a, b in zip(x, y) if np.isfinite(b)], [b for b in y if np.isfinite(b)])
        axis.scatter([a for a, b in zip(x, y) if not np.isfinite(b)], [0.947 for b in y if not np.isfinite(b)], marker="x", color="tab:red", label="no crossover by 0.92")
        axis.set(xlabel=descriptor, ylabel="First positive-HR density", title=f"Crossover vs {descriptor}")
        axis.set_ylim(0.85, 0.952)
        axis.grid(alpha=0.25)
        if any(not np.isfinite(value) for value in y): axis.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(outdir / "fast_heating_crossover_by_descriptor.png", dpi=150); plt.close(fig)

    matrix = np.full((len(points), len(TARGETS)), np.nan)
    for i, point in enumerate(points):
        for j, target in enumerate(TARGETS):
            row = next(row for row in rows if row["design_id"] == point.design_id and row["target_density"] == target)
            if row["eligible_target"]:
                matrix[i, j] = np.mean([row[f"{name}_reached"] for name in ("slow", "fast", "high", "two_step")])
    fig, axis = plt.subplots(figsize=(9, 10))
    image = axis.imshow(matrix, aspect="auto", cmap="magma", interpolation="nearest")
    axis.set_xticks(range(len(TARGETS)), [f"{target:.2f}" for target in TARGETS])
    axis.set_yticks(range(len(points)), [point.design_id for point in points], fontsize=7)
    axis.set(xlabel="Target density", title="Canonical protocol attainability fraction")
    fig.colorbar(image, ax=axis); fig.tight_layout()
    fig.savefig(outdir / "attainability_map.png", dpi=150); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/initial_condition_sensitivity")
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    base = model.Params(memory_model="pore_bin_redistribution")
    oat, factorial = oat_design(), factorial_design()
    points = oat + factorial
    assert_shared_mechanism_parameters(points, base)
    rows = run_design(points, base)
    oat_rows = [row for row in rows if row["design_type"] == "one_at_a_time"]
    factorial_rows = [row for row in rows if row["design_type"] == "corner_factorial"]
    unattainable = [row for row in rows if row["eligible_target"] and (not row["HR_scored"] or not row["TS_scored"])]
    summary = sensitivity_summary(rows, points)
    write_csv(outdir / "factorial_results.csv", factorial_rows)
    write_csv(outdir / "one_at_a_time_results.csv", oat_rows)
    write_csv(outdir / "sensitivity_summary.csv", summary)
    write_csv(outdir / "unattainable_cases.csv", unattainable)
    make_plots(outdir, rows, points)
    print(f"conditions={len(points)} rows={len(rows)} unattainable={len(unattainable)}")


if __name__ == "__main__":
    main()
