#!/usr/bin/env python3
"""Deterministic held-out schedule and ablation audit for topology memory."""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model

HEATING_RATES = (0.2, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0)
TARGETS = (0.88, 0.90, 0.92)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def at_target(result: dict, target: float) -> tuple[float, bool]:
    grain, reached = model.value_at_density(result, target)
    return grain * 1e9, reached


def finite_median(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.median(values[np.isfinite(values)])) if np.any(np.isfinite(values)) else np.nan


def heating_rate_audit(params: model.Params, target: float) -> tuple[list[dict], dict[float, dict]]:
    rows, runs = [], {}
    for rate in HEATING_RATES:
        result = model.run(params, model.RampHold(rate), target)
        runs[rate] = result
        grain_nm, reached = at_target(result, target)
        rows.append({
            "heating_rate_C_min": rate,
            "target_density": target,
            "reached_target": reached,
            "final_density": result["rho"][-1],
            "G_at_target_nm": grain_nm,
            "topology_damage_final": result["topology_damage"][-1],
            "damage_rate_max_s": np.max(result["topology_damage_rate"]),
            "median_E_G": finite_median(result["E_G"]),
            "HR_pct_vs_0p2": np.nan,
        })
    reference = rows[0]["G_at_target_nm"]
    for row in rows:
        row["HR_pct_vs_0p2"] = model.percent_gain(reference, row["G_at_target_nm"])
    return rows, runs


def two_step_audit(params: model.Params, target: float) -> tuple[list[dict], dict]:
    rows, runs = [], {}
    high_runs = {}
    for T1 in (1300.0, 1350.0, 1400.0):
        high_runs[T1] = model.run(replace(params, G0=75e-9), model.Iso(T1), target)
        high_G, high_reached = at_target(high_runs[T1], target)
        for T2 in (1200.0, 1250.0, 1300.0):
            for switch in (0.82, 0.85, 0.88):
                two = model.run(replace(params, G0=75e-9), model.TwoStep(T1, T2, switch), target)
                runs[(T1, T2, switch)] = two
                two_G, two_reached = at_target(two, target)
                success = high_reached and two_reached
                rows.append({
                    "T1_C": T1, "T2_C": T2, "switch_density": switch, "target_density": target,
                    "high_reached": high_reached, "two_step_reached": two_reached,
                    "high_final_density": high_runs[T1]["rho"][-1], "two_step_final_density": two["rho"][-1],
                    "G_high_nm": high_G, "G_two_step_nm": two_G,
                    "TS_pct": model.percent_gain(high_G, two_G) if success else np.nan,
                    "benefit_positive": bool(success and model.percent_gain(high_G, two_G) > 0),
                    "high_damage_final": high_runs[T1]["topology_damage"][-1],
                    "two_step_damage_final": two["topology_damage"][-1],
                })
    return rows, {"high": high_runs, "two_step": runs}


def target_density_audit(params: model.Params) -> list[dict]:
    rows = []
    for target in TARGETS:
        cases = {
            "slow": model.run(params, model.RampHold(0.2), target),
            "fast": model.run(params, model.RampHold(20.0), target),
            "high": model.run(replace(params, G0=75e-9), model.Iso(1350), target),
            "two_step": model.run(replace(params, G0=75e-9), model.TwoStep(), target),
        }
        values = {name: at_target(result, target) for name, result in cases.items()}
        rows.append({
            "target_density": target,
            **{f"{name}_reached": reached for name, (_, reached) in values.items()},
            **{f"{name}_final_density": cases[name]["rho"][-1] for name in cases},
            **{f"G_{name}_nm": grain for name, (grain, _) in values.items()},
            "HR_pct": model.percent_gain(values["slow"][0], values["fast"][0]) if values["slow"][1] and values["fast"][1] else np.nan,
            "TS_pct": model.percent_gain(values["high"][0], values["two_step"][0]) if values["high"][1] and values["two_step"][1] else np.nan,
        })
    return rows


def ablation_parameters(base: model.Params) -> Iterable[tuple[str, model.Params]]:
    yield "memory_disabled", replace(base, enable_topology_memory=False)
    yield "memory_enabled", base
    yield "coverage_only", replace(base, damage_isolation_strength=0.0)
    yield "isolation_only", replace(base, damage_coverage_strength=0.0)
    yield "window_lower", replace(base, surface_damage_T_mid_C=850.0)
    yield "window_higher", replace(base, surface_damage_T_mid_C=1200.0)
    yield "damage_rate_half", replace(base, surface_damage_rate_s=0.5 * base.surface_damage_rate_s)
    yield "damage_rate_double", replace(base, surface_damage_rate_s=2.0 * base.surface_damage_rate_s)


def ablation_audit(base: model.Params, target: float) -> list[dict]:
    rows = []
    for name, params in ablation_parameters(base):
        runs = {
            "slow": model.run(params, model.RampHold(0.2), target),
            "fast": model.run(params, model.RampHold(20.0), target),
            "high": model.run(replace(params, G0=75e-9), model.Iso(1350), target),
            "two_step": model.run(replace(params, G0=75e-9), model.TwoStep(), target),
        }
        values = {key: at_target(value, target) for key, value in runs.items()}
        reached = all(flag for _, flag in values.values())
        rows.append({
            "ablation": name, "all_reached": reached,
            "slow_reached": values["slow"][1], "fast_reached": values["fast"][1],
            "high_reached": values["high"][1], "two_step_reached": values["two_step"][1],
            "HR_pct": model.percent_gain(values["slow"][0], values["fast"][0]) if values["slow"][1] and values["fast"][1] else np.nan,
            "TS_pct": model.percent_gain(values["high"][0], values["two_step"][0]) if values["high"][1] and values["two_step"][1] else np.nan,
            "slow_damage": runs["slow"]["topology_damage"][-1], "fast_damage": runs["fast"]["topology_damage"][-1],
            "high_damage": runs["high"]["topology_damage"][-1], "two_step_damage": runs["two_step"]["topology_damage"][-1],
        })
    return rows


def plot_heating_diagnostics(outdir: Path, runs: dict[float, dict]) -> None:
    fig, axes = plt.subplots(4, 2, figsize=(13, 17))
    for rate, result in runs.items():
        label = f"{rate:g} C/min"
        axes[0, 0].plot(result["T_C"], result["topology_damage"], label=label)
        axes[0, 1].plot(result["rho"], result["topology_damage"], label=label)
        axes[1, 0].plot(result["T_C"], result["topology_damage_rate"], label=label)
        axes[1, 1].plot(result["rho"], result["f_pore"], label=label)
        axes[2, 0].plot(result["rho"], result["isolated_pore_fraction"], label=label)
        axes[2, 1].plot(result["rho"], result["E_G"], label=label)
        axes[3, 0].plot(result["rho"], result["pore_mean_radius"] * 1e9, label=label)
        axes[3, 1].plot(result["rho"], result["large_pore_fraction"], label=label)
    labels = [
        ("Temperature [C]", "Topology damage", "Damage vs temperature"),
        ("Density", "Topology damage", "Damage vs density"),
        ("Temperature [C]", "Damage rate [1/s]", "Damage rate vs temperature"),
        ("Density", "f_pore", "Removable pore coverage"),
        ("Density", "Isolated pore fraction", "Pore isolation"),
        ("Density", "E_G", "Trajectory efficiency"),
        ("Density", "Mean pore radius [nm]", "Pore-bin mean radius"),
        ("Density", "Large-pore fraction", "Large-pore fraction"),
    ]
    for axis, (xlabel, ylabel, title) in zip(axes.flat, labels):
        axis.set(xlabel=xlabel, ylabel=ylabel, title=title)
        axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "held_out_heating_diagnostics.png", dpi=150)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8, 6))
    for rate, result in runs.items():
        axis.plot(result["G"] * 1e9, result["rho"], label=f"{rate:g} C/min")
    axis.set(xlabel="Grain size [nm]", ylabel="Density", title="Held-out heating-rate trajectories")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "held_out_heating_grain_density.png", dpi=150)
    plt.close(fig)


def plot_two_step_window(outdir: Path, rows: list[dict], runs: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
    for axis, T1 in zip(axes, (1300.0, 1350.0, 1400.0)):
        subset = [row for row in rows if row["T1_C"] == T1]
        for switch in (0.82, 0.85, 0.88):
            group = [row for row in subset if row["switch_density"] == switch]
            axis.plot([row["T2_C"] for row in group], [row["TS_pct"] for row in group], marker="o", label=f"switch {switch:.2f}")
        axis.axhline(0, color="black", linewidth=1, linestyle="--")
        axis.set(xlabel="T2 [C]", title=f"T1 = {T1:g} C")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("TS_pct")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "held_out_two_step_window.png", dpi=150)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(9, 6))
    colors = {1300.0: "tab:blue", 1350.0: "tab:orange", 1400.0: "tab:green"}
    for T1, result in runs["high"].items():
        axis.plot(result["G"] * 1e9, result["rho"], color=colors[T1], linewidth=2, linestyle="--", label=f"{T1:g} C high")
    for (T1, T2, switch), result in runs["two_step"].items():
        axis.plot(result["G"] * 1e9, result["rho"], color=colors[T1], alpha=0.28, linewidth=1)
    axis.set(xlabel="Grain size [nm]", ylabel="Density", title="All held-out high-T and two-step trajectories")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "held_out_two_step_grain_density.png", dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/topology_memory_stress")
    parser.add_argument("--target", type=float, default=0.90)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    params = model.Params()

    heating_rows, heating_runs = heating_rate_audit(params, args.target)
    two_step_rows, two_step_runs = two_step_audit(params, args.target)
    target_rows = target_density_audit(params)
    ablation_rows = ablation_audit(params, args.target)
    write_csv(outdir / "held_out_heating_rates.csv", heating_rows)
    write_csv(outdir / "held_out_two_step_grid.csv", two_step_rows)
    write_csv(outdir / "target_density_sweep.csv", target_rows)
    write_csv(outdir / "topology_memory_ablations.csv", ablation_rows)
    plot_heating_diagnostics(outdir, heating_runs)
    plot_two_step_window(outdir, two_step_rows, two_step_runs)
    positive = sum(row["benefit_positive"] for row in two_step_rows)
    reached = sum(row["high_reached"] and row["two_step_reached"] for row in two_step_rows)
    print(f"held-out two-step positive/reached/total: {positive}/{reached}/{len(two_step_rows)}")


if __name__ == "__main__":
    main()
