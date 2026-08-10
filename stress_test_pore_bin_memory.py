#!/usr/bin/env python3
"""Compare no, empirical, and observable pore-bin memory models."""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from stress_test_topology_memory import (
    at_target,
    heating_rate_audit,
    target_density_audit,
    two_step_audit,
    write_csv,
)

MODES = ("none", "empirical_topology_damage", "pore_bin_redistribution")


def reference_runs(params: model.Params, target: float) -> dict[str, dict]:
    return {
        "slow": model.run(params, model.RampHold(0.2), target),
        "fast": model.run(params, model.RampHold(20.0), target),
        "high": model.run(replace(params, G0=75e-9), model.Iso(1350), target),
        "two_step": model.run(replace(params, G0=75e-9), model.TwoStep(), target),
    }


def mode_summary(target: float) -> tuple[list[dict], dict[str, dict[str, dict]]]:
    rows, all_runs = [], {}
    for mode in MODES:
        runs = reference_runs(model.Params(memory_model=mode), target)
        all_runs[mode] = runs
        values = {name: at_target(result, target) for name, result in runs.items()}
        rows.append({
            "memory_model": mode,
            "all_reached": all(flag for _, flag in values.values()),
            "HR_pct": model.percent_gain(values["slow"][0], values["fast"][0]) if values["slow"][1] and values["fast"][1] else np.nan,
            "TS_pct": model.percent_gain(values["high"][0], values["two_step"][0]) if values["high"][1] and values["two_step"][1] else np.nan,
            "slow_mean_radius_nm": runs["slow"]["pore_mean_radius"][-1] * 1e9,
            "fast_mean_radius_nm": runs["fast"]["pore_mean_radius"][-1] * 1e9,
            "slow_large_pore_fraction": runs["slow"]["large_pore_fraction"][-1],
            "fast_large_pore_fraction": runs["fast"]["large_pore_fraction"][-1],
            "slow_removable_fine_fraction": runs["slow"]["removable_fine_pore_fraction"][-1],
            "fast_removable_fine_fraction": runs["fast"]["removable_fine_pore_fraction"][-1],
            "slow_cumulative_redistributed_volume": runs["slow"]["cumulative_redistributed_pore_volume"][-1],
            "fast_cumulative_redistributed_volume": runs["fast"]["cumulative_redistributed_pore_volume"][-1],
        })
    return rows, all_runs


def held_out_tables(target: float) -> tuple[list[dict], list[dict], list[dict]]:
    heating, two_step, targets = [], [], []
    for mode in MODES:
        params = model.Params(memory_model=mode)
        heat_rows, _ = heating_rate_audit(params, target)
        grid_rows, _ = two_step_audit(params, target)
        target_rows = target_density_audit(params)
        heating.extend({"memory_model": mode, **row} for row in heat_rows)
        two_step.extend({"memory_model": mode, **row} for row in grid_rows)
        targets.extend({"memory_model": mode, **row} for row in target_rows)
    return heating, two_step, targets


def plot_mode_comparison(outdir: Path, all_runs: dict[str, dict[str, dict]]) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    colors = {"slow": "tab:blue", "fast": "tab:orange"}
    for mode in MODES:
        for schedule in ("slow", "fast"):
            result = all_runs[mode][schedule]
            label = f"{mode}: {schedule}"
            axes[0, 0].plot(result["rho"], result["pore_mean_radius"] * 1e9, label=label, color=colors[schedule], linestyle={"none":":","empirical_topology_damage":"--","pore_bin_redistribution":"-"}[mode])
            axes[0, 1].plot(result["rho"], result["large_pore_fraction"], label=label, color=colors[schedule], linestyle={"none":":","empirical_topology_damage":"--","pore_bin_redistribution":"-"}[mode])
            axes[0, 2].plot(result["rho"], result["removable_fine_pore_fraction"], label=label, color=colors[schedule], linestyle={"none":":","empirical_topology_damage":"--","pore_bin_redistribution":"-"}[mode])
            axes[1, 0].plot(result["T_C"], result["cumulative_redistributed_pore_volume"], label=label, color=colors[schedule], linestyle={"none":":","empirical_topology_damage":"--","pore_bin_redistribution":"-"}[mode])
            axes[1, 1].plot(result["rho"], result["f_pore"], label=label, color=colors[schedule], linestyle={"none":":","empirical_topology_damage":"--","pore_bin_redistribution":"-"}[mode])
            axes[1, 2].plot(result["G"] * 1e9, result["rho"], label=label, color=colors[schedule], linestyle={"none":":","empirical_topology_damage":"--","pore_bin_redistribution":"-"}[mode])
    settings = [
        ("Density", "Mean pore radius [nm]", "Observable mean-radius memory"),
        ("Density", "Large-pore fraction", "Large-pore memory"),
        ("Density", "Removable fine-pore fraction", "Fine-pore depletion"),
        ("Temperature [C]", "Cumulative redistributed volume", "Redistribution vs temperature"),
        ("Density", "Pore-boundary coverage", "Coverage response"),
        ("Grain size [nm]", "Density", "Matched-density trajectory"),
    ]
    for axis, (xlabel, ylabel, title) in zip(axes.flat, settings):
        axis.set(xlabel=xlabel, ylabel=ylabel, title=title)
        axis.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(outdir / "memory_mode_comparison.png", dpi=150)
    plt.close(fig)

    radii = model.initial_state(model.Params()).pore_radii * 1e9
    positions = np.arange(len(radii))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for axis, schedule in zip(axes, ("slow", "fast")):
        result = all_runs["pore_bin_redistribution"][schedule]
        initial = result["pore_phi"][0] / result["pore_phi"][0].sum()
        final = result["pore_phi"][-1] / result["pore_phi"][-1].sum()
        axis.bar(positions - 0.2, initial, width=0.4, label="initial")
        axis.bar(positions + 0.2, final, width=0.4, label="at rho=0.90")
        axis.set_xticks(positions, [f"{radius:.0f}" for radius in radii], rotation=35)
        axis.set(xlabel="Pore radius [nm]", title=f"{schedule} heating")
        axis.grid(alpha=0.25)
    axes[0].set_ylabel("Pore-volume fraction of distribution")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(outdir / "pore_bin_distributions_slow_fast.png", dpi=150)
    plt.close(fig)


def plot_redistribution_flux(outdir: Path, runs: dict[str, dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for schedule, result in runs.items():
        if schedule not in ("slow", "fast"):
            continue
        axes[0].plot(result["T_C"], np.sum(result["redistribution_flux_by_bin"], axis=1), label=schedule)
        axes[1].plot(result["rho"], result["E_G"], label=schedule)
    axes[0].set(xlabel="Temperature [C]", ylabel="Redistribution flux [pore fraction/s]", title="Conservative adjacent-bin flux")
    axes[1].set(xlabel="Density", ylabel="E_G", title="Trajectory efficiency")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    fig.tight_layout()
    fig.savefig(outdir / "redistribution_flux_and_efficiency.png", dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/pore_bin_memory_stress")
    parser.add_argument("--target", type=float, default=0.90)
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summaries, runs = mode_summary(args.target)
    heating, two_step, targets = held_out_tables(args.target)
    write_csv(outdir / "memory_mode_summary.csv", summaries)
    write_csv(outdir / "held_out_heating_by_mode.csv", heating)
    write_csv(outdir / "held_out_two_step_by_mode.csv", two_step)
    write_csv(outdir / "target_density_by_mode.csv", targets)
    plot_mode_comparison(outdir, runs)
    plot_redistribution_flux(outdir, runs["pore_bin_redistribution"])
    print(*summaries, sep="\n")


if __name__ == "__main__":
    main()
