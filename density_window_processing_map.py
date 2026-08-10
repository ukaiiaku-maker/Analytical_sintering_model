#!/usr/bin/env python3
"""Ensemble processing map over initial microstructure and density window."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model

TARGETS = (0.80, 0.85, 0.88, 0.90, 0.92, 0.95, 0.98)
HEATING_RATES = (0.2, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0)
T1_VALUES = (1300.0, 1350.0, 1400.0)
T2_VALUES = (1150.0, 1200.0, 1250.0, 1300.0)
SWITCH_DENSITIES = (0.78, 0.82, 0.85, 0.88, 0.90)
WINDOWS = (
    ("early", 0.70, 0.80),
    ("intermediate", 0.80, 0.90),
    ("late", 0.90, 0.95),
    ("near_final", 0.95, 0.99),
)


@dataclass(frozen=True)
class InitialClass:
    name: str
    rho0: float
    G0_nm: float
    pore_radius0_nm: float
    pore_ln_sigma: float

    def params(self) -> model.Params:
        return model.Params(
            memory_model="pore_bin_redistribution",
            rho0=self.rho0,
            G0=self.G0_nm * 1e-9,
            pore_radius0=self.pore_radius0_nm * 1e-9,
            pore_ln_sigma=self.pore_ln_sigma,
        )


INITIAL_CLASSES = (
    InitialClass("loose_fine", 0.70, 100.0, 16.0, 0.85),
    InitialClass("baseline_intermediate", 0.75, 150.0, 22.0, 0.65),
    InitialClass("predensified_partially_isolated", 0.82, 250.0, 35.0, 0.45),
)


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sample_at_density(result: dict, target: float, rho0: float) -> dict:
    """Sample only attainable densities; never extrapolate a failed target."""
    maximum = float(np.max(result["rho"]))
    eligible = target > rho0 + 1e-12
    reached = maximum >= target - 1e-12
    if not reached:
        return {
            "eligible_target": eligible, "reached_target": False,
            "time_to_target_h": np.nan, "G_at_target_nm": np.nan,
            "pore_mean_radius_nm": np.nan, "large_pore_fraction": np.nan,
            "removable_fine_pore_fraction": np.nan, "f_pore": np.nan,
            "isolated_pore_fraction": np.nan, "cumulative_redistributed_pore_volume": np.nan,
            "E_G": np.nan,
        }
    rho = np.asarray(result["rho"], dtype=float)
    sample = lambda key: float(np.interp(target, rho, np.asarray(result[key], dtype=float)))
    return {
        "eligible_target": eligible, "reached_target": True,
        "time_to_target_h": sample("t") / 3600.0,
        "G_at_target_nm": sample("G") * 1e9,
        "pore_mean_radius_nm": sample("pore_mean_radius") * 1e9,
        "large_pore_fraction": sample("large_pore_fraction"),
        "removable_fine_pore_fraction": sample("removable_fine_pore_fraction"),
        "f_pore": sample("f_pore"),
        "isolated_pore_fraction": sample("isolated_pore_fraction"),
        "cumulative_redistributed_pore_volume": sample("cumulative_redistributed_pore_volume"),
        "E_G": sample("E_G"),
    }


def protocol_row(initial: InitialClass, kind: str, protocol_id: str, result: dict, target: float, budget_h: float, **metadata) -> dict:
    sampled = sample_at_density(result, target, initial.rho0)
    return {
        "initial_class": initial.name,
        "rho0": initial.rho0,
        "G0_nm": initial.G0_nm,
        "pore_radius0_nm": initial.pore_radius0_nm,
        "pore_ln_sigma": initial.pore_ln_sigma,
        "protocol_kind": kind,
        "protocol_id": protocol_id,
        "target_density": target,
        "time_budget_h": budget_h,
        "final_density": float(result["rho"][-1]),
        **sampled,
        "heating_rate_C_min": metadata.get("heating_rate_C_min", np.nan),
        "T1_C": metadata.get("T1_C", np.nan),
        "T2_C": metadata.get("T2_C", np.nan),
        "switch_density": metadata.get("switch_density", np.nan),
    }


def run_ensemble() -> tuple[list[dict], dict]:
    rows, cache = [], {}
    total_runs = len(INITIAL_CLASSES) * (len(HEATING_RATES) + len(T1_VALUES) + len(T1_VALUES) * len(T2_VALUES) * len(SWITCH_DENSITIES))
    completed = 0
    for initial in INITIAL_CLASSES:
        params = initial.params()
        for rate in HEATING_RATES:
            protocol = model.RampHold(rate)
            result = model.run(params, protocol)
            cache[(initial.name, "heating", rate)] = result
            budget_h = min(protocol.t_end, params.t_max_s) / 3600.0
            for target in TARGETS:
                rows.append(protocol_row(initial, "heating", f"ramp_{rate:g}", result, target, budget_h, heating_rate_C_min=rate))
            completed += 1
            print(f"[{completed}/{total_runs}] {initial.name} heating {rate:g} C/min", flush=True)
        for T1 in T1_VALUES:
            protocol = model.Iso(T1)
            result = model.run(replace(params, G0=params.G0), protocol)
            cache[(initial.name, "high", T1)] = result
            budget_h = min(protocol.t_end, params.t_max_s) / 3600.0
            for target in TARGETS:
                rows.append(protocol_row(initial, "high_T", f"iso_{T1:g}", result, target, budget_h, T1_C=T1))
            completed += 1
            print(f"[{completed}/{total_runs}] {initial.name} high {T1:g} C", flush=True)
        for T1 in T1_VALUES:
            for T2 in T2_VALUES:
                for switch in SWITCH_DENSITIES:
                    protocol = model.TwoStep(T1, T2, switch)
                    result = model.run(replace(params, G0=params.G0), protocol)
                    cache[(initial.name, "two_step", T1, T2, switch)] = result
                    budget_h = min(protocol.t_end, params.t_max_s) / 3600.0
                    for target in TARGETS:
                        rows.append(protocol_row(initial, "two_step", f"two_{T1:g}_{T2:g}_{switch:.2f}", result, target, budget_h, T1_C=T1, T2_C=T2, switch_density=switch))
                    completed += 1
                    print(f"[{completed}/{total_runs}] {initial.name} two {T1:g}/{T2:g}/{switch:.2f}", flush=True)
    return rows, cache


def heating_curves(cache: dict) -> list[dict]:
    rows = []
    for initial in INITIAL_CLASSES:
        reference = cache[(initial.name, "heating", 0.2)]
        for target in TARGETS:
            slow = sample_at_density(reference, target, initial.rho0)
            for rate in HEATING_RATES:
                current = sample_at_density(cache[(initial.name, "heating", rate)], target, initial.rho0)
                scored = slow["eligible_target"] and current["eligible_target"] and slow["reached_target"] and current["reached_target"]
                hr = model.percent_gain(slow["G_at_target_nm"], current["G_at_target_nm"]) if scored else np.nan
                rows.append({
                    "initial_class": initial.name, "rho0": initial.rho0, "target_density": target,
                    "heating_rate_C_min": rate, "slow_reached": slow["reached_target"], "current_reached": current["reached_target"],
                    "eligible_target": slow["eligible_target"] and current["eligible_target"],
                    "HR_pct_vs_0p2": hr, **{f"current_{key}": value for key, value in current.items() if key not in ("eligible_target", "reached_target")},
                })
    return rows


def two_step_grid(cache: dict) -> list[dict]:
    rows = []
    for initial in INITIAL_CLASSES:
        for T1 in T1_VALUES:
            high_result = cache[(initial.name, "high", T1)]
            for T2 in T2_VALUES:
                for switch in SWITCH_DENSITIES:
                    two_result = cache[(initial.name, "two_step", T1, T2, switch)]
                    for target in TARGETS:
                        high = sample_at_density(high_result, target, initial.rho0)
                        two = sample_at_density(two_result, target, initial.rho0)
                        scored = high["eligible_target"] and two["eligible_target"] and high["reached_target"] and two["reached_target"]
                        ts = model.percent_gain(high["G_at_target_nm"], two["G_at_target_nm"]) if scored else np.nan
                        rows.append({
                            "initial_class": initial.name, "rho0": initial.rho0, "T1_C": T1, "T2_C": T2,
                            "switch_density": switch, "target_density": target, "eligible_target": high["eligible_target"] and two["eligible_target"],
                            "high_step_executed": switch > initial.rho0 + 1e-12,
                            "target_after_switch": target > switch + 1e-12,
                            "high_reached": high["reached_target"], "two_step_reached": two["reached_target"], "TS_pct": ts,
                            "high_time_h": high["time_to_target_h"], "two_step_time_h": two["time_to_target_h"],
                            "G_high_nm": high["G_at_target_nm"], "G_two_step_nm": two["G_at_target_nm"],
                            "high_redistributed_volume": high["cumulative_redistributed_pore_volume"],
                            "two_step_redistributed_volume": two["cumulative_redistributed_pore_volume"],
                        })
    return rows


def canonical_density_map(cache: dict) -> list[dict]:
    """Canonical 0.2/20 and 1350/1250/switch-0.85 comparisons."""
    rows = []
    for initial in INITIAL_CLASSES:
        protocols = {
            "slow": cache[(initial.name, "heating", 0.2)],
            "fast": cache[(initial.name, "heating", 20.0)],
            "high": cache[(initial.name, "high", 1350.0)],
            "two_step": cache[(initial.name, "two_step", 1350.0, 1250.0, 0.85)],
        }
        for target in TARGETS:
            values = {name: sample_at_density(result, target, initial.rho0) for name, result in protocols.items()}
            hr_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("slow", "fast"))
            ts_scored = all(values[name]["eligible_target"] and values[name]["reached_target"] for name in ("high", "two_step"))
            hr = model.percent_gain(values["slow"]["G_at_target_nm"], values["fast"]["G_at_target_nm"]) if hr_scored else np.nan
            ts = model.percent_gain(values["high"]["G_at_target_nm"], values["two_step"]["G_at_target_nm"]) if ts_scored else np.nan
            if not values["slow"]["eligible_target"]:
                classification = "neutral"
            elif not hr_scored or not ts_scored:
                classification = "requires late-stage physics" if target >= 0.95 else "unattainable"
            elif hr > 0.5 and ts > 0.5:
                classification = "beneficial"
            elif abs(hr) <= 0.5 and abs(ts) <= 0.5:
                classification = "neutral"
            else:
                classification = "harmful"
            rows.append({
                "initial_class": initial.name, "rho0": initial.rho0, "target_density": target,
                "HR_pct": hr, "TS_pct": ts, "HR_scored": hr_scored, "TS_scored": ts_scored,
                "slow_reached": values["slow"]["reached_target"], "fast_reached": values["fast"]["reached_target"],
                "high_reached": values["high"]["reached_target"], "two_step_reached": values["two_step"]["reached_target"],
                "classification": classification,
                **{f"{name}_{metric}": values[name][metric] for name in values for metric in ("time_to_target_h", "G_at_target_nm", "pore_mean_radius_nm", "large_pore_fraction", "removable_fine_pore_fraction", "f_pore", "isolated_pore_fraction", "cumulative_redistributed_pore_volume", "E_G")},
            })
    return rows


def initial_summary() -> list[dict]:
    rows = []
    for initial in INITIAL_CLASSES:
        state = model.initial_state(initial.params())
        diag = model.pore_distribution_diagnostics(state.pore_phi, state.pore_radii, initial.params())
        rows.append({
            "initial_class": initial.name, "rho0": initial.rho0, "G0_nm": initial.G0_nm,
            "pore_radius0_nm": initial.pore_radius0_nm, "pore_ln_sigma": initial.pore_ln_sigma,
            "f_pore_initial": state.topology.f_pore, "connectivity_initial": state.topology.connectivity,
            "isolated_fraction_initial": state.topology.isolated_pore_fraction,
            "pore_mean_radius_initial_nm": diag["pore_mean_radius"] * 1e9,
            "large_pore_fraction_initial": diag["large_pore_fraction"],
            "removable_fine_fraction_initial": diag["removable_fine_pore_fraction"],
        })
    return rows


def window_for_density(target: float) -> str:
    for name, low, high in WINDOWS:
        if low <= target <= high:
            return name
    raise ValueError(target)


def window_summary(canonical: list[dict]) -> list[dict]:
    rows = []
    for initial in INITIAL_CLASSES:
        subset = [row for row in canonical if row["initial_class"] == initial.name and row["target_density"] > initial.rho0]
        for window, _, _ in WINDOWS:
            group = [row for row in subset if window_for_density(row["target_density"]) == window]
            hr = [row["HR_pct"] for row in group if row["HR_scored"]]
            ts = [row["TS_pct"] for row in group if row["TS_scored"]]
            attainable = [row for row in group if row["HR_scored"] and row["TS_scored"]]
            times = [row[f"{name}_time_to_target_h"] for row in attainable for name in ("slow", "fast", "high", "two_step")]
            grains = [row[f"{name}_G_at_target_nm"] for row in attainable for name in ("slow", "fast", "high", "two_step")]
            rows.append({
                "initial_class": initial.name, "density_window": window, "eligible_points": len(group),
                "fraction_attainable": len(attainable) / len(group) if group else np.nan,
                "fraction_positive_HR": np.mean(np.asarray(hr) > 0) if hr else np.nan,
                "fraction_positive_TS": np.mean(np.asarray(ts) > 0) if ts else np.nan,
                "mean_HR_pct": np.mean(hr) if hr else np.nan, "mean_TS_pct": np.mean(ts) if ts else np.nan,
                "median_time_to_target_h": np.median(times) if times else np.nan,
                "median_grain_size_nm": np.median(grains) if grains else np.nan,
            })
    return rows


def plot_metric_curves(outdir: Path, canonical: list[dict], metric: str, filename: str, ylabel: str) -> None:
    fig, axis = plt.subplots(figsize=(8, 5.5))
    for initial in INITIAL_CLASSES:
        rows = [row for row in canonical if row["initial_class"] == initial.name]
        axis.plot([row["target_density"] for row in rows], [row[metric] for row in rows], marker="o", label=initial.name)
    if metric in ("HR_pct", "TS_pct"):
        axis.axhline(0, color="black", linestyle="--", linewidth=1)
    axis.set(xlabel="Target density", ylabel=ylabel)
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=150)
    plt.close(fig)


def plot_heatmap(outdir: Path, canonical: list[dict], value, filename: str, title: str, cmap="viridis") -> None:
    matrix = np.full((len(INITIAL_CLASSES), len(TARGETS)), np.nan)
    for i, initial in enumerate(INITIAL_CLASSES):
        for j, target in enumerate(TARGETS):
            row = next(row for row in canonical if row["initial_class"] == initial.name and row["target_density"] == target)
            matrix[i, j] = value(row)
    fig, axis = plt.subplots(figsize=(9, 4.2))
    image = axis.imshow(matrix, aspect="auto", cmap=cmap, interpolation="nearest")
    axis.set_xticks(range(len(TARGETS)), [f"{target:.2f}" for target in TARGETS])
    axis.set_yticks(range(len(INITIAL_CLASSES)), [item.name for item in INITIAL_CLASSES])
    axis.set(xlabel="Target density", title=title)
    fig.colorbar(image, ax=axis)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=150)
    plt.close(fig)


def finite_median(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=float)
    return float(np.median(array[np.isfinite(array)])) if np.any(np.isfinite(array)) else np.nan


def make_plots(outdir: Path, canonical: list[dict]) -> None:
    plot_metric_curves(outdir, canonical, "HR_pct", "HR_pct_vs_density.png", "HR_pct (0.2 vs 20 C/min)")
    plot_metric_curves(outdir, canonical, "TS_pct", "TS_pct_vs_density.png", "TS_pct (1350 vs 1350/1250)")
    attainable = lambda row: np.mean([row[f"{name}_reached"] for name in ("slow", "fast", "high", "two_step")])
    plot_heatmap(outdir, canonical, attainable, "attainable_density_map.png", "Fraction of canonical protocols reaching target", "magma")
    plot_heatmap(outdir, canonical, lambda row: finite_median(row[f"{name}_time_to_target_h"] for name in ("slow", "fast", "high", "two_step")), "time_to_target_map.png", "Median time to target [h]")
    for metric, filename, ylabel in (
        ("G_at_target_nm", "grain_size_vs_density.png", "Grain size [nm]"),
        ("pore_mean_radius_nm", "pore_mean_radius_vs_density.png", "Pore mean radius [nm]"),
        ("large_pore_fraction", "large_pore_fraction_vs_density.png", "Large-pore fraction"),
        ("removable_fine_pore_fraction", "removable_fine_fraction_vs_density.png", "Removable fine-pore fraction"),
        ("cumulative_redistributed_pore_volume", "cumulative_redistribution_vs_density.png", "Cumulative redistributed pore volume"),
    ):
        fig, axis = plt.subplots(figsize=(9, 6))
        for initial in INITIAL_CLASSES:
            rows = [row for row in canonical if row["initial_class"] == initial.name]
            for schedule, linestyle in (("slow", "-"), ("fast", "--"), ("high", ":"), ("two_step", "-.")):
                axis.plot([row["target_density"] for row in rows], [row[f"{schedule}_{metric}"] for row in rows], marker="o", linestyle=linestyle, label=f"{initial.name}: {schedule}")
        axis.set(xlabel="Target density", ylabel=ylabel)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(outdir / filename, dpi=150)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/density_window_processing_map")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_rows, cache = run_ensemble()
    heating = heating_curves(cache)
    two_step = two_step_grid(cache)
    canonical = canonical_density_map(cache)
    initial = initial_summary()
    windows = window_summary(canonical)
    unattainable = [row for row in all_rows if row["eligible_target"] and not row["reached_target"]]

    write_csv(outdir / "all_protocol_targets.csv", all_rows)
    write_csv(outdir / "heating_rate_density_curves.csv", heating)
    write_csv(outdir / "two_step_density_grid.csv", two_step)
    write_csv(outdir / "initial_microstructure_summary.csv", initial)
    write_csv(outdir / "window_summary.csv", windows)
    write_csv(outdir / "unattainable_cases.csv", unattainable)
    write_csv(outdir / "canonical_density_map.csv", canonical)
    make_plots(outdir, canonical)
    print(f"wrote {len(all_rows)} protocol-target rows; unattainable={len(unattainable)}")


if __name__ == "__main__":
    main()
