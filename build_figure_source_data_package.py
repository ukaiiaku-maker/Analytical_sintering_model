#!/usr/bin/env python3
"""Build a figure-ready data package from existing exact/diagnostic outputs.

This script performs no simulations and imports no model implementation.  It only
reads frozen CSV/JSON/ZIP evidence, reshapes it, adds provenance, and writes a
compressed source-data package for external figure construction.
"""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
import zipfile
from collections import defaultdict
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "figure_source_data_package"
ARCHIVE = ROOT / "results" / "1_Backup_of_prior_runs.zip"
ARCHIVE_FAST_MEMBER = (
    "visual_inspection_candidate_plots_v2/histories/dense_fast_histories.csv"
)
COMMIT = subprocess.run(
    ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
).stdout.strip()
BRANCH = subprocess.run(
    ["git", "branch", "--show-current"],
    cwd=ROOT,
    check=True,
    capture_output=True,
    text=True,
).stdout.strip()

META_COLUMNS = [
    "figure_id",
    "panel_id",
    "candidate_id",
    "evidence_level",
    "model_layer",
    "source_file",
    "source_table",
    "source_branch_or_commit",
    "units_reference",
    "human_readable_label_reference",
    "derivation_notes",
]

DATASETS: list[dict[str, object]] = []
MISSING: list[dict[str, str]] = []


def read_csv(path: str | Path, **kwargs: object) -> pd.DataFrame:
    return pd.read_csv(ROOT / path, low_memory=False, **kwargs)


def read_archive_csv(member: str) -> pd.DataFrame:
    with zipfile.ZipFile(ARCHIVE) as archive:
        with archive.open(member) as handle:
            return pd.read_csv(handle, low_memory=False)


def safe_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def ensure_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column not in out:
            out[column] = np.nan
    return out


def label_for(column: str) -> str:
    special = {
        "rho": "Relative density",
        "rho1": "Relative density after first step",
        "rho2": "Final relative density",
        "T_C": "Temperature",
        "T1_C": "First-step temperature",
        "T2_C": "Second-step temperature",
        "G0_nm": "Initial mean grain size",
        "G1_nm": "Mean grain size after first step",
        "G2_nm": "Final mean grain size",
        "G_mean_nm": "Mean grain size",
        "G50_nm": "Median grain size proxy",
        "G90_nm": "90th-percentile grain size proxy",
        "D50_nm": "Median pore-diameter proxy",
        "D90_nm": "90th-percentile pore-diameter proxy",
        "R_fast": "Matched-density reference/fast grain-size ratio",
        "R_TS": "Matched-density high-temperature/two-step grain-size ratio",
        "reduction_TS": "Two-step grain-size reduction fraction",
        "reduction_fast": "Fast-firing grain-size reduction fraction",
        "rho_dot": "Densification rate",
        "G_dot": "Grain-growth rate",
        "physical_time_s": "Continuous physical time",
        "physical_time_h": "Continuous physical time",
        "closed_fraction": "Closed-pore fraction of total pore volume",
        "phi_closed": "Closed-pore volume fraction",
        "phi_open": "Open-pore volume fraction",
    }
    return special.get(column, column.replace("_", " ").strip().title())


def unit_for(column: str) -> str:
    low = column.lower()
    if column in {"figure_id", "panel_id", "candidate_id", "property_id"}:
        return "identifier"
    if any(token in low for token in ("label", "classification", "status", "role", "reason", "note", "source", "mode", "stage", "mechanism", "measurement", "interpretation", "formula")):
        return "categorical/text"
    if low.endswith("_c") or "temperature" in low:
        return "degC"
    if low.endswith("_nm") or "d50" in low or "d90" in low:
        return "nm"
    if low.endswith("_kj") or low.endswith("_kj_mol") or "q_" in low and "delta" in low:
        return "kJ mol^-1"
    if low.endswith("_s") or low.startswith("tau_"):
        return "s"
    if low.endswith("_h") or "hold_time_h" in low or "time_budget_h" in low:
        return "h"
    if "rate_c_min" in low:
        return "degC min^-1"
    if "rho_dot" in low or "shrinkage_flux" in low:
        return "s^-1"
    if "g_dot" in low:
        return "nm s^-1"
    if "window_width" in low or low.endswith("boundary_c"):
        return "degC"
    if low.startswith("n_") or low.endswith("_count") or low == "count":
        return "count or model number-density proxy"
    if low.startswith("phi_") or any(
        token in low
        for token in (
            "fraction",
            "ratio",
            "reduction",
            "activity",
            "eligibility",
            "theta",
            "pi_",
            "gamma",
            "share",
            "span",
            "factor",
            "flag",
            "attained",
            "complete",
            "boundary_present",
            "success",
        )
    ):
        return "dimensionless/proxy"
    if low in {"rho", "rho0", "rho1", "rho2", "rho_switch", "rho_target", "target_density", "final_density"}:
        return "dimensionless"
    if low.startswith("k_") or "prefactor" in low:
        return "relative factor unless source specifies otherwise"
    return "model-native or categorical; see definition/provenance"


def export_table(
    df: pd.DataFrame,
    relative_path: str,
    *,
    figure_id: str,
    panel_id: str,
    evidence_level: str,
    model_layer: str,
    source_file: str,
    source_table: str,
    notes: str,
    candidate_id: str | int | None = None,
) -> Path:
    out = df.copy()
    if "candidate_id" not in out:
        out["candidate_id"] = candidate_id if candidate_id is not None else ""
    values = {
        "figure_id": figure_id,
        "panel_id": panel_id,
        "evidence_level": evidence_level,
        "model_layer": model_layer,
        "source_file": source_file,
        "source_table": source_table,
        "source_branch_or_commit": f"{BRANCH}@{COMMIT[:12]}",
        "units_reference": "results/figure_source_data_package/all_columns_dictionary.csv",
        "human_readable_label_reference": "results/figure_source_data_package/all_columns_dictionary.csv",
        "derivation_notes": notes,
    }
    for key, value in values.items():
        out[key] = value
    ordered = [column for column in META_COLUMNS if column in out] + [
        column for column in out if column not in META_COLUMNS
    ]
    out = out[ordered]
    path = OUT / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    out.to_csv(path, index=False, compression=compression)
    DATASETS.append(
        {
            "figure_id": figure_id,
            "panel_id": panel_id,
            "relative_path": str(path.relative_to(ROOT)),
            "rows": len(out),
            "columns": len(out.columns),
            "bytes": path.stat().st_size,
            "evidence_level": evidence_level,
            "model_layer": model_layer,
            "source_file": source_file,
            "source_table": source_table,
            "data_status": "available",
            "notes": notes,
        }
    )
    return path


def matched_density_curves(
    histories: pd.DataFrame,
    reference_filter: pd.Series,
    comparison_filter: pd.Series,
    *,
    reference_name: str,
    comparison_name: str,
    rho_min: float | None = None,
    rho_max: float | None = None,
    step: float = 0.001,
) -> pd.DataFrame:
    ref = histories.loc[reference_filter, ["rho", "G_nm"]].dropna().sort_values("rho")
    cmp = histories.loc[comparison_filter, ["rho", "G_nm"]].dropna().sort_values("rho")
    ref = ref.groupby("rho", as_index=False).last()
    cmp = cmp.groupby("rho", as_index=False).last()
    lo = max(ref.rho.min(), cmp.rho.min(), rho_min if rho_min is not None else -np.inf)
    hi = min(ref.rho.max(), cmp.rho.max(), rho_max if rho_max is not None else np.inf)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return pd.DataFrame()
    grid = np.arange(math.ceil(lo / step) * step, hi + step / 2, step)
    g_ref = np.interp(grid, ref.rho, ref.G_nm)
    g_cmp = np.interp(grid, cmp.rho, cmp.G_nm)
    return pd.DataFrame(
        {
            "rho": grid,
            "reference_path": reference_name,
            "comparison_path": comparison_name,
            "G_ref_nm": g_ref,
            "G_fast_nm": g_cmp,
            "R_fast": g_ref / np.maximum(g_cmp, 1e-30),
            "reduction_fast": 1.0 - g_cmp / np.maximum(g_ref, 1e-30),
            "ratio_threshold_1p2": g_ref / np.maximum(g_cmp, 1e-30) >= 1.2,
            "ratio_threshold_1p5": g_ref / np.maximum(g_cmp, 1e-30) >= 1.5,
            "ratio_threshold_2p0": g_ref / np.maximum(g_cmp, 1e-30) >= 2.0,
            "both_paths_attained": True,
        }
    )


def candidate_matched_curves(histories: pd.DataFrame) -> pd.DataFrame:
    renamed = histories.copy()
    path_names = {
        "highT_reference": "highT_reference",
        "success": "two_step_success",
        "lower_failure": "lower_T2_density_failure",
        "upper_failure": "upper_T2_growth_failure",
    }
    renamed["path_label"] = renamed.path_label.map(path_names).fillna(renamed.path_label)
    high = renamed[renamed.path_label == "highT_reference"]
    frames: list[pd.DataFrame] = []
    for path in ["two_step_success", "lower_T2_density_failure", "upper_T2_growth_failure"]:
        other = renamed[renamed.path_label == path]
        h = high[["rho", "G_mean_nm"]].dropna().sort_values("rho").groupby("rho", as_index=False).last()
        o = other[["rho", "G_mean_nm"]].dropna().sort_values("rho").groupby("rho", as_index=False).last()
        lo = max(0.90, h.rho.min(), o.rho.min())
        hi = min(0.99, h.rho.max(), o.rho.max())
        if hi <= lo:
            continue
        grid = np.arange(math.ceil(lo * 1000) / 1000, hi + 0.0005, 0.001)
        gh = np.interp(grid, h.rho, h.G_mean_nm)
        gt = np.interp(grid, o.rho, o.G_mean_nm)
        frames.append(
            pd.DataFrame(
                {
                    "candidate_id": 693168,
                    "path_pair_id": f"highT_reference__{path}",
                    "rho": grid,
                    "G_highT_nm": gh,
                    "G_two_step_nm": gt,
                    "R_TS": gh / np.maximum(gt, 1e-30),
                    "reduction_TS": 1.0 - gt / np.maximum(gh, 1e-30),
                    "highT_attained": True,
                    "two_step_attained": True,
                    "interpolation_support_flag": True,
                    "density_window_label": np.where(
                        (grid >= 0.95) & (grid <= 0.98), "highlight_0p95_0p98", "context"
                    ),
                }
            )
        )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def build_mechanism_chains() -> None:
    nodes = [
        ("FF01", "fast_firing", "Slow ramp", "Long intermediate-temperature residence", "T(t); low_activity_exposure", "negative-control comparator"),
        ("FF02", "fast_firing", "Low-activity interval", "Growth/redistribution occurs before efficient density gain", "activity; E_G", "final attribution"),
        ("FF03", "fast_firing", "Nucleation waiting", "Disconnection/source renewal remains activation limited", "tau_nuc; Theta_nuc", "final attribution"),
        ("FF04", "fast_firing", "Non-densifying coarsening / redistribution", "Slow exposure increases grain size at matched density", "G; PR_exposure", "diagnostic consequence"),
        ("FF05", "fast_firing", "Delayed densification", "Density gain starts after additional growth", "rho_dot; activity", "final attribution"),
        ("FF06", "fast_firing", "Fast ramp", "Rapid crossing of the low-activity interval", "heating_rate_C_min", "final comparator"),
        ("FF07", "fast_firing", "Reduced low-activity exposure", "Less pre-densification growth", "low_activity_exposure", "final attribution"),
        ("FF08", "fast_firing", "Finer G(rho)", "Fast path retains smaller grains at matched density", "R_fast", "observable final evidence"),
        ("TS01", "two_step", "First-step PR preparation", "First step prepares a closed-pore/topology state", "PR_damage_memory; H_PR", "conditional final attribution"),
        ("TS02", "two_step", "Closed-pore store / accommodation state", "Bounded accommodation controls late density gain", "phi_closed; closed_accommodation_factor", "conditional calibration target"),
        ("TS03", "two_step", "Low-T2 density exhaustion", "Insufficient accommodation/shrinkage prevents target attainment", "rho2; closed_shrinkage_flux", "lower boundary"),
        ("TS04", "two_step", "Intermediate-T2 success", "Target density is reached with bounded growth", "classification; R_TS", "conditional Tier B"),
        ("TS05", "two_step", "High-T2 grain-growth activation", "Migration activates above the success band", "G_dot; growth_fraction", "upper boundary"),
        ("TS06", "two_step", "Complete Chen window", "Finite band bracketed by both failures", "window_width_C", "conditional final evidence"),
        ("SH01", "shared", "Source/sink renewal", "Serial event supply", "tau_nuc; activity", "shared framework"),
        ("SH02", "shared", "Point-defect exchange", "Exchange completion timescale", "tau_exchange", "shared framework"),
        ("SH03", "shared", "Point-defect transport", "Transport completion timescale", "tau_transport", "shared framework"),
        ("SH04", "shared", "Pore-connected eligibility", "Only eligible pore topology supports named shrinkage", "phi_connected; connected_eligibility", "shared framework"),
        ("SH05", "shared", "Matched-density scoring", "Compare grain size only where both paths attain density", "rho; R_fast; R_TS", "guardrail"),
    ]
    node_df = pd.DataFrame(nodes, columns=["node_id", "mechanism_group", "label", "short_description", "variable_links", "final_interpretation_role"])
    node_df["supporting_table_or_report"] = "docs/FINAL_FAST_FIRING_AND_TWO_STEP_MECHANISM_SYNTHESIS.md"
    export_table(node_df, "figure_01_mechanism_chains/mechanism_nodes.csv", figure_id="Figure 1", panel_id="nodes", evidence_level="schematic", model_layer="final_synthesis", source_file="docs/FINAL_FAST_FIRING_AND_TWO_STEP_MECHANISM_SYNTHESIS.md", source_table="mechanism synthesis text", notes="Schematic nodes transcribed from the exact-only final synthesis.")
    edges = [
        ("FF01", "FF02", "extends exposure", "final_attribution"), ("FF02", "FF03", "reveals waiting", "final_attribution"),
        ("FF03", "FF04", "permits growth before renewal", "final_attribution"), ("FF04", "FF05", "precedes", "diagnostic"),
        ("FF06", "FF07", "shortens exposure", "final_attribution"), ("FF07", "FF08", "preserves", "final_attribution"),
        ("TS01", "TS02", "prepares", "final_attribution"), ("TS02", "TS03", "exhausts at low T2", "final_attribution"),
        ("TS02", "TS04", "supports within finite range", "final_attribution"), ("TS04", "TS05", "loses bounded growth above range", "final_attribution"),
        ("TS03", "TS06", "defines lower boundary", "final_attribution"), ("TS05", "TS06", "defines upper boundary", "final_attribution"),
        ("SH01", "FF03", "sets renewal waiting", "final_attribution"), ("SH02", "SH05", "completion required", "diagnostic"),
        ("SH03", "SH05", "completion required", "diagnostic"), ("SH04", "TS02", "sets eligible pore inventory", "hypothesis"),
    ]
    edge_df = pd.DataFrame(edges, columns=["source_node", "target_node", "edge_label", "causal_status"])
    edge_df["supporting_evidence"] = "Exact synthesis, causal ablations, and negative-control progression"
    export_table(edge_df, "figure_01_mechanism_chains/mechanism_edges.csv", figure_id="Figure 1", panel_id="edges", evidence_level="schematic", model_layer="final_synthesis", source_file="docs/FINAL_MECHANISM_SYNTHESIS_CAPTIONS.md", source_table="Figure 1 caption", notes="Causal-status labels preserve final-attribution versus hypothesis distinctions.")
    annotations = pd.DataFrame([
        {"annotation_id": "A1", "text": "Fast firing is nucleation-limited and may persist with PR disabled.", "scope": "fast_firing"},
        {"annotation_id": "A2", "text": "Two-step behavior is conditional Tier B and depends on the bounded closed-accommodation trajectory.", "scope": "two_step"},
        {"annotation_id": "A3", "text": "Surrogate screening is promotion-only evidence.", "scope": "shared"},
    ])
    export_table(annotations, "figure_01_mechanism_chains/mechanism_chain_annotations.csv", figure_id="Figure 1", panel_id="annotations", evidence_level="schematic", model_layer="final_synthesis", source_file="docs/NON_CLAIMS_AND_GUARDRAILS.md", source_table="non-claims", notes="Manuscript guardrail annotations.")


def build_exact_property_phase_map() -> tuple[pd.DataFrame, pd.DataFrame]:
    score = read_csv("results/relative_material_property_window_attribution/source_tables/material_property_window_scorecard.csv")
    exact = read_csv("results/relative_material_property_window_attribution/source_tables/material_property_window_exact_promotions.csv")
    merged = exact.merge(score, on="property_id", how="left", suffixes=("_exact_table", ""))
    rename = {
        "classification_exact": "exact_behavior_class",
        "closed_fraction_switch_exact": "closed_fraction_at_switch",
        "closed_accommodation_fraction_exact": "closed_accommodation_fraction",
        "I_low_slow": "low_activity_exposure",
        "I_low_PR_slow": "PR_preparation_memory",
        "Gamma_mig": "growth_activation_number",
    }
    merged = merged.rename(columns=rename)
    merged["closed_store_positive_exact"] = (
        merged["closed_fraction_at_switch"].fillna(0) > 0
    ) & (merged["closed_accommodation_fraction"].fillna(0) > 0)
    merged["surrogate_both_flag"] = merged["classification_screen"].eq("both_pass")
    merged["exact_both_flag"] = merged["exact_behavior_class"].eq("both_pass")
    required = [
        "property_id", "design_stage", "candidate_id", "exact_behavior_class",
        "fast_firing_pass_exact", "two_step_pass_exact", "exact_fast_attained", "exact_two_attained",
        "R_fast_exact", "span_fast_1p5_exact", "reduction_TS_exact", "span_TS_20_exact",
        "Chen_window_width_C_exact", "lower_boundary_present_exact", "upper_boundary_present_exact",
        "closed_store_positive_exact", "surrogate_both_flag", "exact_both_flag",
        "Q_nuc_delta_kJ", "Q_exchange_delta_kJ", "Q_transport_delta_kJ", "Q_growth_delta_kJ",
        "Q_PR_delta_kJ", "Q_closed_delta_kJ", "k_nuc_factor", "k_exchange_factor",
        "k_transport_factor", "k_growth_factor", "k_PR_factor", "k_closed_factor",
        "Theta_nuc", "S_closed_growth", "Pi_PR", "closed_fraction_at_switch",
        "closed_accommodation_fraction", "low_activity_exposure", "PR_preparation_memory",
        "growth_activation_number",
    ]
    merged = ensure_columns(merged, required)[required]
    export_table(merged, "figure_02_exact_property_phase_map/exact_property_phase_map.csv", figure_id="Figure 2", panel_id="phase_map", evidence_level="exact", model_layer="final_synthesis", source_file="results/relative_material_property_window_attribution/source_tables/material_property_window_exact_promotions.csv", source_table="exact promoted union joined to saved screen descriptors", notes="All 1,903 exact-promoted cases; screen descriptors are axes only and do not control exact classification.")
    count_rows = [
        ("screened_rows", 50655, "surrogate_screen"), ("unique_exact_cases", 1903, "exact"),
        ("fast_only", 485, "exact"), ("two_step_only", 119, "exact"),
        ("both_pass", 73, "exact"), ("neither", 1226, "exact"),
        ("surrogate_both", 19880, "surrogate_screen"),
    ]
    counts = pd.DataFrame(count_rows, columns=["behavior_class_or_total", "count", "count_evidence"])
    export_table(counts, "figure_02_exact_property_phase_map/exact_behavior_classification_counts.csv", figure_id="Figure 2", panel_id="counts", evidence_level="exact", model_layer="final_synthesis", source_file="results/final_mechanism_synthesis_and_property_windows/source_tables/exact_behavior_classification_counts.csv", source_table="exact counts plus saved screening totals", notes="Exact counts reproduce the final synthesis; surrogate both is labeled screening evidence.")
    axis = pd.DataFrame([
        ("Theta_nuc", "tau_nuc/(tau_exchange+tau_transport)", "nucleation dominance", "dimensionless/proxy"),
        ("S_closed_growth", "saved closed-shrinkage/growth selectivity group", "two-step selectivity", "dimensionless/proxy"),
        ("Pi_PR", "saved PR preparation group", "PR preparation memory", "dimensionless/proxy"),
        ("growth_activation_number", "saved Gamma_mig", "migration activation", "dimensionless/proxy"),
    ], columns=["axis_variable", "formula_or_source", "human_readable_axis", "units"])
    export_table(axis, "figure_02_exact_property_phase_map/dimensionless_group_axis_options.csv", figure_id="Figure 2", panel_id="axis_options", evidence_level="diagnostic", model_layer="final_synthesis", source_file="results/relative_material_property_window_attribution/source_tables/dimensionless_groups_by_path.csv", source_table="saved dimensionless groups", notes="Axis options only; exact behavior classes remain the evidence.")
    return score, merged


def build_fast_firing() -> pd.DataFrame:
    fast = read_archive_csv(ARCHIVE_FAST_MEMBER)
    fast["ablation"] = fast.ablation.replace({"no_PR_redistribution": "PR-off", "no_nucleation_limitation": "nucleation-facile"})
    fast = fast[fast.ablation.isin(["full", "PR-off", "nucleation-facile", "transport_only", "exchange_limited_variant"])]
    curves = []
    for (material, ablation), group in fast.groupby(["material_id", "ablation"]):
        rates = sorted(group.rate_C_min.unique())
        if 1.0 not in rates:
            continue
        fast_rate = 20.0 if 20.0 in rates else rates[-1]
        curve = matched_density_curves(
            group,
            group.rate_C_min.eq(1.0),
            group.rate_C_min.eq(fast_rate),
            reference_name="1 C/min",
            comparison_name=f"{fast_rate:g} C/min",
            rho_min=0.70,
            rho_max=0.92,
        )
        if curve.empty:
            continue
        curve["material_id"] = material
        curve["property_id"] = "BASE" if ablation == "full" else ablation
        curve["schedule_id"] = f"{material}_{ablation}_1_vs_{fast_rate:g}"
        curve["ablation"] = ablation
        curve["reference_rate_C_min"] = 1.0
        curve["fast_rate_C_min"] = fast_rate
        curve["peak_T_C"] = 1500.0
        curve["hold_time_h"] = 0.25
        curves.append(curve)
    matched = pd.concat(curves, ignore_index=True)
    export_table(matched[matched.ablation.eq("full")], "figure_03_fast_firing_property_window/fast_firing_matched_density_curves.csv", figure_id="Figure 3", panel_id="matched_density", evidence_level="exact", model_layer="fast_firing_material", source_file=f"results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="frozen dense E0021/E0142 histories", notes="Matched-density interpolation of saved exact histories; no simulation or parameter change.")
    export_table(matched, "figure_03_fast_firing_property_window/fast_firing_ablation_curves.csv", figure_id="Figure 3", panel_id="ablations", evidence_level="exact", model_layer="fast_firing_material", source_file=f"results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="frozen dense fast-firing ablations", notes="Includes full, PR-off, nucleation-facile, and diagnostic transport/exchange controls where archived.")
    oat = read_csv("results/final_mechanism_synthesis_and_property_windows/source_tables/fast_firing_OAT_window.csv")
    oat = oat[(oat.perturbed_parameter.eq("Q_nuc_delta_kJ")) | oat.property_id.eq("BASE")]
    export_table(oat, "figure_03_fast_firing_property_window/fast_firing_OAT_Qnuc_window.csv", figure_id="Figure 3", panel_id="Qnuc_OAT", evidence_level="exact", model_layer="fast_firing_material", source_file="results/final_mechanism_synthesis_and_property_windows/source_tables/fast_firing_OAT_window.csv", source_table="exact OAT Q_nuc rows", notes="Contains the 0 and +50 kJ/mol survival edge and adjacent -25/+75 kJ/mol failures.")
    times = fast[(fast.ablation.eq("full")) & fast.rate_C_min.isin([1.0, 20.0])].copy()
    times = times.rename(columns={"tau_nuc": "tau_nuc_s", "tau_exchange": "tau_exchange_s", "tau_transport": "tau_transport_s", "rate_C_min": "heating_rate_C_min"})
    denom = times.tau_exchange_s + times.tau_transport_s
    times["tau_cycle_s"] = times.tau_nuc_s + denom
    times["Theta_nuc"] = times.tau_nuc_s / np.maximum(denom, 1e-30)
    inv = pd.DataFrame({name: 1.0 / np.maximum(times[name], 1e-30) for name in ["tau_nuc_s", "tau_exchange_s", "tau_transport_s"]})
    inv_sum = inv.sum(axis=1)
    times["f_nuc"] = inv.tau_nuc_s / inv_sum
    times["f_exchange"] = inv.tau_exchange_s / inv_sum
    times["f_transport"] = inv.tau_transport_s / inv_sum
    times["low_activity_exposure"] = (1.0 - times.activity.clip(0, 1)) * times.physical_time_s.diff().fillna(0).clip(lower=0)
    times["H_PR"] = times.PR_exposure
    times["H_dens"] = times.groupby(["material_id", "heating_rate_C_min"]).rho.transform(lambda s: s - s.iloc[0])
    keep = ["material_id", "heating_rate_C_min", "rho", "T_C", "physical_time_h", "tau_nuc_s", "tau_exchange_s", "tau_transport_s", "tau_cycle_s", "activity", "Theta_nuc", "f_nuc", "f_exchange", "f_transport", "low_activity_exposure", "PR_exposure", "H_PR", "H_dens"]
    export_table(times[keep], "figure_03_fast_firing_property_window/fast_firing_timescale_groups.csv", figure_id="Figure 3", panel_id="timescales", evidence_level="diagnostic", model_layer="fast_firing_material", source_file=f"results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="frozen dense histories with derived serial timescale groups", notes="Derived algebraically from saved instantaneous timescales; diagnostic, not a new mechanism.")

    clean = fast[fast.ablation.isin(["full", "PR-off", "nucleation-facile"])].copy()
    clean = clean.rename(columns={"rate_C_min": "heating_rate_C_min", "tau_nuc": "tau_nuc_s", "tau_exchange": "tau_exchange_s", "tau_transport": "tau_transport_s"})
    clean["reference_rate_C_min"] = 1.0
    clean["peak_T_C"] = 1500.0
    clean["hold_h"] = 0.25
    clean["pass_fast_firing_rule"] = np.nan
    export_table(clean, "fast_firing/clean_fast_firing_T_rho_G_time.csv", figure_id="Figure 3", panel_id="clean_time_histories", evidence_level="exact", model_layer="fast_firing_material", source_file=f"results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="frozen dense fast histories", notes="Clean external-plotting history package for E0021/E0142 full, PR-off, and nucleation-facile modes.")
    export_table(clean[[c for c in clean.columns if c in {"material_id", "ablation", "heating_rate_C_min", "rho", "G_nm", "activity", "tau_nuc_s", "tau_exchange_s", "tau_transport_s", "PR_exposure", "reference_rate_C_min", "peak_T_C", "hold_h"}]], "fast_firing/clean_fast_firing_G_vs_rho.csv", figure_id="Figure 3", panel_id="clean_G_rho", evidence_level="exact", model_layer="fast_firing_material", source_file=f"results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="frozen dense fast histories", notes="Direct G-rho rows; no matched-density interpolation.")
    export_table(matched, "fast_firing/clean_fast_firing_ratio_vs_rho.csv", figure_id="Figure 3", panel_id="clean_ratio", evidence_level="exact", model_layer="fast_firing_material", source_file=f"results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="matched-density interpolation", notes="All matched-density ratios are scored only over jointly attained density ranges.")
    rate_rows = []
    for (material, ablation, rate), group in fast.groupby(["material_id", "ablation", "rate_C_min"]):
        if rate == 1.0:
            continue
        curve = matched_density_curves(group=pd.DataFrame()) if False else matched_density_curves(fast[(fast.material_id.eq(material)) & fast.ablation.eq(ablation)], fast[(fast.material_id.eq(material)) & fast.ablation.eq(ablation)].rate_C_min.eq(1.0), fast[(fast.material_id.eq(material)) & fast.ablation.eq(ablation)].rate_C_min.eq(rate), reference_name="1 C/min", comparison_name=f"{rate:g} C/min", rho_min=0.70, rho_max=0.92)
        if curve.empty:
            continue
        rate_rows.append({"material_id": material, "ablation": ablation, "heating_rate_C_min": rate, "reference_rate_C_min": 1.0, "max_R_fast": curve.R_fast.max(), "median_R_fast": curve.R_fast.median(), "span_fast_1p5": curve.loc[curve.R_fast.ge(1.5), "rho"].max() - curve.loc[curve.R_fast.ge(1.5), "rho"].min() if curve.R_fast.ge(1.5).sum() > 1 else 0.0, "pass_fast_firing_rule": bool(curve.R_fast.ge(1.5).any() and (curve.loc[curve.R_fast.ge(1.5), "rho"].max() - curve.loc[curve.R_fast.ge(1.5), "rho"].min() >= 0.03))})
    rate_map = pd.DataFrame(rate_rows)
    export_table(rate_map, "fast_firing/clean_fast_firing_heating_rate_map.csv", figure_id="Figure 3", panel_id="heating_rate_map", evidence_level="exact", model_layer="fast_firing_material", source_file=f"results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="matched-density rate summary", notes="Frozen 1 C/min reference compared with each archived faster ramp.")
    preserve = read_csv("results/local_region_decoder_corrected_dynamic_search/fast_firing_preservation.csv")
    export_table(preserve, "fast_firing/clean_fast_firing_ablation_summary.csv", figure_id="Figure 3", panel_id="ablation_summary", evidence_level="exact", model_layer="fast_firing_material", source_file="results/local_region_decoder_corrected_dynamic_search/fast_firing_preservation.csv", source_table="frozen causal preservation audit", notes="Confirms full and PR-off retention and nucleation-facile loss for E0021/E0142.")
    return clean


def build_two_step() -> pd.DataFrame:
    hist = read_csv("results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_histories.csv")
    curves = candidate_matched_curves(hist)
    repro = read_csv("results/audit_candidate_693168_closed_accommodation/candidate_693168_reproduction_summary.csv").iloc[0]
    curves["T1_C"] = repro.T_at_switch_C
    curves["T2_C"] = curves.path_pair_id.map({"highT_reference__lower_T2_density_failure": 900.0, "highT_reference__two_step_success": 1100.0, "highT_reference__upper_T2_growth_failure": 1220.0})
    curves["rho_switch"] = repro.exact_switch_density
    curves["G0_nm"] = repro.G1_nm / (1 + repro.first_step_growth_fraction)
    curves["G1_nm"] = repro.G1_nm
    curves["first_step_growth_fraction"] = repro.first_step_growth_fraction
    export_table(curves, "figure_04_two_step_property_window/two_step_matched_density_curves_693168.csv", figure_id="Figure 4", panel_id="matched_density", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_histories.csv", source_table="frozen exact dense histories", notes="Matched-density interpolation for success, lower-failure, and upper-failure paths over available overlap; 0.95-0.98 highlighted where attained.", candidate_id=693168)
    oat = read_csv("results/final_mechanism_synthesis_and_property_windows/source_tables/two_step_OAT_window.csv")
    for parameter, filename, panel in [
        ("Q_closed_delta_kJ", "two_step_OAT_Qclosed_window.csv", "Qclosed_OAT"),
        ("k_PR_factor", "two_step_PR_prefactor_threshold.csv", "PR_threshold"),
        ("k_growth_factor", "two_step_growth_prefactor_threshold.csv", "growth_threshold"),
    ]:
        subset = oat[(oat.perturbed_parameter.eq(parameter)) | oat.property_id.eq("BASE")]
        export_table(subset, f"figure_04_two_step_property_window/{filename}", figure_id="Figure 4", panel_id=panel, evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/final_mechanism_synthesis_and_property_windows/source_tables/two_step_OAT_window.csv", source_table=f"exact OAT {parameter}", notes="Frozen exact OAT; no target, budget, or parameter changes.", candidate_id=693168)
    summaries = []
    for path, group in hist.groupby("path_label"):
        group = group.sort_values("physical_time_s")
        dt = group.physical_time_s.diff().fillna(0).clip(lower=0)
        summaries.append({
            "candidate_id": 693168, "path_label": path,
            "closed_fraction_at_switch": group.loc[group.stage.eq("step2"), "closed_fraction"].iloc[0] if group.stage.eq("step2").any() else np.nan,
            "closed_fraction_at_target": group.closed_fraction.iloc[-1],
            "closed_accommodation_fraction": group.closed_accommodation_factor.iloc[-1],
            "PR_damage_memory": group.PR_damage_memory.iloc[-1],
            "closed_shrinkage_share": float((group.closed_shrinkage_flux * dt).sum() / max(((group.closed_shrinkage_flux + group.open_shrinkage_flux) * dt).sum(), 1e-30)),
            "open_shrinkage_share": float((group.open_shrinkage_flux * dt).sum() / max(((group.closed_shrinkage_flux + group.open_shrinkage_flux) * dt).sum(), 1e-30)),
            "Chen_window_width_C": repro.window_width_C,
            "S_closed_growth": np.nan,
        })
    selectivity = pd.DataFrame(summaries)
    export_table(selectivity, "figure_04_two_step_property_window/two_step_selectivity_groups.csv", figure_id="Figure 4", panel_id="selectivity", evidence_level="diagnostic", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_histories.csv", source_table="integrated saved open/closed flux histories", notes="Pathwise diagnostic groups; unavailable S_closed_growth is explicitly NaN.", candidate_id=693168)
    return hist


def build_chen_and_candidate(hist: pd.DataFrame, clean_fast: pd.DataFrame) -> None:
    fine = read_csv("results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_classification_points_fine.csv")
    repro = read_csv("results/audit_candidate_693168_closed_accommodation/candidate_693168_reproduction_summary.csv").iloc[0]
    fine["T1_C"] = repro.T_at_switch_C
    fine["rho_switch"] = repro.exact_switch_density
    fine["G0_nm"] = repro.G1_nm / (1 + repro.first_step_growth_fraction)
    fine["G1_nm"] = repro.G1_nm
    fine["first_step_growth_fraction"] = repro.first_step_growth_fraction
    fine["target_density"] = 0.98
    fine["lower_boundary_flag"] = safe_bool(fine.lower_boundary_marker)
    fine["upper_boundary_flag"] = safe_bool(fine.upper_boundary_marker)
    fine["success_band_flag"] = fine.classification.eq("SUCCESS")
    fine["practical_T2_less_than_T1"] = fine.T2_C < fine.T1_C
    fine["window_width_C"] = repro.window_width_C
    rename = {"closed_accommodation_state": "closed_accommodation_used", "PR_damage_state": "PR_damage_memory"}
    fine = fine.rename(columns=rename)
    requested = ["candidate_id", "T1_C", "T2_C", "rho_switch", "G0_nm", "G1_nm", "first_step_growth_fraction", "rho2", "G2_nm", "growth_fraction", "target_density", "target_attained", "classification", "lower_boundary_flag", "upper_boundary_flag", "success_band_flag", "practical_T2_less_than_T1", "window_width_C", "closed_fraction", "phi_closed", "phi_open", "phi_connected", "phi_iso", "closed_accommodation_capacity", "closed_accommodation_used", "closed_accommodation_factor", "open_shrinkage_flux", "closed_shrinkage_flux", "open_shrinkage_contribution", "closed_shrinkage_contribution", "PR_damage_memory", "H_PR", "H_dens", "w_PR", "w_dens", "migration_factor", "G_dot", "rho_dot"]
    fine = ensure_columns(fine, requested)[requested]
    export_table(fine, "figure_05_chen_window_boundaries/candidate_693168_fine_T2_classification.csv", figure_id="Figure 5", panel_id="fine_classification", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_classification_points_fine.csv", source_table="101-point exact fine T2 scan", notes="Unavailable per-point internal channels remain NaN and are listed in the missing-channel table.", candidate_id=693168)
    bounds = read_csv("results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_window_boundaries_fine.csv")
    export_table(bounds, "figure_05_chen_window_boundaries/candidate_693168_fine_Chen_boundaries.csv", figure_id="Figure 5", panel_id="boundaries", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_window_boundaries_fine.csv", source_table="fine exact boundary summary", notes="Complete lower and upper boundary brackets for the frozen candidate.", candidate_id=693168)
    diag = read_csv("results/publication_style_sintering_figures_693168/source_tables/T2_diagnostic_curves.csv")
    export_table(diag, "figure_05_chen_window_boundaries/candidate_693168_T2_diagnostic_curves.csv", figure_id="Figure 5", panel_id="T2_diagnostics", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/publication_style_sintering_figures_693168/source_tables/T2_diagnostic_curves.csv", source_table="publication-style exact T2 diagnostics", notes="Final density, growth, pore state, and open/closed integrated shrinkage versus T2.", candidate_id=693168)
    tiles = fine[["candidate_id", "T1_C", "T2_C", "rho_switch", "G0_nm", "G1_nm", "rho2", "G2_nm", "growth_fraction", "classification", "success_band_flag", "lower_boundary_flag", "upper_boundary_flag", "practical_T2_less_than_T1", "window_width_C"]]
    export_table(tiles, "figure_05_chen_window_boundaries/candidate_693168_filled_map_tiles.csv", figure_id="Figure 5", panel_id="filled_tiles", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_classification_points_fine.csv", source_table="fine exact classification tiles", notes="Contains failures and successes; not a sparse success-only table.", candidate_id=693168)

    ab = read_csv("results/local_region_decoder_corrected_dynamic_search/ablation_summary.csv")
    full_ab = ab[ab.ablation.eq("full_model")][
        ["candidate_id", "lower_bracketed", "upper_bracketed"]
    ].rename(
        columns={
            "lower_bracketed": "lower_boundary_present",
            "upper_bracketed": "upper_boundary_present",
        }
    )
    family = read_csv(
        "results/reframe_tierB_experimental_plausibility/tierB_candidate_reinterpretation.csv"
    ).merge(full_ab, on="candidate_id", how="left")
    family = family.rename(columns={"span_above_20pct": "span20", "ablations_that_destroy_result": "destructive_ablations"})
    family["tier_status"] = "conditional Tier B"
    family["validation_status"] = "not validation"
    keep = ["candidate_id", "G0_nm", "G1_nm", "first_step_growth_fraction", "median_reduction", "minimum_reduction", "maximum_reduction", "span20", "window_width_C", "T2_success_range", "lower_boundary_present", "upper_boundary_present", "closed_fraction_at_switch", "closed_fraction_at_target", "closed_shrinkage_share", "open_shrinkage_share", "PR_damage_contribution", "destructive_ablations", "robustness_count", "interpretation", "artifact_flag", "tier_status", "validation_status"]
    family = family.rename(columns={"minimum_reduction": "min_reduction", "maximum_reduction": "max_reduction"})
    keep = [c.replace("minimum_reduction", "min_reduction").replace("maximum_reduction", "max_reduction") for c in keep]
    family = ensure_columns(family, keep)[keep]
    export_table(family, "figure_06_six_tierB_family/six_TierB_candidate_summary.csv", figure_id="Figure 6", panel_id="candidate_summary", evidence_level="exact", model_layer="local_region_two_step", source_file="results/reframe_tierB_experimental_plausibility/tierB_candidate_reinterpretation.csv", source_table="six exact candidate reinterpretation", notes="All six exact Tier-B candidates; 693168 remains conditional and unvalidated.")
    ab["retained_TierB"] = ab.tier.eq("Tier_B")
    ab["high_density_attainment"] = (ab.rho_high_final >= 0.98) & (ab.rho_two_final >= 0.98)
    ab["mechanism_loss_category"] = np.select(
        [~ab.high_density_attainment, ~safe_bool(ab.lower_bracketed), ~safe_bool(ab.upper_bracketed), ab.median_reduction.lt(0.2)],
        ["loses_high_density_attainment", "loses_lower_boundary", "loses_upper_boundary", "loses_trajectory_reduction"],
        default="noncausal",
    )
    export_table(ab, "figure_06_six_tierB_family/six_TierB_ablation_matrix.csv", figure_id="Figure 6", panel_id="ablation_matrix", evidence_level="exact", model_layer="local_region_two_step", source_file="results/local_region_decoder_corrected_dynamic_search/ablation_summary.csv", source_table="six-candidate exact ablation matrix", notes="Mechanism-loss category is derived from saved attainment, boundary, and reduction fields.")
    robust = read_csv("results/local_region_decoder_corrected_dynamic_search/robustness_summary.csv")
    export_table(robust, "figure_06_six_tierB_family/six_TierB_robustness_summary.csv", figure_id="Figure 6", panel_id="robustness", evidence_level="exact", model_layer="local_region_two_step", source_file="results/local_region_decoder_corrected_dynamic_search/robustness_summary.csv", source_table="bounded exact initial-condition robustness", notes="No selective budget extension or new cases.")
    corr = family[["candidate_id", "closed_fraction_at_switch", "median_reduction", "min_reduction", "max_reduction", "window_width_C", "tier_status", "validation_status"]]
    export_table(corr, "figure_06_six_tierB_family/six_TierB_closed_fraction_vs_reduction.csv", figure_id="Figure 6", panel_id="closed_fraction_relation", evidence_level="exact", model_layer="local_region_two_step", source_file="results/reframe_tierB_experimental_plausibility/tierB_candidate_reinterpretation.csv", source_table="six-candidate exact summary", notes="Descriptive correlation source; not a universal material law.")

    path_map = {"success": "two_step_success", "lower_failure": "lower_T2_density_failure", "upper_failure": "upper_T2_growth_failure"}
    candidate_hist = hist.copy()
    candidate_hist["path_label"] = candidate_hist.path_label.replace(path_map)
    candidate_hist["target_crossing_flag"] = candidate_hist.groupby("path_label").target_attained.transform(lambda s: safe_bool(s).diff().fillna(False).astype(bool))
    candidate_hist["rho_target_crossed"] = safe_bool(candidate_hist.target_attained)
    candidate_hist["switch_flag"] = candidate_hist.stage.eq("step2") & ~candidate_hist.groupby("path_label").stage.shift().eq("step2")
    candidate_hist["T1_reached_flag"] = candidate_hist.T_C.ge(repro.T_at_switch_C - 1e-6)
    required_candidate = ["candidate_id", "path_label", "stage", "physical_time_s", "physical_time_h", "T_C", "rho", "G_mean_nm", "G50_nm", "G90_nm", "rho_dot", "G_dot", "target_crossing_flag", "rho_target_crossed", "switch_flag", "T1_reached_flag", "phi_open", "phi_connected", "phi_GBseg", "phi_TJ", "phi_iso", "phi_closed", "N_open", "N_connected", "N_GBseg", "N_TJ", "N_iso", "N_closed", "closed_fraction", "isolated_fraction", "connected_fine_pore_fraction", "removable_pore_fraction", "D50_nm", "D90_nm", "large_pore_fraction", "closed_accommodation_capacity", "closed_accommodation_used", "closed_accommodation_factor", "closed_shrinkage_flux", "open_shrinkage_flux", "closed_pore_contribution_to_rho_dot", "open_pore_contribution_to_rho_dot", "PR_damage_memory", "H_PR", "H_dens", "w_PR", "w_dens", "cumulative_PR_surface_energy_loss", "cumulative_densifying_work", "cumulative_non_densifying_work", "migration_factor", "X_J", "Lambda_TJ", "K_TJ", "Lambda_over_K_TJ", "P_comp_TJ", "residual_stress"]
    candidate_hist = ensure_columns(candidate_hist, required_candidate)[required_candidate]
    fast_rows = clean_fast[
        clean_fast.material_id.eq("E0021")
        & clean_fast.ablation.eq("full")
        & clean_fast.heating_rate_C_min.isin([1.0, 20.0])
    ].copy()
    fast_pack = pd.DataFrame(index=fast_rows.index)
    fast_pack["candidate_id"] = ""
    fast_pack["path_label"] = np.where(
        fast_rows.heating_rate_C_min.eq(1.0),
        "fast_firing_reference",
        "fast_firing_fast",
    )
    fast_pack["stage"] = fast_rows.stage
    for target, source in [
        ("physical_time_s", "physical_time_s"),
        ("physical_time_h", "physical_time_h"),
        ("T_C", "T_C"),
        ("rho", "rho"),
        ("G_mean_nm", "G_nm"),
        ("rho_dot", "rho_dot"),
        ("G_dot", "G_dot_nm_s"),
    ]:
        fast_pack[target] = fast_rows[source].to_numpy()
    fast_pack = ensure_columns(fast_pack, required_candidate)[required_candidate]
    combined_hist = pd.concat([candidate_hist, fast_pack], ignore_index=True)
    export_table(combined_hist, "candidate_693168/dense_time_histories.csv", figure_id="Candidate 693168", panel_id="dense_histories", evidence_level="exact", model_layer="candidate_693168_audit", source_file=f"results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_histories.csv; results/1_Backup_of_prior_runs.zip:{ARCHIVE_FAST_MEMBER}", source_table="frozen exact candidate and E0021 fast-firing histories", notes="Continuous time for highT, success, lower-failure, upper-failure, fast-reference, and fast paths. TJ-specific and residual-stress channels unavailable here remain NaN.", candidate_id=693168)
    candidate_curves = candidate_matched_curves(hist)
    export_table(candidate_curves, "candidate_693168/matched_density_ratio_curves.csv", figure_id="Candidate 693168", panel_id="matched_density", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_histories.csv", source_table="matched-density interpolation", notes="Jointly attained density ranges only.", candidate_id=693168)
    export_table(fine, "candidate_693168/fine_T2_scan.csv", figure_id="Candidate 693168", panel_id="fine_T2", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_classification_points_fine.csv", source_table="101-point exact T2 scan", notes="Complete fine scan including failures.", candidate_id=693168)
    export_table(bounds, "candidate_693168/Chen_boundaries_fine.csv", figure_id="Candidate 693168", panel_id="Chen_boundaries", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_window_boundaries_fine.csv", source_table="fine boundary summary", notes="Both boundaries bracketed.", candidate_id=693168)
    for filename, cols, panel in [
        ("closed_accommodation_histories.csv", [c for c in hist.columns if c in {"candidate_id", "path_label", "stage", "physical_time_s", "physical_time_h", "T_C", "rho", "phi_closed", "closed_fraction", "closed_accommodation_capacity", "closed_accommodation_available", "closed_accommodation_used", "closed_accommodation_factor", "closed_shrinkage_flux", "open_shrinkage_flux"}], "closed_accommodation"),
        ("pore_store_histories.csv", [c for c in hist.columns if c.startswith("phi_") or c.startswith("N_") or c in {"candidate_id", "path_label", "stage", "physical_time_s", "physical_time_h", "T_C", "rho", "D50_nm", "D90_nm", "large_pore_fraction", "connected_fine_pore_fraction", "removable_pore_fraction", "closed_fraction", "isolated_fraction"}], "pore_store"),
        ("PR_energy_partition_histories.csv", [c for c in hist.columns if c in {"candidate_id", "path_label", "stage", "physical_time_s", "physical_time_h", "T_C", "rho", "PR_damage_memory", "PR_redistribution_rate", "H_PR", "H_dens", "w_PR", "w_dens", "cumulative_PR_surface_energy_loss", "cumulative_densifying_work", "cumulative_non_densifying_work"}], "PR_partition"),
    ]:
        export_table(hist[cols], f"candidate_693168/{filename}", figure_id="Candidate 693168", panel_id=panel, evidence_level="diagnostic", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_histories.csv", source_table="named saved diagnostic channels", notes="Named model diagnostics; work channels are model proxies, not absolute energies.", candidate_id=693168)
    ab693 = read_csv("results/audit_candidate_693168_closed_accommodation/candidate_693168_ablation_audit.csv")
    export_table(ab693, "candidate_693168/ablation_summary.csv", figure_id="Candidate 693168", panel_id="ablations", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/candidate_693168_ablation_audit.csv", source_table="exact ablation audit", notes="Frozen causal audit; conditional Tier B, not validation.", candidate_id=693168)
    robust693 = read_csv("results/audit_candidate_693168_closed_accommodation/candidate_693168_extended_robustness.csv")
    export_table(robust693, "candidate_693168/extended_robustness.csv", figure_id="Candidate 693168", panel_id="robustness", evidence_level="exact", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/candidate_693168_extended_robustness.csv", source_table="bounded robustness audit", notes="Existing bounded cases only.", candidate_id=693168)
    preserve = read_csv("results/audit_candidate_693168_closed_accommodation/fast_firing_preservation_audit.csv")
    export_table(preserve, "candidate_693168/fast_firing_preservation.csv", figure_id="Candidate 693168", panel_id="fast_preservation", evidence_level="exact", model_layer="fast_firing_material", source_file="results/audit_candidate_693168_closed_accommodation/fast_firing_preservation_audit.csv", source_table="frozen E0021/E0142 preservation audit", notes="Two-step topology layer does not modify shared-state material densification.", candidate_id=693168)
    missing_channels = ["Lambda_TJ", "K_TJ", "Lambda_over_K_TJ", "P_comp_TJ", "migration_factor", "X_J", "residual_stress", "trapped_gas_pressure", "absolute_energy_channels"]
    miss_df = pd.DataFrame({"channel": missing_channels, "availability": "unavailable in frozen candidate-693168 exact history", "handling": "column retained as NaN where requested; no value invented"})
    export_table(miss_df, "candidate_693168/missing_candidate_693168_channels.csv", figure_id="Candidate 693168", panel_id="missing_channels", evidence_level="diagnostic", model_layer="candidate_693168_audit", source_file="results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_histories.csv", source_table="schema comparison", notes="Closed-accommodation is an implemented bounded proxy; no gas law or hidden closed-pore Lambda/K law is asserted.", candidate_id=693168)


def build_chen_maps() -> None:
    specs = [
        ("chen_map_T1_T2_classification.csv", "T1_T2"),
        ("chen_map_G1_T2_classification.csv", "G1_T2"),
        ("chen_map_switch_density_T2_classification.csv", "switch_density_T2"),
    ]
    for filename, map_type in specs:
        source = f"results/publication_style_sintering_figures_693168/source_tables/{filename}"
        df = read_csv(source)
        df["map_type"] = map_type
        df["rho_target"] = df.get("target_density", 0.98)
        df["rho2"] = df.get("rho_final")
        df["G2_nm"] = df.get("G_final_nm")
        df["success_flag"] = df.classification.eq("SUCCESS")
        df["lower_boundary_region_flag"] = df.classification.eq("DENSIFICATION_EXHAUSTION_FAILURE")
        df["upper_boundary_region_flag"] = df.classification.eq("GRAIN_GROWTH_FAILURE")
        df["unattainable_flag"] = df.classification.str.contains("UNATTAINABLE|INELIGIBLE", na=False)
        df["censored_flag"] = df.classification.str.contains("CENSORED", na=False) | safe_bool(df.get("numerical_censor", pd.Series(False, index=df.index)))
        df["practical_T2_less_than_T1"] = df.T2_C < df.T1_C
        export_table(df, f"chen_maps/{filename}", figure_id="Figure 5", panel_id=map_type, evidence_level="exact", model_layer="candidate_693168_audit", source_file=source, source_table="filled exact publication map", notes="Complete grid including all saved failure/success/censored/ineligible rows; no sparse filtering.", candidate_id=693168)
    dictionary = pd.DataFrame([
        ("classification", "Categorical outcome at fixed map point", "categorical"),
        ("rho2", "Final second-step density", "dimensionless"),
        ("G2_nm", "Final second-step mean grain size", "nm"),
        ("growth_fraction", "Second-step fractional grain growth", "dimensionless"),
        ("success_flag", "Density attained with bounded growth", "boolean"),
        ("lower_boundary_region_flag", "Density-exhaustion region", "boolean"),
        ("upper_boundary_region_flag", "Grain-growth region", "boolean"),
    ], columns=["column", "definition", "units"])
    export_table(dictionary, "chen_maps/chen_map_plot_dictionary.csv", figure_id="Figure 5", panel_id="map_dictionary", evidence_level="diagnostic", model_layer="candidate_693168_audit", source_file="docs/PUBLICATION_STYLE_SINTERING_FIGURES_693168.md", source_table="Chen map definitions", notes="External plotting dictionary.", candidate_id=693168)


def build_surrogate_warning(score: pd.DataFrame, exact: pd.DataFrame) -> None:
    counts = pd.DataFrame([
        ("surrogate_both", 19880, "surrogate_screen", "promotion tool only"),
        ("exact_both", 73, "exact", "final classification"),
        ("overprediction_factor", 19880 / 73, "derived", "surrogate_both/exact_both"),
    ], columns=["metric", "value", "count_evidence", "interpretation"])
    export_table(counts, "figure_07_surrogate_vs_exact/surrogate_vs_exact_counts.csv", figure_id="Figure 7", panel_id="counts", evidence_level="diagnostic", model_layer="final_synthesis", source_file="results/final_mechanism_synthesis_and_property_windows/source_tables/surrogate_vs_exact_comparison.csv", source_table="saved screening/exact counts", notes="Surrogate rows are not final evidence.")
    joined = score[["property_id", "classification_screen"]].merge(exact[["property_id", "exact_behavior_class", "exact_fast_attained", "exact_two_attained", "fast_firing_pass_exact", "two_step_pass_exact", "lower_boundary_present_exact", "upper_boundary_present_exact"]], on="property_id", how="inner")
    confusion = joined.groupby(["classification_screen", "exact_behavior_class"], dropna=False).size().reset_index(name="count").rename(columns={"classification_screen": "surrogate_predicted_class", "exact_behavior_class": "exact_class"})
    confusion["pass_fail_reason"] = np.where(confusion.surrogate_predicted_class.eq(confusion.exact_class), "agreement", "surrogate/exact disagreement")
    confusion["false_positive_reason"] = np.where((confusion.surrogate_predicted_class.eq("both_pass")) & ~confusion.exact_class.eq("both_pass"), "one or both exact trajectory rules fail", "not_applicable")
    export_table(confusion, "figure_07_surrogate_vs_exact/surrogate_exact_confusion.csv", figure_id="Figure 7", panel_id="confusion", evidence_level="diagnostic", model_layer="final_synthesis", source_file="results/relative_material_property_window_attribution/source_tables/material_property_window_scorecard.csv; material_property_window_exact_promotions.csv", source_table="exact-union confusion matrix", notes="Only promoted rows have an exact class; screening is not treated as evidence.")
    phase = exact.copy()
    phase["surrogate_predicted_class"] = phase.get("classification_screen")
    phase["exact_class"] = phase.exact_behavior_class
    export_table(phase, "figure_07_surrogate_vs_exact/surrogate_exact_phase_points.csv", figure_id="Figure 7", panel_id="phase_points", evidence_level="diagnostic", model_layer="final_synthesis", source_file="results/relative_material_property_window_attribution/source_tables/material_property_window_exact_promotions.csv", source_table="exact promoted union with saved screening axes", notes="Exact-promoted points only; surrogate values are labels/axes, not final evidence.")


def build_falsification() -> None:
    rows = [
        ("Matched-density G_mean, G50, G90", "fast-firing and two-step grain trajectory", "Fast and successful two-step paths retain smaller distributions", "No attained matched-density separation", "0.80-0.98", "fast ramps and two-step", "highest", "interrupted microscopy/statistical grain sizing", "G_mean_nm; G50_nm; G90_nm", "trajectory"),
        ("Densification-onset time during ramps", "nucleation-limited onset", "Fast ramp crosses low-activity interval before substantial growth", "Onset insensitive to nucleation-facile perturbation", "0.65-0.90", "1-100 C/min ramps", "highest", "dilatometry with interruption", "activity; tau_nuc_s; rho_dot", "fast_firing"),
        ("Exchange and transport relaxation", "serial completion", "Exchange/transport do not alone reproduce nucleation-limited ratio", "Transport-only control reproduces full effect", "0.65-0.90", "ramp relaxation", "medium", "tracer/relaxation measurement", "tau_exchange_s; tau_transport_s", "fast_firing"),
        ("3D open/closed pore fraction", "closed-pore transition", "Closed fraction rises along successful preparation trajectory", "No schedule-dependent closed fraction", "0.88-0.98", "interrupted first/second step", "highest", "FIB-SEM or X-ray nanotomography", "phi_open; phi_closed; closed_fraction", "two_step"),
        ("Pore D50/D90 and large-pore tail", "PR redistribution", "First-step preparation shifts pore distribution", "No observable redistribution despite predicted PR memory", "0.75-0.95", "interrupted first step", "high", "3D tomography/TEM", "D50_nm; D90_nm; large_pore_fraction", "two_step"),
        ("Connected fine-pore fraction/percolation", "open-pore eligibility", "Connectivity persists through the shrinkage-relevant interval", "Connectivity vanishes before predicted open shrinkage", "0.75-0.95", "interrupted paths", "high", "3D segmentation/percolation analysis", "connected_fine_pore_fraction; phi_connected", "two_step"),
        ("Trapped-gas pressure or accommodation proxy", "bounded closed accommodation", "Finite accommodation is consumed during successful T2 hold", "Unlimited shrinkage or absent accommodation signature", "0.90-0.98", "T2 holds", "highest", "gas analysis, pressure proxy, pore-volume evolution", "closed_accommodation_capacity; closed_accommodation_used", "two_step"),
        ("Grain-growth mobility versus T2", "upper Chen boundary", "Mobility activates above the success band", "No upper mobility boundary", "0.88-0.98", "fine T2 scan", "highest", "interrupted grain-size kinetics", "G_dot; growth_fraction", "Chen_map"),
        ("Interrupted first/second-step tomography", "topology persistence", "Prepared topology persists into successful second step", "Immediate collapse to isothermal topology", "0.75-0.98", "paired interruptions", "highest", "correlative 3D tomography", "pore stores; PR_damage_memory", "two_step"),
    ]
    columns = ["measurement", "mechanism_constrained", "expected_signature_if_model_correct", "falsification_signature", "density_range", "schedule", "priority", "possible_experimental_method", "model_variable_constrained", "figure_panel_group"]
    targets = pd.DataFrame(rows, columns=columns)
    export_table(targets, "figure_08_experimental_falsification/falsification_targets.csv", figure_id="Figure 8", panel_id="targets", evidence_level="schematic", model_layer="final_synthesis", source_file="docs/EXPERIMENTAL_CALIBRATION_PLAN.md", source_table="calibration priorities", notes="Experimental design guidance; not simulated evidence.")
    export_table(targets[["measurement", "mechanism_constrained", "model_variable_constrained", "expected_signature_if_model_correct", "falsification_signature"]], "figure_08_experimental_falsification/measurement_to_mechanism_map.csv", figure_id="Figure 8", panel_id="measurement_map", evidence_level="schematic", model_layer="final_synthesis", source_file="docs/EXPERIMENTAL_CALIBRATION_PLAN.md", source_table="measurement-mechanism mapping", notes="Direct measurement-to-model mapping.")
    priority = targets.sort_values("priority", key=lambda s: s.map({"highest": 0, "high": 1, "medium": 2})).copy()
    priority["calibration_sequence"] = np.arange(1, len(priority) + 1)
    export_table(priority, "figure_08_experimental_falsification/calibration_priority_table.csv", figure_id="Figure 8", panel_id="priorities", evidence_level="schematic", model_layer="final_synthesis", source_file="docs/EXPERIMENTAL_CALIBRATION_PLAN.md", source_table="ordered calibration plan", notes="Priority ordering only; no validation claim.")


def build_equation_maps() -> None:
    eq = read_csv("results/equation_functional_form_audit/equation_registry.csv")
    variables = read_csv("results/equation_functional_form_audit/variable_definitions.csv")
    params = read_csv("results/equation_functional_form_audit/parameter_definitions.csv").rename(columns={"parameter": "variable", "evidence_scope": "scope"})
    defs = pd.concat([variables, params[["variable", "definition", "units", "scope"]]], ignore_index=True).drop_duplicates("variable")
    defs["symbol"] = defs.variable
    defs["OMML_or_LaTeX_label"] = defs.variable
    defs["plain_text_key"] = defs.variable
    defs["human_readable_name"] = defs.variable.map(label_for)
    defs["relevant_equation_ids"] = defs.variable.map(lambda symbol: ";".join(eq.loc[eq.variable_definitions.astype(str).str.contains(re.escape(str(symbol)), regex=True, na=False), "equation_id"].astype(str).tolist()))
    defs["source_file/function"] = "results/equation_functional_form_audit/equation_registry.csv"
    defs["recommended_axis_label"] = defs.human_readable_name + " [" + defs.units.fillna("dimensionless/proxy") + "]"
    keep = ["symbol", "OMML_or_LaTeX_label", "plain_text_key", "human_readable_name", "units", "definition", "relevant_equation_ids", "source_file/function", "recommended_axis_label"]
    export_table(defs[keep], "supplement/equation_label_map.csv", figure_id="Supplement", panel_id="equation_labels", evidence_level="diagnostic", model_layer="equation_audit", source_file="results/equation_functional_form_audit/equation_registry.csv", source_table="QC equation/variable registry", notes="Labels derive from audited definitions; proxy quantities retain proxy units/status.")
    export_table(defs[keep], "supplement/variable_axis_label_map.csv", figure_id="Supplement", panel_id="axis_labels", evidence_level="diagnostic", model_layer="equation_audit", source_file="results/equation_functional_form_audit/variable_definitions.csv; parameter_definitions.csv", source_table="QC definitions", notes="Consistent external plotting labels.")
    export_table(defs[["symbol", "units", "definition"]], "supplement/units_map.csv", figure_id="Supplement", panel_id="units", evidence_level="diagnostic", model_layer="equation_audit", source_file="results/equation_functional_form_audit/variable_definitions.csv; parameter_definitions.csv", source_table="QC units", notes="Dimensionless/proxy status preserved.")


def build_citation_and_panel_maps() -> None:
    captions = {
        1: "Fast firing is attributed to nucleation-limited onset; two-step sintering to PR-prepared bounded closed-pore accommodation memory.",
        2: "The 1,903 exact-promoted cases contain 485 fast-only, 119 two-step-only, 73 both-pass, and 1,226 neither cases.",
        3: "The exact fast-firing Q_nuc window is finite and nucleation-facile ablation removes the effect while PR-off may preserve it.",
        4: "Candidate 693168 shows a conditional two-step property window; it is Tier B and not validation.",
        5: "The fine T2 map contains density-exhaustion, success, and grain-growth regimes.",
        6: "Six exact Tier-B candidates demonstrate a family of conditional mechanisms and sensitivities.",
        7: "Surrogate screening predicts 19,880 both-pass rows versus 73 exact both-pass cases and is promotion-only evidence.",
        8: "Interrupted grain, pore, connectivity, and accommodation measurements provide falsification targets.",
    }
    primary_paths = {
        1: "figure_01_mechanism_chains/mechanism_nodes.csv",
        2: "figure_02_exact_property_phase_map/exact_property_phase_map.csv",
        3: "figure_03_fast_firing_property_window/fast_firing_matched_density_curves.csv",
        4: "figure_04_two_step_property_window/two_step_matched_density_curves_693168.csv",
        5: "figure_05_chen_window_boundaries/candidate_693168_fine_T2_classification.csv",
        6: "figure_06_six_tierB_family/six_TierB_candidate_summary.csv",
        7: "figure_07_surrogate_vs_exact/surrogate_vs_exact_counts.csv",
        8: "figure_08_experimental_falsification/falsification_targets.csv",
    }
    rows = []
    for i in range(1, 9):
        rows.append({"figure_id": f"Figure {i}", "figure_title": captions[i].split(".")[0], "source_document": "docs/FINAL_MECHANISM_SYNTHESIS_CAPTIONS.md", "caption_text": captions[i], "cited_in_section": f"Final synthesis Figure {i}", "required_panels": ";".join(sorted({str(d["panel_id"]) for d in DATASETS if d["figure_id"] == f"Figure {i}"})), "source_data_file": f"results/figure_source_data_package/{primary_paths[i]}", "data_status": "available", "notes": "Controlling exact-only synthesis reference."})
    pub = read_csv("results/publication_style_sintering_figures_693168/source_tables/publication_style_figure_inventory.csv")
    broad = {
        "fast_firing": "fast_firing/clean_fast_firing_T_rho_G_time.csv",
        "two_step": "candidate_693168/dense_time_histories.csv",
        "chen_map": "chen_maps/chen_map_T1_T2_classification.csv",
        "supplement": "candidate_693168/dense_time_histories.csv",
        "main": "candidate_693168/dense_time_histories.csv",
    }
    for row in pub.itertuples():
        number = int(str(row.figure_id).lstrip("F"))
        if number <= 5:
            source_path = broad["fast_firing"]
        elif 12 <= number <= 16:
            source_path = broad["chen_map"]
        else:
            source_path = broad.get(row.category, broad["main"])
        rows.append({"figure_id": f"Candidate-{row.figure_id}", "figure_title": row.title, "source_document": "docs/PUBLICATION_STYLE_FIGURE_CAPTIONS_693168.md", "caption_text": f"{row.title}. Candidate 693168 is conditional Tier B, not validation.", "cited_in_section": row.category, "required_panels": row.source_tables, "source_data_file": f"results/figure_source_data_package/{source_path}", "data_status": "available", "notes": row.candidate_status})
    citation = pd.DataFrame(rows)
    citation.to_csv(OUT / "figure_citation_map.csv", index=False)
    panel = pd.DataFrame(DATASETS)[["figure_id", "panel_id", "relative_path", "evidence_level", "model_layer", "source_file", "source_table", "notes"]]
    panel.to_csv(OUT / "figure_panel_map.csv", index=False)


def build_dictionaries_and_docs() -> None:
    manifest = pd.DataFrame(DATASETS)
    manifest.to_csv(OUT / "figure_data_manifest.csv", index=False)
    missing = pd.DataFrame(MISSING, columns=["figure_id", "panel_id", "requirement", "status", "reason", "handling"])
    if missing.empty:
        missing = pd.DataFrame([{"figure_id": "package", "panel_id": "all", "requirement": "No unresolved file-level source requirement", "status": "complete_with_channel_NaNs", "reason": "Unavailable diagnostic channels are separately enumerated.", "handling": "No values invented."}])
    missing.to_csv(OUT / "missing_figure_data_requirements.csv", index=False)
    archive_manifest = pd.DataFrame([
        {"input_file": str(ARCHIVE.relative_to(ROOT)), "archive_member": ARCHIVE_FAST_MEMBER, "access_mode": "direct_archive_read", "extracted_path": "", "purpose": "Frozen E0021/E0142 fast-firing dense histories and causal ablations", "regenerated": False, "source_commit": COMMIT},
    ])
    archive_dir = OUT / "archive_manifest"
    archive_dir.mkdir(parents=True, exist_ok=True)
    archive_manifest.to_csv(archive_dir / "selective_extraction_manifest.csv", index=False)
    archive_manifest.to_csv(archive_dir / "raw_outputs_manifest.csv", index=False)

    all_rows = []
    for item in DATASETS:
        path = ROOT / str(item["relative_path"])
        try:
            frame = pd.read_csv(path, nrows=50, low_memory=False)
        except Exception:
            continue
        for column in frame.columns:
            all_rows.append({"column": column, "human_readable_label": label_for(column), "units": unit_for(column), "definition": label_for(column), "proxy_or_status": "dimensionless/proxy" if "proxy" in unit_for(column) else "measured/model-output/metadata", "appears_in": str(item["relative_path"])})
    dictionary = pd.DataFrame(all_rows).groupby(["column", "human_readable_label", "units", "definition", "proxy_or_status"], as_index=False).appears_in.agg(lambda values: ";".join(sorted(set(values))))
    dictionary.to_csv(OUT / "all_columns_dictionary.csv", index=False)
    recommendations = pd.DataFrame([
        ("Figure 1", "nodes", "node_id", "mechanism_group", "causal_status", "mechanism_group", "categorical", "Mechanism node", "Mechanism group", "schematic", "Use external layout; no aesthetics encoded."),
        ("Figure 2", "phase_map", "Theta_nuc", "S_closed_growth", "exact_behavior_class", "design_stage", "dimensionless", "Nucleation dominance", "Closed-shrinkage/growth selectivity", "scatter", "Use exact classes only."),
        ("Figure 3", "matched_density", "rho", "R_fast", "material_id", "ablation", "dimensionless", "Relative density", "G_ref/G_fast", "line", "Shade jointly attained threshold intervals."),
        ("Figure 4", "matched_density", "rho", "R_TS", "path_pair_id", "density_window_label", "dimensionless", "Relative density", "G_highT/G_two-step", "line", "Highlight 0.95-0.98 only where attained."),
        ("Figure 5", "fine_classification", "G1_nm", "T2_C", "classification", "candidate_id", "mixed", "Prepared grain size [nm]", "T2 [degC]", "filled_categorical_map", "Retain all failure categories."),
        ("Figure 6", "candidate_summary", "closed_fraction_at_switch", "median_reduction", "candidate_id", "tier_status", "dimensionless", "Closed fraction at switch", "Median reduction", "scatter", "Conditional Tier-B family."),
        ("Figure 7", "confusion", "surrogate_predicted_class", "exact_class", "count", "false_positive_reason", "count", "Surrogate class", "Exact class", "heatmap", "Surrogate is screening only."),
        ("Figure 8", "targets", "measurement", "priority", "mechanism_constrained", "figure_panel_group", "categorical", "Measurement", "Priority", "bar", "Planning schematic, not evidence."),
    ], columns=["figure_id", "panel_id", "x_column", "y_column", "hue_column", "style_column", "units", "suggested_axis_label_x", "suggested_axis_label_y", "recommended_plot_type", "notes"])
    recommendations["suggested_scale"] = recommendations.recommended_plot_type.map({"line": "linear", "scatter": "linear", "heatmap": "categorical", "filled_categorical_map": "categorical", "bar": "categorical", "schematic": "categorical"})
    recommendations.to_csv(OUT / "plotting_recommendations.csv", index=False)

    readme = f"""# Figure source data package

Built from frozen evidence at `{BRANCH}@{COMMIT}`. No simulation, search, parameter
tuning, classification change, or model-physics import is performed by the builder.

The package contains {len(DATASETS)} figure/source tables plus citation, panel,
provenance, dictionary, missing-channel, and QC metadata. Final synthesis claims use
only `evidence_level=exact`. Surrogate rows are explicitly labeled
`surrogate_screen` and are promotion evidence only.

Candidate 693168 remains conditional Tier B, not validation. Its bounded
closed-pore accommodation trajectory is the primary calibration target. Large
high-temperature/two-step grain-size separation is not inherently unphysical, but
the magnitude remains uncalibrated. No closed-pore Poisson/gas-transport law and no
hidden closed-pore Lambda/K law are claimed.

## Loading in Python

```python
import pandas as pd
df = pd.read_csv("results/figure_source_data_package/figure_02_exact_property_phase_map/exact_property_phase_map.csv")
dense = pd.read_csv("results/figure_source_data_package/candidate_693168/dense_time_histories.csv")
```

Use `all_columns_dictionary.csv` for units/labels and `plotting_recommendations.csv`
for panel mappings. No journal aesthetics are hard-coded.
"""
    (ROOT / "docs" / "FIGURE_SOURCE_DATA_PACKAGE_README.md").write_text(readme)
    dictionary_text = "# Figure source data dictionary\n\nThe machine-readable dictionary is `results/figure_source_data_package/all_columns_dictionary.csv`. Every exported scientific table contains provenance, evidence level, model layer, and references to that units/label dictionary. Dimensionless proxies are labeled explicitly; unavailable channels remain NaN.\n"
    (ROOT / "docs" / "FIGURE_SOURCE_DATA_DICTIONARY.md").write_text(dictionary_text)
    (OUT / "FIGURE_SOURCE_DATA_DICTIONARY.md").write_text(dictionary_text)
    (ROOT / "docs" / "FIGURE_SOURCE_DATA_PROVENANCE.md").write_text(f"# Figure source data provenance\n\nSource branch/commit: `{BRANCH}@{COMMIT}`.\n\nAll inputs are existing compact exact, diagnostic, or screening tables. The sole archived input is `{ARCHIVE_FAST_MEMBER}`, read directly from `{ARCHIVE.relative_to(ROOT)}`. No archive directory was restored. No deterministic model regeneration was required. See `archive_manifest/selective_extraction_manifest.csv` and `figure_data_manifest.csv`.\n")
    (ROOT / "docs" / "FIGURE_SOURCE_DATA_MISSING_ITEMS.md").write_text("# Figure source data missing items\n\nNo required source file is missing. Candidate-693168 histories do not expose TJ-specific Lambda/K, P_comp_TJ, X_J, a separately named migration factor, residual stress, trapped-gas pressure, or absolute energy. Requested columns are retained as NaN where applicable and listed in `candidate_693168/missing_candidate_693168_channels.csv`. Closed accommodation is an implemented bounded proxy—not a derived Poisson or gas-transport law—and no missing quantity was invented.\n")


def build_archive() -> Path:
    archive_path = ROOT / "results" / f"figure_source_data_package_{date.today():%Y%m%d}.zip"
    if archive_path.exists():
        archive_path.unlink()
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(OUT.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(ROOT / "results"))
        for name in ["FIGURE_SOURCE_DATA_PACKAGE_README.md", "FIGURE_SOURCE_DATA_DICTIONARY.md", "FIGURE_SOURCE_DATA_PROVENANCE.md", "FIGURE_SOURCE_DATA_MISSING_ITEMS.md"]:
            path = ROOT / "docs" / name
            archive.write(path, Path("figure_source_data_package") / "docs" / name)
    return archive_path


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    for directory in [
        "figure_01_mechanism_chains", "figure_02_exact_property_phase_map",
        "figure_03_fast_firing_property_window", "figure_04_two_step_property_window",
        "figure_05_chen_window_boundaries", "figure_06_six_tierB_family",
        "figure_07_surrogate_vs_exact", "figure_08_experimental_falsification",
        "candidate_693168", "fast_firing", "chen_maps", "supplement", "qc",
        "archive_manifest",
    ]:
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    build_mechanism_chains()
    score, exact = build_exact_property_phase_map()
    clean_fast = build_fast_firing()
    hist = build_two_step()
    build_chen_and_candidate(hist, clean_fast)
    build_chen_maps()
    build_surrogate_warning(score, exact)
    build_falsification()
    build_equation_maps()
    build_citation_and_panel_maps()
    build_dictionaries_and_docs()
    archive = build_archive()
    print(f"datasets={len(DATASETS)}")
    print(f"archive={archive.relative_to(ROOT)}")
    print(f"archive_bytes={archive.stat().st_size}")
    print(f"source_commit={COMMIT}")


if __name__ == "__main__":
    main()
