#!/usr/bin/env python3
"""Synthesize existing exact attribution results without rerunning model physics."""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd


ATTR = Path("results/relative_material_property_window_attribution/source_tables")
REFRAME = Path("results/reframe_tierB_experimental_plausibility")
AUDIT = Path("results/audit_candidate_693168_closed_accommodation/final_tables")
OUT = Path("results/final_mechanism_synthesis_and_property_windows")
SRC = OUT / "source_tables"

Q_PARAMS = ["Q_nuc_delta_kJ", "Q_exchange_delta_kJ", "Q_transport_delta_kJ",
            "Q_growth_delta_kJ", "Q_PR_delta_kJ", "Q_closed_delta_kJ"]
K_PARAMS = ["k_nuc_factor", "k_exchange_factor", "k_transport_factor",
            "k_growth_factor", "k_PR_factor", "k_closed_factor"]


def as_bool(series: pd.Series) -> np.ndarray:
    return series.astype("boolean").fillna(False).to_numpy(bool)


def perturbed_parameter(row: pd.Series) -> tuple[str, float]:
    changed = [(p, float(row[p])) for p in Q_PARAMS if abs(float(row[p])) > 1e-12]
    changed += [(p, float(row[p])) for p in K_PARAMS if abs(float(row[p]) - 1.0) > 1e-12]
    return changed[0] if len(changed) == 1 else ("multiple_or_none", np.nan)


def oat_tables(merged: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    oat = merged[merged.design_stage.eq("OAT")].copy()
    pv = oat.apply(perturbed_parameter, axis=1, result_type="expand")
    oat["perturbed_parameter"] = pv[0]
    oat["perturbed_value"] = pv[1]
    fast = oat[["property_id", "perturbed_parameter", "perturbed_value",
                "R_fast_exact", "span_fast_1p5_exact", "exact_fast_attained",
                "fast_firing_pass_exact", "nucleation_facile_pass_exact", "PR_off_pass_exact"]].copy()
    fast["evidence_level"] = "exact OAT"
    fast["boundary_role"] = "interior"
    fast.loc[(fast.perturbed_parameter == "Q_nuc_delta_kJ") & fast.perturbed_value.isin([-25, 75]), "boundary_role"] = "adjacent_failure"
    fast.loc[(fast.perturbed_parameter == "Q_nuc_delta_kJ") & fast.perturbed_value.isin([0, 50]), "boundary_role"] = "survival_edge"
    base = merged[merged.property_id.eq("BASE")].iloc[0]
    fast_base = {"property_id": "BASE", "perturbed_parameter": "Q_nuc_delta_kJ", "perturbed_value": 0.0,
                 "R_fast_exact": base.R_fast_exact, "span_fast_1p5_exact": base.span_fast_1p5_exact,
                 "exact_fast_attained": base.exact_fast_attained, "fast_firing_pass_exact": base.fast_firing_pass_exact,
                 "nucleation_facile_pass_exact": base.nucleation_facile_pass_exact, "PR_off_pass_exact": base.PR_off_pass_exact,
                 "evidence_level": "exact base/OAT center", "boundary_role": "survival_edge"}
    fast = pd.concat([fast, pd.DataFrame([fast_base])], ignore_index=True).drop_duplicates(["property_id", "perturbed_parameter"])

    two = oat[["property_id", "perturbed_parameter", "perturbed_value",
               "reduction_TS_exact", "span_TS_20_exact", "Chen_window_width_C_exact",
               "exact_two_attained", "lower_boundary_present_exact", "upper_boundary_present_exact",
               "two_step_pass_exact"]].copy()
    two["evidence_level"] = "exact OAT"
    two["boundary_role"] = "interior"
    two.loc[(two.perturbed_parameter == "Q_closed_delta_kJ") & (two.perturbed_value == -50), "boundary_role"] = "lower_boundary_lost"
    two.loc[(two.perturbed_parameter == "Q_closed_delta_kJ") & two.perturbed_value.isin([-25, 100]), "boundary_role"] = "survival_edge"
    two.loc[(two.perturbed_parameter == "k_PR_factor") & (two.perturbed_value == .1), "boundary_role"] = "below_threshold_failure"
    two.loc[(two.perturbed_parameter == "k_PR_factor") & (two.perturbed_value == .3), "boundary_role"] = "threshold_pass"
    two.loc[(two.perturbed_parameter == "k_growth_factor") & (two.perturbed_value == .03), "boundary_role"] = "upper_boundary_lost"
    two.loc[(two.perturbed_parameter == "k_growth_factor") & (two.perturbed_value == .1), "boundary_role"] = "threshold_pass"
    two_base = {"property_id": "BASE", "perturbed_parameter": "Q_closed_delta_kJ", "perturbed_value": 0.0,
                "reduction_TS_exact": base.reduction_TS_exact, "span_TS_20_exact": base.span_TS_20_exact,
                "Chen_window_width_C_exact": base.Chen_window_width_C_exact, "exact_two_attained": base.exact_two_attained,
                "lower_boundary_present_exact": base.lower_boundary_present_exact,
                "upper_boundary_present_exact": base.upper_boundary_present_exact,
                "two_step_pass_exact": base.two_step_pass_exact, "evidence_level": "exact base/OAT center",
                "boundary_role": "interior"}
    two = pd.concat([two, pd.DataFrame([two_base])], ignore_index=True).drop_duplicates(["property_id", "perturbed_parameter"])
    return fast, two


def family_table() -> pd.DataFrame:
    family = pd.read_csv(REFRAME / "tierB_candidate_reinterpretation.csv")
    compact = family[["candidate_id", "G0_nm", "G1_nm", "first_step_growth_fraction",
                      "median_reduction", "window_width_C", "T2_success_range",
                      "closed_fraction_at_switch", "closed_fraction_at_target",
                      "closed_shrinkage_share", "ablations_that_destroy_result",
                      "robustness_count", "interpretation", "artifact_flag"]].copy()
    compact["tier_status"] = "conditional Tier B"
    compact["destructive_ablation_count"] = compact.ablations_that_destroy_result.fillna("").str.split(";").map(len)
    compact["representative_role"] = np.where(compact.candidate_id.eq(693168),
                                                "best strong-separation conditional comparator",
                                                "Tier-B family comparator")
    compact["validation_status"] = "not validation"
    compact["evidence_level"] = "exact candidate base results; material robustness reduced transfer only"
    return compact


def falsification_targets() -> pd.DataFrame:
    rows = [
        ("Matched-density G_mean, G50, G90", "fast and two-step trajectory separation", "fast path smaller after onset; two-step distribution remains finer"),
        ("Densification-onset time during ramps", "nucleation-limited waiting", "nucleation-facile material weakens heating-rate separation"),
        ("Exchange and transport relaxation", "post-nucleation attainment", "completion remains fast after the onset threshold"),
        ("3D open/closed pore fraction", "PR-prepared closed store", "first step increases a persistent closed inventory"),
        ("Pore D50, D90, and large-pore tail", "PR/surface redistribution", "schedule-dependent distribution at matched density"),
        ("Connected fine-pore fraction/percolation", "removable-pore topology", "connectivity changes across first-step preparation"),
        ("Trapped-gas pressure or accommodation proxy", "finite closed accommodation", "low-T2 exhaustion accompanies accommodation depletion"),
        ("Grain-growth mobility versus T2", "upper Chen boundary", "migration rises across the upper success boundary"),
        ("Interrupted first/second-step tomography", "memory persistence", "prepared topology persists long enough to affect G(rho)"),
    ]
    return pd.DataFrame(rows, columns=["measurement", "mechanism_constrained", "falsifiable_expected_signature"])


def main() -> None:
    SRC.mkdir(parents=True, exist_ok=True)
    score = pd.read_csv(ATTR / "material_property_window_scorecard.csv", low_memory=False)
    exact = pd.read_csv(ATTR / "material_property_window_exact_promotions.csv")
    merged = exact.merge(score, on="property_id", how="left", suffixes=("", "_screen"))
    counts = exact.classification_exact.value_counts().reindex(["fast_only", "two_step_only", "both_pass", "neither"], fill_value=0)
    count_table = pd.DataFrame({"classification": counts.index, "exact_count": counts.values})
    count_table["evidence_level"] = "exact promoted union"
    count_table.to_csv(SRC / "exact_behavior_classification_counts.csv", index=False)

    fast_oat, two_oat = oat_tables(merged)
    fast_oat.to_csv(SRC / "fast_firing_OAT_window.csv", index=False)
    two_oat.to_csv(SRC / "two_step_OAT_window.csv", index=False)

    thresholds = pd.read_csv(ATTR / "dimensionless_thresholds.csv")
    thresholds["final_claim_eligible"] = True
    thresholds["scope_note"] = "exact promoted subset; coverage-limited, not universal bounds"
    thresholds.to_csv(SRC / "dimensionless_group_thresholds_final.csv", index=False)

    family = family_table()
    family.to_csv(SRC / "six_TierB_family_mechanism_summary.csv", index=False)

    surrogate_counts = score.classification_screen.value_counts()
    comparison = pd.DataFrame([
        {"classification": c, "surrogate_count": int(surrogate_counts.get(c, 0)),
         "exact_count": int(counts.get(c, 0)), "evidence_note": "surrogate is screening only; exact is final"}
        for c in ["fast_only", "two_step_only", "both_pass", "neither"]
    ])
    comparison["surrogate_to_exact_ratio"] = comparison.surrogate_count / comparison.exact_count.replace(0, np.nan)
    comparison.to_csv(SRC / "surrogate_vs_exact_comparison.csv", index=False)

    targets = falsification_targets()
    targets.to_csv(SRC / "experimental_falsification_targets.csv", index=False)

    summary_rows = [
        ("campaign", "perturbations_screened", len(score), "rows", "surrogate screen", "not final evidence"),
        ("campaign", "exact_fast_promotions", int(exact.fast_firing_pass_exact.notna().sum()), "rows", "exact", "final"),
        ("campaign", "exact_two_step_promotions", int(exact.two_step_pass_exact.notna().sum()), "rows", "exact", "final"),
        ("campaign", "unique_exact_cases", len(exact), "rows", "exact union", "final"),
        ("classification", "fast_only", int(counts.fast_only), "cases", "exact", "final"),
        ("classification", "two_step_only", int(counts.two_step_only), "cases", "exact", "final"),
        ("classification", "both_pass", int(counts.both_pass), "cases", "exact", "final"),
        ("classification", "neither", int(counts.neither), "cases", "exact", "final"),
        ("screen_warning", "surrogate_both", int(surrogate_counts.get("both_pass", 0)), "cases", "surrogate", "screen overprediction"),
        ("screen_warning", "exact_both", int(counts.both_pass), "cases", "exact", "final"),
        ("fast_OAT", "Delta_Q_nuc_survival", "0 to +50", "kJ/mol", "exact OAT", "-25 and +75 fail"),
        ("two_step_OAT", "Delta_Q_closed_survival", "-25 to +100", "kJ/mol", "exact OAT", "lower boundary absent at -50"),
        ("two_step_OAT", "Delta_Q_growth_tested_viable", "-100 to +100", "kJ/mol", "exact OAT", "limit not found"),
        ("two_step_OAT", "PR_prefactor_threshold", 0.3, "x base", "exact OAT", "0.1x fails"),
        ("two_step_OAT", "growth_prefactor_threshold", 0.1, "x base", "exact OAT", "0.03x loses upper boundary"),
        ("family", "best_TierB_candidate", 693168, "candidate_id", "exact base", "conditional Tier B; not validation"),
        ("family", "TierB_family_size", len(family), "candidates", "exact base", "reduced material transfer only"),
        ("guardrail", "model_physics_changed", False, "boolean", "synthesis metadata", "no"),
        ("guardrail", "topology_parameters_changed", False, "boolean", "synthesis metadata", "no"),
        ("guardrail", "validation_claim", False, "boolean", "synthesis metadata", "conditional attribution only"),
    ]
    summary = pd.DataFrame(summary_rows, columns=["category", "metric", "value", "units", "evidence_level", "interpretation"])
    summary.to_csv(SRC / "final_property_window_summary.csv", index=False)

    state = {
        "status": "complete", "new_search_run": False, "model_physics_changed": False,
        "topology_parameters_changed": False, "material_parameters_retuned": False,
        "screened_rows": int(len(score)), "unique_exact_cases": int(len(exact)),
        "exact_fast_only": int(counts.fast_only), "exact_two_step_only": int(counts.two_step_only),
        "exact_both_pass": int(counts.both_pass), "exact_neither": int(counts.neither),
        "surrogate_both": int(surrogate_counts.get("both_pass", 0)),
        "candidate_693168_status": "conditional Tier B; not validation",
    }
    (OUT / "run_state.json").write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state, indent=2))


if __name__ == "__main__":
    main()
