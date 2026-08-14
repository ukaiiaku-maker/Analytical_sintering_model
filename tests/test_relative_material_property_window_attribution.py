"""Guardrails for the frozen-topology material-property attribution campaign."""
from pathlib import Path

import numpy as np
import pandas as pd

from mechanism_dimensionless_groups import artifact_reasons, longest_span


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/relative_material_property_window_attribution"
SRC = OUT / "source_tables"


def tables():
    score = pd.read_csv(SRC / "material_property_window_scorecard.csv", low_memory=False)
    exact = pd.read_csv(SRC / "material_property_window_exact_promotions.csv")
    return score, exact.merge(score, on="property_id", how="left", suffixes=("", "_screen"))


def b(series):
    return series.astype("boolean").fillna(False).to_numpy(bool)


def test_primary_scan_freezes_topology_parameters():
    score, exact = tables()
    primary = score[~score.diagnostic_only.astype(bool)]
    assert len(primary) >= 50_000
    assert not primary.topology_parameters_modified.astype(bool).any()
    assert not b(exact.topology_parameters_modified).any()


def test_dimensionless_groups_are_persisted_and_finite_per_scored_behavior():
    _, exact = tables()
    groups = pd.read_csv(SRC / "dimensionless_groups_by_path.csv")
    assert len(groups) == 2 * len(exact)
    assert groups.groupby("property_id").behavior.nunique().eq(2).all()
    fast = groups[groups.behavior.eq("fast_firing")]
    two = groups[groups.behavior.eq("two_step")]
    fast_cols = ["Theta_nuc", "f_nuc", "f_exchange", "f_transport", "I_low_slow", "I_low_PR_slow", "Pi_PR"]
    two_cols = ["S_closed_growth", "A_closed_fraction", "M_PR_closed", "Gamma_mig"]
    assert np.isfinite(fast[fast_cols].to_numpy(float)).all()
    assert np.isfinite(two[two_cols].to_numpy(float)).all()


def test_fast_pass_requires_ratio_span_attainment_and_causal_ablation():
    _, exact = tables()
    passed = exact[b(exact.fast_firing_pass_exact)]
    assert len(passed) > 0
    assert b(passed.exact_fast_attained).all()
    assert (passed.R_fast_exact >= 1.5).all()
    assert (passed.span_fast_1p5_exact >= 0.03 - 1e-12).all()
    assert not b(passed.nucleation_facile_pass_exact).any()
    # PR-off survival is admissible, and occurs in the exact base case.
    base = passed[passed.property_id.eq("BASE")].iloc[0]
    assert bool(base.PR_off_pass_exact)


def test_two_step_pass_requires_finite_attained_high_density_window():
    _, exact = tables()
    passed = exact[b(exact.two_step_pass_exact)]
    assert len(passed) > 0
    assert b(passed.exact_two_attained).all()
    assert (passed.reduction_TS_exact >= 0.20).all()
    assert (passed.span_TS_20_exact >= 0.02 - 1e-12).all()
    assert (passed.Chen_window_width_C_exact >= 25.0).all()
    assert b(passed.lower_boundary_present_exact).all()
    assert b(passed.upper_boundary_present_exact).all()
    # The exact solver, not the reduced promotion screen, is authoritative.
    # Active high-density support is represented by a nonzero closed store.
    assert (passed.closed_fraction_switch_exact > 0).all()


def test_large_attained_reduction_is_not_an_artifact():
    reasons = artifact_reasons(attained=True, bounded=True, interpolation_supported=True)
    assert reasons == []
    _, exact = tables()
    base = exact[exact.property_id.eq("BASE")].iloc[0]
    assert base.reduction_TS_exact > 0.8
    assert pd.isna(base.two_artifact_reasons) or not str(base.two_artifact_reasons).strip()


def test_span_helper_requires_a_continuous_interval():
    rho = np.array([0.95, 0.96, 0.97, 0.98])
    assert np.isclose(longest_span(rho, [1.6, 1.6, 1.2, 1.6], 1.5), 0.01)


def test_all_six_tierB_candidates_and_proxy_evidence_are_visible():
    family = pd.read_csv(SRC / "tierB_family_material_window_comparison.csv")
    assert set(family.candidate_id.astype(int)) == {693168, 822940, 581668, 295003, 366094, 85161}
    assert family.evidence_level.str.contains("not exact family validation", case=False).all()


def test_reports_retain_conditional_nonvalidation_language():
    reports = [
        "MATERIAL_PROPERTY_WINDOW_MECHANISM_ATTRIBUTION.md",
        "FAST_FIRING_MECHANISM_ATTRIBUTION.md",
        "TWO_STEP_MECHANISM_ATTRIBUTION.md",
        "RELATIVE_ACTIVATION_ENERGY_WINDOW.md",
        "TIERB_FAMILY_MATERIAL_PROPERTY_ROBUSTNESS.md",
        "GENERAL_PARTICLE_SYSTEM_TRENDS.md",
        "MATERIAL_PROPERTY_WINDOW_FIGURE_MANIFEST.md",
        "MATERIAL_PROPERTY_WINDOW_CAPTION_DRAFTS.md",
    ]
    for name in reports:
        text = (ROOT / "docs" / name).read_text().lower()
        assert "not validation" in text or "not a validated" in text or "not universal" in text


def test_required_figures_exist_in_vector_and_raster_forms():
    stems = [
        "mechanism_ingredient_summary", "fast_firing_mechanism_attribution",
        "two_step_mechanism_attribution", "relative_barrier_window_fast_firing",
        "relative_barrier_window_two_step", "joint_behavior_phase_map",
        "sensitivity_tornado_fast_firing", "sensitivity_tornado_two_step",
        "closed_accommodation_plausibility_window", "six_TierB_material_window_comparison",
        "common_particle_system_trends",
    ]
    for stem in stems:
        for suffix in (".pdf", ".png"):
            path = OUT / "figures" / (stem + suffix)
            assert path.exists() and path.stat().st_size > 10_000
