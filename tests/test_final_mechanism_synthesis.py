"""Guardrails for the no-new-search final mechanism synthesis."""
from pathlib import Path
import subprocess

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/final_mechanism_synthesis_and_property_windows"
SRC = OUT / "source_tables"


def truth(value) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def test_final_counts_come_from_exact_promoted_union():
    exact = pd.read_csv(ROOT / "results/relative_material_property_window_attribution/source_tables/material_property_window_exact_promotions.csv")
    expected = exact.classification_exact.value_counts()
    final = pd.read_csv(SRC / "exact_behavior_classification_counts.csv").set_index("classification")
    for classification in ("fast_only", "two_step_only", "both_pass", "neither"):
        assert int(final.loc[classification, "exact_count"]) == int(expected[classification])
    assert final.evidence_level.str.contains("exact", case=False).all()


def test_surrogate_rows_are_not_final_both_pass_evidence():
    comparison = pd.read_csv(SRC / "surrogate_vs_exact_comparison.csv").set_index("classification")
    assert int(comparison.loc["both_pass", "surrogate_count"]) == 19_880
    assert int(comparison.loc["both_pass", "exact_count"]) == 73
    assert comparison.loc["both_pass", "surrogate_to_exact_ratio"] > 270
    summary = pd.read_csv(SRC / "final_property_window_summary.csv")
    exact_both = summary[summary.metric.eq("both_pass")].iloc[0]
    assert exact_both.evidence_level == "exact"


def test_fast_firing_OAT_window_has_survival_and_failure_edges():
    d = pd.read_csv(SRC / "fast_firing_OAT_window.csv")
    q = d[d.perturbed_parameter.eq("Q_nuc_delta_kJ")].set_index("perturbed_value")
    assert truth(q.loc[0.0, "fast_firing_pass_exact"])
    assert truth(q.loc[50.0, "fast_firing_pass_exact"])
    assert not truth(q.loc[-25.0, "fast_firing_pass_exact"])
    assert not truth(q.loc[75.0, "fast_firing_pass_exact"])
    assert q.loc[-25.0, "boundary_role"] == "adjacent_failure"
    assert q.loc[50.0, "boundary_role"] == "survival_edge"


def test_two_step_OAT_window_has_closed_survival_and_lower_bound_loss():
    d = pd.read_csv(SRC / "two_step_OAT_window.csv")
    q = d[d.perturbed_parameter.eq("Q_closed_delta_kJ")].set_index("perturbed_value")
    assert truth(q.loc[-25.0, "two_step_pass_exact"])
    assert truth(q.loc[100.0, "two_step_pass_exact"])
    assert not truth(q.loc[-50.0, "two_step_pass_exact"])
    assert not truth(q.loc[-50.0, "lower_boundary_present_exact"])
    assert q.loc[-50.0, "boundary_role"] == "lower_boundary_lost"
    pr = d[d.perturbed_parameter.eq("k_PR_factor")].set_index("perturbed_value")
    assert not truth(pr.loc[0.1, "two_step_pass_exact"])
    assert truth(pr.loc[0.3, "two_step_pass_exact"])


def test_693168_and_all_six_candidates_remain_conditional_tierB():
    family = pd.read_csv(SRC / "six_TierB_family_mechanism_summary.csv")
    assert set(family.candidate_id.astype(int)) == {693168, 822940, 581668, 295003, 366094, 85161}
    assert family.tier_status.eq("conditional Tier B").all()
    assert family.validation_status.eq("not validation").all()
    assert not family.tier_status.str.contains("Tier A", case=False).any()
    best = family[family.candidate_id.eq(693168)].iloc[0]
    assert "best" in best.representative_role


def test_all_required_synthesis_figures_exist_in_both_formats():
    stems = [
        "mechanism_chain_fast_vs_twostep", "relative_property_phase_map_exact",
        "fast_firing_property_window", "two_step_property_window",
        "chen_window_mechanism_boundaries", "six_TierB_family_property_summary",
        "surrogate_vs_exact_warning", "experimental_falsification_targets",
    ]
    for stem in stems:
        for suffix in (".pdf", ".png"):
            p = OUT / "figures" / (stem + suffix)
            assert p.exists() and p.stat().st_size > 10_000
    inventory = pd.read_csv(SRC / "final_figure_inventory.csv")
    assert set(inventory.figure) == set(stems)


def test_reports_retain_no_validation_and_surrogate_warning_language():
    main = (ROOT / "docs/FINAL_FAST_FIRING_AND_TWO_STEP_MECHANISM_SYNTHESIS.md").read_text().lower()
    captions = (ROOT / "docs/FINAL_MECHANISM_SYNTHESIS_CAPTIONS.md").read_text().lower()
    assert "not validation" in main
    assert "not validation" in captions
    assert "surrogate" in main and "not final evidence" in main
    assert "surrogate screens overpredict exact overlap" in captions
    assert "conditional tier b" in main


def test_synthesis_changed_no_model_or_topology_physics_files():
    allowed = {
        "final_mechanism_synthesis_and_property_windows.py",
        "final_mechanism_synthesis_plots.py",
        "tests/test_final_mechanism_synthesis.py",
        "equation_functional_form_audit.py",
        "tests/test_equation_audit_outputs.py",
        "conftest.py",
    }
    changed = subprocess.check_output(
        ["git", "diff", "--name-only", "fd3be2d", "--", "*.py"], cwd=ROOT, text=True
    ).splitlines()
    untracked = subprocess.check_output(
        ["git", "ls-files", "--others", "--exclude-standard", "--", "*.py"], cwd=ROOT, text=True
    ).splitlines()
    assert set(changed + untracked).issubset(allowed)
    state = pd.read_json(OUT / "run_state.json", typ="series")
    assert not bool(state.model_physics_changed)
    assert not bool(state.topology_parameters_changed)
    assert not bool(state.material_parameters_retuned)
    assert not bool(state.new_search_run)
