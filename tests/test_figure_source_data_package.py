from pathlib import Path
import json

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/figure_source_data_package"


def read(relative):
    return pd.read_csv(OUT / relative, low_memory=False)


def test_all_cited_figures_have_nonempty_source_data():
    citations = read("figure_citation_map.csv")
    assert len(citations) >= 33
    for source in citations.source_data_file:
        path = ROOT / source
        assert path.exists() and path.stat().st_size > 0


def test_exact_counts_and_surrogate_warning_are_reproduced():
    counts = read("figure_02_exact_property_phase_map/exact_behavior_classification_counts.csv")
    observed = dict(zip(counts.behavior_class_or_total, counts["count"]))
    assert observed == {
        "screened_rows": 50655,
        "unique_exact_cases": 1903,
        "fast_only": 485,
        "two_step_only": 119,
        "both_pass": 73,
        "neither": 1226,
        "surrogate_both": 19880,
    }
    warning = read("figure_07_surrogate_vs_exact/surrogate_vs_exact_counts.csv")
    values = dict(zip(warning.metric, warning.value))
    assert values["surrogate_both"] == 19880
    assert values["exact_both"] == 73


def test_candidate_693168_package_has_dense_curves_maps_and_ablations():
    required = [
        "candidate_693168/dense_time_histories.csv",
        "candidate_693168/matched_density_ratio_curves.csv",
        "candidate_693168/fine_T2_scan.csv",
        "candidate_693168/Chen_boundaries_fine.csv",
        "candidate_693168/ablation_summary.csv",
    ]
    assert all((OUT / name).exists() and (OUT / name).stat().st_size > 0 for name in required)
    curves = read(required[1])
    success = curves[curves.path_pair_id.str.contains("two_step_success")]
    assert success.rho.min() <= 0.95
    assert success.rho.max() >= 0.98


def test_fast_firing_package_contains_full_PR_off_and_nucleation_facile():
    histories = read("fast_firing/clean_fast_firing_T_rho_G_time.csv")
    assert {"E0021", "E0142"} <= set(histories.material_id)
    assert {"full", "PR-off", "nucleation-facile"} <= set(histories.ablation)


def test_complete_Chen_maps_retain_all_three_outcome_regions():
    required_classes = {
        "DENSIFICATION_EXHAUSTION_FAILURE",
        "SUCCESS",
        "GRAIN_GROWTH_FAILURE",
    }
    for filename in [
        "chen_maps/chen_map_T1_T2_classification.csv",
        "chen_maps/chen_map_G1_T2_classification.csv",
        "chen_maps/chen_map_switch_density_T2_classification.csv",
    ]:
        frame = read(filename)
        assert required_classes <= set(frame.classification)
        assert (~frame.success_flag.astype(bool)).any()


def test_package_QC_passes_and_guardrails_are_explicit():
    checks = read("qc/figure_source_data_audit.csv")
    assert checks.passed.astype(bool).all()
    summary = json.loads((OUT / "qc/figure_source_data_package_summary.json").read_text())
    assert summary["passed"]
    text = " ".join(
        [
            (ROOT / "docs/FIGURE_SOURCE_DATA_PACKAGE_README.md").read_text(),
            (ROOT / "docs/FIGURE_SOURCE_DATA_MISSING_ITEMS.md").read_text(),
        ]
    ).lower()
    assert "conditional tier b" in text
    assert "not validation" in text
    assert "bounded proxy" in text
    assert "hidden closed-pore lambda/k law" in text
