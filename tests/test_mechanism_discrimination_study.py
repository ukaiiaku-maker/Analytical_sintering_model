import csv
import math
from pathlib import Path

import mechanism_discrimination_study as study

RESULTS=Path(__file__).parents[1]/"results"/"mechanism_discrimination"


def test_matched_density_and_fast_firing_designs_are_complete():
    matched=list(csv.DictReader((RESULTS/"matched_density_results.csv").open()))
    fast=list(csv.DictReader((RESULTS/"fast_firing_results.csv").open()))
    assert len(matched)==len(study.MODEL_STYLES)*len(study.G0_MATCHED)*len(study.MATCH_DENSITIES)*len(study.HISTORY_RATES)
    assert len(fast)==2*len(study.FAST_G0)*len(study.FAST_RHO0)*len(study.FAST_RATES)
    assert all(r["history_reached"]=="True" for r in matched)


def test_no_failed_fast_target_is_scored():
    rows=list(csv.DictReader((RESULTS/"fast_firing_results.csv").open()))
    for row in rows:
        score=float(row["HR_pct_vs_slow"])
        if not math.isnan(score):assert row["reached_target"]=="True"


def test_summary_matrix_has_required_model_discrimination_rows():
    rows=list(csv.DictReader((RESULTS/"robustness_matrix.csv").open()))
    expected={"finite two-step window","nanoscale two-step robustness","two-step size boundary","matched-density topology identifiability",
              "high-heating-rate advantage","heating-rate robustness across particle size","heating-rate robustness across initial density",
              "grain-growth failure reproduced","densification-exhaustion failure reproduced","unattainable-first-step region reproduced"}
    assert {r["criterion"] for r in rows}==expected
    assert set(rows[0])=={"criterion",*study.MODEL_STYLES}
