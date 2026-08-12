from pathlib import Path
import pandas as pd

ROOT=Path("results/reduced_model_falsification_summary")

def test_required_summary_csvs_exist():
    assert (ROOT/"mechanism_scorecard.csv").is_file();assert (ROOT/"failure_mode_matrix.csv").is_file()

def test_every_major_mechanism_branch_is_present():
    q=pd.read_csv(ROOT/"mechanism_scorecard.csv");text=" ".join(q.mechanism.str.lower())
    required=("aggregate","pore-placement","pore-location","persistent junction","preparation","production pr","observable","heterogeneous","persistent defect","local connected","late-stage")
    assert all(x in text for x in required)

def test_no_mechanism_is_labeled_meaningful_without_ratio_and_span():
    q=pd.read_csv(ROOT/"mechanism_scorecard.csv");m=q[q.trajectory_meaningful==True]
    assert ((m.max_ratio>=1.5)&(m.span_ge_1p5>=.03)).all()

def test_high_density_support_requires_attainment():
    q=pd.read_csv(ROOT/"mechanism_scorecard.csv")
    assert not ((q.rho_ref_max<.95)|(q.rho_fast_max<.95))[q.density_range.str.contains("0.98|0.99",regex=True)].empty
    # No scorecard row claims meaningful high-density response.
    assert not q.trajectory_meaningful.any()

def test_readme_contains_corrected_status():
    text=Path("README.md").read_text();assert "Current scientific status" in text;assert "has **not** been achieved" in text;assert "spatial/network" in text
