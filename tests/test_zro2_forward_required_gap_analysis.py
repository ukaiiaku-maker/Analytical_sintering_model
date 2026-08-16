from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

ROOT=Path(__file__).parents[1];SRC=ROOT/"results/zro2_forward_950C_sensitivity_chen_failure_audit";OUT=ROOT/"results/zro2_forward_required_chen_physics_gap_analysis"

def test_source_commit_tables_are_present_and_complete():
    assert (SRC/"run_state.json").exists();x=pd.read_csv(SRC/"chen_failure_decomposition_full.csv");assert len(x)==195
    assert subprocess.check_output(["git","merge-base","--is-ancestor","acf6692","HEAD"],cwd=ROOT).decode()==""

def test_gap_only_exists_with_both_boundaries():
    b=pd.read_csv(OUT/"boundary_gap_by_state_ranked.csv");both=b.T_lower_density_C.notna()&b.T_upper_growth_C.notna();assert b.loc[~both,"gap_C"].isna().all();assert np.allclose(b.loc[both,"gap_C"],b.loc[both,"T_upper_growth_C"]-b.loc[both,"T_lower_density_C"])

def test_negative_gaps_are_not_windows_and_require_positive_shift():
    b=pd.read_csv(OUT/"boundary_gap_by_state_ranked.csv");n=b[b.gap_C<0];assert len(n)>0 and not n.finite_window_present.any();assert n.required_shift_C.gt(0).all();assert np.allclose(n.required_shift_C,-n.gap_C+25)

def test_strict_result_and_relaxed_diagnostic_tags_are_preserved():
    t=pd.read_csv(OUT/"threshold_relaxation_transition_table.csv");strict=t.query("density_target==.976 and grain_threshold_um==.29");assert len(strict)==1 and strict.success_count.iloc[0]==0 and strict.finite_window_count.iloc[0]==0 and not strict.diagnostic_only.iloc[0];assert t.drop(strict.index).diagnostic_only.all()

def test_OAT_cannot_overwrite_strict_outcome():
    o=pd.read_csv(OUT/"boundary_shift_required_by_parameter.csv");assert o.not_a_tuning_recommendation.all();tri=pd.read_csv(SRC/"candidate_triage_950C_sensitivity.csv");assert not tri.strict_Chen_window.any() and not tri.pathway_consistent.any()

def test_common_state_analysis_has_no_method_specific_initialization():
    p=pd.read_csv(OUT/"common_state_pathway_consistency_summary.csv");assert p.method_specific_initialization_absent.all();assert not p.all_pathway_gates_pass.any();assert not p.boundary_gap_identifiable.any()

def test_no_model_source_files_changed_on_analysis_branch():
    names=subprocess.check_output(["git","diff","--name-only","acf6692"],cwd=ROOT,text=True).splitlines();assert not any(n.startswith("zro2_forward/") for n in names)

def test_reports_explicitly_reject_validation_claim():
    reports=["ZRO2_FORWARD_REQUIRED_CHEN_PHYSICS_GAP_ANALYSIS.md","ZRO2_FORWARD_BOUNDARY_ORDERING_INTERPRETATION.md","ZRO2_FORWARD_COMMON_STATE_SENSITIVITY_INTERPRETATION.md","ZRO2_FORWARD_EXPERIMENTAL_DATA_NEEDED_FOR_NEXT_MODEL.md","ZRO2_FORWARD_DECISION_AFTER_950C_AUDIT.md"]
    for name in reports:
        text=(ROOT/"docs"/name).read_text().lower();assert "not validat" in text or "no validation claim" in text
