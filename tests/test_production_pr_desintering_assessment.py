import csv,inspect
from dataclasses import replace
from pathlib import Path
import numpy as np

import agentic_mechanism_model as discovery
import pr_desintering_memory_model as memory
import production_pr_desintering_assessment as production

ROOT=Path(__file__).parents[1]/"results"/"production_pr_desintering_assessment"


def rows(name):
    with (ROOT/name).open(newline="") as f:return list(csv.DictReader(f))


def test_production_candidates_are_frozen_and_q0_q1_visible():
    c=production.candidates();assert set(c)=={"mech_009","mech_019","mech_009_q0","mech_019_q0"}
    assert all(p.early_memory_mode=="PR_plus_connected_fine_attrition" and p.k_PR_ref_s==2e-4 for p in c.values())
    assert {p.base.q_TJ for p in c.values()}=={0,1}


def test_persisted_production_score_is_joint_positive_and_nonuniversal():
    score=rows("production_joint_scorecard.csv");enabled=[r for r in score if r["candidate_id"]!="disabled_control"]
    assert len(enabled)==4 and all(r["joint_positive"]=="True" and r["universal_fast"]=="False" and int(r["beneficial_fast_count"])>0 for r in enabled)
    disabled=[r for r in score if r["candidate_id"]=="disabled_control"][0];assert disabled["joint_positive"]=="False" and disabled["beneficial_fast_count"]=="0"


def test_complete_practical_windows_require_T2_below_T1_and_both_bounds():
    successful=rows("successful_practical_windows.csv");assert successful
    for r in successful:
        assert r["map_type"]=="practical" and r["boundary_status"]=="COMPLETE_WINDOW"
        assert r["lower_bracketed"]==r["upper_bracketed"]=="True"
        assert float(r["T_last_success_C"])<float(r["T1_C"])
        assert float(r["first_step_growth_fraction"])<=float(r["prep_growth_tolerance"])+1e-12


def test_fast_successes_are_scored_only_when_both_paths_attain_and_failures_persist():
    success=rows("successful_fast_firing_cases.csv");assert success and all(r["comparison_attained"]=="True" and float(r["HR_pct"])>1 for r in success)
    score=rows("production_joint_scorecard.csv");enabled=[r for r in score if r["candidate_id"]!="disabled_control"]
    assert all(int(r["unattainable_count"])>0 for r in enabled)


def test_TJ_constraint_modes_modify_migration_not_densification():
    p0=production.candidates()["mech_009"];state=memory.initial_state(p0);rates=[]
    for mode in discovery.TJ_CONSTRAINT_MODES:
        p=replace(p0,base=replace(p0.base,TJ_constraint_mode=mode));d=memory.local_competition(state,1200,p);rates.append(d["rho_dot"])
        assert d["P_TJ_assisted_densification"]==d["P_TJ_dens"]
    assert np.allclose(rates,rates[0],rtol=0,atol=0)


def test_current_all_TJ_mode_is_exact_default_and_diagnostics_are_separate():
    p=production.candidates()["mech_009"];a=memory.local_competition(memory.initial_state(p),1200,p);explicit=replace(p,base=replace(p.base,TJ_constraint_mode="current_all_TJ_multihit"));b=memory.local_competition(memory.initial_state(explicit),1200,explicit)
    for key in ("rho_dot","G_dot","Lambda_TJ","K_TJ","P_comp_TJ","P_TJ_multihit"):assert a[key]==b[key]
    assert a["C_TJ_pore"]==a["C_TJ_constraint"] and a["P_TJ_pore_drag"]==0


def test_TJ_ablation_is_nonuniversal_and_keeps_power_channels_separate():
    table=rows("TJ_constraint_ablation.csv");assert {r["TJ_constraint_mode"] for r in table}==set(discovery.TJ_CONSTRAINT_MODES)
    assert all(int(r["n_fast_beneficial"])<int(r["n_fast_attained"]) for r in table)
    pinned=[r for r in table if r["TJ_constraint_mode"]=="pore_pinned_TJ_drag"]
    assert all(float(r["P_TJ_pore_drag_median"])>0 and r["joint_positive"]=="True" for r in pinned)


def test_local_TJ_and_PR_laws_have_no_schedule_leakage():
    source=(inspect.getsource(discovery.local_mechanism)+inspect.getsource(memory.local_competition)).lower()
    for word in ("protocol","schedule","ramp_rate","slow","fast","rho_target"):assert word not in source
