import csv,gzip,math
from dataclasses import replace
from pathlib import Path

import expanded_phase_space_exploration as study

RESULTS=Path(__file__).parents[1]/"results"/"expanded_phase_space"


def test_unattainable_and_attained_rows_have_uniform_schema():
    unattained=study.two_step_group_task(("test","density",150.,.70,1150.,.98,(1000.,)))[0]
    attained=study.two_step_group_task(("test","density",150.,.70,1300.,.75,(1200.,)))[0]
    assert set(unattained)==set(attained)
    assert not unattained["first_step_attained"] and math.isnan(unattained["rho2"])


def test_already_reached_target_is_ineligible_not_success():
    row={"first_step_attained":True,"rho1":.925,"rho2":.94,"growth_fraction":0.,"second_connected_coverage":.5,"second_isolation":0.,"second_final_rho_dot":1e-6}
    assert study.expanded_classify(row,.90,.05)=="INELIGIBLE_TARGET_ALREADY_REACHED"


def test_manifest_runtime_and_refinement_counts_are_persisted():
    manifest=list(csv.DictReader((RESULTS/"sweep_manifest.csv").open()))
    runtime=list(csv.DictReader((RESULTS/"runtime_summary.csv").open()))[0]
    assert len(manifest)==992
    assert int(runtime["coarse_trajectories"])==18848
    assert int(runtime["refined_trajectories"])==8134
    assert int(runtime["fast_trajectories"])==3780
    refined=list(csv.DictReader((RESULTS/"refined_two_step_trajectories.csv").open()))
    assert any(float(r["T2_C"])%25!=0 for r in refined)


def test_complete_classification_table_never_scores_ineligible_target():
    with gzip.open(RESULTS/"all_two_step_classifications.csv.gz","rt") as stream:
        rows=csv.DictReader(stream)
        count=0
        for row in rows:
            count+=1
            if row["eligible_second_step_target"]=="False":
                assert row["classification"] in ("INELIGIBLE_TARGET_ALREADY_REACHED","UNATTAINABLE_FIRST_STEP")
    assert count==30534*len(study.TARGETS)*len(study.TOLERANCES)


def test_fixed_parameters_and_budgets_across_styles():
    base=study.base_params();params=[replace(base,smoothing_gate_mode=s,G0=g*1e-9,rho0=r) for s in study.STYLES for g in (25.,2000.) for r in (.60,.80)]
    study.assert_fixed(params,base)
    with (RESULTS/"two_step_trajectories.csv").open() as stream:
        for i,row in enumerate(csv.DictReader(stream)):
            assert float(row["first_step_budget_h"])==float(row["second_step_budget_h"])==study.STEP_BUDGET_S/3600
            if i>1000:break


def test_failed_fast_targets_have_no_hr_score():
    with (RESULTS/"fast_firing_surface.csv").open() as stream:
        for row in csv.DictReader(stream):
            if row["reached_0p90"]=="False":
                assert math.isnan(float(row["HR_pct_vs_0p2"]))
