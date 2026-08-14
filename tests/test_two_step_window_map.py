import csv
from dataclasses import replace
from pathlib import Path

import numpy as np

import topology_constrained_sintering as model
import two_step_window_map as window

ROOT=Path(__file__).parents[1]
RESULTS=ROOT/"results"/"two_step_window_map"


def test_classification_truth_table_is_mutually_exclusive():
    cases=((True,.91,.01,.90,.05,"SUCCESS"),(True,.91,.08,.90,.05,"GRAIN_GROWTH_FAILURE"),
           (True,.89,.01,.90,.05,"DENSIFICATION_EXHAUSTION_FAILURE"),(True,.89,.08,.90,.05,"MIXED_FAILURE"),
           (False,np.nan,np.nan,.90,.05,"UNATTAINABLE_FIRST_STEP"))
    for attained,rho,growth,target,tolerance,expected in cases:
        result=window.classify(attained,rho,growth,target,tolerance)
        assert result==expected and result in window.CLASSIFICATIONS


def test_success_requires_density_and_growth_and_high_targets_are_not_hidden():
    rows=list(csv.DictReader((RESULTS/"two_step_window_points.csv").open()))
    assert {float(r["rho_target"]) for r in rows}==set(window.RHO_TARGETS)
    assert {float(r["growth_tolerance"]) for r in rows}==set(window.GROWTH_TOLERANCES)
    for row in rows:
        if row["classification"]=="SUCCESS":
            assert row["first_step_attained"]=="True"
            assert float(row["rho_final"])>=float(row["rho_target"])-1e-12
            assert float(row["growth_fraction"])<=float(row["growth_tolerance"])+1e-12
        if row["classification"]=="UNATTAINABLE_FIRST_STEP":
            assert row["first_step_attained"]=="False"
    assert not any(r["classification"]=="SUCCESS" for r in rows if float(r["rho_target"]) in (.95,.98))


def test_time_budgets_and_material_parameters_are_fixed():
    rows=list(csv.DictReader((RESULTS/"two_step_window_points.csv").open()))
    assert {float(r["first_step_budget_h"]) for r in rows}=={window.FIRST_STEP_BUDGET_S/3600}
    assert {float(r["second_step_budget_h"]) for r in rows}=={window.SECOND_STEP_BUDGET_S/3600}
    base=window.base_params();params=[replace(base,G0=g*1e-9,smoothing_gate_mode=mode) for g in window.G0_NM for mode in window.GATE_MODES]
    window.assert_fixed_parameters(params,base)


def test_restart_uses_real_state_and_preserves_pores():
    params=replace(window.base_params(),t_max_s=3600)
    first=model.run(params,model.Iso(1300,1800),stop_at_rho=.72)
    state=model.state_from_result(first,params)
    assert np.isclose(state.rho,first["rho"][-1])
    second=model.run(params,model.Iso(1200,1800),initial=state)
    assert np.all(second["pore_phi"]>=0) and np.all(second["pore_N"]>=0)
    assert np.allclose(second["rho"],1-second["pore_phi"].sum(axis=1),atol=1e-12)


def test_size_summary_does_not_misclassify_attained_float_sizes():
    rows=list(csv.DictReader((RESULTS/"size_window_summary.csv").open()))
    assert all(int(float(r["n_first_step_states"]))==12 for r in rows)
    assert not any(r["size_window_classification"]=="FIRST_STEP_UNATTAINABLE" for r in rows)
