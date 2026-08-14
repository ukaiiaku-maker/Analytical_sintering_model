import csv
import inspect
import math
from dataclasses import replace
from pathlib import Path

import numpy as np

import growth_mechanism_sensitivity as study
import topology_constrained_sintering as model

RESULTS=Path(__file__).parents[1]/"results"/"growth_mechanism_sensitivity"


def test_baseline_is_exact_identity_and_reproduces_negative_control():
    p=study.base_params();s=model.initial_state(p)
    d=model.growth_mobility_diagnostics(s,1100.,p)
    assert d["growth_mobility_factor"]==1. and d["junction_time_ratio"]==0.
    rows=list(csv.DictReader((RESULTS/"baseline_recovery.csv").open()))
    assert max(float(r["rho2_abs_difference"]) for r in rows)<1e-12
    assert max(float(r["G2_nm_abs_difference"]) for r in rows)<1e-9


def test_junction_mobility_is_local_bounded_and_thermally_activated():
    p=replace(study.base_params(),growth_mode="junction_limited")
    small=model.initial_state(replace(p,G0=50e-9));large=model.initial_state(replace(p,G0=300e-9))
    low=model.growth_mobility_diagnostics(small,1000.,p)["growth_mobility_factor"]
    hot=model.growth_mobility_diagnostics(small,1300.,p)["growth_mobility_factor"]
    coarse=model.growth_mobility_diagnostics(large,1000.,p)["growth_mobility_factor"]
    assert 0.<low<hot<=1. and low<coarse<=1.
    source=inspect.getsource(model.growth_mobility_diagnostics).lower()
    assert not any(label in source for label in ("protocol","schedule","ramp_rate","slow","fast","target"))


def test_growth_mode_does_not_change_instantaneous_densification_propensity():
    base=study.base_params();state=model.initial_state(base)
    a=model.clean_GB_migration(state,1200.,base)
    b=model.clean_GB_migration(state,1200.,replace(base,growth_mode="junction_limited"))
    assert a.rho_dot==b.rho_dot==0.
    assert math.isclose(a.power,b.power,rel_tol=0.,abs_tol=0.)
    assert b.G_dot<a.G_dot


def test_map_uses_fixed_parameters_and_uniform_time_budgets():
    base=study.base_params()
    params=[replace(base,growth_mode=m,G0=g*1e-9) for m in study.GROWTH_MODES for g in (25.,600.)]
    study.assert_only_growth_mode_and_size_vary(params,base)
    with (RESULTS/"growth_mode_trajectories.csv").open() as stream:
        for i,row in enumerate(csv.DictReader(stream)):
            assert float(row["first_step_budget_h"])==float(row["second_step_budget_h"])==study.STEP_BUDGET_S/3600
            if i>1000:break


def test_classifications_are_exclusive_and_success_requires_both_conditions():
    with (RESULTS/"growth_mode_classifications.csv").open() as stream:
        for row in csv.DictReader(stream):
            category=row["classification"]
            assert category in {"SUCCESS","GRAIN_GROWTH_FAILURE","DENSIFICATION_EXHAUSTION_FAILURE","MIXED_FAILURE","UNATTAINABLE_FIRST_STEP"}
            if category=="SUCCESS":
                assert row["first_step_attained"]=="True"
                matching=next(r for r in csv.DictReader((RESULTS/"growth_mode_trajectories.csv").open())
                              if all(r[k]==row[k] for k in ("growth_mode","G0_nm","T1_C","rho_switch","T2_C")))
                assert float(matching["rho2"])>=study.TARGET-1e-12
                assert float(matching["growth_fraction"])<=float(row["growth_tolerance"])+1e-12


def test_map_retains_both_physical_failure_boundaries():
    rows=list(csv.DictReader((RESULTS/"growth_mode_classifications.csv").open()))
    junction={r["classification"] for r in rows if r["growth_mode"]=="junction_limited"}
    assert "DENSIFICATION_EXHAUSTION_FAILURE" in junction
    assert "GRAIN_GROWTH_FAILURE" in junction


def test_all_growth_modes_preserve_pore_conservation_and_nonnegative_bins():
    for mode in study.GROWTH_MODES:
        p=replace(study.base_params(),growth_mode=mode,G0=100e-9)
        h=model.run(p,model.Iso(1200.,1800.))
        assert np.min(h["pore_phi"])>=-1e-15 and np.min(h["pore_N"])>=-1e-15
        assert np.max(np.abs(h["rho"]-(1.-np.sum(h["pore_phi"],axis=1))))<1e-12
