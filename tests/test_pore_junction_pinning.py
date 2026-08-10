import csv
import inspect
from dataclasses import replace
from pathlib import Path

import numpy as np

import growth_mechanism_sensitivity as common
import pore_junction_pinning_sensitivity as study
import topology_constrained_sintering as model

RESULTS=Path(__file__).parents[1]/"results"/"pore_junction_pinning"


def combined_rates(state,params,T=1200.):
    _,mechanisms=model.evaluate_mechanisms(state,T,params)
    weights=model.solve_dissipation_partition(state,state.topology,mechanisms,params)
    return model.combine(mechanisms,weights)


def test_pinning_is_local_observable_and_schedule_free():
    source=inspect.getsource(model.growth_mobility_diagnostics).lower()
    assert not any(word in source for word in ("protocol","schedule","ramp_rate","slow","fast","target"))
    p=replace(common.base_params(),growth_mode="pore_junction_pinning",G0=150e-9)
    s=model.initial_state(p);base=model.growth_mobility_diagnostics(s,1200.,p)
    s.pore_phi[:-1]*=.2;s.pore_phi[-1]+=.8*float(np.sum(model.initial_state(p).pore_phi[:-1]))
    s.rho=1-float(np.sum(s.pore_phi));s.pore_N=model.pore_number(s.pore_phi,s.pore_radii)
    s.topology=model.infer_topology(s.rho,s.G,s.pore_radii,s.pore_phi,p)
    changed=model.growth_mobility_diagnostics(s,1200.,p)
    assert 0.<base["growth_mobility_factor"]<=1.
    assert changed["pore_line_density_m2"]!=base["pore_line_density_m2"]
    assert changed["growth_mobility_factor"]!=base["growth_mobility_factor"]


def test_pinning_changes_migration_but_not_instantaneous_densification():
    base=common.base_params();state=model.initial_state(base)
    reference=combined_rates(state,base)
    pinned=combined_rates(state,replace(base,growth_mode="pore_junction_pinning"))
    assert pinned.rho_dot==reference.rho_dot
    assert pinned.G_dot<reference.G_dot


def test_baseline_recovery_is_exact_and_budgets_are_uniform():
    rows=list(csv.DictReader((RESULTS/"baseline_recovery.csv").open()))
    assert max(float(r["rho2_abs_difference"]) for r in rows)<1e-12
    assert max(float(r["G2_nm_abs_difference"]) for r in rows)<1e-9
    with (RESULTS/"pinning_trajectories.csv").open() as stream:
        for i,row in enumerate(csv.DictReader(stream)):
            assert float(row["first_step_budget_h"])==float(row["second_step_budget_h"])==common.STEP_BUDGET_S/3600
            if i>1000:break


def test_pinning_retains_both_failure_boundaries_and_does_not_fake_nanoscale_success():
    rows=list(csv.DictReader((RESULTS/"pinning_classifications.csv").open()))
    pin=[r for r in rows if r["growth_mode"]=="pore_junction_pinning"]
    categories={r["classification"] for r in pin}
    assert "DENSIFICATION_EXHAUSTION_FAILURE" in categories
    assert "GRAIN_GROWTH_FAILURE" in categories
    boundaries=list(csv.DictReader((RESULTS/"pinning_boundaries.csv").open()))
    nanoscale=[r for r in boundaries if r["growth_mode"]=="pore_junction_pinning" and float(r["G0_nm"])<=300]
    assert all(r["window_width_C"] in ("","nan") for r in nanoscale)


def test_pinning_mode_preserves_pore_conservation_and_nonnegative_bins():
    p=replace(common.base_params(),growth_mode="pore_junction_pinning",G0=100e-9)
    h=model.run(p,model.Iso(1200.,3600.))
    assert np.min(h["pore_phi"])>=-1e-15 and np.min(h["pore_N"])>=-1e-15
    assert np.max(np.abs(h["rho"]-(1.-np.sum(h["pore_phi"],axis=1))))<1e-12
    assert np.all((h["growth_mobility_factor"]>0)&(h["growth_mobility_factor"]<=1))
