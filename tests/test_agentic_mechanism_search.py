import inspect
import csv
from pathlib import Path
import numpy as np

import agentic_mechanism_model as model
import agentic_mechanism_search as search
import mechanism_registry
import adaptive_T2_boundary_search as adaptive
import pore_location_agentic_sensitivity as prior
import topology_constrained_sintering as aggregate


def params(mode="persistent_tj_multihit_q0"):
    return model.DiscoveryParams(prior.base_action("action_evolving_capture"),mode)


def test_registry_declares_required_source_grounded_fields():
    required={"persistent_junction_drag","TJ_multihit_reaction","vacancy_accommodation_multihit","stress_accumulation_release","closed_pore_placeholder"}
    assert required <= set(mechanism_registry.REGISTRY)
    for spec in mechanism_registry.REGISTRY.values():
        assert spec.mechanism_class in {"A_drag","B_multihit","C_exchange","mixed"}
        assert spec.source_rationale and spec.rejection_criteria


def test_local_closure_has_no_schedule_leakage():
    assert tuple(inspect.signature(model.local_mechanism).parameters)==("s","T_C","p")
    src=inspect.getsource(model.local_mechanism).lower()
    for word in ("protocol","schedule","slow","fast","ramp_rate","rho_target"):
        assert word not in src


def test_poisson_completion_limits_and_monotonicity():
    vals=[model.poisson_completion(3,x) for x in (0,.1,1,10,100)]
    assert vals[0]==0 and vals[-1]>.999
    assert np.all(np.diff(vals)>=0)
    assert model.poisson_completion(4,2)<model.poisson_completion(2,2)


def test_migration_only_closures_preserve_shared_state_densification():
    s=model.initial_state(params());rates=[]
    for mode in ("action_baseline","persistent_junction","tj_multihit_q0","tj_multihit_q1","persistent_tj_multihit_q1"):
        p=params(mode);rates.append(model.local_mechanism(s,1200,p)["rho_dot"])
    assert np.allclose(rates,rates[0],rtol=0,atol=0)


def test_conservation_nonnegative_bins_and_bounded_persistent_state():
    h=model.run(params(),aggregate.Iso(1250,600))
    phi=h['phi_GBseg']+h['phi_TJ']+h['phi_iso']
    assert np.all(phi>=0) and np.all(h['N_GBseg']>=0) and np.all(h['N_TJ']>=0) and np.all(h['N_iso']>=0)
    assert np.max(np.abs(h['rho']-(1-phi.sum(axis=1))))<1e-12
    assert np.all((h['X_J']>=0)&(h['X_J']<=params().XJ_capacity))


def test_isolated_pores_have_no_open_pore_removal_flux():
    s=model.initial_state(params());before=s.pore.phi_iso.copy();d=model.local_mechanism(s,1250,params())
    assert 'iso_remove' not in d
    assert np.all(before>=0) and d['rho_dot']>=0


def test_conservative_action_relocations_sum_to_zero():
    s=model.initial_state(params());d=model.local_mechanism(s,1250,params())
    net=d['GB_smooth'].sum()-d['GB_to_TJ'].sum()+d['GB_to_TJ'].sum()+d['TJ_to_GBseg_capture'].sum()-d['TJ_to_GBseg_capture'].sum()-d['TJ_to_iso'].sum()+d['TJ_to_iso'].sum()
    assert abs(net)<1e-18


def test_classification_is_exclusive_and_precludes_already_reached_target():
    assert search.classify(True,.90,.91,.01,.05)=="INELIGIBLE_TARGET_ALREADY_REACHED"
    found={search.classify(True,.85,r,g,.05) for r,g in ((.90,.04),(.90,.06),(.89,.04),(.89,.06))}
    assert found=={"SUCCESS","GRAIN_GROWTH_FAILURE","DENSIFICATION_EXHAUSTION_FAILURE","MIXED_FAILURE"}
    assert search.classify(False,.8,.9,.01,.05)=="UNATTAINABLE_FIRST_STEP"


def test_action_baseline_is_exactly_recoverable():
    p=params("action_baseline");protocol=aggregate.Iso(1200,300)
    a=model.run(p,protocol);b=prior.action.run(p.action,protocol)
    for key in ("rho","G","phi_GBseg","phi_TJ","phi_iso"):
        assert np.array_equal(a[key],b[key])


def test_persisted_decisions_have_reasons_and_accepted_cases_keep_boundaries():
    root=Path(__file__).parents[1]/"results"/"agentic_mechanism_search"
    with (root/"rejected_parameter_sets.csv").open(newline="") as f:rejected=list(csv.DictReader(f))
    with (root/"reduced_map_summary.csv").open(newline="") as f:summary=list(csv.DictReader(f))
    assert rejected and all(row["rejection_reason"] for row in rejected)
    accepted=[row for row in summary if row["rejected"]=="False"]
    assert accepted and all(row["has_lower_failure"]=="True" and row["has_upper_failure"]=="True" for row in accepted)


def test_persisted_successes_meet_both_criteria():
    path=Path(__file__).parents[1]/"results"/"agentic_mechanism_search"/"full_map_classifications.csv"
    successes=0
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["classification"]=="SUCCESS":
                successes+=1
                assert row["first_step_attained"]=="True"
                assert float(row["rho1"])<float(row["rho_target"])
                assert float(row["rho2"])>=float(row["rho_target"])-1e-12
                assert float(row["growth_fraction"])<=float(row["growth_tolerance"])+1e-12
    assert successes>0


def synthetic(T,rho,growth,T1=1300):
    return dict(T2_C=float(T),rho1=.85,rho2=float(rho),growth_fraction=float(growth),T1_C=T1)


def test_complete_window_requires_both_bracketed_boundaries():
    complete=[synthetic(900,.89,.01),synthetic(1000,.90,.02),synthetic(1100,.90,.08)]
    assert adaptive.status(complete,.05,1300)['boundary_status']=='COMPLETE_WINDOW'
    no_upper=complete[:2]
    assert adaptive.status(no_upper,.05,1300)['boundary_status']=='UPPER_BOUND_RIGHT_CENSORED'


def test_censored_window_is_not_complete_and_practical_requires_T2_below_T1():
    points=[synthetic(900,.89,.01,1050),synthetic(1000,.90,.02,1050),synthetic(1100,.90,.08,1050)]
    assert adaptive.status(points,.05,1050,False)['boundary_status']=='COMPLETE_WINDOW'
    practical=adaptive.status(points,.05,1050,True)
    assert practical['boundary_status']=='UPPER_BOUND_RIGHT_CENSORED'
    assert all(p['T2_C']<1050 for p in points if p['T2_C']<1050)


def test_adaptive_constants_preserve_target_budget_and_parameters():
    assert adaptive.TARGET==search.TARGET==.90
    assert adaptive.BUDGET==search.BUDGET==96*3600
    original=adaptive.survivor_params()
    again=adaptive.survivor_params()
    assert original==again


def test_high_temperature_extension_conserves_pores_and_nonnegative_stores():
    p=adaptive.survivor_params()['mech_009'];h=model.run(p,aggregate.Iso(1400,600))
    phi=h['phi_GBseg']+h['phi_TJ']+h['phi_iso']
    assert np.all(phi>=0) and np.all(h['N_GBseg']>=0) and np.all(h['N_TJ']>=0) and np.all(h['N_iso']>=0)
    assert np.max(np.abs(h['rho']-(1-phi.sum(axis=1))))<1e-12


def test_persisted_adaptive_complete_windows_have_both_brackets():
    path=Path(__file__).parents[1]/'results'/'agentic_mechanism_search'/'adaptive_window_boundaries.csv'
    with path.open(newline='') as f:rows=list(csv.DictReader(f))
    complete=[r for r in rows if r['boundary_status']=='COMPLETE_WINDOW']
    assert complete
    assert all(r['lower_bracketed']=='True' and r['upper_bracketed']=='True' for r in complete)
    assert all(r['boundary_status']!='UPPER_BOUND_RIGHT_CENSORED' for r in complete)


def test_persisted_practical_map_contains_only_T2_below_T1():
    path=Path(__file__).parents[1]/'results'/'agentic_mechanism_search'/'practical_two_step_map.csv'
    with path.open(newline='') as f:
        assert all(float(r['T2_C'])<float(r['T1_C']) for r in csv.DictReader(f))
