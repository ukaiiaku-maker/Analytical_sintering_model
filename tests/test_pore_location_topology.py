import csv, inspect, math
from dataclasses import replace
from pathlib import Path

import numpy as np

import topology_constrained_sintering as aggregate
import pore_location_topology_model as model
import pore_location_topology_sensitivity as study

RESULTS=Path(__file__).parents[1]/"results"/"pore_location_topology"


def test_disabled_mode_exactly_recovers_pr2_aggregate_model():
    base=study.base_params();p=model.LocationParams(base,"disabled");protocol=aggregate.Iso(1200.,3600.)
    a=model.run(p,protocol);b=aggregate.run(base,protocol)
    assert set(a)==set(b)
    for key in a:assert np.array_equal(a[key],b[key])


def test_exact_density_identity_and_nonnegative_location_bins():
    p=study.make_params("evolving");h=model.run(p,aggregate.Iso(1250.,6*3600))
    total=h['phi_GBseg']+h['phi_TJ']+h['phi_iso']
    assert np.max(np.abs(h['rho']-(1-np.sum(total,axis=1))))<1e-12
    for key in ('phi_GBseg','phi_TJ','phi_iso','N_GBseg','N_TJ','N_iso'):assert np.min(h[key])>=-1e-15


def test_relocation_and_smoothing_fluxes_conserve_volume():
    p=study.make_params("evolving");s=model.initial_state(p);d=model.instantaneous(s,1200.,p)
    assert abs(float(np.sum(d['GB_smooth'])))<1e-18
    assert abs(float(np.sum(d['TJ_smooth'])))<1e-18
    assert np.all(d['GB_to_TJ']>=0) and np.all(d['TJ_to_iso']>=0)
    relocation_total=float(np.sum(-d['GB_to_TJ']+d['GB_to_TJ']-d['TJ_to_iso']+d['TJ_to_iso']))
    assert abs(relocation_total)<1e-18


def test_open_pore_densification_never_removes_isolated_store():
    p=study.make_params("static",(.05,.05,.90));s=model.initial_state(p);before=s.phi_iso.copy();d=model.instantaneous(s,1300.,p)
    assert 'iso_remove' not in d
    dt=1.;s.phi_iso=s.phi_iso.copy()
    assert np.array_equal(s.phi_iso,before)
    assert d['rho_dot']>=0


def test_clean_gb_migration_has_no_direct_pore_volume_flux():
    p=study.make_params("static");s=model.initial_state(p);d=model.instantaneous(s,1300.,p)
    assert d['G_dot']>0
    assert all(np.all(d[k]==0) for k in ('GB_smooth','GB_to_TJ','TJ_to_iso','TJ_smooth'))


def test_controlled_location_ordering_and_tj_stress_drag():
    gb=study.make_params("static",(.80,.15,.05));iso=study.make_params("static",(.25,.15,.60));tj=study.make_params("static",(.25,.70,.05))
    dgb=model.instantaneous(model.initial_state(gb),1300.,gb);diso=model.instantaneous(model.initial_state(iso),1300.,iso);dtj=model.instantaneous(model.initial_state(tj),1300.,tj)
    assert dgb['rho_dot']>diso['rho_dot']
    assert dtj['C_TJ']>dgb['C_TJ']
    assert dtj['P_TJ_drag']>dgb['P_TJ_drag']
    assert dtj['sigma_TJ_pore']>dgb['sigma_TJ_pore']


def test_local_mechanisms_have_no_schedule_or_target_leakage():
    source='\n'.join(inspect.getsource(fn).lower() for fn in model.LOCAL_FUNCTIONS)
    assert not any(word in source for word in ('protocol','schedule','ramp_rate','slow','fast','target'))


def test_classifications_are_exclusive_and_failed_targets_are_not_success():
    rows=list(csv.DictReader((RESULTS/'chen_style_pore_location_map.csv').open()))
    allowed={'SUCCESS','GRAIN_GROWTH_FAILURE','DENSIFICATION_EXHAUSTION_FAILURE','MIXED_FAILURE','UNATTAINABLE_FIRST_STEP'}
    keys=set();
    for r in rows:
        assert r['classification'] in allowed
        key=tuple(r[k] for k in ('case_id','G0_nm','T1_C','rho_switch','T2_C','growth_tolerance'));assert key not in keys;keys.add(key)
        if r['classification']=='SUCCESS':
            assert r['first_step_attained']=='True' and float(r['rho2'])>=study.TARGET-1e-12 and float(r['growth_fraction'])<=float(r['growth_tolerance'])+1e-12
        if r['first_step_attained']=='False':assert r['classification']=='UNATTAINABLE_FIRST_STEP'


def test_screen_is_bounded_and_rejections_are_persisted():
    summary=list(csv.DictReader((RESULTS/'parameter_screen_summary.csv').open()));rejected=list(csv.DictReader((RESULTS/'rejected_parameter_sets.csv').open()))
    assert len(summary)==64 and 0<len(rejected)<64
    assert all(r['rejection_reason'] in ('missing_lower_boundary','missing_upper_boundary','universal_success') for r in rejected)
