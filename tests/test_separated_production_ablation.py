from dataclasses import replace
import inspect,re
import separated_fast_chen_model as m

def state():
    p=m.MaterialKinetics();s=m.initial_state(p);return p,s

def test_no_pr_changes_only_pr_channel_at_shared_state():
    p,s=state();a=m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,p);b=m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,replace(p,ablation_mode='no_PR_redistribution'))
    assert b['PR_propensity']==0
    for key in ('rho_dot','tau_nuc','tau_exchange','tau_transport','growth_base'):assert a[key]==b[key]

def test_no_nucleation_changes_nucleation_activity_not_transport():
    p,s=state();a=m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,p);b=m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,replace(p,ablation_mode='no_nucleation_limitation'))
    assert b['tau_nuc']<a['tau_nuc'];assert b['tau_exchange']==a['tau_exchange'];assert b['tau_transport']==a['tau_transport'];assert b['activity']!=a['activity']

def test_local_ablation_has_no_schedule_identity():
    src=inspect.getsource(m.material_rates).lower()
    assert not any(re.search(rf'\b{x}\b',src) for x in ('protocol','schedule','ramp_rate','slow','fast','target'))

def test_extreme_promotion_rule_is_explicit():
    src=inspect.getsource(__import__('separated_mechanism_production_search').metrics)
    assert 'promotion_blocked' in src and '>100' in src
