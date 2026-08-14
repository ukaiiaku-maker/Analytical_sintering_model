import inspect,re
import numpy as np
import latent_topology_memory_model as m
import latent_topology_memory_search as s
import latent_topology_objectives as o

def test_three_dynamic_families_are_explicit():
    x=m.LatentState();assert hasattr(x,'removable_pore_memory') and hasattr(x,'X_J') and hasattr(x,'stress_memory')

def test_dynamic_advance_bounds_stores():
    x=m.LatentState();p=s.decode(np.mean(np.array(s.BOUNDS),axis=1));rates={'PR_propensity':1e-4,'activity':.2,'rho_dot':1e-6,'growth_base':1e-12};m.advance(x,rates,1300,p,10)
    assert 0<=x.removable_pore_memory<=1 and 0<=x.X_J<=p['XJ_capacity'] and 0<=x.stress_memory<=p['stress_cap']

def test_migration_factor_does_not_change_rho_dot():
    import separated_fast_chen_model as b
    mat=b.MaterialKinetics();st=b.initial_state(mat);p=s.decode(np.mean(np.array(s.BOUNDS),axis=1));rate=m.shared_rates(mat,st,1300)['rho_dot'];m.migration_factor(m.LatentState(X_J=.2),1300,st['G'],p);assert m.shared_rates(mat,st,1300)['rho_dot']==rate

def test_projection_only_is_not_registered_and_tiers_reject_prep_growth():
    assert o.components(1,1,50,300,.2,.01)['tier']=='Tier_C'

def test_local_law_has_no_schedule_words():
    for fn in (m.derivatives,m.advance,m.migration_factor):
        src=inspect.getsource(fn).lower();assert not any(re.search(rf'\b{x}\b',src) for x in ('protocol','schedule','ramp_rate','slow','fast','target'))
