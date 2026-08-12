import inspect,re
from dataclasses import replace
import dynamic_chen_topology_search as s
import separated_fast_chen_model as m

def test_isolated_success_is_not_window():
    pts=[dict(T2_C=t,T1_C=1300,rho1=.85,rho2=.90 if t==1100 else (.89 if t<1100 else .91),growth_fraction=.01 if t<=1100 else .2,prep_growth_fraction=.01) for t in (1050,1100,1150)]
    assert s.summarize(pts)['rejection_reason']=='isolated_success_not_window'

def test_finite_window_requires_both_boundaries():
    pts=[dict(T2_C=t,T1_C=1300,rho1=.85,rho2=.89 if t==1000 else .91,growth_fraction=.01 if t<1150 else .2,prep_growth_fraction=.01) for t in (1000,1050,1100,1150)]
    x=s.summarize(pts);assert x['outcome']=='finite_window' and x['window_width_C']==50

def test_q0_q1_and_rejection_reasons_visible():
    d=s.design(64);assert {p.q_TJ for _,_,p in d}=={0,1}
    assert s.summarize([])['rejection_reason']=='no_success'

def test_topology_migration_only_and_material_frozen():
    mat=m.MaterialKinetics();x=m.initial_state(mat);rate=m.material_rates(x['rho'],x['G'],x['phi'],x['radii'],1300,mat)['rho_dot']
    for _,_,top in s.design(16):
        m.topology_growth_factor({'G':x['G'],'X_J':0,'connected_coverage':.5},1300,top)
        assert m.material_rates(x['rho'],x['G'],x['phi'],x['radii'],1300,mat)['rho_dot']==rate
    assert mat==replace(mat)

def test_local_laws_have_no_schedule_labels():
    for fn in m.LOCAL_FUNCTIONS:
        src=inspect.getsource(fn).lower();assert not any(re.search(rf'\b{x}\b',src) for x in ('protocol','schedule','ramp_rate','slow','fast','target'))
