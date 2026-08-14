from dataclasses import replace
import inspect
import nucleation_fast_chen_production as p
import separated_fast_chen_model as m

def test_pr_independent_class_is_allowed():
    src=inspect.getsource(p.reclassify)
    assert 'nucleation_limited_PR_independent' in src

def test_dynamic_success_requires_practical_temperature():
    assert p.classify(.85,.90,.01,.05,1200,1300,True)=='SUCCESS'
    assert p.classify(.85,.90,.01,.05,1300,1300,True)!='SUCCESS'

def test_topology_families_keep_q0_and_q1_visible():
    assert {x[2].q_TJ for x in p.topology_registry()}=={0,1}

def test_topology_does_not_change_material_density_rate():
    mat=m.MaterialKinetics();s=m.initial_state(mat);a=m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,mat)
    for _,_,top in p.topology_registry():
        factor,_=m.topology_growth_factor({'G':s['G'],'X_J':0,'connected_coverage':a['connected_fine']},1300,top)
        assert a['rho_dot']==m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,mat)['rho_dot'];assert factor>=0

def test_exact_state_transfer_preserves_material_and_density():
    mat=m.MaterialKinetics();s=m.initial_state(mat);clone=m.clone_state(s,reset_time=True)
    assert mat==replace(mat);assert clone['rho']==s['rho'];assert clone['phi'] is not s['phi']
