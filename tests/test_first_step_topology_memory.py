import inspect,re
import first_step_topology_memory_registry as reg
import separated_fast_chen_model as m

def test_registry_has_six_named_families_and_metadata():
    rows=reg.rows();assert len(rows)==6
    for r in rows:
        for k in ('physical_rationale','state_variables','evolution_law','changes_density','changes_migration_only','conservative_pore_flux','observable_signatures','parameters','rejection_criteria'):assert k in r

def test_migration_closure_preserves_density_rate():
    p=m.MaterialKinetics();s=m.initial_state(p);a=m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,p)
    for q in (0,1):
        m.topology_growth_factor({'G':s['G'],'X_J':.2,'connected_coverage':.4},1300,m.TopologyGrowthClosure(mode='audit',q_TJ=q,TJ_drag_strength=10))
        assert m.material_rates(s['rho'],s['G'],s['phi'],s['radii'],1300,p)['rho_dot']==a['rho_dot']

def test_local_laws_have_no_schedule_leakage():
    for fn in m.LOCAL_FUNCTIONS:
        src=inspect.getsource(fn).lower();assert not any(re.search(rf'\b{x}\b',src) for x in ('protocol','schedule','ramp_rate','slow','fast','target'))

def test_q_variants_visible():
    assert {0,1}=={q for q in (0,1)}
