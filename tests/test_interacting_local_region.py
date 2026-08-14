import inspect,numpy as np
import interacting_local_region_model as m
import interacting_local_region_objectives as o
def total(s):return float(s.weights@(s.phi_GBseg+s.phi_TJ+s.phi_iso+s.phi_closed))
def test_global_density_identity():
 s=m.initial();assert np.isclose(m.global_observables(s)['rho_global'],1-total(s))
def test_conservative_transfer():
 s=m.initial();p=m.defaults();p.update(k_open=0.,k_closed=0.);v=total(s);m.advance(s,1000,p,1000);assert np.isclose(v,total(s))
def test_closed_not_open_and_zero_flux():
 s=m.initial();p=m.defaults();assert np.all(m.local_fluxes(s,1200,p)['rho_dot_closed']==0)
def test_locality():
 src=inspect.getsource(m.local_fluxes).lower();assert not any(x in src for x in ('schedule','protocol','slow','fast','two_step','target'))
def test_exact_required_for_tier():
 s=dict(attained=True,span20=.03,median_reduction=.3,max_reduction=.4);w=dict(complete=True);assert o.assign_tier(s,w,False)=='unscored'
def test_migration_only_parameters_do_not_change_density_flux():
 s=m.initial();p=m.defaults();base=m.local_fluxes(s,1250,p)
 p.update(attached_drag=1000,junction_drag=1000,lambda_TJ=.01,K_TJ=50,q_TJ=2,stress_migration=100)
 changed=m.local_fluxes(s,1250,p)
 assert np.allclose(base['rho_dot_open'],changed['rho_dot_open'])
 assert np.allclose(base['rho_dot_closed'],changed['rho_dot_closed'])
def test_closed_shrinkage_requires_closed_volume():
 s=m.initial();p=m.defaults();s.phi_closed[:]=0
 assert np.all(m.local_fluxes(s,1250,p)['rho_dot_closed']==0)
def test_nonnegative_pore_stores_under_extreme_decoder():
 import interacting_local_region_decoder as d
 p={**m.defaults(),**d.decode(np.full(len(d.NAMES),.999))};s=m.initial(p['N_regions'],p=p)
 for T in (900,1200,1450):m.advance(s,T,p,1800,m.network_adjacency(p['N_regions'],p))
 for name in ('phi_GBseg','phi_TJ','phi_iso','phi_closed'):assert np.all(getattr(s,name)>=0)
 assert np.allclose(s.rho,1-(s.phi_GBseg+s.phi_TJ+s.phi_iso+s.phi_closed))
def test_closed_transition_does_not_change_density_without_densification():
 p=m.defaults();p.update(k_open=0.,k_closed=0.,closed_transition=1e-3)
 s=m.initial();rho=m.global_observables(s)['rho_global'];m.advance(s,1250,p,1800)
 assert np.isclose(rho,m.global_observables(s)['rho_global'])
