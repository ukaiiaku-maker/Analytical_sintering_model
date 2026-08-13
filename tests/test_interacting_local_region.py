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
