import inspect,numpy as np
import grain_growth_pore_coalescence_model as m
import massive_latent_topology_objectives as o
def volume(s):return s.phi_connected_fine+s.phi_large_attached+s.phi_large_TJ+s.phi_isolated+s.phi_closed
def test_coalescence_conserves_volume_and_reduces_number():
 p=m.defaults();p.update(k_open=0.,closed_pore_shrinkage_prefactor=0.,k_drag_detach=0.,k_recapture=0.,k_closed_transition=0.);s=m.initial_state();v=volume(s);n=m.diagnostics(s)['pore_number_reduction_factor'];m.advance(s,1400,p,1e5);assert np.isclose(v,volume(s));assert m.diagnostics(s)['pore_number_reduction_factor']<=n
def test_open_flux_does_not_remove_isolated_closed():
 p=m.defaults();p.update(k_sweep_coalesce=0.,k_drag_detach=0.,k_recapture=0.,k_closed_transition=0.,closed_pore_shrinkage_prefactor=0.);s=m.initial_state();s.phi_closed=.01;a=(s.phi_isolated,s.phi_closed);m.advance(s,1400,p,1e4);assert a==(s.phi_isolated,s.phi_closed)
def test_density_identity_and_nonnegative_stores():
 s=m.initial_state();m.advance(s,1500,m.defaults(),1e6);assert np.isclose(s.rho,1-volume(s));assert min(s.phi_connected_fine,s.phi_large_attached,s.phi_large_TJ,s.phi_isolated,s.phi_closed)>=0
def test_locality():
 src=inspect.getsource(m.local_rates).lower();assert not any(x in src for x in ('schedule','protocol','slow','fast','two_step','target','success'))
def test_tier_requires_attainment_and_finite_span():
 h={'rho':np.array([.7,.98]),'G_nm':np.array([100.,200.])};assert o.trajectory_score(h,h)['tier'] not in ('Tier_A','Tier_B')
