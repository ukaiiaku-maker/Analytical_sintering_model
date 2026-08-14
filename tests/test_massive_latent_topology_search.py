import inspect
import numpy as np
import massive_latent_topology_models as m
import massive_latent_topology_objectives as o

def test_conservative_transfer_and_density_identity():
 p=m.default_parameters();s=m.initial_state(.70,100);before=1-s.rho
 p.update(k_open=0.,k_closed=0.);m.step(s,1300,p,1000)
 assert np.isclose(before,s.phi_open+s.phi_connected+s.phi_isolated+s.phi_closed)
 assert np.isclose(s.rho,1-(s.phi_open+s.phi_connected+s.phi_isolated+s.phi_closed))
 assert min(s.phi_open,s.phi_connected,s.phi_isolated,s.phi_closed)>=0

def test_closed_store_not_removed_by_open_channel():
 p=m.default_parameters();p.update(k_closed=0.,closure_rate=0.,detachment_rate=0.)
 s=m.initial_state(.90,100);s.phi_closed=.02;s.phi_open-=.02;v=s.phi_closed;m.step(s,1400,p,1000);assert s.phi_closed==v

def test_projection_cannot_be_a_tier_and_unattained_is_rejected():
 h={'rho':np.array([.7,.94]),'G_nm':np.array([100.,110.])};assert o.trajectory_score(h,h)['tier']=='reject'

def test_local_law_has_no_schedule_leakage():
 src=inspect.getsource(m.rates).lower()
 assert not any(x in src for x in ('protocol','schedule','ramp_rate','slow','fast','target'))

def test_success_window_requires_both_boundaries():
 pts=[{'T2_C':900,'classification':'density_exhaustion'},{'T2_C':1000,'classification':'success'},{'T2_C':1050,'classification':'success'},{'T2_C':1100,'classification':'grain_growth'}]
 assert o.chen_window(pts)['complete'];assert not o.chen_window(pts[:-1])['complete']

def test_extreme_local_step_never_creates_negative_density_or_store():
 p=m.default_parameters();p.update(k_open=1e-2,connected_loss=1e-2,connected_recovery=1e-2)
 s=m.initial_state(.70,100);m.step(s,1500,p,1e6)
 assert 0<=s.rho<=1 and min(s.phi_open,s.phi_connected,s.phi_isolated,s.phi_closed)>=0
