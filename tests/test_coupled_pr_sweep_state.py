import inspect,numpy as np
import coupled_pr_sweep_state_model as m
def test_PR_conservative_and_no_direct_density():
 p=m.defaults();p.update(k_open=0.,k_closed=0.,k_sweep=0.,closed_transition=0.);s=m.State();v=m.pore_volume(s);m.advance(s,1000,p,1e4);assert np.isclose(v,m.pore_volume(s)) and np.isclose(s.rho,1-v)
def test_sweep_uses_damaged_store_and_reduces_number():
 p=m.defaults();p.update(k_PR=0.,k_open=0.,k_closed=0.,closed_transition=0.);s=m.State(phi_PR_damaged_connected=.05,N_PR_damaged_connected=.2,phi_connected_fine=.13);n=m.pore_number(s);m.advance(s,1400,p,1e5);assert s.sweep_coalescence_memory>0 and m.pore_number(s)<=n
def test_closed_named_zero_without_store():assert m.local_rates(m.State(),1200,m.defaults())['rho_dot_closed']==0
def test_capacity_bounded_nonnegative():
 s=m.State();m.advance(s,1400,m.defaults(),1e6);assert s.closed_accommodation_capacity>=0
def test_locality():
 src=inspect.getsource(m.local_rates).lower();assert not any(x in src for x in ('schedule','protocol','slow','fast','target','success'))
