import inspect,numpy as np
import pr_lower_bound_coalescence_model as m
def vol(q):return q.phi_connected_fine+q.phi_large_attached+q.phi_large_TJ+q.phi_isolated+q.phi_closed
def test_PR_transfer_conservative_and_observable():
 p=m.defaults();p.update(k_open=0.,closed_pore_shrinkage_prefactor=0.);s=m.PRState(m.base.initial_state());v=vol(s.pore);m.advance(s,1000,p,1e4);assert np.isclose(v,vol(s.pore));assert s.cumulative_PR_surface_energy_loss>=0 and s.connected_fine_pore_loss_from_PR>=0
def test_closed_flux_zero_without_closed_store():
 s=m.PRState(m.base.initial_state());p=m.defaults();assert s.pore.phi_closed==0 and m.local_rates(s,1000,p)['rho_dot_closed']==0
def test_drive_only_is_diagnostic():assert 'PR_drive_loss_no_pore_topology' in m.MODES
def test_local_gate_has_no_labels():
 src=inspect.getsource(m.local_rates).lower();assert not any(x in src for x in ('schedule','protocol','slow','fast','two_step','target','success'))
def test_open_removal_does_not_touch_isolated_closed():
 p=m.defaults();p.update(mode='baseline_coalescence',k_sweep_coalesce=0.,k_TJ_coalescence=0.,k_drag_detach=0.,k_recapture=0.,k_closed_transition=0.,closed_pore_shrinkage_prefactor=0.);s=m.PRState(m.base.initial_state());a=(s.pore.phi_isolated,s.pore.phi_closed);m.advance(s,1400,p,1e3);assert a==(s.pore.phi_isolated,s.pore.phi_closed)
