import numpy as np
import topology_constrained_sintering as m
def test_larger_pores_reduce_coverage_at_fixed_volume():
    p=m.Params(n_bins=3);phi=np.array([.1,0,0]);r=np.array([20e-9,40e-9,80e-9]);assert m.infer_topology(.9,200e-9,4*r,phi,p).f_pore<m.infer_topology(.9,200e-9,r,phi,p).f_pore
def test_partition_nonnegative_and_conservative():
    p=m.Params();s=m.initial_state(p);k,mech=m.evaluate_mechanisms(s,1300,p);w=m.solve_dissipation_partition(s,s.topology,mech,p);assert all(x>=0 for x in w.values());assert np.isclose(sum(w.values()),1);assert k['tau_event']>=k['tau_exchange']+k['tau_transport']
def test_pore_conservation_and_nonnegative_bins():
    r=m.run(m.Params(rho0=.83,G0=100e-9,t_max_s=3600),m.Iso(1300,3600));assert np.all(r['pore_phi']>=0);assert np.all(r['pore_N']>=0);assert np.allclose(r['rho'],1-r['pore_phi'].sum(axis=1),atol=1e-12)
def test_required_diagnostics():
    r=m.run(m.Params(t_max_s=60),m.Iso(1300,60));required={'rho','G','pore_phi','f_pore','f_clean','f_PR','f_TL','sigma_base','sigma_concentration','sigma_local','r_nuc','tau_exchange','tau_transport','tau_TL','activity','rho_dot','dGdt','E_G','power_renewal_densification'};assert required<=r.keys()
