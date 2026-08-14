import numpy as np
from zro2_forward.pore_population import initial_population, transfer_fluxes, removal_weights

def test_initial_state_and_identity():
    p=initial_population(); assert p.phi_closed.sum()==0 and p.phi_iso.sum()==0
    assert np.isclose(1-p.total,2.88/5.95)

def test_redistribution_conserves_volume_and_does_not_change_topology():
    p=initial_population(); a,b,c,_=transfer_fluxes(p,.7,1e-20,.1,1.,C_iso=0,C_close=0)
    assert np.isclose(a.sum(),0,atol=1e-25) and np.all(b==0) and np.all(c==0)

def test_removal_scaling():
    p=initial_population(bins=2); p.radii_m=np.array([1.,2.]); p.phi_open=np.ones(2)
    w=removal_weights(p); assert np.isclose(w[0]/w[1],16.)

def test_transitions_conserve_volume():
    p=initial_population(); a,b,c,_=transfer_fluxes(p,.95,1e-20,.1,0.,C_PR=0)
    assert np.isclose((a+b+c).sum(),0,atol=1e-25)
