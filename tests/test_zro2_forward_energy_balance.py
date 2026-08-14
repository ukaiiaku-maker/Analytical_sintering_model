from zro2_forward.energy_balance import solve_effective_stress
def test_bounded_and_excess_reported():
    x=solve_effective_stress(.5,-1e9,1.,lambda s:1e-20,sigma_max=2.5e8)
    assert x.sigma_eff_Pa==2.5e8 and x.stress_bound_hit and x.P_excess_W_m3>0
