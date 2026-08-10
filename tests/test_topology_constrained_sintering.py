import ast
import csv
import inspect
from pathlib import Path

import numpy as np
import topology_constrained_sintering as m
def test_larger_pores_reduce_coverage_at_fixed_volume():
    p=m.Params(n_bins=3);phi=np.array([.1,0,0]);r=np.array([20e-9,40e-9,80e-9]);assert m.infer_topology(.9,200e-9,4*r,phi,p).f_pore<m.infer_topology(.9,200e-9,r,phi,p).f_pore
def test_partition_nonnegative_and_conservative():
    p=m.Params();s=m.initial_state(p);k,mech=m.evaluate_mechanisms(s,1300,p);w=m.solve_dissipation_partition(s,s.topology,mech,p);assert all(x>=0 for x in w.values());assert np.isclose(sum(w.values()),1);assert k['tau_event']>=k['tau_exchange']+k['tau_transport']
def test_pore_conservation_and_nonnegative_bins():
    r=m.run(m.Params(rho0=.83,G0=100e-9,t_max_s=3600),m.Iso(1300,3600));assert np.all(r['pore_phi']>=0);assert np.all(r['pore_N']>=0);assert np.allclose(r['rho'],1-r['pore_phi'].sum(axis=1),atol=1e-12)
def test_required_diagnostics():
    r=m.run(m.Params(t_max_s=60),m.Iso(1300,60));required={'rho','G','pore_phi','f_pore','f_clean','f_PR','f_TL','topology_damage','topology_damage_rate','sigma_base','sigma_concentration','sigma_local','r_nuc','tau_exchange','tau_transport','tau_TL','activity','rho_dot','dGdt','E_G','power_renewal_densification'};assert required<=r.keys()
def test_topology_memory_flips_heating_rate_sign_and_preserves_two_step():
    from dataclasses import replace
    target=.90;p=m.Params(memory_model='empirical_topology_damage');off=replace(p,enable_topology_memory=False)
    def protocols(q):return [m.run(q,m.RampHold(.2),target),m.run(q,m.RampHold(20),target),m.run(replace(q,G0=75e-9),m.Iso(1350),target),m.run(replace(q,G0=75e-9),m.TwoStep(),target)]
    enabled=protocols(p);disabled=protocols(off)
    assert all(m.value_at_density(r,target)[1] for r in enabled+disabled)
    ge=[m.value_at_density(r,target)[0] for r in enabled];gd=[m.value_at_density(r,target)[0] for r in disabled]
    assert m.percent_gain(ge[0],ge[1])>0 and m.percent_gain(ge[2],ge[3])>0
    assert enabled[0]['topology_damage'][-1]>enabled[1]['topology_damage'][-1]
    disabled_hr=m.percent_gain(gd[0],gd[1])
    assert np.isclose(disabled_hr,-6.450343467465089,rtol=1e-10,atol=1e-10)
    for r in enabled+disabled:
        assert np.all((r['topology_damage']>=0)&(r['topology_damage']<=1))
        assert np.all(r['pore_phi']>=0) and np.all(r['pore_N']>=0)
        assert np.allclose(r['rho'],1-r['pore_phi'].sum(axis=1),atol=1e-12)
        assert np.max(r['sigma_local'])<=p.stress_cap
def test_damage_rate_has_no_schedule_inputs_or_labels():
    signature=inspect.signature(m.topology_damage_rate)
    assert tuple(signature.parameters)==('s','T_C','p','k')
    source=inspect.getsource(m.topology_damage_rate)
    names={node.id for node in ast.walk(ast.parse(source)) if isinstance(node,ast.Name)}
    assert not names&{'protocol','schedule','ramp_rate','heating_rate','slow','fast'}
def test_reported_successes_really_reach_target():
    root=Path(__file__).parents[1]/'results'/'topology_memory_stress'
    with (root/'held_out_heating_rates.csv').open() as stream:
        for row in csv.DictReader(stream):
            if row['reached_target']=='True':assert float(row['final_density'])>=float(row['target_density'])-1e-12
    with (root/'held_out_two_step_grid.csv').open() as stream:
        for row in csv.DictReader(stream):
            target=float(row['target_density'])
            if row['high_reached']=='True':assert float(row['high_final_density'])>=target-1e-12
            if row['two_step_reached']=='True':assert float(row['two_step_final_density'])>=target-1e-12
    pore_root=Path(__file__).parents[1]/'results'/'pore_bin_memory_stress'
    with (pore_root/'held_out_heating_by_mode.csv').open() as stream:
        for row in csv.DictReader(stream):
            if row['reached_target']=='True':assert float(row['final_density'])>=float(row['target_density'])-1e-12
    with (pore_root/'held_out_two_step_by_mode.csv').open() as stream:
        for row in csv.DictReader(stream):
            target=float(row['target_density'])
            if row['high_reached']=='True':assert float(row['high_final_density'])>=target-1e-12
            if row['two_step_reached']=='True':assert float(row['two_step_final_density'])>=target-1e-12
def test_surface_redistribution_is_local_and_volume_conservative():
    signature=inspect.signature(m.surface_smoothing_redistribution)
    assert tuple(signature.parameters)==('s','T','p','k')
    source=inspect.getsource(m.surface_smoothing_redistribution)
    names={node.id for node in ast.walk(ast.parse(source)) if isinstance(node,ast.Name)}
    assert not names&{'protocol','schedule','ramp_rate','heating_rate','slow','fast'}
    p=m.Params(memory_model='pore_bin_redistribution');s=m.initial_state(p)
    k=m.kinetic_diagnostics(s,1025,p);flux=m.surface_smoothing_redistribution(s,1025,p,k)
    assert np.isclose(np.sum(flux.pore_phi_dot),0,atol=1e-18)
    assert flux.rho_dot==0
    assert np.all(s.pore_phi+flux.pore_phi_dot*p.dt_max_s>=0)
def test_pore_bin_memory_produces_observable_memory_and_correct_signs():
    from dataclasses import replace
    target=.90;p=m.Params(memory_model='pore_bin_redistribution')
    runs=[m.run(p,m.RampHold(.2),target),m.run(p,m.RampHold(20),target),m.run(replace(p,G0=75e-9),m.Iso(1350),target),m.run(replace(p,G0=75e-9),m.TwoStep(),target)]
    values=[m.value_at_density(r,target) for r in runs]
    assert all(ok for _,ok in values)
    assert m.percent_gain(values[0][0],values[1][0])>0
    assert m.percent_gain(values[2][0],values[3][0])>0
    assert runs[0]['pore_mean_radius'][-1]>runs[1]['pore_mean_radius'][-1]
    assert runs[0]['large_pore_fraction'][-1]>runs[1]['large_pore_fraction'][-1]
    assert runs[0]['removable_fine_pore_fraction'][-1]<runs[1]['removable_fine_pore_fraction'][-1]
    assert runs[0]['cumulative_redistributed_pore_volume'][-1]>runs[1]['cumulative_redistributed_pore_volume'][-1]
    for r in runs:
        assert np.all(r['pore_phi']>=0) and np.all(r['pore_N']>=0)
        assert np.allclose(r['rho'],1-r['pore_phi'].sum(axis=1),atol=1e-12)
        assert np.max(r['sigma_local'])<=p.stress_cap
def test_memory_none_recovers_old_negative_heating_rate_result():
    p=m.Params(memory_model='none');slow=m.run(p,m.RampHold(.2),.9);fast=m.run(p,m.RampHold(20),.9)
    gs=m.value_at_density(slow,.9)[0];gf=m.value_at_density(fast,.9)[0]
    assert np.isclose(m.percent_gain(gs,gf),-6.450343467465089,rtol=1e-10,atol=1e-10)
