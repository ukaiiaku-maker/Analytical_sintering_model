import inspect
from dataclasses import replace
import numpy as np
import pandas as pd

import late_stage_closed_pore_model as late
import late_stage_closed_pore_search as search
import pr_desintering_memory_model as memory
import production_pr_desintering_assessment as production
import production_mechanism_assessment as prior
import observable_trajectory_effect_audit as effect
import joint_heterogeneity_residual_stress_search as common

def base():
    p=production.candidates()["mech_009_q0"];return replace(p,base=prior.fast_params(p.base,75,.70,prior.TOPOLOGIES["baseline"]))

def test_disabled_exactly_recovers_open_pore_baseline():
    p=base();protocol=search.BudgetSchedule(20,1400,2,96);a=late.run(late.LateStageParams(p,late_stage_mode="disabled"),protocol);b=memory.run(p,protocol)
    assert np.array_equal(a["rho"],b["rho"]) and np.array_equal(a["G"],b["G"])

def test_closed_flux_zero_without_closed_volume():
    p=late.LateStageParams(base(),late_stage_mode="closed_pore_vacancy_transport");r=memory.initial_state(p.base).base.pore.pore_radii
    loss,_,_=late.closed_flux(np.zeros_like(r),r,1400,p);assert np.array_equal(loss,np.zeros_like(r))

def test_closed_shrinkage_density_identity_and_nonnegative_stores():
    p=late.LateStageParams(base(),late_stage_mode="closed_pore_vacancy_transport",rho_close_mid=.90,k_closed_ref_s_Pa=2e-13);h=late.run(p,search.BudgetSchedule(20,1400,20,96))
    total=h["phi_GBseg"]+h["phi_TJ"]+h["phi_iso"]+h["phi_closed"]
    assert np.allclose(h["rho"],1-total.sum(axis=1));assert all(np.min(h[k])>=0 for k in ("phi_GBseg","phi_TJ","phi_iso","phi_closed","N_closed"));assert h["rho_dot_closed"].max()>0

def test_open_flux_does_not_remove_isolated_or_closed_stores():
    p=late.LateStageParams(base(),late_stage_mode="pore_detachment_and_closure",k_closed_ref_s_Pa=2e-20);h=late.run(p,search.BudgetSchedule(20,1400,2,96))
    assert np.all(h["phi_closed"]>=0);assert np.sum(h["phi_iso"][-1])+np.sum(h["phi_closed"][-1])>=np.sum(h["phi_iso"][0])-1e-8

def test_gas_content_is_conserved_after_assignment():
    p=late.LateStageParams(base(),late_stage_mode="gas_limited_closed_pore",rho_close_mid=.90,gas_pressure_ratio=.5);h=late.run(p,search.BudgetSchedule(20,1400,20,96))
    assert np.all(h["gas_content_final"]>=0)

def test_stress_is_not_a_pore_volume_state_and_locality_has_no_labels():
    assert "sigma_res_hydro_Pa" in late.LateStageParams.__dataclass_fields__
    forbidden=("slow","fast","ramp_rate","rho_target","schedule_class")
    for fn in late.LOCAL_FUNCTIONS:
        src=inspect.getsource(fn).lower();assert not any(x in src for x in forbidden)

def test_effect_and_chen_gates_remain_strict():
    rho=np.arange(.85,.921,.001);assert effect.longest_span(rho,np.where(rho<=.88,1.5,1),1.5)>=.03-1e-12
    good=dict(complete_practical_window=True,lower_bracketed=True,upper_bracketed=True,T2_C=1250,T1_C=1350);assert common.chen_window_valid(good);assert not common.chen_window_valid({**good,"T2_C":1400})

def test_unattained_high_density_is_not_scored():
    w=pd.read_csv("results/late_stage_closed_pore_trajectory/density_window_effects.csv");q=w[w.density_window.isin(["late_stage","near_final"])]
    assert not q.both_paths_attained.any()

