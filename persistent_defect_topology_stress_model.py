#!/usr/bin/env python3
"""Persistent PR-created defect/topology-stress memory below rho=0.92."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

import pr_desintering_memory_model as memory
import residual_stress_memory_model as residual

MODES=("disabled","persistent_defect_memory")


@dataclass
class PersistentParams:
    mode:str="disabled";defect_generation_per_flux:float=12.;defect_decay_time_s:float=8e5
    sink_relaxation_strength:float=1.5;eligibility_damage_strength:float=1.8
    defect_coarsening_strength:float=3.;stored_work_scale:float=8e5
    stress_storage_fraction:float=.12;stress_release_strength:float=2.


@dataclass
class PersistentState:
    f_defect_large_pore:float=0.;f_crack_like_pore:float=0.;defect_D90:float=0.
    defect_connectedness:float=1.;stored_PR_work:float=0.;stored_shear_coupled_stress:float=0.


def local_persistent_coupling(d,state:PersistentState,T_C,p:PersistentParams):
    """Instantaneous state-only coupling; schedule identity is absent."""
    if p.mode=="disabled":return {**d,"persistent_eligibility":1.,"persistent_growth_factor":1.}
    damage=float(np.clip(state.f_defect_large_pore+state.f_crack_like_pore,0,1))
    stress_gate=1+state.stored_shear_coupled_stress/max(p.stored_work_scale,1e-30)
    eligibility=math.exp(-p.eligibility_damage_strength*damage*stress_gate)
    d={**d};d["GBseg_remove"]=d["GBseg_remove"]*eligibility;d["TJ_remove"]=d["TJ_remove"]*eligibility;d["rho_dot"]=float(np.sum(d["GBseg_remove"]+d["TJ_remove"]))
    thermal=1/(1+math.exp(-(T_C-1150)/80));growth=1+p.defect_coarsening_strength*damage*thermal
    d["G_dot"]=d["G_dot"]*growth
    return {**d,"persistent_eligibility":eligibility,"persistent_growth_factor":growth}


def derivatives(state,d,p:PersistentParams):
    if p.mode=="disabled":return PersistentState()
    pr=max(d.get("PR_desintering_flux",0),0);dens=max(d.get("rho_dot",0),0);shear=max(d.get("P_clean_GB",0)+d.get("P_TJ_multihit",0),0)
    generation=p.defect_generation_per_flux*pr*(1-state.f_defect_large_pore)
    relax=(1/p.defect_decay_time_s+p.sink_relaxation_strength*dens*state.defect_connectedness)
    df=generation-relax*state.f_defect_large_pore;dc=.15*generation-.5*relax*state.f_crack_like_pore
    dwork=max(d.get("PR_work_dot",0),0)-relax*state.stored_PR_work
    dstress=p.stress_storage_fraction*(max(d.get("PR_work_dot",0),0)+1e-6*shear)-p.stress_release_strength*dens*state.stored_shear_coupled_stress-state.stored_shear_coupled_stress/p.defect_decay_time_s
    return PersistentState(df,dc,generation*max(d.get("pore_mean_radius",0),0)-relax*state.defect_D90,-.2*generation*state.defect_connectedness+.1*dens*(1-state.defect_connectedness),dwork,dstress)


def update(s,ds,dt):
    s.f_defect_large_pore=float(np.clip(s.f_defect_large_pore+ds.f_defect_large_pore*dt,0,1));s.f_crack_like_pore=float(np.clip(s.f_crack_like_pore+ds.f_crack_like_pore*dt,0,1-s.f_defect_large_pore));s.defect_D90=max(s.defect_D90+ds.defect_D90*dt,0);s.defect_connectedness=float(np.clip(s.defect_connectedness+ds.defect_connectedness*dt,0,1));s.stored_PR_work=max(s.stored_PR_work+ds.stored_PR_work*dt,0);s.stored_shear_coupled_stress=max(s.stored_shear_coupled_stress+ds.stored_shear_coupled_stress*dt,0)


def run(base_params,protocol,residual_params,p:PersistentParams,cohort_sign=0.):
    if p.mode=="disabled":return residual.run(base_params,protocol,residual_params,cohort_sign)
    s=memory.initial_state(base_params);rs=residual.initial_state(residual_params,cohort_sign);ps=PersistentState(defect_D90=base_params.base.action.location.base.pore_radius0*4)
    scalar="t T_C rho G connected_fine_pore_fraction pore_mean_radius large_pore_fraction cumulative_PR_desintering_work sigma_res_GBseg sigma_res_TJ sigma_res_large_pore sigma_res_crack_like f_defect_large_pore f_crack_like_pore defect_D90 defect_connectedness stored_PR_work stored_shear_coupled_stress persistent_eligibility persistent_growth_factor PR_work_dot rho_dot G_dot residual_defect_flux".split();h={k:[] for k in scalar};h.update(phi_GBseg=[],phi_TJ=[],phi_iso=[],N_GBseg=[],N_TJ=[],N_iso=[]);lp=memory.action.effective_location_params(base_params.base.action)
    while s.base.pore.t<min(protocol.t_end,lp.base.t_max_s) and s.base.pore.rho<lp.base.rho_cap:
        pore=s.base.pore;T_C=protocol.T(pore.t,pore.rho);d=residual.local_residual_coupling(memory.local_competition(s,T_C,base_params),rs,T_C,residual_params);d=local_persistent_coupling(d,ps,T_C,p)
        vals={**d,"t":pore.t,"T_C":T_C,"rho":pore.rho,"G":pore.G,"cumulative_PR_desintering_work":s.cumulative_PR_desintering_work,**rs.__dict__,**ps.__dict__}
        for k in scalar:h[k].append(vals[k])
        for k in ("phi_GBseg","phi_TJ","phi_iso","N_GBseg","N_TJ","N_iso"):h[k].append(getattr(pore,k).copy())
        dt=min(lp.base.dt_max_s,protocol.t_end-pore.t);dT=abs(protocol.T(pore.t+1,pore.rho)-T_C)
        if dT:dt=min(dt,lp.base.dT_max_C/dT)
        if d["rho_dot"]>0:dt=min(dt,lp.base.drho_max/d["rho_dot"])
        if d["G_dot"]>0:dt=min(dt,lp.base.dG_fraction_max*pore.G/d["G_dot"])
        outgb=-d["GBseg_remove"]+d["GB_smooth"]-d["GB_to_TJ"]+d["TJ_to_GBseg_capture"];outtj=-d["TJ_remove"]+d["GB_to_TJ"]-d["TJ_to_GBseg_capture"]-d["TJ_to_iso"];outiso=d["TJ_to_iso"]
        loss=max(float(np.max(np.maximum(-outgb,0)/np.maximum(pore.phi_GBseg,1e-300))),float(np.max(np.maximum(-outtj,0)/np.maximum(pore.phi_TJ,1e-300))))
        if loss>0:dt=min(dt,.2/loss)
        dt=max(lp.base.dt_min_s,dt);pore.phi_GBseg=np.maximum(pore.phi_GBseg+outgb*dt,0);pore.phi_TJ=np.maximum(pore.phi_TJ+outtj*dt,0);pore.phi_iso=np.maximum(pore.phi_iso+outiso*dt,0);pore.rho=1-float(np.sum(pore.phi_total));pore.N_GBseg=memory.location._number(pore.phi_GBseg,pore.pore_radii);pore.N_TJ=memory.location._number(pore.phi_TJ,pore.pore_radii);pore.N_iso=memory.location._number(pore.phi_iso,pore.pore_radii);pore.G=max(pore.G+d["G_dot"]*dt,1e-9);s.base.X_J=float(np.clip(s.base.X_J+d["X_J_dot"]*dt,0,base_params.base.XJ_capacity));s.cumulative_PR_desintering_work+=d["PR_work_dot"]*dt;residual.update(rs,residual.derivatives(rs,d,T_C,residual_params),dt);update(ps,derivatives(ps,d,p),dt);pore.t+=dt
    out={k:np.asarray(v,float) for k,v in h.items()};out["pore_radii"]=s.base.pore.pore_radii.copy();return out


LOCAL_FUNCTIONS=(local_persistent_coupling,derivatives)
