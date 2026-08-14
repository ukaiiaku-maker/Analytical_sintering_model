#!/usr/bin/env python3
"""Ablatable cohort-level residual-stress coupling to local PR competition."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

import pr_desintering_memory_model as memory

MODES=("disabled","initial_only","evolving_PR_generated","evolving_shear_generated","mixed_evolving")


@dataclass
class ResidualStressParams:
    mode:str="disabled";sigma_res_scale:float=0.;stress_sign:str="mixed"
    sigma_ref_Pa:float=5e7;alpha_comp:float=.35;alpha_tensile_penalty:float=.25
    beta_PR:float=1.;beta_crack:float=1.5;A_PR_stress:float=2e-3
    A_shear_stress:float=2e-7;k_relax_dens:float=8.;tau_res_ref_s:float=2e4
    Q_res_J_mol:float=180e3;T_ref_C:float=1100.;large_defect_fraction:float=.02;max_steps:int=1000
    large_defect_radius_factor:float=8.;crack_like_aspect_proxy:float=4.;defect_stress_concentration:float=2.


@dataclass
class ResidualStressState:
    sigma_res_GBseg:float=0.;sigma_res_TJ:float=0.;sigma_res_large_pore:float=0.;sigma_res_crack_like:float=0.


def initial_state(p:ResidualStressParams,cohort_sign=0.):
    if p.mode not in MODES:raise ValueError("invalid residual_stress_mode")
    amp=p.sigma_res_scale*p.sigma_ref_Pa;sign=cohort_sign
    if p.stress_sign=="compressive":return ResidualStressState(amp,.5*amp,0,0)
    if p.stress_sign=="tensile":return ResidualStressState(0,0,amp,p.defect_stress_concentration*amp)
    return ResidualStressState(.5*amp,.25*amp,max(sign,.5)*amp,max(sign,.5)*p.defect_stress_concentration*amp)


def local_residual_coupling(d,state:ResidualStressState,T_C,p:ResidualStressParams):
    """Modify hazards/flux allocation locally; schedule identity is absent."""
    if p.mode=="disabled":return {**d,"residual_densification_factor":1.,"residual_PR_factor":1.,"residual_defect_flux":0.}
    comp=max(state.sigma_res_GBseg+state.sigma_res_TJ,0)/p.sigma_ref_Pa
    tensile=max(state.sigma_res_large_pore+state.sigma_res_crack_like,0)/p.sigma_ref_Pa
    dens=max(0.,1+p.alpha_comp*comp-p.alpha_tensile_penalty*tensile)
    pr=math.exp(float(np.clip(p.beta_PR*tensile,-10,10)))
    # Densifying removal stays a densification channel. Residual stress itself
    # neither adds nor removes pore volume.
    d={**d};d["GBseg_remove"]=d["GBseg_remove"]*dens;d["TJ_remove"]=d["TJ_remove"]*dens;d["rho_dot"]=float(np.sum(d["GBseg_remove"]+d["TJ_remove"]))
    # Conservative rare-defect transfer from connected large bins toward the
    # next bin; zero sum by construction.
    smooth=d["GB_smooth"].copy();source=max(p.large_defect_fraction*(pr-1),0)*1e-5
    defect=np.zeros_like(smooth);move=source*np.maximum(smooth[:-1]*0+1,0)
    if len(move):move*=max(d.get("large_pore_fraction",0),0)/len(move);defect[:-1]-=move;defect[1:]+=move
    d["GB_smooth"]=smooth+defect
    return {**d,"residual_densification_factor":dens,"residual_PR_factor":pr,"residual_defect_flux":float(np.sum(np.maximum(defect,0)))}


def derivatives(state,d,T_C,p):
    if p.mode in ("disabled","initial_only"):gen_pr=gen_sh=0.
    else:
        gen_pr=p.A_PR_stress*d.get("PR_work_dot",0) if p.mode in ("evolving_PR_generated","mixed_evolving") else 0.
        gen_sh=p.A_shear_stress*(d.get("P_clean_GB",0)+d.get("P_TJ_multihit",0)) if p.mode in ("evolving_shear_generated","mixed_evolving") else 0.
    T=T_C+273.15;Tr=p.T_ref_C+273.15;tau=p.tau_res_ref_s*math.exp(float(np.clip(p.Q_res_J_mol/memory.R*(1/T-1/Tr),-30,30)))
    relax=p.k_relax_dens*max(d.get("rho_dot",0),0)+1/max(tau,1e-30)
    return ResidualStressState(gen_sh-relax*state.sigma_res_GBseg,gen_sh-relax*state.sigma_res_TJ,gen_pr-relax*state.sigma_res_large_pore,p.defect_stress_concentration*gen_pr-relax*state.sigma_res_crack_like)


def update(state,deriv,dt):
    for name in state.__dataclass_fields__:setattr(state,name,float(np.clip(getattr(state,name)+getattr(deriv,name)*dt,-1e10,1e10)))


def run(base_params:memory.PRMemoryParams,protocol,p:ResidualStressParams,cohort_sign=0.):
    """Integrate the existing PR model with only the declared local coupling."""
    if p.mode=="disabled":return memory.run(base_params,protocol)
    s=memory.initial_state(base_params);rs=initial_state(p,cohort_sign)
    scalar="t T_C rho G connected_fine_pore_fraction pore_mean_radius large_pore_fraction cumulative_PR_desintering_work residual_densification_factor residual_PR_factor residual_defect_flux sigma_res_GBseg sigma_res_TJ sigma_res_large_pore sigma_res_crack_like PR_work_dot rho_dot G_dot".split()
    h={k:[] for k in scalar};h.update(phi_GBseg=[],phi_TJ=[],phi_iso=[],N_GBseg=[],N_TJ=[],N_iso=[])
    lp=memory.action.effective_location_params(base_params.base.action)
    steps=0
    while s.base.pore.t<min(protocol.t_end,lp.base.t_max_s) and s.base.pore.rho<lp.base.rho_cap and steps<p.max_steps:
        steps+=1
        pore=s.base.pore;T_C=protocol.T(pore.t,pore.rho);raw=memory.local_competition(s,T_C,base_params);d=local_residual_coupling(raw,rs,T_C,p)
        values={**d,"t":pore.t,"T_C":T_C,"rho":pore.rho,"G":pore.G,"cumulative_PR_desintering_work":s.cumulative_PR_desintering_work,
                "sigma_res_GBseg":rs.sigma_res_GBseg,"sigma_res_TJ":rs.sigma_res_TJ,"sigma_res_large_pore":rs.sigma_res_large_pore,"sigma_res_crack_like":rs.sigma_res_crack_like}
        for k in scalar:h[k].append(values[k])
        for k in ("phi_GBseg","phi_TJ","phi_iso","N_GBseg","N_TJ","N_iso"):h[k].append(getattr(pore,k).copy())
        dt=min(lp.base.dt_max_s,protocol.t_end-pore.t);dT=abs(protocol.T(pore.t+1,pore.rho)-T_C)
        if dT:dt=min(dt,lp.base.dT_max_C/dT)
        if d["rho_dot"]>0:dt=min(dt,lp.base.drho_max/d["rho_dot"])
        if d["G_dot"]>0:dt=min(dt,lp.base.dG_fraction_max*pore.G/d["G_dot"])
        outgb=-d["GBseg_remove"]+d["GB_smooth"]-d["GB_to_TJ"]+d["TJ_to_GBseg_capture"]
        outtj=-d["TJ_remove"]+d["GB_to_TJ"]-d["TJ_to_GBseg_capture"]-d["TJ_to_iso"];outiso=d["TJ_to_iso"]
        loss=max(float(np.max(np.maximum(-outgb,0)/np.maximum(pore.phi_GBseg,1e-300))),float(np.max(np.maximum(-outtj,0)/np.maximum(pore.phi_TJ,1e-300))))
        if loss>0:dt=min(dt,.2/loss)
        dt=max(lp.base.dt_min_s,dt)
        pore.phi_GBseg=np.maximum(pore.phi_GBseg+outgb*dt,0);pore.phi_TJ=np.maximum(pore.phi_TJ+outtj*dt,0);pore.phi_iso=np.maximum(pore.phi_iso+outiso*dt,0)
        pore.rho=1-float(np.sum(pore.phi_total));pore.N_GBseg=memory.location._number(pore.phi_GBseg,pore.pore_radii);pore.N_TJ=memory.location._number(pore.phi_TJ,pore.pore_radii);pore.N_iso=memory.location._number(pore.phi_iso,pore.pore_radii)
        pore.G=max(pore.G+d["G_dot"]*dt,1e-9);s.base.X_J=float(np.clip(s.base.X_J+d["X_J_dot"]*dt,0,base_params.base.XJ_capacity));s.cumulative_PR_desintering_work+=d["PR_work_dot"]*dt
        update(rs,derivatives(rs,d,T_C,p),dt);pore.t+=dt
    out={k:np.asarray(v,float) for k,v in h.items()};out["pore_radii"]=s.base.pore.pore_radii.copy();out["numerical_censored"]=steps>=p.max_steps and s.base.pore.t<min(protocol.t_end,lp.base.t_max_s);return out


LOCAL_FUNCTIONS=(local_residual_coupling,derivatives)
