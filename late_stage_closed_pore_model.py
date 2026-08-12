#!/usr/bin/env python3
"""Explicit late-stage closed-pore extension of the frozen PR open-pore model."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

import pr_desintering_memory_model as memory

MODES=("disabled","closed_pore_vacancy_transport","gas_limited_closed_pore","pore_detachment_and_closure","combined_late_stage")
COUPLINGS=("disabled","assist_shrinkage","oppose_or_crack_like","mixed_signed")

@dataclass
class LateStageParams:
    base:memory.PRMemoryParams;late_stage_mode:str="disabled"
    rho_close_mid:float=.92;rho_close_width:float=.02;closure_rate_s:float=2e-5
    k_closed_ref_s_Pa:float=2e-14;Q_closed_J_mol:float=400e3
    T_closed_ref_C:float=1400.;q_closed:float=2.;gas_pressure_ratio:float=.25
    residual_closed_pore_coupling:str="disabled";alpha_sigma:float=.25
    sigma_res_hydro_Pa:float=0.;pore_detachment_rate_s:float=2e-6
    max_steps:int=5000

def validate(p):
    if p.late_stage_mode not in MODES:raise ValueError("invalid late_stage_mode")
    if p.residual_closed_pore_coupling not in COUPLINGS:raise ValueError("invalid residual coupling")
    if min(p.rho_close_width,p.k_closed_ref_s_Pa,p.closure_rate_s)<=0:raise ValueError("positive rates/width required")

def sigmoid(x):return float(1/(1+math.exp(-float(np.clip(x,-60,60)))))

def closed_flux(phi_closed,radii,T_C,p,gas_content=None,sigma_res=None):
    """Named closed-pore shrinkage; zero without closed pore volume."""
    validate(p);phi=np.maximum(np.asarray(phi_closed,float),0);r=np.asarray(radii,float)
    zero=np.zeros_like(phi)
    if p.late_stage_mode=="disabled" or not np.any(phi):return zero,zero,zero
    T=T_C+273.15;Tr=p.T_closed_ref_C+273.15;k=p.k_closed_ref_s_Pa*math.exp(float(np.clip(-p.Q_closed_J_mol/memory.R*(1/T-1/Tr),-40,40)))
    cap=2*p.base.base.action.location.base.gamma_s/np.maximum(r,1e-30)
    gas=np.zeros_like(phi) if gas_content is None else np.asarray(gas_content)/np.maximum(phi,1e-300)
    sig=0 if sigma_res is None else np.asarray(sigma_res)
    if p.residual_closed_pore_coupling=="assist_shrinkage":stress=np.maximum(sig,0)*p.alpha_sigma
    elif p.residual_closed_pore_coupling=="oppose_or_crack_like":stress=-np.abs(sig)*p.alpha_sigma
    elif p.residual_closed_pore_coupling=="mixed_signed":stress=sig*p.alpha_sigma
    else:stress=0
    drive=np.maximum(cap-gas+stress,0);transport=(p.base.base.action.location.base.pore_radius0/np.maximum(r,1e-30))**p.q_closed
    loss=k*phi*drive*transport
    return loss,gas,cap

def run(p,protocol):
    validate(p)
    if p.late_stage_mode=="disabled":
        out=memory.run(p.base,protocol);out["phi_closed"]=np.zeros_like(out["phi_iso"]);out["N_closed"]=np.zeros_like(out["N_iso"]);out["rho_dot_closed"]=np.zeros(len(out["rho"]));out["rho_dot_open"]=out["rho_dot"].copy();out["numerical_censored"]=False;return out
    s=memory.initial_state(p.base);pore=s.base.pore;r=pore.pore_radii;closed=np.zeros_like(r);gas_content=np.zeros_like(r);sigma=np.full_like(r,p.sigma_res_hydro_Pa)
    keys="t T_C rho G rho_dot rho_dot_open rho_dot_closed closure_gate closed_pore_fraction closed_mean_radius closed_D90 gas_capillary_ratio sigma_res_hydro cumulative_PR_desintering_work".split();h={k:[] for k in keys};h.update(phi_GBseg=[],phi_TJ=[],phi_iso=[],phi_closed=[],N_GBseg=[],N_TJ=[],N_iso=[],N_closed=[])
    lp=memory.action.effective_location_params(p.base.base.action);steps=0
    while pore.t<min(protocol.t_end,lp.base.t_max_s) and pore.rho<lp.base.rho_cap and steps<p.max_steps:
        steps+=1;T_C=protocol.T(pore.t,pore.rho);d=memory.local_competition(s,T_C,p.base);total_pore=float(np.sum(pore.phi_total)+np.sum(closed));iso_frac=float(np.sum(pore.phi_iso))/max(total_pore,1e-300);connected=float(np.sum(pore.phi_GBseg+pore.phi_TJ))/max(total_pore,1e-300);gate=sigmoid((pore.rho-p.rho_close_mid)/p.rho_close_width)*iso_frac*(1-connected)
        loss,pgas,pcap=closed_flux(closed,r,T_C,p,gas_content,sigma);rho_closed=float(np.sum(loss));rho_open=float(d["rho_dot"])
        z=max(float(np.sum(closed)),1e-300);cdf=np.cumsum(closed)/z;d90=float(r[min(np.searchsorted(cdf,.9),len(r)-1)]) if np.sum(closed)>0 else 0
        vals=dict(t=pore.t,T_C=T_C,rho=pore.rho,G=pore.G,rho_dot=rho_open+rho_closed,rho_dot_open=rho_open,rho_dot_closed=rho_closed,closure_gate=gate,closed_pore_fraction=float(np.sum(closed))/max(total_pore,1e-300),closed_mean_radius=float(np.sum(closed*r))/z if np.sum(closed)>0 else 0,closed_D90=d90,gas_capillary_ratio=float(np.sum(closed*pgas)/max(np.sum(closed*pcap),1e-300)),sigma_res_hydro=float(np.mean(sigma)),cumulative_PR_desintering_work=s.cumulative_PR_desintering_work)
        for k in keys:h[k].append(vals[k])
        for k in ("phi_GBseg","phi_TJ","phi_iso","N_GBseg","N_TJ","N_iso"):h[k].append(getattr(pore,k).copy())
        h["phi_closed"].append(closed.copy());h["N_closed"].append(memory.location._number(closed,r))
        dt=min(lp.base.dt_max_s,protocol.t_end-pore.t);dT=abs(protocol.T(pore.t+1,pore.rho)-T_C)
        if dT:dt=min(dt,lp.base.dT_max_C/dT)
        if rho_open+rho_closed>0:dt=min(dt,lp.base.drho_max/(rho_open+rho_closed))
        if d["G_dot"]>0:dt=min(dt,lp.base.dG_fraction_max*pore.G/d["G_dot"])
        # Open removal never acts on phi_iso/phi_closed. Closure and detachment
        # are conservative transfers; only `loss` changes total pore volume.
        outgb=-d["GBseg_remove"]+d["GB_smooth"]-d["GB_to_TJ"]+d["TJ_to_GBseg_capture"]
        outtj=-d["TJ_remove"]+d["GB_to_TJ"]-d["TJ_to_GBseg_capture"]-d["TJ_to_iso"]
        outiso=d["TJ_to_iso"]
        close=p.closure_rate_s*gate*pore.phi_iso
        if p.late_stage_mode in ("pore_detachment_and_closure","combined_late_stage"):
            detach=p.pore_detachment_rate_s*gate*(pore.phi_GBseg+pore.phi_TJ);outgb-=detach*pore.phi_GBseg/np.maximum(pore.phi_GBseg+pore.phi_TJ,1e-300);outtj-=detach*pore.phi_TJ/np.maximum(pore.phi_GBseg+pore.phi_TJ,1e-300);close+=detach
        maxrate=max(float(np.max(close/np.maximum(pore.phi_iso+close,1e-300))),float(np.max(loss/np.maximum(closed,1e-300))))
        if maxrate>0:dt=min(dt,.2/maxrate)
        dt=max(lp.base.dt_min_s,dt);transfer=np.minimum(close*dt,pore.phi_iso+np.maximum(close-p.closure_rate_s*gate*pore.phi_iso,0)*dt);newclosed=np.maximum(closed+transfer-loss*dt,0)
        if not np.any(gas_content) and p.gas_pressure_ratio>0:gas_content+=transfer*p.gas_pressure_ratio*(2*p.base.base.action.location.base.gamma_s/np.maximum(r,1e-30))
        pore.phi_GBseg=np.maximum(pore.phi_GBseg+outgb*dt,0);pore.phi_TJ=np.maximum(pore.phi_TJ+outtj*dt,0);pore.phi_iso=np.maximum(pore.phi_iso+outiso*dt-transfer,0);closed=newclosed;pore.rho=1-float(np.sum(pore.phi_total)+np.sum(closed));pore.N_GBseg=memory.location._number(pore.phi_GBseg,r);pore.N_TJ=memory.location._number(pore.phi_TJ,r);pore.N_iso=memory.location._number(pore.phi_iso,r);pore.G=max(pore.G+d["G_dot"]*dt,1e-9);s.base.X_J=float(np.clip(s.base.X_J+d["X_J_dot"]*dt,0,p.base.base.XJ_capacity));s.cumulative_PR_desintering_work+=d["PR_work_dot"]*dt;sigma*=math.exp(-dt/2e5);pore.t+=dt
    out={k:np.asarray(v,float) for k,v in h.items()};out["pore_radii"]=r.copy();out["gas_content_final"]=gas_content.copy();out["numerical_censored"]=steps>=p.max_steps and pore.t<min(protocol.t_end,lp.base.t_max_s);return out

LOCAL_FUNCTIONS=(closed_flux,)
