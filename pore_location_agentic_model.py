#!/usr/bin/env python3
"""Local constrained topology-action layer for the pore-placement model."""
from __future__ import annotations
from dataclasses import dataclass,replace
import math
from typing import Optional

import numpy as np

import pore_location_topology_model as fixed
import topology_constrained_sintering as aggregate

R=aggregate.R


@dataclass
class TopologyAction:
    name:str;source:str;destination:str;conservative:bool;affects_density:bool;affects_growth:bool
    propensity:float;driving_power:float;resistance:float;diagnostics:dict


@dataclass
class ActionParams:
    location:fixed.LocationParams
    action_mode:str='fixed_flux_baseline'
    w_GBseg_shrink:float=1.
    w_TJ_shrink:float=.25
    w_GBseg_smooth:float=.5
    w_GBseg_to_TJ:float=.7
    w_TJ_capture:float=.8
    w_TJ_iso:float=.7
    w_clean_GB_growth:float=1.
    k_capture_0_s:float=8e4
    Q_capture_J_mol:float=325e3
    capture_fine_exp:float=.5
    high_TJ_drag_factor:float=2.0


ACTION_NAMES=('GBseg_shrink','TJ_shrink','GBseg_smooth','GBseg_to_TJ','TJ_to_GBseg_capture','TJ_to_iso','clean_GB_growth','GBseg_pin_drag','TJ_pin_drag')


def validate(p:ActionParams):
    allowed={'fixed_flux_baseline','action_static','action_evolving_no_capture','action_evolving_capture','action_evolving_capture_high_TJ_drag'}
    if p.action_mode not in allowed:raise ValueError(f'invalid action_mode {p.action_mode!r}')
    for name in ('w_GBseg_shrink','w_TJ_shrink','w_GBseg_smooth','w_GBseg_to_TJ','w_TJ_capture','w_TJ_iso','w_clean_GB_growth'):
        if getattr(p,name)<0:raise ValueError(f'{name} must be nonnegative')


def effective_location_params(p:ActionParams):
    if p.action_mode=='fixed_flux_baseline':return replace(p.location,pore_location_mode='evolving')
    if p.action_mode=='action_static':return replace(p.location,pore_location_mode='static')
    q=replace(p.location,pore_location_mode='evolving')
    if p.action_mode=='action_evolving_capture_high_TJ_drag':q=replace(q,A_TJ_pin=q.A_TJ_pin*p.high_TJ_drag_factor)
    return q


def _thermal_ratio(Q,T_C,Tref_C=1250.):
    T=T_C+273.15;Tr=Tref_C+273.15;return math.exp(float(np.clip(-Q/R*(1/T-1/Tr),-30,30)))


def score_actions(s:fixed.PoreLocationState,T_C:float,p:ActionParams)->tuple[dict[str,TopologyAction],dict]:
    """Score admissible instantaneous actions; no path identity enters."""
    validate(p);lp=effective_location_params(p);d=fixed.instantaneous(s,T_C,lp);top=fixed.topology_diagnostics(s,lp);total=max(float(np.sum(s.phi_total)),1e-300)
    fgb=float(np.sum(s.phi_GBseg))/total;ftj=float(np.sum(s.phi_TJ))/total;fiso=float(np.sum(s.phi_iso))/total
    capillary=4*lp.base.gamma_s/max(s.G,1e-30);kin_res=max(d['tau_event'],1e-300);mig_res=1/max(d['growth_mobility_factor'],1e-12)
    capture_allowed=p.action_mode in ('action_evolving_capture','action_evolving_capture_high_TJ_drag')
    capture_gate=top['C_TJ']*top['f_clean_GB']*d['growth_mobility_factor'];thermal_capture=_thermal_ratio(p.Q_capture_J_mol,T_C)
    iso_gate=aggregate.sig((s.rho-.86)/.025)*(1-top['C_GBseg']+.1)
    scores={
      'GBseg_shrink':p.w_GBseg_shrink*fgb*top['C_GBseg']*d['activity']/kin_res,
      'TJ_shrink':p.w_TJ_shrink*ftj*top['C_TJ']*d['activity']/kin_res,
      'GBseg_smooth':p.w_GBseg_smooth*fgb*(1-d['activity'])*_thermal_ratio(lp.Q_GB_smooth_J_mol,T_C),
      'GBseg_to_TJ':p.w_GBseg_to_TJ*fgb*top['f_clean_GB']*(1-d['activity'])*_thermal_ratio(lp.Q_GB_to_TJ_J_mol,T_C),
      'TJ_to_GBseg_capture':p.w_TJ_capture*ftj*capture_gate*thermal_capture if capture_allowed else 0.,
      'TJ_to_iso':p.w_TJ_iso*ftj*iso_gate*_thermal_ratio(lp.Q_TJ_iso_J_mol,T_C),
      'clean_GB_growth':p.w_clean_GB_growth*top['f_clean_GB']*d['growth_mobility_factor'],
      'GBseg_pin_drag':top['C_GBseg']*d['R_GBseg_pin'],
      'TJ_pin_drag':top['C_TJ']*d['R_TJ_pin'],
    }
    resist={'GBseg_shrink':kin_res,'TJ_shrink':kin_res/max(lp.eta_TJ_dens,1e-12),'GBseg_smooth':1/max(_thermal_ratio(lp.Q_GB_smooth_J_mol,T_C),1e-12),'GBseg_to_TJ':1/max(_thermal_ratio(lp.Q_GB_to_TJ_J_mol,T_C),1e-12),'TJ_to_GBseg_capture':1/max(thermal_capture,1e-12),'TJ_to_iso':1/max(_thermal_ratio(lp.Q_TJ_iso_J_mol,T_C),1e-12),'clean_GB_growth':mig_res,'GBseg_pin_drag':1+ d['R_GBseg_pin'],'TJ_pin_drag':1+d['R_TJ_pin']}
    meta={'GBseg_shrink':('GBseg','removed',False,True,False),'TJ_shrink':('TJ','removed',False,True,False),'GBseg_smooth':('GBseg_i','GBseg_i+1',True,False,False),'GBseg_to_TJ':('GBseg','TJ',True,False,False),'TJ_to_GBseg_capture':('TJ','GBseg',True,False,False),'TJ_to_iso':('TJ','iso',True,False,False),'clean_GB_growth':('clean_GB','larger_G',False,False,True),'GBseg_pin_drag':('GBseg','GBseg',True,False,False),'TJ_pin_drag':('TJ','TJ',True,False,False)}
    actions={}
    for name in ACTION_NAMES:
        source,dest,cons,den,grow=meta[name];actions[name]=TopologyAction(name,source,dest,cons,den,grow,max(scores[name],0),max(scores[name],0)*capillary,resist[name],{'availability':{'GBseg':fgb,'TJ':ftj,'iso':fiso}.get(source.split('_')[0],1.),'capillary_drive':capillary})
    total_score=max(sum(a.propensity for a in actions.values()),1e-300);diag={f'action_propensity_{n}':a.propensity for n,a in actions.items()};diag.update({f'action_weight_{n}':a.propensity/total_score for n,a in actions.items()})
    return actions,{**diag,**d}


def allocated_fluxes(s:fixed.PoreLocationState,T_C:float,p:ActionParams)->dict:
    actions,d=score_actions(s,T_C,p);lp=effective_location_params(p);r=s.pore_radii;rref=lp.base.pore_radius0
    def pair(a,b):
        x=actions[a].propensity;y=actions[b].propensity;z=x+y
        return (x/z,y/z) if z>0 else (0.,0.)
    shrink_gb,shrink_tj=pair('GBseg_shrink','TJ_shrink');rho_dot=d['rho_dot'];gbw=s.phi_GBseg*(rref/r)**lp.q_GBseg;tjw=s.phi_TJ*(rref/r)**lp.q_TJ
    gb_remove=rho_dot*shrink_gb*gbw/max(float(np.sum(gbw)),1e-300);tj_remove=rho_dot*shrink_tj*tjw/max(float(np.sum(tjw)),1e-300)
    gb_s_share,gb_tj_share=pair('GBseg_smooth','GBseg_to_TJ');gb_capacity=float(np.sum(np.maximum(d['GB_to_TJ'],0)))+float(np.sum(np.maximum(-d['GB_smooth'],0)))
    gb_to_tj=gb_capacity*gb_tj_share*s.phi_GBseg/max(float(np.sum(s.phi_GBseg)),1e-300)
    smooth_transfer=s.phi_GBseg[:-1]*(rref/r[:-1])**2;smooth_transfer*=gb_capacity*gb_s_share/max(float(np.sum(smooth_transfer)),1e-300) if gb_capacity*gb_s_share>0 else 0
    gb_smooth=np.zeros_like(r);gb_smooth[:-1]-=smooth_transfer;gb_smooth[1:]+=smooth_transfer
    capture_share,iso_share=pair('TJ_to_GBseg_capture','TJ_to_iso')
    capture_nominal=p.k_capture_0_s*math.exp(float(np.clip(-p.Q_capture_J_mol/(R*(T_C+273.15)),-700,700)))*d['C_TJ']*d['f_clean_GB']*d['growth_mobility_factor']*s.phi_TJ*(rref/r)**p.capture_fine_exp
    tj_capacity=float(np.sum(capture_nominal))+float(np.sum(np.maximum(d['TJ_to_iso'],0)))
    capture=s.phi_TJ*(rref/r)**p.capture_fine_exp;capture*=tj_capacity*capture_share/max(float(np.sum(capture)),1e-300) if tj_capacity*capture_share>0 else 0
    isolate=s.phi_TJ*(r/rref)**.5;isolate*=tj_capacity*iso_share/max(float(np.sum(isolate)),1e-300) if tj_capacity*iso_share>0 else 0
    powers={'P_GBseg_dens':d['sigma_act_total']*rho_dot*shrink_gb,'P_TJ_dens':d['sigma_act_total']*rho_dot*shrink_tj,'P_clean_GB':d['P_clean_GB'],'P_GBseg_drag':d['P_GBseg_drag'],'P_TJ_drag':d['P_TJ_drag'],'P_GBseg_to_TJ':lp.base.gamma_s*float(np.sum(gb_to_tj/r)),'P_TJ_to_GBseg_capture':lp.base.gamma_s*float(np.sum(capture/r)),'P_TJ_iso':lp.base.gamma_s*float(np.sum(isolate/r))}
    flux={'GBseg_remove':gb_remove,'TJ_remove':tj_remove,'GB_smooth':gb_smooth,'GB_to_TJ':gb_to_tj,'TJ_to_GBseg_capture':capture,'TJ_to_iso':isolate}
    diagnostics={f'action_flux_{k}':float(np.sum(np.abs(v))) for k,v in flux.items()};return {**d,**powers,**flux,**diagnostics,'rho_dot':rho_dot*(shrink_gb+shrink_tj)}


def run(p:ActionParams,protocol,stop_at_rho:Optional[float]=None,initial:Optional[fixed.PoreLocationState]=None):
    validate(p);lp=effective_location_params(p)
    if p.action_mode in ('fixed_flux_baseline','action_static'):return fixed.run(lp,protocol,stop_at_rho=stop_at_rho,initial=initial)
    s=fixed.initial_state(lp) if initial is None else fixed.clone_state(initial);action_keys=[f'action_{kind}_{n}' for kind in ('propensity','weight') for n in ACTION_NAMES];flux_keys=[f'action_flux_{n}' for n in ('GBseg_remove','TJ_remove','GB_smooth','GB_to_TJ','TJ_to_GBseg_capture','TJ_to_iso')]
    scalar='t T_C rho G C_GBseg C_TJ f_clean_GB f_iso sigma_GBseg_pore sigma_TJ_pore sigma_clean_GB sigma_iso sigma_act_total activity rho_dot G_dot E_G growth_mobility_factor'.split();power='P_GBseg_dens P_TJ_dens P_clean_GB P_GBseg_drag P_TJ_drag P_GBseg_to_TJ P_TJ_to_GBseg_capture P_TJ_iso'.split();h={k:[] for k in scalar+power+action_keys+flux_keys};h.update(phi_GBseg=[],phi_TJ=[],phi_iso=[],N_GBseg=[],N_TJ=[],N_iso=[])
    while s.t<min(protocol.t_end,lp.base.t_max_s) and s.rho<lp.base.rho_cap:
        T=protocol.T(s.t,s.rho);actions,d0=score_actions(s,T,p);d=allocated_fluxes(s,T,p);vals={'t':s.t,'T_C':T,'rho':s.rho,'G':s.G,'E_G':d['rho_dot']/(d['G_dot']/max(s.G,1e-30)+1e-30),**d0,**d}
        for k in scalar+power+action_keys+flux_keys:h[k].append(vals[k])
        for k in ('phi_GBseg','phi_TJ','phi_iso','N_GBseg','N_TJ','N_iso'):h[k].append(getattr(s,k).copy())
        if stop_at_rho is not None and s.rho>=stop_at_rho:break
        dt=min(lp.base.dt_max_s,protocol.t_end-s.t);dT=abs(protocol.T(s.t+1,s.rho)-T)
        if dT:dt=min(dt,lp.base.dT_max_C/dT)
        if d['rho_dot']>0:dt=min(dt,lp.base.drho_max/d['rho_dot'])
        if d['G_dot']>0:dt=min(dt,lp.base.dG_fraction_max*s.G/d['G_dot'])
        out_gb=-d['GBseg_remove']+d['GB_smooth']-d['GB_to_TJ']+d['TJ_to_GBseg_capture'];out_tj=-d['TJ_remove']+d['GB_to_TJ']-d['TJ_to_GBseg_capture']-d['TJ_to_iso'];out_iso=d['TJ_to_iso']
        loss=max(float(np.max(np.maximum(-out_gb,0)/np.maximum(s.phi_GBseg,1e-300))),float(np.max(np.maximum(-out_tj,0)/np.maximum(s.phi_TJ,1e-300))),0)
        if loss>0:dt=min(dt,.2/loss)
        dt=max(lp.base.dt_min_s,dt);s.phi_GBseg=np.maximum(s.phi_GBseg+out_gb*dt,0);s.phi_TJ=np.maximum(s.phi_TJ+out_tj*dt,0);s.phi_iso=np.maximum(s.phi_iso+out_iso*dt,0);s.rho=1-float(np.sum(s.phi_total));s.N_GBseg=fixed._number(s.phi_GBseg,s.pore_radii);s.N_TJ=fixed._number(s.phi_TJ,s.pore_radii);s.N_iso=fixed._number(s.phi_iso,s.pore_radii);s.G=max(s.G+d['G_dot']*dt,1e-9);s.t+=dt
    out={k:np.asarray(v,float) for k,v in h.items()};out['pore_radii']=s.pore_radii.copy();return out


LOCAL_ACTION_FUNCTIONS=(score_actions,allocated_fluxes)
