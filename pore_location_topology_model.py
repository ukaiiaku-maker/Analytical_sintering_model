#!/usr/bin/env python3
"""Reduced pore-size x pore-location topology model.

The disabled mode delegates exactly to the PR #2 aggregate model.  Static and
evolving modes resolve GB-segment, triple-junction, and isolated pore stores.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import inspect, math
from typing import Optional

import numpy as np
import topology_constrained_sintering as aggregate

R=aggregate.R


@dataclass
class LocationParams:
    base: aggregate.Params
    pore_location_mode: str="disabled"
    f_GBseg_init: float=.60
    f_TJ_init: float=.30
    f_iso_init: float=.10
    location_size_bias: str="neutral"
    chi_GB: float=1.8
    chi_TJ: float=1.4
    eta_TJ_dens: float=.25
    q_GBseg: float=2.
    q_TJ: float=2.
    k_GB_smooth_0_s: float=1.0e5
    Q_GB_smooth_J_mol: float=300e3
    k_GB_to_TJ_0_s: float=4.0e4
    Q_GB_to_TJ_J_mol: float=315e3
    k_TJ_iso_0_s: float=2.0e5
    Q_TJ_iso_J_mol: float=345e3
    k_TJ_smooth_0_s: float=5.0e4
    Q_TJ_smooth_J_mol: float=320e3
    A_GBseg_pin: float=12.
    A_TJ_pin: float=28.
    pin_Q_GB_J_mol: float=470e3
    pin_Q_TJ_J_mol: float=525e3
    pin_T_ref_C: float=1250.
    pin_G_ref: float=150e-9
    pin_size_exp_GB: float=1.
    pin_size_exp_TJ: float=2.
    K_GBseg: float=.45
    K_TJ: float=1.1
    K_clean: float=.08
    stress_cap: float=5e8


@dataclass
class PoreLocationState:
    rho: float
    G: float
    pore_radii: np.ndarray
    phi_GBseg: np.ndarray
    phi_TJ: np.ndarray
    phi_iso: np.ndarray
    N_GBseg: np.ndarray
    N_TJ: np.ndarray
    N_iso: np.ndarray
    t: float=0.

    @property
    def phi_total(self): return self.phi_GBseg+self.phi_TJ+self.phi_iso


def _number(phi,r): return aggregate.pore_number(phi,r)
def _arr(A,Q,T_C): return A*math.exp(float(np.clip(-Q/(R*(T_C+273.15)),-700,700)))
def _mean_radius(phi,r):
    total=float(np.sum(phi));return float(np.sum(phi*r)/total) if total>0 else math.nan


def validate_params(p:LocationParams)->None:
    if p.pore_location_mode not in {"disabled","static","evolving"}:raise ValueError("invalid pore_location_mode")
    fractions=np.array([p.f_GBseg_init,p.f_TJ_init,p.f_iso_init])
    if np.any(fractions<0) or abs(float(np.sum(fractions))-1)>1e-12:raise ValueError("initial location fractions must be nonnegative and sum to one")
    if not 0<=p.eta_TJ_dens<=1:raise ValueError("eta_TJ_dens must lie in [0,1]")


def initial_state(p:LocationParams)->PoreLocationState:
    validate_params(p);base=aggregate.initial_state(p.base);r=base.pore_radii;total=base.pore_phi
    if p.location_size_bias=="neutral":gbw=tjw=isow=np.ones_like(r)
    elif p.location_size_bias=="separated":
        x=r/r[0];gbw=x**-1.5;tjw=np.ones_like(r);isow=x**1.5
    elif p.location_size_bias=="TJ_large":
        x=r/r[0];gbw=x**-1.;tjw=x;isow=x**1.5
    else:raise ValueError("invalid location_size_bias")
    def allocate(frac,w):
        raw=total*w;return frac*float(np.sum(total))*raw/max(float(np.sum(raw)),1e-300)
    gb,tj,iso=(allocate(p.f_GBseg_init,gbw),allocate(p.f_TJ_init,tjw),allocate(p.f_iso_init,isow))
    return PoreLocationState(1-float(np.sum(gb+tj+iso)),p.base.G0,r,gb,tj,iso,_number(gb,r),_number(tj,r),_number(iso,r))


def topology_diagnostics(s:PoreLocationState,p:LocationParams)->dict:
    A_GB=2/max(s.G,1e-30);L_TJ=6/max(s.G**2,1e-30)
    area_GB=float(np.sum(s.N_GBseg*math.pi*s.pore_radii**2))
    line_TJ=float(np.sum(s.N_TJ*2*math.pi*s.pore_radii))
    arg_GB=p.chi_GB*area_GB/max(A_GB,1e-300);arg_TJ=p.chi_TJ*line_TJ/max(L_TJ,1e-300)
    C_GB=1-math.exp(-max(arg_GB,0));C_TJ=1-math.exp(-max(arg_TJ,0))
    total=max(float(np.sum(s.phi_total)),1e-300);fiso=float(np.sum(s.phi_iso))/total
    return dict(C_GBseg=C_GB,C_TJ=C_TJ,f_clean_GB=1-C_GB,f_iso=fiso,
        pore_bearing_GB_fraction=C_GB,pore_free_GB_fraction=1-C_GB,
        GBseg_pore_mean_radius=_mean_radius(s.phi_GBseg,s.pore_radii),TJ_pore_mean_radius=_mean_radius(s.phi_TJ,s.pore_radii),
        iso_pore_mean_radius=_mean_radius(s.phi_iso,s.pore_radii),A_GB=A_GB,L_TJ=L_TJ,
        GBseg_projected_area_density=area_GB,TJ_pore_line_density=line_TJ,
        coverage_argument_GB=arg_GB,occupancy_argument_TJ=arg_TJ)


def stress_diagnostics(s:PoreLocationState,p:LocationParams,top:dict)->dict:
    base=4*p.base.gamma_s/max(s.G,1e-30)
    sigma_GB=base*top['C_GBseg'];sigma_TJ=base*top['C_TJ'];sigma_clean=base*top['f_clean_GB'];sigma_iso=base*top['f_iso']
    total=min(base+p.K_GBseg*sigma_GB+p.K_TJ*sigma_TJ+p.K_clean*sigma_clean,p.stress_cap)
    return dict(sigma_base=base,sigma_GBseg_pore=sigma_GB,sigma_TJ_pore=sigma_TJ,sigma_clean_GB=sigma_clean,sigma_iso=sigma_iso,sigma_act_total=total)


def location_fluxes(s:PoreLocationState,T_C:float,activity:float,p:LocationParams,top:dict)->dict:
    """Conservative placement and within-location smoothing fluxes."""
    z=np.zeros_like(s.phi_GBseg)
    if p.pore_location_mode!='evolving':return {k:z.copy() for k in ('GB_smooth','GB_to_TJ','TJ_to_iso','TJ_smooth')}
    low=(1-activity);r=s.pore_radii;rref=p.base.pore_radius0
    gb_s=_arr(p.k_GB_smooth_0_s,p.Q_GB_smooth_J_mol,T_C)*low*s.phi_GBseg[:-1]*(rref/r[:-1])**2
    tj_s=_arr(p.k_TJ_smooth_0_s,p.Q_TJ_smooth_J_mol,T_C)*low*s.phi_TJ[:-1]*(rref/r[:-1])**2
    relocation=_arr(p.k_GB_to_TJ_0_s,p.Q_GB_to_TJ_J_mol,T_C)*low*(top['f_clean_GB']+.2)*s.phi_GBseg
    isolate_gate=aggregate.sig((s.rho-.86)/.025)*(.25+.75*top['C_TJ'])*(1-top['C_GBseg']+.1)
    isolation=_arr(p.k_TJ_iso_0_s,p.Q_TJ_iso_J_mol,T_C)*isolate_gate*s.phi_TJ*(r/rref)**.5
    def upward(j):
        out=np.zeros_like(r);out[:-1]-=j;out[1:]+=j;return out
    return dict(GB_smooth=upward(gb_s),GB_to_TJ=relocation,TJ_to_iso=isolation,TJ_smooth=upward(tj_s))


def instantaneous(s:PoreLocationState,T_C:float,p:LocationParams)->dict:
    top=topology_diagnostics(s,p);stress=stress_diagnostics(s,p,top)
    fake=aggregate.State(s.rho,s.G,s.pore_radii,s.phi_total,_number(s.phi_total,s.pore_radii),
        aggregate.TopologyState(top['C_GBseg'],top['f_clean_GB'],0,top['C_TJ'],1-top['f_iso'],top['f_iso']),
        aggregate.StressState(stress['sigma_base'],stress['sigma_act_total']-stress['sigma_base'],stress['sigma_act_total']))
    kin=aggregate.kinetic_diagnostics(fake,T_C,p.base);event_rate=1/max(kin['tau_event'],1e-300);r=s.pore_radii;rref=p.base.pore_radius0
    gb_weight=s.phi_GBseg*(rref/r)**p.q_GBseg;tj_weight=s.phi_TJ*(rref/r)**p.q_TJ
    gb_elig=top['C_GBseg']*float(np.sum(gb_weight));tj_elig=p.eta_TJ_dens*top['C_TJ']*float(np.sum(tj_weight))
    eligibility=(gb_elig+tj_elig)/max(float(np.sum(s.phi_total)),1e-300)
    rho_dot=min((1-s.rho)*p.base.event_strain*event_rate*eligibility,.05)
    gb_share=gb_elig/max(gb_elig+tj_elig,1e-300);tj_share=1-gb_share
    gb_remove=rho_dot*gb_share*gb_weight/max(float(np.sum(gb_weight)),1e-300)
    tj_remove=rho_dot*tj_share*tj_weight/max(float(np.sum(tj_weight)),1e-300)
    T=T_C+273.15;Tref=p.pin_T_ref_C+273.15;size=p.pin_G_ref/max(s.G,1e-30)
    Rgb=p.A_GBseg_pin*top['C_GBseg']*size**p.pin_size_exp_GB*math.exp(float(np.clip((p.pin_Q_GB_J_mol-p.base.Q_growth)/R*(1/T-1/Tref),-60,60)))
    Rtj=p.A_TJ_pin*top['C_TJ']*(1-top['f_iso'])*size**p.pin_size_exp_TJ*math.exp(float(np.clip((p.pin_Q_TJ_J_mol-p.base.Q_growth)/R*(1/T-1/Tref),-60,60)))
    mobility=1/(1+Rgb+Rtj);free_G=top['f_clean_GB']*aggregate.arr(p.base.growth_prefactor_m2_s,p.base.Q_growth,T)/max(s.G,1e-30)
    G_dot=free_G*mobility+p.base.event_growth_fraction*s.G*event_rate*(1-min(eligibility,1))
    flux=location_fluxes(s,T_C,kin['activity'],p,top)
    powers=dict(P_GBseg_dens=stress['sigma_act_total']*rho_dot*gb_share,P_TJ_dens=stress['sigma_act_total']*rho_dot*tj_share,
        P_clean_GB=p.base.gamma_gb*free_G/max(s.G**2,1e-300),P_GBseg_drag=p.base.gamma_gb*free_G*(1-mobility)*Rgb/max(Rgb+Rtj,1e-300)/max(s.G**2,1e-300),
        P_TJ_drag=p.base.gamma_gb*free_G*(1-mobility)*Rtj/max(Rgb+Rtj,1e-300)/max(s.G**2,1e-300),
        P_GB_to_TJ_relocation=p.base.gamma_s*float(np.sum(flux['GB_to_TJ']/r)),P_TJ_iso_conversion=p.base.gamma_s*float(np.sum(flux['TJ_to_iso']/r)),P_iso=0.,P_persistent_drag=0.)
    return {**top,**stress,**kin,**powers,'rho_dot':rho_dot,'G_dot':G_dot,'GBseg_remove':gb_remove,'TJ_remove':tj_remove,
        'R_GBseg_pin':Rgb,'R_TJ_pin':Rtj,'growth_mobility_factor':mobility,'densification_eligibility':eligibility,**flux}


def _aggregate_disabled(p:LocationParams,protocol,stop_at_rho,initial):
    if initial is not None:raise ValueError("disabled mode accepts aggregate-model initial states through aggregate.run directly")
    return aggregate.run(p.base,protocol,stop_at_rho=stop_at_rho)


def run(p:LocationParams,protocol,stop_at_rho:Optional[float]=None,initial:Optional[PoreLocationState]=None)->dict:
    validate_params(p)
    if p.pore_location_mode=='disabled':return _aggregate_disabled(p,protocol,stop_at_rho,initial)
    s=initial_state(p) if initial is None else clone_state(initial)
    scalar_keys='t T_C rho G C_GBseg C_TJ f_clean_GB f_iso pore_bearing_GB_fraction pore_free_GB_fraction GBseg_pore_mean_radius TJ_pore_mean_radius iso_pore_mean_radius A_GB L_TJ GBseg_projected_area_density TJ_pore_line_density coverage_argument_GB occupancy_argument_TJ sigma_base sigma_GBseg_pore sigma_TJ_pore sigma_clean_GB sigma_iso sigma_act_total activity rho_dot G_dot E_G R_GBseg_pin R_TJ_pin growth_mobility_factor densification_eligibility'.split()
    power_keys='P_GBseg_dens P_TJ_dens P_clean_GB P_GBseg_drag P_TJ_drag P_GB_to_TJ_relocation P_TJ_iso_conversion P_iso P_persistent_drag'.split()
    h={k:[] for k in scalar_keys+power_keys};h.update(phi_GBseg=[],phi_TJ=[],phi_iso=[],N_GBseg=[],N_TJ=[],N_iso=[])
    while s.t<min(protocol.t_end,p.base.t_max_s) and s.rho<p.base.rho_cap:
        T=protocol.T(s.t,s.rho);d=instantaneous(s,T,p)
        vals={'t':s.t,'T_C':T,'rho':s.rho,'G':s.G,'E_G':d['rho_dot']/(d['G_dot']/max(s.G,1e-30)+1e-30),**d}
        for k in scalar_keys+power_keys:h[k].append(vals[k])
        for k in ('phi_GBseg','phi_TJ','phi_iso','N_GBseg','N_TJ','N_iso'):h[k].append(getattr(s,k).copy())
        if stop_at_rho is not None and s.rho>=stop_at_rho:break
        dt=min(p.base.dt_max_s,protocol.t_end-s.t);dT=abs(protocol.T(s.t+1,s.rho)-T)
        if dT:dt=min(dt,p.base.dT_max_C/dT)
        if d['rho_dot']>0:dt=min(dt,p.base.drho_max/d['rho_dot'])
        if d['G_dot']>0:dt=min(dt,p.base.dG_fraction_max*s.G/d['G_dot'])
        out_gb=-d['GBseg_remove']+d['GB_smooth']-d['GB_to_TJ']
        out_tj=-d['TJ_remove']+d['TJ_smooth']+d['GB_to_TJ']-d['TJ_to_iso']
        out_iso=d['TJ_to_iso']
        max_loss=max(float(np.max(np.maximum(-out_gb,0)/np.maximum(s.phi_GBseg,1e-300))),float(np.max(np.maximum(-out_tj,0)/np.maximum(s.phi_TJ,1e-300))),0)
        if max_loss>0:dt=min(dt,.2/max_loss)
        dt=max(p.base.dt_min_s,dt)
        s.phi_GBseg=np.maximum(s.phi_GBseg+out_gb*dt,0);s.phi_TJ=np.maximum(s.phi_TJ+out_tj*dt,0);s.phi_iso=np.maximum(s.phi_iso+out_iso*dt,0)
        s.rho=1-float(np.sum(s.phi_total));s.N_GBseg=_number(s.phi_GBseg,s.pore_radii);s.N_TJ=_number(s.phi_TJ,s.pore_radii);s.N_iso=_number(s.phi_iso,s.pore_radii)
        s.G=max(s.G+d['G_dot']*dt,1e-9);s.t+=dt
    out={k:np.asarray(v,float) for k,v in h.items()};out['pore_radii']=s.pore_radii.copy();return out


def clone_state(s:PoreLocationState,reset_time:bool=False)->PoreLocationState:
    return PoreLocationState(s.rho,s.G,s.pore_radii.copy(),s.phi_GBseg.copy(),s.phi_TJ.copy(),s.phi_iso.copy(),s.N_GBseg.copy(),s.N_TJ.copy(),s.N_iso.copy(),0. if reset_time else s.t)


def state_from_result(h:dict,index:int=-1)->PoreLocationState:
    r=np.asarray(h.get('pore_radii',[]))
    if r.size==0:raise ValueError("result does not carry radii; pass state from run_with_state")
    return PoreLocationState(float(h['rho'][index]),float(h['G'][index]),r,np.asarray(h['phi_GBseg'][index]).copy(),np.asarray(h['phi_TJ'][index]).copy(),np.asarray(h['phi_iso'][index]).copy(),np.asarray(h['N_GBseg'][index]).copy(),np.asarray(h['N_TJ'][index]).copy(),np.asarray(h['N_iso'][index]).copy())


def final_state(h:dict,p:LocationParams,index:int=-1)->PoreLocationState:
    r=aggregate.initial_state(p.base).pore_radii
    return PoreLocationState(float(h['rho'][index]),float(h['G'][index]),r,np.asarray(h['phi_GBseg'][index]).copy(),np.asarray(h['phi_TJ'][index]).copy(),np.asarray(h['phi_iso'][index]).copy(),np.asarray(h['N_GBseg'][index]).copy(),np.asarray(h['N_TJ'][index]).copy(),np.asarray(h['N_iso'][index]).copy())


LOCAL_FUNCTIONS=(topology_diagnostics,stress_diagnostics,location_fluxes,instantaneous)
