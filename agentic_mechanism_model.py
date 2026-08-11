#!/usr/bin/env python3
"""Persistent-junction and Class-B TJ multihit discovery closures."""
from __future__ import annotations
from dataclasses import dataclass,replace
import math
from typing import Optional
import numpy as np

import pore_location_agentic_model as prior
import pore_location_topology_model as fixed
import topology_constrained_sintering as aggregate

R=aggregate.R
MODES=("action_baseline","persistent_junction","tj_multihit_q0","tj_multihit_q1","persistent_tj_multihit_q0","persistent_tj_multihit_q1")
TJ_CONSTRAINT_MODES=("current_all_TJ_multihit","pore_relaxed_TJ_multihit","pore_pinned_TJ_drag","mixed_TJ_constraint")


@dataclass
class DiscoveryParams:
    action: prior.ActionParams
    mechanism_mode: str="action_baseline"
    A_J: float=18.
    XJ_capacity: float=1.
    XJ_prod_TJ: float=2.0
    XJ_prod_relocation: float=8.0
    XJ_prod_sweep: float=.4
    tau_J_ref_s: float=2e5
    Q_relax_J_mol: float=380e3
    J_size_exp: float=1.
    lambda_TJ_ref: float=3.
    Q_TJ_event_J_mol: float=400e3
    K_TJ0: float=3.
    q_TJ: int=0
    packet_G_ref: float=150e-9
    TJ_constraint_mode: str="current_all_TJ_multihit"
    eta_TJ_pore_relax: float=.6
    TJ_structural_background: float=.25
    A_TJ_pore_drag: float=8.
    Q_TJ_unpin_J_mol: float=220e3


@dataclass
class DiscoveryState:
    pore: fixed.PoreLocationState
    X_J: float=0.


def validate(p):
    if p.mechanism_mode not in MODES:raise ValueError("invalid mechanism_mode")
    if p.A_J<0 or p.XJ_capacity<=0 or p.tau_J_ref_s<=0 or p.lambda_TJ_ref<0 or p.K_TJ0<1:raise ValueError("invalid nonnegative mechanism parameter")
    if p.q_TJ not in (0,1):raise ValueError("q_TJ must be 0 or 1")
    if p.TJ_constraint_mode not in TJ_CONSTRAINT_MODES:raise ValueError("invalid TJ_constraint_mode")
    if not 0<=p.eta_TJ_pore_relax<=1 or p.A_TJ_pore_drag<0:raise ValueError("invalid TJ constraint parameters")


def poisson_completion(k,lam):
    """P[N>=ceil(k)] for a Poisson variable, evaluated without SciPy."""
    n=max(1,int(math.ceil(k)));x=max(float(lam),0.)
    if x>100 and x>5*n:return 1.
    term=math.exp(-min(x,745.));cdf=term
    for j in range(1,n):term*=x/j;cdf+=term
    return float(np.clip(1-cdf,0,1))


def local_mechanism(s:DiscoveryState,T_C:float,p:DiscoveryParams):
    """Instantaneous local mechanism closure using state, temperature, parameters."""
    validate(p);d=prior.allocated_fluxes(s.pore,T_C,p.action);top=fixed.topology_diagnostics(s.pore,prior.effective_location_params(p.action));T=T_C+273.15;Tr=1250+273.15
    persistent="persistent" in p.mechanism_mode;multihit="multihit" in p.mechanism_mode
    capture=d['action_flux_TJ_to_GBseg_capture'];reloc=d['action_flux_GB_to_TJ'];sweep=max(d['G_dot'],0)/max(s.pore.G,1e-30)*top['C_TJ']*top['f_clean_GB']
    production=(p.XJ_prod_TJ*top['C_TJ']*d['rho_dot']+p.XJ_prod_relocation*(capture+reloc)+p.XJ_prod_sweep*sweep)*(1-s.X_J/p.XJ_capacity) if persistent else 0.
    tau=p.tau_J_ref_s*math.exp(float(np.clip(p.Q_relax_J_mol/R*(1/T-1/Tr),-50,50)))
    relaxation=s.X_J/max(tau,1e-30) if persistent else 0.;Xdot=production-relaxation
    size=(p.packet_G_ref/max(s.pore.G,1e-30))**p.J_size_exp
    Rpersistent=p.A_J*max(s.X_J,0)*size if persistent else 0.
    q=1 if p.mechanism_mode.endswith('q1') else p.q_TJ
    K=p.K_TJ0*(max(s.pore.G,1e-30)/p.packet_G_ref)**q
    thermal=math.exp(float(np.clip(-p.Q_TJ_event_J_mol/R*(1/T-1/Tr),-50,50)))
    C_pore=float(np.clip(top['C_TJ'],0,1));C_structural=float(np.clip(C_pore+p.TJ_structural_background*(1-C_pore),0,1));mode=p.TJ_constraint_mode
    if mode=='current_all_TJ_multihit':C_constraint=C_pore;C_relaxed=0.;C_pinned=0.
    elif mode=='pore_relaxed_TJ_multihit':C_relaxed=p.eta_TJ_pore_relax*C_pore;C_constraint=max(C_pore-C_relaxed,0);C_pinned=0.
    elif mode=='pore_pinned_TJ_drag':C_relaxed=C_pore;C_constraint=max(C_structural-C_pore,0);C_pinned=C_pore
    else:C_relaxed=p.eta_TJ_pore_relax*C_pore;C_pinned=(1-p.eta_TJ_pore_relax)*C_pore;C_constraint=max(C_structural-C_relaxed,0)
    available=max(C_constraint+.05,0)*max(top['f_clean_GB'],.02)
    Lambda=p.lambda_TJ_ref*available*thermal*(p.packet_G_ref/max(s.pore.G,1e-30))
    gamma=poisson_completion(K,Lambda) if multihit else 1.
    oldG=d['G_dot'];migration_part=max(d['P_clean_GB'],0)
    thermal_unpin=math.exp(float(np.clip(-p.Q_TJ_unpin_J_mol/R*(1/T-1/Tr),-30,30)));R_pore=p.A_TJ_pore_drag*C_pinned*size/max(thermal_unpin,1e-12)
    mobility_multiplier=gamma/(1+Rpersistent+R_pore);Gdot=oldG*mobility_multiplier
    preg=max(d['P_clean_GB'],0)*(1-mobility_multiplier);den=d['rho_dot']
    regime='smooth' if Lambda/max(K,1e-30)>3 else ('stagnant' if Lambda/max(K,1e-30)<.3 else 'intermittent')
    denom=max(Rpersistent+R_pore+(1-gamma),1e-30)
    return {**d,'rho_dot':den,'G_dot':Gdot,'X_J':s.X_J,'X_J_dot':Xdot,'X_J_production':production,'X_J_relaxation':relaxation,'R_J_persistent':Rpersistent,'C_TJ_pore':C_pore,'C_TJ_constraint':C_constraint,'C_TJ_relaxed':C_relaxed,'C_TJ_pinned':C_pinned,'R_TJ_pore_drag':R_pore,'Lambda_TJ':Lambda,'K_TJ':K,'Lambda_over_K_TJ':Lambda/max(K,1e-30),'P_comp_TJ':gamma,'TJ_regime':regime,'growth_mobility_discovery':d['growth_mobility_factor']*mobility_multiplier,'P_persistent_junction_drag':preg*Rpersistent/denom if persistent else 0.,'P_TJ_multihit':preg*(1-gamma)/denom if multihit else 0.,'P_TJ_pore_drag':preg*R_pore/denom,'P_TJ_assisted_densification':d['P_TJ_dens'],'P_clean_GB_discovery':migration_part*mobility_multiplier}


def initial_state(p):return DiscoveryState(fixed.initial_state(prior.effective_location_params(p.action)),0.)
def clone_state(s,reset_time=False):return DiscoveryState(fixed.clone_state(s.pore,reset_time),s.X_J)


def run(p:DiscoveryParams,protocol,stop_at_rho:Optional[float]=None,initial:Optional[DiscoveryState]=None):
    validate(p)
    if p.mechanism_mode=='action_baseline':return prior.run(p.action,protocol,stop_at_rho=stop_at_rho,initial=None if initial is None else initial.pore)
    s=initial_state(p) if initial is None else clone_state(initial);scalars='t T_C rho G C_GBseg C_TJ f_clean_GB f_iso activity rho_dot G_dot E_G growth_mobility_factor growth_mobility_discovery sigma_base sigma_GBseg_pore sigma_TJ_pore sigma_clean_GB sigma_iso sigma_act_total X_J X_J_dot X_J_production X_J_relaxation R_J_persistent C_TJ_pore C_TJ_constraint C_TJ_relaxed C_TJ_pinned R_TJ_pore_drag Lambda_TJ K_TJ Lambda_over_K_TJ P_comp_TJ P_persistent_junction_drag P_TJ_multihit P_TJ_pore_drag P_TJ_assisted_densification P_clean_GB_discovery'.split();powers='P_GBseg_dens P_TJ_dens P_clean_GB P_GBseg_drag P_TJ_drag P_GBseg_to_TJ P_TJ_to_GBseg_capture P_TJ_iso'.split();fluxes='action_flux_GBseg_remove action_flux_TJ_remove action_flux_GB_smooth action_flux_GB_to_TJ action_flux_TJ_to_GBseg_capture action_flux_TJ_to_iso'.split();h={k:[] for k in scalars+powers+fluxes};h.update(TJ_regime=[],phi_GBseg=[],phi_TJ=[],phi_iso=[],N_GBseg=[],N_TJ=[],N_iso=[])
    lp=prior.effective_location_params(p.action)
    while s.pore.t<min(protocol.t_end,lp.base.t_max_s) and s.pore.rho<lp.base.rho_cap:
        T=protocol.T(s.pore.t,s.pore.rho);d=local_mechanism(s,T,p);vals={'t':s.pore.t,'T_C':T,'rho':s.pore.rho,'G':s.pore.G,'E_G':d['rho_dot']/(d['G_dot']/max(s.pore.G,1e-30)+1e-30),**d}
        for k in scalars+powers+fluxes:h[k].append(vals[k])
        h['TJ_regime'].append(d['TJ_regime'])
        for k in ('phi_GBseg','phi_TJ','phi_iso','N_GBseg','N_TJ','N_iso'):h[k].append(getattr(s.pore,k).copy())
        if stop_at_rho is not None and s.pore.rho>=stop_at_rho:break
        dt=min(lp.base.dt_max_s,protocol.t_end-s.pore.t);dT=abs(protocol.T(s.pore.t+1,s.pore.rho)-T)
        if dT:dt=min(dt,lp.base.dT_max_C/dT)
        if d['rho_dot']>0:dt=min(dt,lp.base.drho_max/d['rho_dot'])
        if d['G_dot']>0:dt=min(dt,lp.base.dG_fraction_max*s.pore.G/d['G_dot'])
        outgb=-d['GBseg_remove']+d['GB_smooth']-d['GB_to_TJ']+d['TJ_to_GBseg_capture'];outtj=-d['TJ_remove']+d['GB_to_TJ']-d['TJ_to_GBseg_capture']-d['TJ_to_iso'];outiso=d['TJ_to_iso']
        loss=max(float(np.max(np.maximum(-outgb,0)/np.maximum(s.pore.phi_GBseg,1e-300))),float(np.max(np.maximum(-outtj,0)/np.maximum(s.pore.phi_TJ,1e-300))),max(-d['X_J_dot']/max(s.X_J,1e-300),0))
        if loss>0:dt=min(dt,.2/loss)
        dt=max(lp.base.dt_min_s,dt);s.pore.phi_GBseg=np.maximum(s.pore.phi_GBseg+outgb*dt,0);s.pore.phi_TJ=np.maximum(s.pore.phi_TJ+outtj*dt,0);s.pore.phi_iso=np.maximum(s.pore.phi_iso+outiso*dt,0);s.pore.rho=1-float(np.sum(s.pore.phi_total));s.pore.N_GBseg=fixed._number(s.pore.phi_GBseg,s.pore.pore_radii);s.pore.N_TJ=fixed._number(s.pore.phi_TJ,s.pore.pore_radii);s.pore.N_iso=fixed._number(s.pore.phi_iso,s.pore.pore_radii);s.pore.G=max(s.pore.G+d['G_dot']*dt,1e-9);s.X_J=float(np.clip(s.X_J+d['X_J_dot']*dt,0,p.XJ_capacity));s.pore.t+=dt
    out={k:np.asarray(v,float if k!='TJ_regime' else object) for k,v in h.items()};out['pore_radii']=s.pore.pore_radii.copy();return out


def final_state(h,p,index=-1):
    return DiscoveryState(fixed.final_state(h,prior.effective_location_params(p.action),index),float(h['X_J'][index]) if 'X_J' in h else 0.)


LOCAL_FUNCTIONS=(local_mechanism,poisson_completion)
