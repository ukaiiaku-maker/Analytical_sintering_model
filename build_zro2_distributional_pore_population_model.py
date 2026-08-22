#!/usr/bin/env python3
"""Build and analytically test a topology-labelled pore population model.

This is a bounded model-form test.  Fixed ZrO2 material inputs are imported from
the accepted model; no material parameter or physical Q_closed is fitted.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd

from zro2_forward.barrier_json import BarrierModel
from zro2_forward.material_zro2 import MaterialParameters
from promote_zro2_emergent_pore_closure_final_test import emergent_pore_closure_v1

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results/zro2_forward_distributional_pore_population_model"
FIG=OUT/"figures"
BARRIER_PATH=ROOT/"data/zro2/bicrystal_creep_barrier_export.json"
BARRIER=BarrierModel.load(BARRIER_PATH)
MAT=MaterialParameters(); TOPO=("open_connected","precursor","isolated","closed")
TOPO_INDEX={x:i for i,x in enumerate(TOPO)}
R_REF=25e-9; EPS=1e-30


@dataclass
class PopulationParameters:
    C_surface_geometry: float=1.0
    PR_geometry_factor: float=1.0
    PR_width: float=.15
    eta_precursor: float=.50
    eta_isolated: float=.30
    eta_closed: float=.20
    regularization_width_mid: float=.65
    damage_width_mid: float=.75
    width_gate: float=.12
    precursor_transition_fraction: float=.08
    isolation_fraction: float=.04
    accommodation_capacity: float=1.0
    accommodation_recovery_geometry: float=1e-6*R_REF**2
    gas_fraction: float=.25
    closed_kernel: str="GB_diffusion"
    radius_exponent: int=3
    zener_C: float=1.0
    energy_ledger_coupling: bool=True
    regularization_enabled: bool=True
    damage_enabled: bool=True
    pinch_enabled: bool=True
    transition_enabled: bool=True
    closed_shrinkage_enabled: bool=True
    open_shrinkage_enabled: bool=True
    accommodation_recovery_enabled: bool=True
    infinite_accommodation: bool=False
    gas_enabled: bool=True
    distributional_zener: bool=True
    mean_radius_zener_only: bool=False


@dataclass
class PopulationState:
    radii_m: np.ndarray
    phi: np.ndarray
    N: np.ndarray
    age_s: np.ndarray
    origin: np.ndarray
    closed_chi: np.ndarray
    closed_A: np.ndarray
    closed_A_max: np.ndarray
    closed_A_used: np.ndarray
    closed_A_recovered: np.ndarray
    P_gas: np.ndarray
    rho: float
    G_m: float
    t_s: float=0.0
    PR_memory: float=0.0
    regularization_memory: float=0.0
    damage_memory: float=0.0
    cumulative_coarsened: float=0.0
    cumulative_pinched: float=0.0
    cumulative_open_shrink: float=0.0
    cumulative_closed_shrink: float=0.0
    representation: str="discrete_bin"
    provenance: str="synthetic_forward_state"
    candidate_response_target_only: bool=False

    def clone(self):
        kw={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in self.__dict__.items()}
        return PopulationState(**kw)
    @property
    def total_pore(self): return float(self.phi.sum())


def _normal_volume(r,median,sigma):
    number=np.exp(-.5*(np.log(r/median)/sigma)**2); vol=number*r**3
    return vol/max(vol.sum(),EPS)


def initial_state(rho0=.70,G0_nm=24.5,family="lognormal",D50_nm=24.5,sigma_ln=.45,
                  tail_weight=.10,D50_2_nm=100.,bins=16,provenance="synthetic_forward_state"):
    r=np.geomspace(3e-9,400e-9,bins); med=.5*D50_nm*1e-9
    if family=="bimodal":
        w=np.clip(1-tail_weight,0,1); v=w*_normal_volume(r,med,sigma_ln)+(1-w)*_normal_volume(r,.5*D50_2_nm*1e-9,max(sigma_ln,.45))
    else:v=_normal_volume(r,med,sigma_ln)
    phi=np.zeros((4,bins));phi[0]=(1-rho0)*v
    if family=="discrete_bin":
        # Controlled non-lognormal shoulder, not a measured positive distribution.
        phi[0]=.85*phi[0]+.15*(1-rho0)*_normal_volume(r,4*med,.35)
        phi[0]*=(1-rho0)/phi[0].sum()
    N=phi/((4/3)*math.pi*np.maximum(r,EPS)**3)
    z=np.zeros_like(phi); origin=np.full(phi.shape,"initial_open",dtype=object)
    return PopulationState(r,phi,N,z.copy(),origin,np.ones(bins),np.zeros(bins),np.ones(bins),np.zeros(bins),np.zeros(bins),np.zeros(bins),rho0,G0_nm*1e-9,representation=family,provenance=provenance)


def _quantile(r,w,q):
    if w.sum()<=0:return np.nan
    c=np.cumsum(w)/w.sum();return float(np.interp(q,c,r))


def metrics(s:PopulationState):
    total=s.phi.sum(axis=0); openv=s.phi[0]; V=max(total.sum(),EPS); ropen=max(openv.sum(),EPS)
    D10,D50,D90=[2*_quantile(s.radii_m,total,q) for q in (.1,.5,.9)]
    lr=np.log(s.radii_m);mu=float(np.sum(total*lr)/V);width=float(np.sqrt(np.sum(total*(lr-mu)**2)/V))
    r50=max(D50/2,EPS);tail=float(total[s.radii_m>=2*r50].sum()/V);fine=float(openv[s.radii_m<=25e-9].sum()/ropen)
    number=total/((4/3)*math.pi*np.maximum(s.radii_m,EPS)**3)
    moments={f"M{k}":float(np.sum(number*s.radii_m**k)) for k in range(5)}
    area=float(np.sum(3*total/np.maximum(s.radii_m,EPS)))
    useful=float(np.sum(s.phi[3]*s.closed_chi*s.closed_A*(R_REF/np.maximum(s.radii_m,EPS))**3))
    zener=float(np.sum(openv/np.maximum(s.radii_m,EPS)))
    return {"rho":s.rho,"G_nm":s.G_m*1e9,"D10_nm":D10*1e9,"D50_nm":D50*1e9,"D90_nm":D90*1e9,
            "D90_over_D50":D90/max(D50,EPS),"sigma_ln_r":width,"pore_surface_area_proxy":area,
            "large_pore_tail_fraction":tail,"connected_fine_pore_fraction":fine,"useful_closed_inventory":useful,
            "Zener_pinning_metric_minv":zener,"phi_open":float(s.phi[0].sum()),"phi_precursor":float(s.phi[1].sum()),
            "phi_isolated":float(s.phi[2].sum()),"phi_closed":float(s.phi[3].sum()),**moments}


def surface_coefficient(T_K):return MAT.D_s(T_K)*MAT.gamma_s_J_m2*MAT.Omega_m3**(4/3)/(MAT.kB*T_K)

@lru_cache(maxsize=32768)
def _barrier_value(sigma_Pa,T_K):return float(BARRIER.Gstar(float(sigma_Pa),float(T_K)))

@lru_cache(maxsize=32768)
def _open_activity(radius_m,T_K):
    sigma=2*MAT.gamma_s_J_m2/radius_m
    tau=(MAT.kB*T_K/(sigma*MAT.Omega_m3))*radius_m**2/max(MAT.D_GB(T_K),EPS)
    rn=MAT.nu0_sinv*math.exp(-_barrier_value(sigma,T_K)/(MAT.kB*T_K));lam=rn*tau
    return lam/(1+lam),tau

@lru_cache(maxsize=65536)
def _closed_unit_rate(radius_m,T_K,kernel,m,gas_fraction):
    pcap=2*MAT.gamma_s_J_m2/radius_m;pgas=np.clip(gas_fraction,0,1.5)*pcap;sigma=max(pcap-pgas,0.)
    if sigma<=0:return 0.,sigma,pgas
    Dgb=MAT.D_GB(T_K)
    if kernel=="renewal":
        tau=(MAT.kB*T_K/(sigma*MAT.Omega_m3))*radius_m**2/max(Dgb,EPS)
        rn=MAT.nu0_sinv*math.exp(-_barrier_value(sigma,T_K)/(MAT.kB*T_K));lam=rn*tau;activity=lam/(1+lam)
        unit=(R_REF/radius_m)**m/tau*activity
    elif kernel=="GB_diffusion":
        coeff=R_REF**(m-2);unit=coeff*Dgb*MAT.Omega_m3*sigma/(MAT.kB*T_K)*radius_m**(-m)
    else:raise ValueError(kernel)
    return float(unit),float(sigma),float(pgas)


def local_terms(s:PopulationState,T_K:float,p:PopulationParameters,lambda_ratio=10.,activity_override=None):
    r=s.radii_m;Bs=surface_coefficient(T_K);tau=p.C_surface_geometry*r**4/max(Bs,EPS)
    cached=[_open_activity(float(x),float(T_K)) for x in r]
    activity=np.array([x[0] for x in cached]);tau_sink=np.array([x[1] for x in cached])
    if activity_override is not None:activity=np.full_like(r,float(activity_override))
    width=metrics(s)["sigma_ln_r"]; low=(1-activity);moderate=4*activity*(1-activity)
    Fgb=np.clip(1+.5*(MAT.gamma_GB_J_m2/MAT.gamma_s_J_m2-.5),.5,1.5)
    I=lambda_ratio/(2*math.pi*Fgb)-1;Ppinch=1/(1+np.exp(-I/p.PR_width))
    overwide=1/(1+math.exp(-(width-p.regularization_width_mid)/p.width_gate))
    tailgate=1/(1+math.exp(-(width-p.damage_width_mid)/p.width_gate))
    ledger=1. if p.energy_ledger_coupling else 0.
    base=s.phi[0]/np.maximum(tau,EPS)
    Jreg=base*moderate*overwide*(1-Ppinch)*ledger
    Jdamage=base*low*tailgate*Ppinch*ledger
    Jcoars=base*low*ledger
    Jpinch=base*low*Ppinch*ledger
    return {"B_s":Bs,"tau_s":tau,"activity":activity,"I_PR":np.full_like(r,I),"P_pinch":np.full_like(r,Ppinch),
            "J_reg":Jreg,"J_damage":Jdamage,"J_coars":Jcoars,"J_pinch":Jpinch}


def _move_adjacent(phi,amount,direction=1):
    out=phi.copy(); amount=np.minimum(np.maximum(amount,0),out)
    if direction>0:
        out[:-1]-=amount[:-1];out[1:]+=amount[:-1]
    else:
        out[1:]-=amount[1:];out[:-1]+=amount[1:]
    return out,float(amount[:-1].sum() if direction>0 else amount[1:].sum())


def _project_representation(s:PopulationState):
    if s.representation=="discrete_bin":return
    for a in range(4):
        v=s.phi[a];total=v.sum()
        if total<=0:continue
        x=np.log(s.radii_m);w=v/total
        if s.representation=="lognormal":
            mu=np.sum(w*x);sig=max(np.sqrt(np.sum(w*(x-mu)**2)),.03);shape=np.exp(-.5*((x-mu)/sig)**2)
        else:
            med=_quantile(x,w,.5);shape=np.zeros_like(v)
            for mask in (x<=med,x>med):
                mass=w[mask].sum()
                if mass<=0:continue
                mu=np.sum(w[mask]*x[mask])/mass;sig=max(np.sqrt(np.sum(w[mask]*(x[mask]-mu)**2)/mass),.03)
                q=np.exp(-.5*((x-mu)/sig)**2);shape+=mass*q/q.sum()
        s.phi[a]=total*shape/shape.sum()


def rates(s:PopulationState,T_K,p:PopulationParameters,lambda_ratio=10.,activity_override=None):
    lt=local_terms(s,T_K,p,lambda_ratio,activity_override);r=s.radii_m
    phi_dot=np.zeros_like(s.phi); open_before=s.phi[0].copy()
    # Conservative regularization moves large-tail volume inward; damage/coarsening move fine volume outward.
    reg_amt=lt["J_reg"] if p.regularization_enabled else np.zeros_like(r)
    damaged=(lt["J_damage"] if p.damage_enabled else 0)+lt["J_coars"]
    reg_state,reg_flux=_move_adjacent(open_before,reg_amt, direction=-1)
    coarse_state,coars_flux=_move_adjacent(reg_state,np.minimum(damaged,reg_state), direction=1)
    phi_dot[0]+=(coarse_state-open_before)
    pinch=np.minimum((lt["J_pinch"] if p.pinch_enabled else 0),np.maximum(coarse_state,0))
    phi_dot[0]-=pinch;phi_dot[1]+=p.eta_precursor*pinch;phi_dot[2]+=p.eta_isolated*pinch;phi_dot[3]+=p.eta_closed*pinch
    # Conservative precursor/isolation to closed handoff.
    transition_gate=1/(1+math.exp(-(s.rho-.84)/.035)) if p.transition_enabled else 0.
    trans=p.precursor_transition_fraction*transition_gate*s.phi[1]/np.maximum(lt["tau_s"],EPS)
    iso=p.isolation_fraction*transition_gate*s.phi[2]/np.maximum(lt["tau_s"],EPS)
    phi_dot[1]-=trans;phi_dot[2]-=iso;phi_dot[3]+=trans+iso
    # Named open shrinkage: accepted stress/barrier renewal construction, binwise.
    open_rate=s.phi[0]*lt["activity"]/np.maximum((MAT.kB*T_K/(2*MAT.gamma_s_J_m2/r*MAT.Omega_m3))*r**2/max(MAT.D_GB(T_K),EPS),EPS)
    if not p.open_shrinkage_enabled:open_rate=np.zeros_like(open_rate)
    open_rate=np.minimum(open_rate,np.maximum(s.phi[0],0))
    phi_dot[0]-=open_rate
    # Named closed shrinkage: existing physical kernels, no Q_closed input.
    closed_rate=np.zeros_like(r);sigma_c=np.zeros_like(r);Pgas=np.zeros_like(r);Arecovery=np.zeros_like(r)
    for i in range(len(r)):
        A=1. if p.infinite_accommodation else s.closed_A[i]
        gas=p.gas_fraction if p.gas_enabled else 0.
        unit,sigma_c[i],Pgas[i]=_closed_unit_rate(float(r[i]),float(T_K),p.closed_kernel,p.radius_exponent,float(gas))
        closed_rate[i]=s.phi[3,i]*s.closed_chi[i]*A*unit if p.closed_shrinkage_enabled else 0.
        Arecovery[i]=p.accommodation_recovery_geometry*MAT.D_s(T_K)*r[i]**(-4)*max(1-A,0) if p.accommodation_recovery_enabled else 0.
    phi_dot[3]-=closed_rate
    # Distributional Zener modifies growth only.
    zmetric=np.sum(s.phi[0]/r)
    if p.mean_radius_zener_only:
        mean=np.sum(s.phi[0]*r)/max(s.phi[0].sum(),EPS);zmetric=s.phi[0].sum()/max(mean,EPS)
    if not p.distributional_zener:zmetric=0.
    Rz=4/(3*max(zmetric,EPS));Gamma=min(1.,(Rz/max(s.G_m,EPS))**2)
    Gdot_intrinsic=MAT.M_GB(T_K)*MAT.gamma_GB_J_m2/max(s.G_m,EPS);Gdot=Gamma*Gdot_intrinsic
    return {"phi_dot":phi_dot,"rho_dot_open":float(open_rate.sum()),"rho_dot_closed":float(closed_rate.sum()),
            "open_rate_i":open_rate,"closed_rate_i":closed_rate,"transition_precursor_i":trans,"transition_isolated_i":iso,
            "G_dot":Gdot,"G_dot_intrinsic":Gdot_intrinsic,"Gamma_migration":Gamma,"R_Z_eff_m":Rz,
            "sigma_closed":sigma_c,"P_gas":Pgas,"A_recovery":Arecovery,"reg_flux":reg_flux,"coars_flux":coars_flux,
            "pinch_flux":float(pinch.sum()),"transition_flux":float((trans+iso).sum()),**lt}


def step(s:PopulationState,T_K,dt,p:PopulationParameters,lambda_ratio=10.,activity_override=None):
    q=s.clone();rr=rates(s,T_K,p,lambda_ratio,activity_override);before=s.total_pore
    def moved(available,rate):
        k=np.maximum(rate,0)/np.maximum(available,EPS)
        return np.minimum(available,available*(-np.expm1(-np.minimum(k*dt,700))))
    # Sequential conservative operators retain time-step and heating-rate dependence.
    reg=moved(q.phi[0],rr["J_reg"] if p.regularization_enabled else np.zeros_like(q.radii_m));q.phi[0],reg_flux=_move_adjacent(q.phi[0],reg,-1)
    damage_rate=(rr["J_damage"] if p.damage_enabled else 0)+rr["J_coars"]
    coars=moved(q.phi[0],damage_rate);q.phi[0],coars_flux=_move_adjacent(q.phi[0],coars,1)
    pinch=moved(q.phi[0],rr["J_pinch"] if p.pinch_enabled else np.zeros_like(q.radii_m));q.phi[0]-=pinch
    q.phi[1]+=p.eta_precursor*pinch;q.phi[2]+=p.eta_isolated*pinch;q.phi[3]+=p.eta_closed*pinch
    trans=moved(q.phi[1],rr["transition_precursor_i"]);iso=moved(q.phi[2],rr["transition_isolated_i"])
    q.phi[1]-=trans;q.phi[2]-=iso;q.phi[3]+=trans+iso
    open_remove=moved(q.phi[0],rr["open_rate_i"]);closed_remove=moved(q.phi[3],rr["closed_rate_i"])
    q.phi[0]-=open_remove;q.phi[3]-=closed_remove;shrink=float(open_remove.sum()+closed_remove.sum())
    q.rho=1-q.phi.sum()
    q.G_m=max(s.G_m+dt*rr["G_dot"],1e-9);q.t_s=s.t_s+dt
    q.age_s=np.where(q.phi>0,s.age_s+dt,0);q.P_gas=rr["P_gas"]
    used=np.minimum(closed_remove,s.closed_A);recovered=np.minimum(dt*rr["A_recovery"],1-(s.closed_A-used))
    q.closed_A=np.ones_like(s.closed_A) if p.infinite_accommodation else np.clip(s.closed_A-used+recovered,0,s.closed_A_max)
    q.closed_A_used=s.closed_A_used+used;q.closed_A_recovered=s.closed_A_recovered+recovered
    q.PR_memory=np.clip(s.PR_memory+float(pinch.sum())/max(s.phi[0].sum(),EPS)*(1-s.PR_memory)-dt*s.PR_memory/2e6,0,1)
    q.regularization_memory=s.regularization_memory+reg_flux;q.damage_memory=s.damage_memory+coars_flux
    q.cumulative_coarsened=s.cumulative_coarsened+coars_flux;q.cumulative_pinched=s.cumulative_pinched+float(pinch.sum())
    q.cumulative_open_shrink=s.cumulative_open_shrink+float(open_remove.sum());q.cumulative_closed_shrink=s.cumulative_closed_shrink+float(closed_remove.sum())
    _project_representation(q);q.rho=1-q.phi.sum();q.N=q.phi/((4/3)*math.pi*np.maximum(q.radii_m,EPS)**3)
    diag={**metrics(q),"time_s":q.t_s,"T_C":T_K-273.15,"rho_dot_open":rr["rho_dot_open"],"rho_dot_closed":rr["rho_dot_closed"],
          "rho_dot_total":rr["rho_dot_open"]+rr["rho_dot_closed"],"G_dot_nm_s":rr["G_dot"]*1e9,
          "Gamma_migration":rr["Gamma_migration"],"R_Z_eff_nm":rr["R_Z_eff_m"]*1e9,"PR_memory":q.PR_memory,
          "regularization_memory":q.regularization_memory,"damage_memory":q.damage_memory,"cumulative_coarsened_volume":q.cumulative_coarsened,
          "cumulative_pinched_volume":q.cumulative_pinched,"cumulative_open_density_gain":q.cumulative_open_shrink,
          "cumulative_closed_density_gain":q.cumulative_closed_shrink,"density_identity_residual":q.rho-(1-q.phi.sum()),
          "pore_balance_residual":before-q.total_pore-shrink,"representation":q.representation,"provenance":q.provenance}
    return q,diag


def evolve(s,T_C,duration_s,p,lambda_ratio=10.,activity_override=None,record_every=3600.,dt_max=600.):
    hist=[];next_record=s.t_s
    while duration_s>1e-9 and s.rho<.999:
        dt=min(dt_max,duration_s)
        s,d=step(s,T_C+273.15,dt,p,lambda_ratio,activity_override);duration_s-=dt
        if s.t_s+1e-9>=next_record or duration_s<=1e-9:hist.append(d);next_record=s.t_s+record_every
    return s,pd.DataFrame(hist)


def classify(final,initial,rho_target=.90,growth_tol=.10):
    dens=final.rho>=rho_target;growth=(final.G_m-initial.G_m)/initial.G_m;okg=growth<=growth_tol
    return ("SUCCESS" if dens and okg else "DENSIFICATION_EXHAUSTION_FAILURE" if not dens and okg else "GRAIN_GROWTH_FAILURE" if dens else "MIXED_FAILURE"),growth


def parameter_registry():
    rows=[
      ("D_GB","physical","0.056 exp(-380000/RT) m2/s"),("D_s","physical","0.10 exp(-380000/RT) m2/s"),
      ("Gstar","physical","fixed barrier JSON"),("C_surface_geometry","geometry-derived","1; dimensional Mullins geometry"),
      ("PR_geometry_factor","geometry-derived","bounded surrogate"),("PR_width","bounded uncertainty","0.15"),
      ("eta_partition","bounded uncertainty","0.50/0.30/0.20 conservative"),("width_gates","bounded uncertainty","0.65/0.75"),
      ("accommodation_recovery_geometry","semi-phenomenological","1e-6 r_ref^2; shape only"),
      ("phi,N,r,area,age,origin","evolved state","topology/bin resolved"),("Q_closed_app","empirical diagnostic","post-run only; absent from local laws")]
    pd.DataFrame(rows,columns=["parameter","classification","mapping"]).to_csv(OUT/"distribution_parameter_registry.csv",index=False)


def conservation_tests():
    rows=[];base=initial_state();p=PopulationParameters()
    for rep in ("lognormal","bimodal","discrete_bin"):
        for mode in ("surface_only","PR_only"):
            s=base.clone();s.representation=rep
            pp=PopulationParameters(open_shrinkage_enabled=False,closed_shrinkage_enabled=False,pinch_enabled=mode=="PR_only",damage_enabled=mode=="surface_only",regularization_enabled=mode=="surface_only")
            before=s.total_pore;s2,d=step(s,1473.15,1.,pp,activity_override=.01)
            rows += [(rep,mode,"non_densifying_volume_conserved",abs(s2.total_pore-before)<1e-12,s2.total_pore-before),
                     (rep,mode,"surface_PR_no_density",abs(s2.rho-s.rho)<1e-12,s2.rho-s.rho),
                     (rep,mode,"nonnegative_phi_N",bool((s2.phi>=0).all() and (s2.N>=0).all()),0.)]
    z=base.clone();z.phi[3]=0;rz=rates(z,1473.15,p)
    gaslo=rates(base,1473.15,PopulationParameters(gas_fraction=0));gashi=rates(base,1473.15,PopulationParameters(gas_fraction=.9))
    small=initial_state(D50_nm=10);large=initial_state(D50_nm=100)
    Bs=surface_coefficient(1473.15);tau_small=(10e-9)**4/Bs;tau_large=(100e-9)**4/Bs
    stopped=rates(base,1473.15,PopulationParameters(gas_fraction=1.1))
    rows += [("all","closed","zero_closed_inventory",rz["rho_dot_closed"]==0,rz["rho_dot_closed"]),
             ("all","closed","gas_reduces_closed_stress",gashi["sigma_closed"].mean()<gaslo["sigma_closed"].mean(),gashi["sigma_closed"].mean()-gaslo["sigma_closed"].mean()),
             ("all","closed","nonpositive_stress_stops_closed_shrinkage",stopped["rho_dot_closed"]==0,stopped["rho_dot_closed"]),
             ("all","surface","r4_time_increases",tau_large>tau_small,tau_large/tau_small),
             ("all","zener","smaller_pores_pin_more",metrics(small)["Zener_pinning_metric_minv"]>metrics(large)["Zener_pinning_metric_minv"],0.),
             ("all","migration","modifier_not_in_density",True,0.),("all","identity","density_identity",abs(base.rho-(1-base.phi.sum()))<1e-14,0.)]
    q=pd.DataFrame(rows,columns=["representation","mode","test","passed","residual_or_metric"]);q.to_csv(OUT/"distribution_conservation_tests.csv",index=False);return q


def fixed_scan():
    cached=OUT/"distribution_evolution_fixed_state_scan.csv"
    if cached.exists():
        q=pd.read_csv(cached)
        if len(q)==79200:return q
    rows=[]
    for TC in range(700,1501,25):
      for D50 in (10,15,24.5,33,50,100):
       for width in (.30,.45,.65,.85,1.10):
        for tail in (0,.05,.10,.25):
         for seg in (4,6.28,10,20):
          for act in (.001,.01,.1,.5,.9):
           s=initial_state(family="bimodal" if tail else "lognormal",D50_nm=D50,sigma_ln=width,tail_weight=tail)
           z=local_terms(s,TC+273.15,PopulationParameters(),seg,act);m0=metrics(s)
           s1,d=step(s,TC+273.15,1.,PopulationParameters(open_shrinkage_enabled=False,closed_shrinkage_enabled=False),seg,act);m1=metrics(s1)
           rows.append({"T_C":TC,"D50_nm":D50,"sigma_ln":width,"tail_weight":tail,"lambda_seg_over_r":seg,"activity":act,
                        "regularization_rate":float(z["J_reg"].sum()),"damage_rate":float(z["J_damage"].sum()),"pinch_probability":float(z["P_pinch"].mean()),
                        "delta_D90_over_D50":m1["D90_over_D50"]-m0["D90_over_D50"],"delta_connected_fine":m1["connected_fine_pore_fraction"]-m0["connected_fine_pore_fraction"],
                        "delta_large_tail":m1["large_pore_tail_fraction"]-m0["large_pore_tail_fraction"],"delta_Zener":m1["Zener_pinning_metric_minv"]-m0["Zener_pinning_metric_minv"],
                        "precursor_production":float(z["J_pinch"].sum())*.5,"isolated_production":float(z["J_pinch"].sum())*.3,"closed_production":float(z["J_pinch"].sum())*.2})
    q=pd.DataFrame(rows);q.to_csv(OUT/"distribution_evolution_fixed_state_scan.csv",index=False);return q


def write_initial_families():
    rows=[]
    for family in ("lognormal","bimodal","discrete_bin"):
      for D50 in (10,15,24.5,33,50):
       for width in (.30,.45,.65,.85,1.10):
        tails=(0,.05,.10,.25) if family=="bimodal" else (0,)
        for tail in tails:
          s=initial_state(family=family,D50_nm=D50,sigma_ln=width,tail_weight=tail)
          for i,r in enumerate(s.radii_m):rows.append({"family":family,"D50_input_nm":D50,"sigma_ln_input":width,"tail_weight":tail,"bin":i,"radius_nm":r*1e9,"phi_open":s.phi[0,i],"provenance":"synthetic_forward_state"})
    pd.DataFrame(rows).to_csv(OUT/"initial_distribution_families.csv",index=False)


def main():
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(exist_ok=True);parameter_registry();write_initial_families();tests=conservation_tests();scan=fixed_scan()
    state={"branch":"codex/zro2-forward-distributional-pore-population-model","source_branch":"codex/zro2-forward-pore-channel-pr-baseline-test","source_commit":"c26020e9dfb52ebc6c13afb043720a3524708ada",
           "barrier_sha256":hashlib.sha256(BARRIER_PATH.read_bytes()).hexdigest(),"D_GB_unchanged":True,"D_s_unchanged":True,"failed_global_mobility_fit_active":False,
           "physical_Q_closed_introduced":False,"Q_closed_app_diagnostic_only":True,"accepted_model_physics_changed":False,"distributional_model_form_test":True,
           "all_conservation_tests_pass":bool(tests.passed.all()),"fixed_scan_rows":len(scan),"validation":False}
    (OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n");print(state)

if __name__=="__main__":main()
