#!/usr/bin/env python3
"""Explicit topology, renewal-event, and dissipation-partition prototype."""
from __future__ import annotations
from dataclasses import dataclass, field
import math
from typing import Optional, Protocol
import numpy as np

KB, R = 1.380649e-23, 8.31446261815324
def arr(A,Q,T): return A*math.exp(float(np.clip(-Q/(R*T),-700,700)))
def sig(x): return float(1/(1+math.exp(-float(np.clip(x,-60,60)))))

@dataclass
class TopologyState:
    f_pore: float; f_clean: float; f_PR: float; f_TL: float
    connectivity: float; isolated_pore_fraction: float
@dataclass
class StressState:
    sigma_base: float; sigma_concentration: float; sigma_local: float
@dataclass
class State:
    rho: float; G: float; pore_radii: np.ndarray; pore_phi: np.ndarray; pore_N: np.ndarray
    topology: TopologyState; stress: StressState; t: float=0.0
@dataclass
class MechanismFlux:
    rho_dot: float=0.; G_dot: float=0.
    pore_phi_dot: np.ndarray=field(default_factory=lambda:np.zeros(0))
    pore_N_dot: np.ndarray=field(default_factory=lambda:np.zeros(0))
    power: float=0.; diagnostics: dict=field(default_factory=dict)
@dataclass
class Params:
    rho0: float=.74; G0: float=150e-9; n_bins: int=9
    pore_radius0: float=22e-9; pore_radius_max_factor: float=18.; pore_ln_sigma: float=.65
    gamma_s: float=1.; gamma_gb: float=.35; event_strain: float=4e-3; event_growth_fraction: float=8e-5
    coverage_chi: float=1.8; connectivity_rho_mid: float=.875; connectivity_rho_width: float=.025
    isolation_rho_mid: float=.90; isolation_rho_width: float=.018; concentration_factor: float=.65; stress_cap: float=5e8
    nucleation_prefactor: float=3e13; Q_nucleation: float=435e3; activation_volume: float=8e-29
    exchange_prefactor_s: float=2e-10; Q_exchange: float=245e3
    transport_prefactor_s_m2: float=1e-3; Q_transport: float=360e3
    tl_prefactor_s: float=2e-8; Q_TL: float=285e3
    growth_prefactor_m2_s: float=5e-5; Q_growth: float=410e3
    pore_drag_resistance: float=3.; tl_drag_resistance: float=2.
    pr_prefactor_s: float=5e-6; Q_PR: float=260e3
    coarsening_prefactor_s: float=2e-5; Q_coarsening: float=300e3
    removal_radius_exp: float=2.; dt_min_s: float=.2; dt_max_s: float=180.
    drho_max: float=1e-3; dG_fraction_max: float=.01; dT_max_C: float=2.; rho_cap: float=.985; t_max_s: float=8e5
    enable_PR: bool=True; enable_TL_drag: bool=True; enable_pore_coarsening: bool=True

class ThermalProtocol(Protocol):
    t_end: float
    def T(self,t:float,rho:float)->float: ...
class RampHold:
    def __init__(self,heating_rate_C_min,target_C=1450.,hold_s=8*3600.,start_C=25.):
        self.rate=heating_rate_C_min/60; self.target_C=target_C; self.start_C=start_C
        self.ramp_s=(target_C-start_C)/self.rate; self.t_end=self.ramp_s+hold_s
    def T(self,t,rho): return min(self.target_C,self.start_C+self.rate*t)
class Iso:
    def __init__(self,T_C,t_end=96*3600): self.T_C,self.t_end=T_C,t_end
    def T(self,t,rho): return self.T_C
class TwoStep:
    def __init__(self,T1_C=1350.,T2_C=1250.,rho_switch=.83,t_end=96*3600):
        self.T1_C,self.T2_C,self.rho_switch,self.t_end=T1_C,T2_C,rho_switch,t_end
    def T(self,t,rho): return self.T2_C if rho>=self.rho_switch else self.T1_C

def pore_number(phi,r): return np.maximum(phi,0)/np.maximum(4*math.pi*r**3/3,1e-300)
def infer_topology(rho,G,radii,phi,p):
    N=pore_number(phi,radii); area=float(np.sum(N*math.pi*radii**2)); gb=2/max(G,1e-30)
    fp=1-math.exp(-max(p.coverage_chi*area/gb,0)); conn=sig((p.connectivity_rho_mid-rho)/p.connectivity_rho_width)
    iso=sig((rho-p.isolation_rho_mid)/p.isolation_rho_width)
    small=float(np.sum(phi[radii<=np.median(radii)])/max(np.sum(phi),1e-300))
    ftl=fp*(1-math.exp(-max(float(np.sum(N*2*math.pi*radii)*G**2),0)))
    return TopologyState(*(float(np.clip(x,0,1)) for x in (fp,1-fp,fp*conn*small,ftl,conn,iso)))
def infer_stress(s,p):
    base=4*p.gamma_s/max(s.G,1e-30); conc=p.concentration_factor*s.topology.f_TL*base
    return StressState(base,conc,min(base+conc,p.stress_cap))
def initial_state(p):
    r=np.geomspace(p.pore_radius0,p.pore_radius0*p.pore_radius_max_factor,p.n_bins)
    w=np.exp(-.5*(np.log(r/p.pore_radius0)/p.pore_ln_sigma)**2); phi=(1-p.rho0)*w/w.sum()
    top=infer_topology(p.rho0,p.G0,r,phi,p); s=State(p.rho0,p.G0,r,phi,pore_number(phi,r),top,StressState(0,0,0))
    s.stress=infer_stress(s,p); return s
def kinetic_diagnostics(s,T_C,p):
    T=T_C+273.15; rn=p.nucleation_prefactor*math.exp(float(np.clip(-p.Q_nucleation/(R*T)+p.activation_volume*s.stress.sigma_local/(KB*T),-700,700)))
    tn=1/max(rn,1e-300); te=p.exchange_prefactor_s*math.exp(float(np.clip(p.Q_exchange/(R*T),-700,700)))
    tt=p.transport_prefactor_s_m2*s.G**2*math.exp(float(np.clip(p.Q_transport/(R*T),-700,700)))
    tl=p.tl_prefactor_s*s.topology.f_TL*math.exp(float(np.clip(p.Q_TL/(R*T),-700,700))) if p.enable_TL_drag else 0
    total=tn+te+tt+tl; completion=1/max(te+tt+tl,1e-300); L=rn/max(completion,1e-300)
    return dict(r_nuc=rn,tau_nuc=tn,tau_exchange=te,tau_transport=tt,tau_TL=tl,tau_event=total,completion_rate=completion,Lambda=L,activity=L/(1+L))
def zeros(s): return np.zeros_like(s.pore_phi),np.zeros_like(s.pore_N)
def renewal_densification(s,T,p,k):
    w=s.pore_phi*(p.pore_radius0/s.pore_radii)**p.removal_radius_exp; w/=max(w.sum(),1e-300)
    rate=1/max(k['tau_event'],1e-300); y=s.topology.f_pore*s.topology.connectivity*(1-s.topology.isolated_pore_fraction)
    rd=min((1-s.rho)*p.event_strain*rate*y,.05); pd=-rd*w; gd=p.event_growth_fraction*s.G*rate*(1-y)
    return MechanismFlux(rd,gd,pd,pd/np.maximum(4*math.pi*s.pore_radii**3/3,1e-300),s.stress.sigma_local*rd,
      {'event_rate':rate,'density_gain_per_event':p.event_strain*y,'grain_growth_per_event':p.event_growth_fraction*s.G*(1-y)})
def clean_GB_migration(s,T,p):
    gd=s.topology.f_clean*arr(p.growth_prefactor_m2_s,p.Q_growth,T+273.15)/max(s.G,1e-30); z=zeros(s)
    return MechanismFlux(G_dot=gd,pore_phi_dot=z[0],pore_N_dot=z[1],power=p.gamma_gb*gd/max(s.G**2,1e-300))
def pore_connected_GB_migration(s,T,p):
    f=clean_GB_migration(s,T,p); f.G_dot*=s.topology.f_pore*s.topology.connectivity/max(s.topology.f_clean,1e-12); f.power=p.gamma_gb*f.G_dot/max(s.G**2,1e-300); return f
def pore_drag(s,T,p):
    gd=clean_GB_migration(s,T,p).G_dot; z=zeros(s); return MechanismFlux(pore_phi_dot=z[0],pore_N_dot=z[1],power=p.pore_drag_resistance*s.topology.f_pore*gd**2/max(s.G,1e-30))
def triple_line_drag(s,T,p,k):
    z=zeros(s); power=p.tl_drag_resistance*s.topology.f_TL/max(k['tau_TL']+k['tau_nuc'],1e-300) if p.enable_TL_drag else 0
    return MechanismFlux(pore_phi_dot=z[0],pore_N_dot=z[1],power=power)
def up_flux(s,rate):
    tr=np.maximum(s.pore_phi[:-1],0)*max(rate,0); pd=np.zeros_like(s.pore_phi); pd[:-1]-=tr; pd[1:]+=tr
    return pd,pd/np.maximum(4*math.pi*s.pore_radii**3/3,1e-300)
def PR_desintering(s,T,p,k):
    rate=arr(p.pr_prefactor_s,p.Q_PR,T+273.15)*s.topology.f_PR*(1-k['activity']) if p.enable_PR else 0; pd,nd=up_flux(s,rate)
    return MechanismFlux(pore_phi_dot=pd,pore_N_dot=nd,power=p.gamma_s*float(np.sum(abs(pd)/s.pore_radii)))
def pore_coarsening(s,T,p):
    rate=arr(p.coarsening_prefactor_s,p.Q_coarsening,T+273.15) if p.enable_pore_coarsening else 0; pd,nd=up_flux(s,rate)
    return MechanismFlux(pore_phi_dot=pd,pore_N_dot=nd,power=p.gamma_s*float(np.sum(abs(pd)/s.pore_radii)))
def exchange_dissipation(s,k):
    z=zeros(s); return MechanismFlux(pore_phi_dot=z[0],pore_N_dot=z[1],power=k['tau_exchange']/k['tau_event'])
def transport_dissipation(s,k):
    z=zeros(s); return MechanismFlux(pore_phi_dot=z[0],pore_N_dot=z[1],power=k['tau_transport']/k['tau_event'])
def evaluate_mechanisms(s,T,p):
    k=kinetic_diagnostics(s,T,p); m={'renewal_densification':renewal_densification(s,T,p,k),'clean_GB_migration':clean_GB_migration(s,T,p),'pore_connected_GB_migration':pore_connected_GB_migration(s,T,p),'pore_drag':pore_drag(s,T,p),'triple_line_drag':triple_line_drag(s,T,p,k),'PR_desintering':PR_desintering(s,T,p,k),'pore_coarsening':pore_coarsening(s,T,p),'exchange_dissipation':exchange_dissipation(s,k),'transport_dissipation':transport_dissipation(s,k)}
    return k,m
def solve_dissipation_partition(state,topology,mechanisms,params):
    compat={'renewal_densification':topology.f_pore*topology.connectivity*(1-topology.isolated_pore_fraction),'clean_GB_migration':topology.f_clean,'pore_connected_GB_migration':topology.f_pore,'pore_drag':topology.f_pore,'triple_line_drag':topology.f_TL,'PR_desintering':topology.f_PR}
    q={n:max(f.power,0)*max(compat.get(n,1),0) for n,f in mechanisms.items()}; total=sum(q.values())
    w={n:(v/total if total>0 else float(n=='renewal_densification')) for n,v in q.items()}
    assert all(v>=0 for v in w.values()) and abs(sum(w.values())-1)<1e-12; return w
def combine(m,w):
    z=np.zeros_like(next(iter(m.values())).pore_phi_dot); out=MechanismFlux(pore_phi_dot=z.copy(),pore_N_dot=z.copy())
    for n,f in m.items(): out.rho_dot+=w[n]*f.rho_dot; out.G_dot+=w[n]*f.G_dot; out.pore_phi_dot+=w[n]*f.pore_phi_dot; out.pore_N_dot+=w[n]*f.pore_N_dot; out.power+=w[n]*f.power
    return out
def run(p,protocol,stop_at_rho:Optional[float]=None):
    s=initial_state(p); keys='t T_C rho G f_pore f_clean f_PR f_TL connectivity isolated_pore_fraction sigma_base sigma_concentration sigma_local r_nuc tau_exchange tau_transport tau_TL activity rho_dot dGdt E_G'.split(); h={k:[] for k in keys}; h.update(pore_phi=[],pore_N=[]); power_names=[]
    while s.t<min(protocol.t_end,p.t_max_s) and s.rho<p.rho_cap:
        T=protocol.T(s.t,s.rho); s.topology=infer_topology(s.rho,s.G,s.pore_radii,s.pore_phi,p); s.stress=infer_stress(s,p); k,m=evaluate_mechanisms(s,T,p); w=solve_dissipation_partition(s,s.topology,m,p); f=combine(m,w)
        vals={'t':s.t,'T_C':T,'rho':s.rho,'G':s.G,**vars(s.topology),**vars(s.stress),**k,'rho_dot':f.rho_dot,'dGdt':f.G_dot,'E_G':f.rho_dot/(f.G_dot/max(s.G,1e-30)+1e-30)}
        for key in keys:h[key].append(vals[key])
        h['pore_phi'].append(s.pore_phi.copy()); h['pore_N'].append(s.pore_N.copy())
        if not power_names:
            power_names=list(m)
            for n in power_names:h['power_'+n]=[];h['weight_'+n]=[]
        for n in power_names:h['power_'+n].append(m[n].power);h['weight_'+n].append(w[n])
        if stop_at_rho is not None and s.rho>=stop_at_rho:break
        dt=min(p.dt_max_s,protocol.t_end-s.t); dT=abs(protocol.T(s.t+1,s.rho)-T)
        if dT:dt=min(dt,p.dT_max_C/dT)
        if f.rho_dot>0:dt=min(dt,p.drho_max/f.rho_dot)
        if f.G_dot>0:dt=min(dt,p.dG_fraction_max*s.G/f.G_dot)
        dt=max(p.dt_min_s,dt); total=max(float(s.pore_phi.sum())-f.rho_dot*dt,0); new=np.maximum(s.pore_phi+f.pore_phi_dot*dt,0)
        if new.sum()>0:new*=total/new.sum()
        s.pore_phi=new;s.rho=1-float(new.sum());s.pore_N=pore_number(new,s.pore_radii);s.G=max(s.G+f.G_dot*dt,1e-9);s.t+=dt
    return {k:np.asarray(v,float) for k,v in h.items()}
def value_at_density(result,target):
    i=np.flatnonzero(result['rho']>=target); return (float(result['G'][i[0]]),True) if i.size else (math.nan,False)
def percent_gain(reference,improved): return 100*(reference-improved)/reference if np.isfinite(reference) and np.isfinite(improved) and reference>0 else math.nan
