#!/usr/bin/env python3
"""Separated material kinetics and migration-only topology closure prototype."""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

R=8.31446261815324;KB=1.380649e-23

@dataclass(frozen=True)
class MaterialKinetics:
    Q_GB_diffusion:float=500e3;D_GB_prefactor:float=2e-5
    Q_surface_diffusion:float=350e3;D_surface_prefactor:float=2e-4
    Q_disconnection_nucleation:float=500e3;v_star:float=8e-29;nu0_nucleation:float=3e13
    Q_exchange:float=245e3;tau_exchange_prefactor:float=2e-10
    Q_transport:float=360e3;tau_transport_prefactor:float=2e-3
    stress_concentration:float=3.;gamma_s:float=1.;gamma_GB:float=.35
    zeta_eta_ratio:float=1.;event_strain:float=4e-3
    pore_radius0:float=22e-9;pore_ln_sigma:float=.65;G0:float=100e-9;rho0:float=.70
    PR_prefactor:float=2e-4;PR_partition:str="balanced";activity_form:str="serial_fraction"

@dataclass(frozen=True)
class TopologyGrowthClosure:
    mode:str="disabled";TJ_drag_strength:float=0.;pore_drag_strength:float=0.
    XJ_capacity:float=.5;XJ_relaxation:float=1e-5;lambda_ref:float=2.
    K_ref:float=2.;q_TJ:int=0;pore_relax_fraction:float=.5;pore_drag_fraction:float=.5

def material_rates(rho,G,phi,radii,T_C,p:MaterialKinetics):
    """Instantaneous material law evaluated from thermodynamic state."""
    T=T_C+273.15;stress=4*p.gamma_s/max(G,1e-30)*(1+p.stress_concentration*np.sum(phi[radii>2*p.pore_radius0])/max(np.sum(phi),1e-300))
    tau_nuc=math.exp(float(np.clip(p.Q_disconnection_nucleation/(R*T)-p.v_star*stress/(KB*T),-50,50)))/p.nu0_nucleation
    tau_exchange=p.tau_exchange_prefactor*math.exp(float(np.clip(p.Q_exchange/(R*T),-50,50)))
    tau_transport=p.tau_transport_prefactor*G**2*math.exp(float(np.clip(p.Q_transport/(R*T),-50,50)))
    cycle=tau_nuc+tau_exchange+tau_transport
    activity=(tau_exchange+tau_transport)/cycle if p.activity_form=="serial_fraction" else 1/(1+tau_nuc/max(tau_exchange+tau_transport,1e-300))
    connected=max(1-(rho-.82)/.18,0);fine=float(np.sum(phi[radii<=2*p.pore_radius0])/max(np.sum(phi),1e-300));geo=connected*fine
    rho_dot=geo*p.event_strain/max(cycle,1e-300)*p.zeta_eta_ratio
    surface=p.D_surface_prefactor*math.exp(float(np.clip(-p.Q_surface_diffusion/(R*T),-50,50)))/p.pore_radius0**2
    pr=p.PR_prefactor*surface*(1-activity)**2*fine
    gb=p.D_GB_prefactor*math.exp(float(np.clip(-p.Q_GB_diffusion/(R*T),-50,50)))
    growth_base=gb*p.gamma_GB/max(G,1e-30)
    return dict(tau_nuc=tau_nuc,tau_exchange=tau_exchange,tau_transport=tau_transport,activity=activity,rho_dot=max(rho_dot,0),PR_propensity=max(pr,0),growth_base=max(growth_base,0),connected_fine=geo,stress=stress)

def topology_growth_factor(state,T_C,p:TopologyGrowthClosure):
    """Migration-only factor; never enters material densification rates."""
    if p.mode=="disabled":return 1.,dict(X_J=0.,Lambda_over_K=0.,pore_drag=0.)
    X=state.get("X_J",0);coverage=state.get("connected_coverage",0);G=state["G"]
    structural=max(X-p.pore_relax_fraction*coverage,0);K=p.K_ref*(G/150e-9)**p.q_TJ;completion=(p.lambda_ref+max(T_C-1200,0)/150)/max(K,1e-30)
    drag=p.TJ_drag_strength*structural/(1+completion)+p.pore_drag_strength*p.pore_drag_fraction*coverage
    return 1/(1+max(drag,0)),dict(X_J=X,Lambda_over_K=completion,pore_drag=drag)

def initial_state(p):
    r=np.geomspace(p.pore_radius0/2,p.pore_radius0*12,7);w=np.exp(-.5*(np.log(r/p.pore_radius0)/p.pore_ln_sigma)**2);phi=(1-p.rho0)*w/w.sum();return dict(t=0.,rho=p.rho0,G=p.G0,radii=r,phi=phi,X_J=0.,PR_exposure=0.)

def run(material,topology,protocol,dt_max=600.,max_steps=6000):
    s=initial_state(material);out={k:[] for k in "t T_C rho G pore_D50 pore_D90 connected_fine large_pore_fraction PR_exposure activity tau_nuc tau_exchange tau_transport rho_dot G_dot X_J Lambda_over_K pore_drag".split()}
    numerical_censored=False
    while s["t"]<protocol.t_end and s["rho"]<.995 and len(out["t"])<max_steps:
        T=protocol.T(s["t"],s["rho"]);d=material_rates(s["rho"],s["G"],s["phi"],s["radii"],T,material);s["connected_coverage"]=d["connected_fine"];gf,td=topology_growth_factor(s,T,topology);Gdot=d["growth_base"]*gf
        z=max(np.sum(s["phi"]),1e-300);cdf=np.cumsum(s["phi"])/z;D50=s["radii"][min(np.searchsorted(cdf,.5),6)];D90=s["radii"][min(np.searchsorted(cdf,.9),6)]
        vals={**d,**td,"t":s["t"],"T_C":T,"rho":s["rho"],"G":s["G"],"pore_D50":D50,"pore_D90":D90,"large_pore_fraction":float(np.sum(s["phi"][s["radii"]>4*material.pore_radius0])/z),"PR_exposure":s["PR_exposure"],"G_dot":Gdot}
        for k in out:out[k].append(vals[k])
        dt=min(dt_max,protocol.t_end-s["t"]);rate=max(d["rho_dot"],0)
        if rate:dt=min(dt,5e-4/rate)
        if Gdot:dt=min(dt,.01*s["G"]/Gdot)
        # Conservative adjacent-bin PR redistribution.
        move=np.minimum(d["PR_propensity"]*s["phi"][:-1]*dt,.2*s["phi"][:-1]);s["phi"][:-1]-=move;s["phi"][1:]+=move
        dr=min(d["rho_dot"]*dt,np.sum(s["phi"])*.2);weights=s["phi"]*(material.pore_radius0/s["radii"])**2;weights/=max(weights.sum(),1e-300);s["phi"]=np.maximum(s["phi"]-dr*weights,0);s["rho"]=1-float(np.sum(s["phi"]));s["G"]+=Gdot*dt;s["PR_exposure"]+=d["PR_propensity"]*dt;s["X_J"]=float(np.clip(s["X_J"]+(d["PR_propensity"]*(topology.XJ_capacity-s["X_J"])-topology.XJ_relaxation*d["activity"]*s["X_J"])*dt,0,topology.XJ_capacity));s["t"]+=dt
    if s["t"]<protocol.t_end and s["rho"]<.995:
        numerical_censored=True
    result={k:np.asarray(v,float) for k,v in out.items()}
    result["numerical_censored"]=numerical_censored
    result["rho_final"]=s["rho"]
    result["G_final"]=s["G"]
    return result

LOCAL_FUNCTIONS=(material_rates,topology_growth_factor)
