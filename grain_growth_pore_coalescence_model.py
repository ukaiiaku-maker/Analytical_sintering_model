#!/usr/bin/env python3
"""Number- and volume-resolved grain-growth-driven pore coalescence memory."""
from dataclasses import dataclass
import numpy as np
from massive_latent_topology_models import R,sigmoid

@dataclass
class PoreSweepState:
 rho:float;G_nm:float
 N_connected_fine_GBseg:float;N_connected_fine_TJ:float;N_large_attached:float;N_large_TJ:float;N_isolated:float;N_closed:float
 phi_connected_fine:float;phi_large_attached:float;phi_large_TJ:float;phi_isolated:float;phi_closed:float
 swept_pore_memory:float=0.;coalesced_pore_fraction:float=0.;dragged_pore_fraction:float=0.;detached_pore_fraction:float=0.;recaptured_pore_fraction:float=0.;X_J:float=0.

def diagnostics(s):
 phi=np.array([s.phi_connected_fine,s.phi_large_attached,s.phi_large_TJ,s.phi_isolated,s.phi_closed]);N=np.maximum(np.array([s.N_connected_fine_GBseg,s.N_large_attached,s.N_large_TJ,s.N_isolated,s.N_closed]),1e-30);r=(phi/N)**(1/3);w=phi/max(phi.sum(),1e-30);cdf=np.cumsum(w);return dict(removable_pore_fraction=s.phi_connected_fine/max(phi.sum(),1e-30),pore_number_reduction_factor=N.sum(),D50_nm=float(r[min(np.searchsorted(cdf,.5),4)]),D90_nm=float(r[min(np.searchsorted(cdf,.9),4)]),large_pore_tail_fraction=float(w[1:3].sum()))

def local_rates(s,T_C,p):
 """Instantaneous pore-sweep and migration laws."""
 T=T_C+273.15;therm=lambda Q:np.exp(-Q/R*(1/T-1/1623.15));base_growth=p['k_growth']*therm(p['Q_growth'])/max(s.G_nm,1);connected=s.phi_connected_fine+s.phi_large_attached+s.phi_large_TJ;total=max(1-s.rho,1e-30)
 completion=(1-np.exp(-p['lambda_TJ']/max(p['K_TJ']*(s.G_nm/150)**p['q_TJ'],1e-12)))
 drag=p['pore_drag_strength']*(s.phi_large_attached+s.phi_large_TJ)/total+p['A_J']*s.X_J
 Gdot=base_growth*completion/(1+drag);sweep=(abs(Gdot)/max(s.G_nm,1))**p['sweep_Gdot_exponent']*therm(p['Q_sweep'])
 coal=p['k_sweep_coalesce']*sweep*s.phi_connected_fine;to_tj=p['k_TJ_coalescence']*sweep*s.phi_connected_fine;detach=p['k_drag_detach']*therm(p['Q_detach'])*sweep*(s.phi_large_attached+s.phi_large_TJ);recap=p['k_recapture']*therm(p['Q_recapture'])*sweep*s.phi_isolated
 close=p['k_closed_transition']*sigmoid((s.rho-p['rho_close_mid'])/p['rho_close_width'])*s.phi_isolated
 removable=s.phi_connected_fine/max(total,1e-30);open_shrink=p['k_open']*therm(p['Q_density'])*s.phi_connected_fine*sigmoid((removable-p['removable_fraction_threshold'])*10)
 closed_shrink=p['closed_pore_shrinkage_prefactor']*therm(p['Q_closed'])*s.phi_closed*(1-p['gas_pressure_ratio'])
 return dict(G_dot=Gdot,sweep_rate=sweep,coalesce=coal,to_TJ=to_tj,detach=detach,recapture=recap,close=close,rho_dot_open=open_shrink,rho_dot_closed=closed_shrink,pore_drag=drag,Lambda_over_K=p['lambda_TJ']/max(p['K_TJ']*(s.G_nm/150)**p['q_TJ'],1e-12),migration_factor=completion/(1+drag))

def advance(s,T,p,dt):
 r=local_rates(s,T,p);cap=lambda x,a:min(max(x*dt,0),max(a,0));co=cap(r['coalesce'],s.phi_connected_fine);tj=cap(r['to_TJ'],s.phi_connected_fine-co);de=cap(r['detach'],s.phi_large_attached+s.phi_large_TJ);re=cap(r['recapture'],s.phi_isolated);cl=cap(r['close'],s.phi_isolated-re);op=cap(r['rho_dot_open'],s.phi_connected_fine-co-tj+re);cs=cap(r['rho_dot_closed'],s.phi_closed+cl)
 la=s.phi_large_attached/max(s.phi_large_attached+s.phi_large_TJ,1e-30);s.phi_connected_fine=max(s.phi_connected_fine-co-tj+re-op,0);s.phi_large_attached=max(s.phi_large_attached+co-de*la,0);s.phi_large_TJ=max(s.phi_large_TJ+tj-de*(1-la),0);s.phi_isolated=max(s.phi_isolated+de-re-cl,0);s.phi_closed=max(s.phi_closed+cl-cs,0)
 # Coalescence transfers volume but merges fine pores into fewer large pores.
 nco=cap(p['number_merge_efficiency']*r['coalesce'],s.N_connected_fine_GBseg);ntj=cap(p['number_merge_efficiency']*r['to_TJ'],s.N_connected_fine_TJ);s.N_connected_fine_GBseg=max(s.N_connected_fine_GBseg-nco,0);s.N_connected_fine_TJ=max(s.N_connected_fine_TJ-ntj,0);s.N_large_attached+=nco/max(p['coalescence_radius_exponent'],1);s.N_large_TJ+=ntj/max(p['coalescence_radius_exponent'],1);s.N_isolated=max(s.N_isolated+de-re-cl,0);s.N_closed=max(s.N_closed+cl-cs,0)
 s.G_nm+=r['G_dot']*dt;s.rho=1-(s.phi_connected_fine+s.phi_large_attached+s.phi_large_TJ+s.phi_isolated+s.phi_closed);s.swept_pore_memory+=r['sweep_rate']*dt;s.coalesced_pore_fraction+=co+tj;s.dragged_pore_fraction+=de;s.detached_pore_fraction+=de;s.recaptured_pore_fraction+=re;s.X_J=np.clip(s.X_J+p['XJ_prod']*(tj+de)-s.X_J/max(p['tau_J'],1)*dt,0,p['XJ_capacity']);return r

def initial_state(rho0=.70,G0=100.):
 v=1-rho0;return PoreSweepState(rho0,G0,1.,.3,.02,.01,.05,0.,.62*v,.08*v,.05*v,.25*v,0.)

def defaults():return dict(k_growth=9e3,Q_growth=500e3,Q_density=475e3,k_open=1.5e-5,k_sweep_coalesce=20.,k_drag_detach=3.,k_recapture=.3,k_TJ_coalescence=5.,k_closed_transition=2e-5,coalescence_radius_exponent=4.,number_merge_efficiency=2.,sweep_Gdot_exponent=1.,removable_fraction_threshold=.08,pore_drag_strength=30.,closed_pore_shrinkage_prefactor=2e-6,Q_sweep=350e3,Q_detach=400e3,Q_recapture=400e3,Q_closed=475e3,q_TJ=1,lambda_TJ=2.,K_TJ=3.,A_J=10.,XJ_prod=4.,tau_J=2e6,XJ_capacity=2.,rho_close_mid=.92,rho_close_width=.02,gas_pressure_ratio=.25)

def simulate(p,T1=1400,T2=None,rho_switch=.88,rho0=.70,G0=100.,hours=500,dt=1800):
 s=initial_state(rho0,G0);keys=list(PoreSweepState.__dataclass_fields__)+['t','T_C','G_dot','rho_dot_open','rho_dot_closed','pore_drag','Lambda_over_K','migration_factor','D50_nm','D90_nm','large_pore_tail_fraction','removable_pore_fraction','pore_number_reduction_factor'];h={k:[] for k in keys};t=0.;sw=False
 while t<=hours*3600 and s.rho<.995:
  T=T2 if sw and T2 is not None else T1;r=local_rates(s,T,p);d=diagnostics(s)
  for k in PoreSweepState.__dataclass_fields__:h[k].append(getattr(s,k))
  h['t'].append(t);h['T_C'].append(T)
  for k in ('G_dot','rho_dot_open','rho_dot_closed','pore_drag','Lambda_over_K','migration_factor'):h[k].append(r[k])
  for k,v in d.items():h[k].append(v)
  if T2 is not None and not sw and s.rho>=rho_switch:sw=True
  advance(s,T,p,min(dt,hours*3600-t));t+=dt
  if t>=hours*3600:break
 out={k:np.asarray(v,float) for k,v in h.items()};out['switched']=sw;return out

LOCAL_FUNCTIONS=(local_rates,)
