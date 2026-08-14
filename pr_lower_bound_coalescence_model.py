#!/usr/bin/env python3
"""Observable PR energy diversion coupled to pore-sweep/coalescence states."""
from dataclasses import dataclass
import numpy as np
import grain_growth_pore_coalescence_model as base
from massive_latent_topology_models import R,sigmoid

MODES=('baseline_coalescence','PR_lower_bound_only','coalescence_plus_PR_lower_bound','coalescence_plus_closed_accommodation','coalescence_plus_PR_plus_closed_accommodation','PR_no_low_activity_gate','PR_conservative_only_no_drive_loss','PR_drive_loss_no_pore_topology')

@dataclass
class PRState:
 pore:base.PoreSweepState
 cumulative_PR_surface_energy_loss:float=0.;cumulative_densifying_work:float=0.;cumulative_non_densifying_work:float=0.
 PR_energy_partition_to_smoothing:float=0.;PR_energy_partition_to_large_tail:float=0.;PR_energy_partition_to_TJ:float=0.;PR_energy_partition_to_isolated:float=0.;PR_energy_partition_to_closed:float=0.
 lost_densification_drive:float=0.;connected_fine_pore_loss_from_PR:float=0.;large_pore_generation_from_PR:float=0.;isolation_generation_from_PR:float=0.;PR_lower_bound_memory:float=0.;effective_densification_drive_after_PR:float=1.

def partitions(p):
 a=np.maximum([p['PR_to_smoothing_fraction'],p['PR_to_large_tail_fraction'],p['PR_to_TJ_fraction'],p['PR_to_isolated_fraction'],p['PR_to_closed_fraction']],0);a/=max(a.sum(),1e-30);return a

def local_rates(s,T_C,p):
 """Local PR competition derived from temperature, activity, and pore state."""
 q=s.pore;r=base.local_rates(q,T_C,p);T=T_C+273.15;activity=sigmoid((T_C-p['activity_T_mid_C'])/p['activity_T_width_C']);gate=1. if p['mode']=='PR_no_low_activity_gate' else sigmoid((p['low_activity_gate_mid']-activity)/p['low_activity_gate_width'])**p['activity_power'];fine=q.phi_connected_fine/max(1-q.rho,1e-30);k=p['k_PR_ref']*np.exp(-p['Q_PR']/R*(1/T-1/(p['T_PR_ref_C']+273.15)));H_PR=k*gate*fine;H_dens=activity*fine*max(r['rho_dot_open'],0);z=H_PR+H_dens+1e-30;wpr=H_PR/z;wd=H_dens/z
 relax=np.exp(-p['Q_PR_damage_relax']/R*(1/T-1/1473.15))/max(p['PR_damage_persistence_tau'],1);memory_dot=H_PR*(1-s.PR_lower_bound_memory)-relax*s.PR_lower_bound_memory
 drive_loss=p['drive_loss_coupling']*s.PR_lower_bound_memory if p['mode']!='PR_conservative_only_no_drive_loss' else 0.;topology=p['mode']!='PR_drive_loss_no_pore_topology';elig=np.exp(-drive_loss)*(1-p['connected_sink_loss_coupling']*s.PR_lower_bound_memory*topology);r['rho_dot_open']*=max(elig,0)
 r.update(activity=activity,H_PR=H_PR,H_dens=H_dens,w_PR=wpr,w_dens=wd,memory_dot=memory_dot,PR_topology_active=topology,effective_drive=max(elig,0));return r

def advance(s,T,p,dt):
 r=local_rates(s,T,p);q=s.pore
 # Advance the base law with its open shrinkage overridden by the named PR eligibility.
 p0={**p,'k_open':p['k_open']*r['effective_drive']};before=1-q.rho;base.advance(q,T,p0,dt);after=1-q.rho
 if p['mode'] not in ('baseline_coalescence','coalescence_plus_closed_accommodation') and r['PR_topology_active']:
  a=partitions(p);move=min(r['H_PR']*dt,q.phi_connected_fine);q.phi_connected_fine-=move;q.phi_large_attached+=move*a[1];q.phi_large_TJ+=move*a[2];q.phi_isolated+=move*a[3];q.phi_closed+=move*a[4];q.phi_connected_fine+=move*a[0]
  # unassigned smoothing fraction stays connected; transfers conserve volume.
  balance=move*(1-a.sum());q.phi_connected_fine+=balance
  nloss=min(q.N_connected_fine_GBseg,p['PR_number_loss']*move);q.N_connected_fine_GBseg-=nloss;q.N_large_attached+=nloss*a[1]/4;q.N_large_TJ+=nloss*a[2]/4;q.N_isolated+=nloss*a[3]/4;q.N_closed+=nloss*a[4]/4
  s.connected_fine_pore_loss_from_PR+=move*(1-a[0]);s.large_pore_generation_from_PR+=move*a[1];s.isolation_generation_from_PR+=move*a[3]
 q.rho=1-(q.phi_connected_fine+q.phi_large_attached+q.phi_large_TJ+q.phi_isolated+q.phi_closed)
 s.PR_lower_bound_memory=float(np.clip(s.PR_lower_bound_memory+r['memory_dot']*dt,0,1));s.lost_densification_drive+=max(0,1-r['effective_drive'])*dt;s.effective_densification_drive_after_PR=r['effective_drive'];s.cumulative_PR_surface_energy_loss+=r['H_PR']*dt;s.cumulative_non_densifying_work+=r['H_PR']*dt;s.cumulative_densifying_work+=r['H_dens']*dt
 a=partitions(p);s.PR_energy_partition_to_smoothing+=r['H_PR']*a[0]*dt;s.PR_energy_partition_to_large_tail+=r['H_PR']*a[1]*dt;s.PR_energy_partition_to_TJ+=r['H_PR']*a[2]*dt;s.PR_energy_partition_to_isolated+=r['H_PR']*a[3]*dt;s.PR_energy_partition_to_closed+=r['H_PR']*a[4]*dt;return r,before,after

def defaults():
 p=base.defaults();p.update(mode='coalescence_plus_PR_lower_bound',k_PR_ref=1e-5,Q_PR=250e3,T_PR_ref_C=1100.,low_activity_gate_mid=.2,low_activity_gate_width=.1,activity_power=2.,activity_T_mid_C=1180.,activity_T_width_C=70.,PR_to_smoothing_fraction=.1,PR_to_large_tail_fraction=.3,PR_to_TJ_fraction=.2,PR_to_isolated_fraction=.3,PR_to_closed_fraction=.1,drive_loss_coupling=1.,connected_sink_loss_coupling=.3,large_tail_penalty_power=2.,PR_damage_persistence_tau=1e6,Q_PR_damage_relax=300e3,PR_number_loss=3.)
 return p

def simulate(p,T1=1400,T2=None,rho_switch=.88,rho0=.70,G0=100.,hours=500,dt=1800):
 s=PRState(base.initial_state(rho0,G0));keys=list(PRState.__dataclass_fields__)[1:]+list(base.PoreSweepState.__dataclass_fields__)+['t','T_C','G_dot','rho_dot_open','rho_dot_closed','pore_drag','Lambda_over_K','migration_factor','activity','H_PR','H_dens','w_PR','w_dens','D50_nm','D90_nm','large_pore_tail_fraction','removable_pore_fraction','pore_number_reduction_factor'];h={k:[] for k in keys};t=0.;sw=False
 while t<=hours*3600 and s.pore.rho<.995:
  T=T2 if sw and T2 is not None else T1;r=local_rates(s,T,p);d=base.diagnostics(s.pore)
  for k in PRState.__dataclass_fields__:
   if k!='pore':h[k].append(getattr(s,k))
  for k in base.PoreSweepState.__dataclass_fields__:h[k].append(getattr(s.pore,k))
  h['t'].append(t);h['T_C'].append(T)
  for k in ('G_dot','rho_dot_open','rho_dot_closed','pore_drag','Lambda_over_K','migration_factor','activity','H_PR','H_dens','w_PR','w_dens'):h[k].append(r[k])
  for k,v in d.items():h[k].append(v)
  if T2 is not None and not sw and s.pore.rho>=rho_switch:sw=True
  advance(s,T,p,min(dt,hours*3600-t));t+=dt
  if t>=hours*3600:break
 out={k:np.asarray(v,float) for k,v in h.items()};out['switched']=sw;return out

LOCAL_FUNCTIONS=(local_rates,)
