#!/usr/bin/env python3
from dataclasses import replace
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import run_zro2_forward_closed_channel_physical_law_comparison as base
from zro2_forward.conditioned_950c import ConditionedTwoStep,run_path,BARRIER
from zro2_forward.resolved_rules import ResolvedRuleModel
from zro2_forward.schedules import RampNoHold,Iso
OUT=Path('results/zro2_forward_closed_pore_evolution_addendum');TARGET_RHO=.976;TARGET_G=.29
MODES={
 'renewal_m4':dict(closed_mapping_mode='mechanistic_renewal_limited',closed_radius_exponent=4),
 'renewal_m3_sensitivity':dict(closed_mapping_mode='mechanistic_renewal_limited',closed_radius_exponent=3),
 'GB_transport_m4':dict(closed_mapping_mode='mechanistic_GB_diffusion',closed_radius_exponent=4),
 'GB_transport_m3':dict(closed_mapping_mode='mechanistic_GB_diffusion',closed_radius_exponent=3),
 'gas_accommodation_m4':dict(closed_mapping_mode='mechanistic_gas_accommodation',closed_radius_exponent=4),
 'surface_accommodation_only':dict(closed_mapping_mode='mechanistic_surface_accommodation',closed_radius_exponent=4),}
def model(mode):return ResolvedRuleModel(parameters=base.params(mechanism_mode='defined_laws_port',**MODES[mode]))
def initial():return base.initial()
def classify(f,s):return base.classify(f,s,True)
def integral(h,col):return float(np.trapezoid(np.maximum(h[col],0),h.t_s)) if len(h)>1 else 0.
def prepare(mode):
 _,h=run_path(model(mode),ConditionedTwoStep(1400,1100,.88,1),initial(),600,600,'prepare');pre=h[h.rho<.88];return base.state_from_row(pre.iloc[-1] if len(pre) else h.iloc[-1])
def injected(seed):return base.state_from_row(seed,rho=.88,G_nm=117,closed_fraction=.649,A=.152,PR=1.)
def detailed(h,run_id,mode,kind):
 x=h.copy();x.insert(0,'run_id',run_id);x.insert(1,'mode',mode);x.insert(2,'state_kind',kind)
 dt=x.t_s.diff().fillna(x.t_s.iloc[0]);x['cumulative_closed_density_gain']=np.cumsum(x.rho_dot_closed_sinv*dt);x['cumulative_open_density_gain']=np.cumsum(x.rho_dot_open_sinv*dt)
 keep=['run_id','mode','state_kind','t_s','T_C','rho','G_um','phi_closed_json','N_closed_i_json','r_closed_i_json','chi_shrink_i_json','C_geom_i_json','P_gas_i_json','sigma_closed_i_json','Gstar_closed_i_json','r_nuc_closed_i_json','tau_transport_closed_i_json','tau_cycle_closed_i_json','A_closed_i_json','A_closed_used_i_json','A_closed_recovered_i_json','rho_dot_closed_i_json','rho_dot_open_i_json','cumulative_closed_density_gain','cumulative_open_density_gain','surface_diffusion_accommodation_rate_i_json','gas_pressure_factor_i_json','radius_exponent_m','closed_law_mode','gas_model']
 keep.insert(7,'closed_inventory')
 keep.extend(['rho_dot_closed_sinv','rho_dot_open_sinv'])
 return x[keep]
def main():
 OUT.mkdir(parents=True,exist_ok=True);fixed=[];scans=[];hist=[];natural={};seed=None
 registry=[]
 for mode,kw in MODES.items():
  registry.append({'mode':mode,'governing_source':'closed-pore evolution law addendum','transport':'D_GB fixed' if 'surface' not in mode else 'D_s fixed; accommodation only','barrier':'fixed Gstar(sigma,T)','radius_exponent_m':kw['closed_radius_exponent'],'dimensional_status':'physical mapping with phenomenological geometry/accommodation' if 'GB_' not in mode else 'semi-phenomenological C_GB_closed','Q_closed_input':False,'validation':False})
  m=model(mode)
  for rate in (5,50):
   f,h=run_path(m,RampNoHold(rate,1500,start_C=950),initial(),300 if rate==5 else 60,600,f'{mode}_{rate}');hist.append(detailed(h,f'{mode}_ramp_{rate}',mode,'PDF_conditioned'))
   fixed.append({'mode':mode,'path':f'{rate}C_min','final_rho':f.rho,'final_G_um':f.G_m*1e6,'Delta_rho_open':integral(h,'rho_dot_open_sinv'),'Delta_rho_closed':min(integral(h,'rho_dot_closed_sinv'),.34),'closed_inventory_formed':max(h.closed_fraction*(1-h.rho)),'final_A_closed':f.A_closed})
  natural[mode]=prepare(mode)
  if seed is None:
   _,ph=run_path(m,ConditionedTwoStep(1400,1100,.88,1),initial(),600,600,'seed');pre=ph[ph.rho<.88];seed=pre.iloc[-1] if len(pre) else ph.iloc[-1]
 pd.DataFrame(registry).to_csv(OUT/'closed_pore_law_mapping_registry.csv',index=False)
 for mode in MODES:
  for kind,state in [('candidate_like_injected',injected(seed)),('naturally_prepared',natural[mode])]:
   for T2 in range(900,1301,25):
    f,h=run_path(model(mode),Iso(T2,40),replace(state,t_s=0.),1800,1800,f'{mode}_{kind}_{T2}');c=classify(f,state);hist.append(detailed(h,f'{mode}_{kind}_{T2}',mode,kind))
    scans.append({'mode':mode,'state_kind':kind,'T2_C':T2,'initial_rho':state.rho,'initial_G_nm':state.G_m*1e9,'initial_closed_fraction':state.pores.phi_closed.sum()/max(state.pores.total,1e-300),'initial_A_closed':state.A_closed,'final_rho':f.rho,'final_G_um':f.G_m*1e6,'Delta_rho_open':integral(h,'rho_dot_open_sinv'),'Delta_rho_closed':min(integral(h,'rho_dot_closed_sinv'),1-state.rho),'classification':c,'strict_success':c=='SUCCESS','diagnostic_injection':kind=='candidate_like_injected'})
 pd.DataFrame(fixed).to_csv(OUT/'fixed_path_summary.csv',index=False);scan=pd.DataFrame(scans);scan.to_csv(OUT/'T2_state_scan.csv',index=False);pd.concat(hist,ignore_index=True).to_csv(OUT/'closed_pore_bin_histories.csv',index=False)
 accept=[]
 for (mode,kind),g in scan.groupby(['mode','state_kind']):
  low=(g[g.T2_C<=950].final_rho<TARGET_RHO).any();mid=((g[g.T2_C.between(975,1175)].final_rho>=TARGET_RHO)&(g[g.T2_C.between(975,1175)].final_G_um<=TARGET_G)).any();high=(g[g.T2_C>=1200].final_G_um>TARGET_G).any();universal=(g.final_rho>=TARGET_RHO).all()
  accept.append({'mode':mode,'state_kind':kind,'low_density_failure':low,'intermediate_success':mid,'high_growth_failure':high,'universal_density_success':universal,'plausible_topology':low and mid and high and not universal,'acceptance':'candidate' if low and mid and high and not universal and kind=='naturally_prepared' else 'diagnostic_or_reject'})
 pd.DataFrame(accept).to_csv(OUT/'boundary_topology_acceptance.csv',index=False)
 # Apparent slopes are post-run diagnostics of state-normalized rate, never inputs.
 slopes=[];allh=pd.concat(hist,ignore_index=True)
 for mode,g in allh[(allh.rho_dot_closed_sinv>0)&(allh.closed_inventory>0)].groupby('mode'):
  k=g.rho_dot_closed_sinv/g.closed_inventory;coef=np.polyfit(1/(g.T_C+273.15),np.log(k),1) if len(g)>2 else [np.nan,np.nan];slopes.append({'mode':mode,'Q_closed_app_kJ_mol':-8.314462618*coef[0]/1000,'diagnostic_only':True,'material_input':False})
 pd.DataFrame(slopes).to_csv(OUT/'apparent_closed_rate_slopes.csv',index=False)
 state={'branch':'codex/zro2-forward-binwise-closed-pore-evolution','source_commit':'b70ac1b40b2b0c34d9791d61bcd160b61a17f92e','barrier_sha256':hashlib.sha256(BARRIER.read_bytes()).hexdigest(),'Q_closed_physical_input':False,'mini_map_run':False,'validation':False};(OUT/'run_state.json').write_text(json.dumps(state,indent=2)+'\n')
 print(pd.DataFrame(accept).to_string(index=False))
if __name__=='__main__':main()
