#!/usr/bin/env python3
from dataclasses import replace
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import run_zro2_forward_closed_channel_physical_law_comparison as base
from zro2_forward.conditioned_950c import ConditionedTwoStep,run_path,matched,BARRIER
from zro2_forward.resolved_rules import ResolvedRuleModel
from zro2_forward.schedules import RampNoHold,Iso
OUT=Path('results/zro2_forward_port_defined_closed_laws');TARGET_RHO=.976;TARGET_G=.29
CONFIGS={
 'resolved_proxy_current':dict(closed_mapping_mode='resolved_proxy_current'),
 'candidate_reduced_transfer':dict(closed_mapping_mode='candidate_reduced_transfer',closed_mapping_rate_factor=26.),
 'mechanistic_renewal_closed_m4':dict(closed_mapping_mode='mechanistic_renewal_closed',closed_radius_exponent=4),
 'mechanistic_renewal_closed_m3':dict(closed_mapping_mode='mechanistic_renewal_closed',closed_radius_exponent=3),
 'mechanistic_GB_diffusion_closed_m4':dict(closed_mapping_mode='mechanistic_GB_diffusion_closed',closed_radius_exponent=4),
 'mechanistic_GB_diffusion_closed_m3':dict(closed_mapping_mode='mechanistic_GB_diffusion_closed',closed_radius_exponent=3),
 'mechanistic_surface_accommodation':dict(closed_mapping_mode='mechanistic_surface_accommodation'),
 'mechanistic_gas_accommodation':dict(closed_mapping_mode='mechanistic_gas_accommodation',closed_radius_exponent=4),
 'empirical_closed_rate_diagnostic':dict(closed_mapping_mode='empirical_closed_rate_diagnostic',closed_mapping_rate_factor=100.),}
def model(config):return ResolvedRuleModel(parameters=base.params(mechanism_mode='defined_closed_laws_port',**CONFIGS[config]))
def initial():return base.initial()
def integ(h,c):return float(np.trapezoid(np.maximum(h[c],0),h.t_s)) if len(h)>1 else 0.
def persisted_history(h,run_id,mode,kind):
 x=h.copy();x.insert(0,'run_id',run_id);x.insert(1,'mode',mode);x.insert(2,'state_kind',kind);dt=x.t_s.diff().fillna(x.t_s.iloc[0]);x['cumulative_open_density_gain']=np.cumsum(x.rho_dot_open_sinv*dt);x['cumulative_closed_density_gain']=np.cumsum(x.rho_dot_closed_sinv*dt)
 aliases={'G':'G_um','phi_open_i':'phi_open_json','phi_iso_i':'phi_iso_json','phi_closed_i':'phi_closed_json','N_closed_i':'N_closed_i_json','r_closed_i':'r_closed_i_json','chi_shrink_i':'chi_shrink_i_json','C_geom_i':'C_geom_i_json','P_gas_i':'P_gas_i_json','sigma_closed_i':'sigma_closed_i_json','Gstar_closed_i':'Gstar_closed_i_json','r_nuc_closed_i':'r_nuc_closed_i_json','tau_transport_closed_i':'tau_transport_closed_i_json','tau_cycle_closed_i':'tau_cycle_closed_i_json','A_closed_i':'A_closed_i_json','A_closed_max_i':'A_closed_max_i_json','A_closed_available_i':'A_closed_available_i_json','A_closed_used_i':'A_closed_used_i_json','A_closed_recovered_i':'A_closed_recovered_i_json','rho_dot_open_i':'rho_dot_open_i_json','rho_dot_closed_i':'rho_dot_closed_i_json','PR_work':'cumulative_PR_work','precursor_closed_i':'precursor_closed_i_json','surface_accommodation_rate':'surface_diffusion_accommodation_rate_i_json','gas_pressure_factor':'gas_pressure_factor_i_json'}
 for dst,src in aliases.items():x[dst]=x[src]
 keep=['run_id','mode','state_kind','T_C','t_s','rho','G',*aliases.keys(),'cumulative_open_density_gain','cumulative_closed_density_gain','PR_memory','closed_law_mode','radius_exponent_m','rho_dot_open_sinv','rho_dot_closed_sinv','closed_inventory']
 return x.loc[:,list(dict.fromkeys(keep))]
def registry():
 rows=[
 ('renewal_open_densification','twostep_renewal_powerchannels.py','renewal_edot lines 226-236','r_nuc; tau_cycle=1/r_nuc+tau_exchange+tau_transport; activity=Lambda/(1+Lambda)','open density',True,False,False,'phi_open,sigma,T,G','Gstar,nu0,D_GB,Omega,b','triple-line geometry','pore size/connectivity','site multiplier,sink factor','','kinetic_state + density_rate','','physical directly implemented'),
 ('capillary_onsager_stress_balance','renewal/Onsager MATLAB export','solve_renewal_onsager_picard lines 2116-2197','W_surf=(1-rho)sigma edot with bounded sigma','stress accounting',False,False,False,'rho,surface area','gamma_s,Gstar','surface area','pore area','stress bounds','','solve_effective_stress','','physical directly implemented'),
 ('conservative_PR_topology_preparation','discussion brief / powerchannels','PR lines 394-438','J_PR=k_PR F_low_activity F_excess F_topology phi_open','topology',False,True,False,'phi_open,PR_memory,activity','D_s','pore bins','pore distribution','C_PR,Q_PR','','conservative_adjacent_PR','power-to-transfer mapping','reduced phenomenological'),
 ('open_to_precursor_to_closed_transition','candidate Tier-B reports','required transition structure','phi_open -> phi_iso -> phi_closed','topology',False,True,False,'phi_open,phi_iso,phi_closed,PR_memory','','pore bins','connectivity/isolation','transition times','','ResolvedRuleModel.rates','transition kinetics','reduced phenomenological'),
 ('finite_accommodation_state','closed-pore addendum','finite accommodation section','A_available=clamp(Amax-Aused+Arecovered); dAused=beta|dphi_removed|','closed accommodation',False,False,False,'A_used,A_recovered,Amax','D_s for shape recovery','bin geometry','accommodation state','beta_A,k_recover','','ClosedPoreEvolution','capacity/recovery measurements','physical state with phenomenological kinetics'),
 ('renewal_limited_closed_shrinkage','closed-pore addendum','renewal closed equation','phi_closed chi A (r_ref/r)^m / tau_cycle_closed','closed density',True,False,False,'phi_closed,N_closed,r,A','Gstar,nu0,D_GB,Omega,gamma_s','C_geom,r,ell','closed size/gas/accommodation','exchange time,transport length','','closed_channel_rates','exchange/length/shrinkability','physical but missing parameters'),
 ('GB_diffusion_closed_shrinkage','closed-pore addendum','GB alternative','C_GB phi chi A D_GB Omega sigma/(kBT) r^-m','closed density',True,False,False,'phi_closed,r,A','D_GB,Omega,gamma_s','C_geom,r','closed size/gas/accommodation','C_GB_closed','','closed_channel_rates','dimensional prefactor derivation','semi-phenomenological'),
 ('surface_diffusion_accommodation_only','closed-pore addendum','surface accommodation','Adot=k_shape D_s r^-4(Amax-A)','shape/accommodation',False,False,False,'A,r','D_s','r','shape/accommodation','k_shape','','closed_channel_rates/advance','k_shape mapping','physical transport; non-densifying'),
 ('gas_accommodation_pressure_limit','closed-pore addendum','gas option','Pgas=P0(V0/V)(T/T0); sigma=max(Cgeom 2gamma/r-Pgas,0)','closed driving stress',False,False,False,'Pgas,V,r,T','gamma_s','C_geom,r','gas content/volume','initial gas fraction','','ClosedPoreEvolution.gas_pressure_Pa','gas content/escape','physical proxy missing measurements'),
 ('Zener_pore_size_pinning','twostep_renewal_powerchannels.py','Zener lines 327-339','R_Z=k_Z r_p/f_v; growth multiplied by pinning','migration',False,False,True,'pore size,volume,G','gamma_GB','r/f_v','pore distribution','k_Z','','growth_state','k_Z transfer','geometry law implemented'),
 ('intrinsic_growth_times_activity','discussion brief','growth rule','Gdot=M_GB gamma_GB/G * Gamma_migration','migration',False,False,True,'G,pore state,PR_memory','M_GB,gamma_GB','','microstructure','growth anchor,drag','','ResolvedRuleModel.rates','mobility uncertainty','physical with bounded uncertainty'),
 ('Chen_window_classification','resolved map utilities','window_rows','lower density failure + success band + upper growth failure','classification',False,False,False,'rho,G,T1,T2','fixed targets','map spacing','trajectory','','','window classifier','','computed diagnostic'),
 ('candidate693168_reduced_closed_channel','candidate 693168 Tier B','required/destructive ablations','PR preparation + transition + closed shrinkage + finite A','reduced comparator',True,True,False,'PR_memory,phi_closed,A','','reduced geometry','candidate state','all candidate rates/capacities','','candidate_reduced_transfer','not transferable as inputs','reduced phenomenological; not validated')]
 cols=['law_id','source_document','source_code_or_section','governing_equation','affected_process','changes_density_directly','conservative_transfer','migration_only','state_variables','physical_inputs','geometry_inputs','measurable_microstructure_inputs','phenomenological_parameters','empirical_diagnostic_parameters','current_forward_mapping','missing_mapping','implementation_status'];return pd.DataFrame(rows,columns=cols)
def parameter_mapping():
 groups={
 'fixed_ZrO2_input':['Gstar_sigma_T','nu0'], 'literature_input':['b','D_GB0','Q_GB','D_s0','Q_s','gamma_s','gamma_GB','Omega'],
 'geometry_derived':['C_TJ','C_GB','rho_TL_area','eps_event','sigma_closed_i','R_Z'],
 'evolved_state_variable':['phi_open_i','phi_iso_i','phi_closed_i','N_closed_i','r_closed_i','P_gas_i','A_closed_i','A_closed_used_i','A_closed_recovered_i','PR_memory','PR_work','precursor_closed_i'],
 'computed_diagnostic':['sigma_eff','tau_nuc','tau_transport','tau_sink','Lambda','activity','rho_dot_open_i','rho_dot_closed_i','Q_closed_app','Gamma_migration','G_dot_intrinsic','G_dot_actual'],
 'global_calibration':['M0_growth'],'bounded_uncertainty':['Q_growth'],
 'reduced_phenomenological':['A_closed_max_i'], 'missing_physical_mapping':['tau_exchange','chi_shrink_i','C_geom_i','k_closed_eff']}
 rows=[]
 for cls,names in groups.items():
  for n in names:rows.append({'parameter_name':n,'parameter_class':cls,'physical_source':'barrier/material input' if cls in ('fixed_ZrO2_input','literature_input') else 'defined law/current state','forward_location':'defined_closed_laws_port','calibration_required':cls in ('global_calibration','reduced_phenomenological','missing_physical_mapping'),'notes':'post-run apparent diagnostic only' if n=='Q_closed_app' else ''})
 rows.append({'parameter_name':'Q_closed_emp,k0_closed_emp','parameter_class':'empirical_diagnostic','physical_source':'none','forward_location':'empirical_closed_rate_diagnostic','calibration_required':False,'notes':'nonphysical diagnostic only'})
 rows.append({'parameter_name':'candidate693168_rate_capacity_set','parameter_class':'reduced_phenomenological','physical_source':'conditional Tier B comparator','forward_location':'candidate_reduced_transfer','calibration_required':True,'notes':'not calibrated ZrO2 inputs'})
 return pd.DataFrame(rows)
def prepare(config):
 _,h=run_path(model(config),ConditionedTwoStep(1400,1100,.88,1),initial(),600,600,'prepare');pre=h[h.rho<.88];row=pre.iloc[-1] if len(pre) else h.iloc[-1];return base.state_from_row(row),row
def inject(row):return base.state_from_row(row,rho=.88,G_nm=117,closed_fraction=.649,A=.152,PR=1.)
def scan(config,state,kind,hist):
 rows=[]
 for T2 in range(900,1301,25):
  f,h=run_path(model(config),Iso(T2,40),replace(state,t_s=0.),1800,1800,config);c=base.classify(f,state,True);hist.append(persisted_history(h,f'{config}_{kind}_{T2}',config,kind))
  rows.append({'mode':config,'state_kind':kind,'T2_C':T2,'initial_rho':state.rho,'initial_G_nm':state.G_m*1e9,'closed_fraction_at_switch':state.pores.phi_closed.sum()/max(state.pores.total,1e-300),'A_closed_at_switch':state.A_closed,'PR_memory_at_switch':state.PR_memory,'final_rho':f.rho,'final_G_um':f.G_m*1e6,'Delta_rho_open':integ(h,'rho_dot_open_sinv'),'Delta_rho_closed':min(integ(h,'rho_dot_closed_sinv'),1-state.rho),'closed_inventory_formed':max(h.closed_fraction*(1-h.rho)),'final_A_closed':f.A_closed,'classification':c,'strict_success':c=='SUCCESS','candidate_state_injected':kind=='candidate_like_injected','diagnostic_only':kind=='candidate_like_injected'})
 return pd.DataFrame(rows)
def main():
 OUT.mkdir(parents=True,exist_ok=True);reg=registry();reg.to_csv(OUT/'defined_closed_law_registry.csv',index=False);pm=parameter_mapping();pm.to_csv(OUT/'defined_closed_parameter_mapping.csv',index=False);reg[['law_id','source_document','source_code_or_section','current_forward_mapping','implementation_status','missing_mapping']].to_csv(OUT/'source_law_traceability.csv',index=False)
 fixed=[];flux=[];hist=[];natural={};seed=None
 for config in CONFIGS:
  hs={}
  for rate in (5,50):
   f,h=run_path(model(config),RampNoHold(rate,1500,start_C=950),initial(),300 if rate==5 else 60,600,config);hs[rate]=h;hist.append(persisted_history(h,f'{config}_ramp_{rate}',config,'PDF_conditioned'))
   fixed.append({'mode':config,'path':f'PDF_conditioned_{rate}C_min','final_rho':f.rho,'final_G_um':f.G_m*1e6,'closed_inventory_formed':max(h.closed_fraction*(1-h.rho)),'final_A_closed':f.A_closed})
   flux.append({'mode':config,'path':f'PDF_conditioned_{rate}C_min','Delta_rho_open':integ(h,'rho_dot_open_sinv'),'Delta_rho_closed':min(integ(h,'rho_dot_closed_sinv'),.34),'PR_topology_transfer':integ(h,'PR_coarsening_flux'),'closed_transition':integ(h,'closure_rate'),'state_density_gain':f.rho-.66})
  ratio=matched(hs[5],hs[50]).G_5_over_G_50.median()
  for r in fixed[-2:]:r['matched_density_G5_over_G50_median']=ratio
  natural[config],row=prepare(config)
  if seed is None:seed=row
 pd.DataFrame(fixed).to_csv(OUT/'fixed_path_mode_summary.csv',index=False);pd.DataFrame(flux).to_csv(OUT/'fixed_path_flux_integrals.csv',index=False)
 cand=pd.concat([scan(c,inject(seed),'candidate_like_injected',hist) for c in CONFIGS],ignore_index=True);nat=pd.concat([scan(c,natural[c],'naturally_prepared',hist) for c in CONFIGS],ignore_index=True);cand.to_csv(OUT/'candidate_state_T2_scan_by_mode.csv',index=False);nat.to_csv(OUT/'natural_state_T2_scan_by_mode.csv',index=False);allh=pd.concat(hist,ignore_index=True);allh.to_csv(OUT/'closed_law_bin_histories.csv.gz',index=False,compression='gzip')
 slopes=[]
 for mode,g in allh[(allh.rho_dot_closed_sinv>0)&(allh.closed_inventory>0)].groupby('mode'):
  k=g.rho_dot_closed_sinv/g.closed_inventory;coef=np.polyfit(1/(g.T_C+273.15),np.log(k),1) if len(g)>2 else [np.nan,np.nan];slopes.append({'mode':mode,'Q_closed_app_kJ_mol':-8.314462618*coef[0]/1000,'diagnostic_only':True,'material_input':False})
 pd.DataFrame(slopes).to_csv(OUT/'apparent_closed_rate_slopes.csv',index=False)
 eligible=[]
 for mode,g in nat.groupby('mode'):
  low=(g[g.T2_C<=950].final_rho<TARGET_RHO).any();mid=((g[g.T2_C.between(975,1175)].final_rho>=TARGET_RHO)&(g[g.T2_C.between(975,1175)].final_G_um<=TARGET_G)).any();high=(g[g.T2_C>=1200].final_G_um>TARGET_G).any()
  if low and mid and high:eligible.append(mode)
 cols=['mode','T1_C','switch_density','T2_C','hold_h','final_rho','final_G_um','classification','strict_success'];mm=pd.DataFrame(columns=cols);mm.to_csv(OUT/'mini_map_classification_points.csv',index=False)
 wcols=['mode','T1_C','switch_density','hold_h','lower_boundary_present','upper_boundary_present','success_count','finite_window','window_width_C','boundary_gap_C'];win=pd.DataFrame(columns=wcols);win.to_csv(OUT/'mini_map_window_boundaries.csv',index=False);pd.DataFrame(columns=['mode','strict_success_count','finite_window_count','lower_boundary_count','upper_boundary_count']).to_csv(OUT/'mini_map_boundary_gap_summary.csv',index=False)
 central=[
 ('PR_damage_rate','candidate range','surface-diffusion/excess-power topology preparation','C_PR,Q_PR,low-activity gate','phenomenological_unmapped','PR/topology-resolved experiment',True,False),
 ('closed_transition_rate','candidate range','precursor-to-closed topology change','closed_transition_tau','phenomenological_unmapped','closed fraction versus first-step history',True,False),
 ('closed_shrinkage_transport','candidate rate scale','GB transport or serial renewal','D_GB,Gstar,r,sigma_closed','mechanistically_derived','closed pore radius/pressure/transport length',True,False),
 ('finite_accommodation_capacity','candidate finite value','available rearrangement capacity','A_max,A_used,A_recovered','measurable_state','shape/coordination evolution',True,False),
 ('accommodation_consumption','candidate coefficient','capacity consumed per removed volume','beta_A','phenomenological_unmapped','in-situ pore shape versus volume',True,False),
 ('accommodation_recovery','candidate rate','surface/topology recovery','k_A_recover,D_s','phenomenological_unmapped','hold-temperature recovery data',True,False),
 ('closed_pore_radius','candidate reduced size','measured/evolved pore radius','r_closed_i,N_closed_i','measurable_state','3D pore size/count data',True,False),
 ('geometry_shrinkability','candidate topology factor','coordination/dihedral factor','C_geom_i,chi_shrink_i','geometry_derived','3D coordination/dihedral data',True,False),
 ('gas_pressure','not explicit in candidate','trapped-gas opposition to capillarity','P_gas_i','measurable_state','gas content and pore volume',True,False),
 ('Q_closed_app','apparent only','post-run composite slope','diagnostic fit','empirical_only','none; never physical input',False,False),
 ('intrinsic_growth_mobility','candidate reduced growth branch','M_GB gamma/G','fixed/bounded growth inputs','direct_physical','high-T grain-growth data',False,True)]
 cm=['parameter_name','reduced_value_or_range','physical_interpretation','forward_model_equivalent','mapping_status','data_needed','controls_lower_boundary','controls_upper_boundary'];pd.DataFrame(central,columns=cm).to_csv(OUT/'phenomenological_to_mechanistic_mapping_status.csv',index=False)
 comp=[]
 for mode,g in cand.groupby('mode'):
  s=g[g.strict_success];comp.append({'model':mode,'switch_density':.88,'switch_G_nm':117,'closed_fraction_at_switch':.649,'A_closed':.152,'PR_memory':1.,'closed_density_contribution':g.Delta_rho_closed.max(),'low_T2_failure':bool((g[g.T2_C<=950].final_rho<TARGET_RHO).any()),'success_band':f'{s.T2_C.min()}-{s.T2_C.max()}' if len(s) else 'none','high_T2_growth_failure':bool((g[g.T2_C>=1200].final_G_um>TARGET_G).any()),'finite_accommodation_role':'bounded binwise state','conditional_comparator':False,'validated':False})
 comp.append({'model':'candidate_693168','switch_density':.88,'switch_G_nm':117,'closed_fraction_at_switch':.649,'A_closed':.152,'PR_memory':1.,'closed_density_contribution':.244,'low_T2_failure':True,'success_band':'conditional inherited interval','high_T2_growth_failure':True,'finite_accommodation_role':'required; infinite accommodation destructive','conditional_comparator':True,'validated':False});pd.DataFrame(comp).to_csv(OUT/'candidate693168_defined_law_comparison.csv',index=False)
 pm[pm.parameter_class.isin(['global_calibration','bounded_uncertainty','reduced_phenomenological','empirical_diagnostic','missing_physical_mapping'])].to_csv(OUT/'unresolved_parameter_mapping.csv',index=False)
 pd.DataFrame([{'decision':'law transfer','action':'already-defined laws ported; no new mechanism'},{'decision':'state preparation','action':'measure PR-to-closed inventory and geometry evolution'},{'decision':'closed mapping','action':'constrain pressure, shrinkability, exchange/transport length, accommodation capacity/recovery'},{'decision':'candidate 693168','action':'conditional comparator only'}]).to_csv(OUT/'next_implementation_decision.csv',index=False)
 state={'branch':'codex/zro2-forward-port-defined-closed-laws','source_commit':'986aabb','prior_binwise_branch_preserved':True,'barrier_sha256':hashlib.sha256(BARRIER.read_bytes()).hexdigest(),'discussion_doc_available':False,'eligible_modes':eligible,'mini_map_run':False,'strict_success_count':0,'finite_window_count':0,'Q_closed_physical_input':False,'validation':False};(OUT/'run_state.json').write_text(json.dumps(state,indent=2)+'\n');print(state)
if __name__=='__main__':main()
