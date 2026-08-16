#!/usr/bin/env python3
from __future__ import annotations
from dataclasses import asdict,replace
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

import run_zro2_forward_closed_channel_physical_law_comparison as prior
from zro2_forward.conditioned_950c import ConditionedTwoStep,run_path,matched,BARRIER
from zro2_forward.resolved_rules import ResolvedRuleModel
from zro2_forward.schedules import RampNoHold,Iso

OUT=Path("results/zro2_forward_defined_law_parameter_mapping");TARGET_RHO=.976;TARGET_G=.29
MODES={
 "baseline_forward_current":dict(mechanism_mode="resolved_rules",closed_mapping_mode="baseline_forward_current"),
 "defined_laws_port":dict(mechanism_mode="defined_laws_port",closed_mapping_mode="defined_laws_port"),
 "reduced_candidate_law_transfer":dict(mechanism_mode="defined_laws_port",closed_mapping_mode="reduced_candidate_law_transfer",closed_mapping_rate_factor=26.),
 "mechanistic_GB_diffusion":dict(mechanism_mode="defined_laws_port",closed_mapping_mode="mechanistic_GB_diffusion"),
 "mechanistic_renewal_limited":dict(mechanism_mode="defined_laws_port",closed_mapping_mode="mechanistic_renewal_limited"),
 "mechanistic_gas_accommodation":dict(mechanism_mode="defined_laws_port",closed_mapping_mode="mechanistic_gas_accommodation"),
 "empirical_rate_scale_diagnostic":dict(mechanism_mode="defined_laws_port",closed_mapping_mode="empirical_rate_scale_diagnostic",closed_mapping_rate_factor=100.),
}
def model(mode):return ResolvedRuleModel(parameters=prior.params(**MODES[mode]))
def initial():return prior.initial()
def integ(h,col):return float(np.trapezoid(np.maximum(h[col],0),h.t_s)) if len(h)>1 else 0.
def classify(f,first,att=True):return prior.classify(f,first,att)
def registry():
 rows=[
 ("serial_renewal_densification","Two-barrier serial renewal","twostep_renewal_powerchannels.py","lines 226-236","tau_event=1/r_nuc+tau_sink; edot=eps_event/tau_event","open densification",True,False,False,"sigma,T,G,phi_open","Gstar,nu0,D_GB,Omega,b","triple-line density,eps_event","site_density_multiplier","sink_time_factor","","kinetic_state/density_rate","physical and directly implemented","","renewal appears once"),
 ("capillary_onsager_stress_balance","Capillary/Onsager work balance","run_sinter_onsager_export_pressureless_memory_3rate_aluminaGG_renewalpicard.m","lines 2116-2197","W_surf=(1-rho)sigma edot with bounded sigma","stress/power accounting",False,False,False,"rho,surface area","gamma_s,Gstar","pore area","","stress bounds","","solve_effective_stress","physical and directly implemented","","power does not set rho_dot directly"),
 ("excess_power_PR_desintering","Excess-power PR topology drive","twostep_renewal_powerchannels.py","lines 394-438","J_PR=k_PR F_low_activity F_excess F_topology phi","topology",False,True,False,"phi_open,activity,PR_memory","D_s","pore bins","","C_PR,Q_PR","","ResolvedRuleModel.rates","reduced phenomenological","absolute excess-power mapping","conservative"),
 ("conservative_pore_population_evolution","Conservative pore stores","twostep_renewal_powerchannels.py","lines 506-517","open loss equals isolated/closed gain except named shrinkage","pore topology",False,True,False,"phi_open,phi_iso,phi_closed","","bin geometry","","transfer rates","","ResolvedRuleModel.rates","physical and directly implemented","","density identity exact"),
 ("removable_to_trapped_or_closed_transition","Removable-to-precursor/closed transfer","discussion5 brief","closed transition rule","phi_open->phi_iso->phi_closed","topology",False,True,False,"phi_open,phi_iso,PR_memory","","pore bins","","transition times","","ResolvedRuleModel.rates","reduced phenomenological","transition kinetics",""),
 ("closed_pore_shrinkage_with_finite_accommodation","Finite-accommodation closed shrinkage","candidate laws / physical-law comparison","closed-law registry","rho_dot_closed=phi_closed A/tau_closed or mechanistic mapping","closed densification",True,False,False,"phi_closed,A_closed,r_i","D_GB,gamma_s,Omega","r_i","","A_max,recovery","k0_closed_emp","closed_channel_rates","physical but missing parameters","accommodation kinetics/pressure","bounded A"),
 ("Zener_pore_size_pinning","Zener pore-size pinning","twostep_renewal_powerchannels.py","lines 327-339","R_Z=k_Z r_p/f_v; Gamma_Z=max(0,1-G/(2R_Z))","grain growth",False,False,True,"pore sizes,volume","gamma_GB","r_p/f_v","zener_strength","","","growth_state","physical and directly implemented","k_Z transfer","smaller pores pin more"),
 ("intrinsic_growth_times_activity","Intrinsic growth times migration activity","discussion5 brief","growth rule","Gdot=M_GB gamma_GB/G * Gamma_migration","grain growth",False,False,True,"G,pore state,PR_memory","M_GB,gamma_GB","","M0 growth anchor","drag factors","","ResolvedRuleModel.rates","physical and directly implemented","mobility uncertainty","never enters density"),
 ("PR_prepared_closed_accommodation_memory","PR-prepared accommodation memory","candidate 693168 Tier B","required mechanism ablations","Adot=prepare(PR_memory)(Amax-A)-rho_dot_closed/capacity","closed preparation",False,False,False,"PR_memory,A_closed","","","","prepare/decay/capacity","","ResolvedRuleModel.step","reduced phenomenological","physical accommodation mapping","finite"),
 ("Chen_window_classification","Strict Chen window topology","resolved-rule map utilities","window_rows","lower density boundary + success band + upper growth boundary","classification",False,False,False,"final rho,G,T1,T2","targets","map spacing","","","","window_rows","computed diagnostic","","not local physics"),
 ("fast_firing_nucleation_onset","Fast-firing nucleation onset","discussion5 brief / barrier JSON","onset rule","r_nuc=nu0 exp[-Gstar(sigma,T)/kBT]","nucleation",False,False,False,"sigma,T","Gstar,nu0","","","","","kinetic_state","physical and directly implemented","sub-950C extrapolation",""),
 ("candidate693168_reduced_closed_channel","Candidate 693168 reduced channel","candidate 693168 Tier B","conditional comparator","PR preparation + transition + finite-A closed shrinkage","reduced comparator",True,True,False,"PR_memory,phi_closed,A_closed","","reduced stores","","all candidate rates/capacities","","reduced transfer mode","reduced phenomenological","not transferable as ZrO2 inputs","not validated"),]
 cols=["law_id","law_name","source_file","source_section_or_line_hint","governing_equation_text","affected_process","changes_density_directly","conservative_transfer","migration_only","state_variables","physical_inputs","geometry_inputs","fitted_or_calibrated_parameters","phenomenological_parameters","empirical_only_parameters","current_forward_function","implementation_status","missing_mapping","notes"]
 return pd.DataFrame(rows,columns=cols)
def mapping():
 fixed={"G*(sigma,T)":"fixed_ZrO2_input","b":"literature_input","nu0":"fixed_ZrO2_input","D_GB0":"literature_input","Q_GB":"literature_input","D_s0":"literature_input","Q_s":"literature_input","M0_growth":"global_calibration","Q_growth":"bounded_uncertainty"}
 geo={x:"geometry_derived" for x in ["C_TJ","C_GB","rho_TL_area","eps_event","Zener_R_Z","Zener_pinning_factor"]}
 state={x:"state_variable" for x in ["phi_open","phi_iso","phi_closed","phi_removable","phi_trapped","PR_memory","PR_work","precursor_closed","A_closed","A_closed_used","A_closed_recovered","r_i"]}
 diag={x:"computed_diagnostic" for x in ["sigma_eff","tau_nuc","tau_exchange","tau_transport","tau_sink","Lambda","activity","pore_D50","pore_D90","tau_remove_i","closed_shrinkage_rate","G_dot_intrinsic","G_dot_actual","Gamma_migration"]}
 special={"A_closed_max":"reduced_phenomenological","Q_closed_app":"empirical_diagnostic","k_closed_eff":"missing_physical_mapping"}
 rows=[]
 for name,cls in {**fixed,**geo,**state,**diag,**special}.items():
  source="MaterialParameters/barrier JSON" if cls in ("fixed_ZrO2_input","literature_input") else "current state/forward diagnostics" if cls in ("state_variable","computed_diagnostic") else "defined reduced/geometry law"
  rows.append({"parameter":name,"parameter_class":cls,"source":source,"used_by":"defined_laws_port","physical_ZrO2_input":cls in ("fixed_ZrO2_input","literature_input"),"requires_calibration":cls in ("global_calibration","reduced_phenomenological","missing_physical_mapping"),"notes":"apparent only; not a material property" if name=="Q_closed_app" else ""})
 for mode in MODES:
  if mode=="reduced_candidate_law_transfer":rows.append({"parameter":"candidate693168_rate_and_capacity_set","parameter_class":"reduced_phenomenological","source":"candidate 693168 conditional Tier B","used_by":mode,"physical_ZrO2_input":False,"requires_calibration":True,"notes":"comparator logic only; values not calibrated ZrO2 inputs"})
  if mode=="empirical_rate_scale_diagnostic":rows.append({"parameter":"Q_closed_emp_and_k0_closed_emp","parameter_class":"empirical_diagnostic","source":"physical-law comparison diagnostic","used_by":mode,"physical_ZrO2_input":False,"requires_calibration":False,"notes":"nonphysical and nonvalidated"})
 return pd.DataFrame(rows)
def injected_state(seed):return prior.state_from_row(seed,rho=.88,G_nm=117,closed_fraction=.649,A=.152,PR=1.)
def prepare(mode):
 m=model(mode);_,h=run_path(m,ConditionedTwoStep(1400,1100,.88,1),initial(),600,600,f"prepare_{mode}");pre=h[h.rho<.88];row=pre.iloc[-1] if len(pre) else h.iloc[-1];return prior.state_from_row(row),row
def scan(mode,state,kind):
 m=model(mode);rows=[]
 for T2 in range(900,1301,25):
  f,h=run_path(m,Iso(T2,40),replace(state,t_s=0.),1800,1800,mode);c=classify(f,state);closed=min(integ(h,"rho_dot_closed_sinv"),1-state.rho)
  rows.append({"mode":mode,"state_kind":kind,"T2_C":T2,"hold_h":40,"initial_rho":state.rho,"initial_G_nm":state.G_m*1e9,"closed_fraction_at_switch":state.pores.phi_closed.sum()/max(state.pores.total,1e-300),"A_closed_at_switch":state.A_closed,"PR_memory_at_switch":state.PR_memory,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"Delta_rho_open":integ(h,"rho_dot_open_sinv"),"Delta_rho_closed":closed,"PR_topology_transfer":integ(h,"PR_coarsening_flux"),"closed_inventory_formed":max(h.closed_fraction*(1-h.rho)),"final_A_closed":f.A_closed,"classification":c,"strict_success":c=="SUCCESS","candidate_state_injected":kind=="candidate_like_injected","diagnostic_only":kind=="candidate_like_injected"})
 return pd.DataFrame(rows)
def boundaries(x):
 rows=[]
 for key,g in x.groupby(["mode","T1_C","switch_density","hold_h"]):
  s=g[g.strict_success];dens=g[g.final_rho.ge(TARGET_RHO)];grain=g[g.final_G_um.le(TARGET_G)];lo=dens.T2_C.min() if len(dens) else np.nan;hi=grain.T2_C.max() if len(grain) else np.nan
  lp=bool(np.isfinite(lo) and (g[g.T2_C<lo].final_rho<TARGET_RHO).any());up=bool(np.isfinite(hi) and (g[g.T2_C>hi].final_G_um>TARGET_G).any());finite=bool(len(s)>=2 and lp and up and s.T2_C.max()>s.T2_C.min() and s.T2_C.max()<key[1])
  rows.append({"mode":key[0],"T1_C":key[1],"switch_density":key[2],"hold_h":key[3],"lower_boundary_present":lp,"upper_boundary_present":up,"success_count":len(s),"success_min_T2_C":s.T2_C.min() if len(s) else np.nan,"success_max_T2_C":s.T2_C.max() if len(s) else np.nan,"finite_window":finite,"window_width_C":s.T2_C.max()-s.T2_C.min() if finite else np.nan,"boundary_gap_C":hi-lo if np.isfinite(lo) and np.isfinite(hi) else np.nan})
 return pd.DataFrame(rows)
def mini_map(eligible):
 cols=["mode","T1_C","switch_density","T2_C","hold_h","final_rho","final_G_um","classification","strict_success"]
 rows=[]
 for mode in eligible:
  m=model(mode)
  for T1 in (1250,1300,1350,1400,1450,1500):
   for switch in (.75,.80,.85,.88,.90):
    _,prep=run_path(m,ConditionedTwoStep(T1,900,switch,.25),initial(),900,900,mode);second=prep[prep.T_C.eq(900)];att=len(second)>0;first=second.iloc[0] if att else prep.iloc[-1];st=prior.state_from_row(first) if att else None
    for hold in (20,40):
     for T2 in range(900,1301,25):
      if T2>=T1:continue
      if att:f,h=run_path(m,Iso(T2,hold),replace(st,t_s=0.),1800,1800,mode)
      else:f=initial()
      c=classify(f,first,att);rows.append({"mode":mode,"T1_C":T1,"switch_density":switch,"T2_C":T2,"hold_h":hold,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"classification":c,"strict_success":c=="SUCCESS"})
 return pd.DataFrame(rows,columns=cols)
def main():
 OUT.mkdir(parents=True,exist_ok=True);reg=registry();reg.to_csv(OUT/"defined_law_registry.csv",index=False);pm=mapping();pm.to_csv(OUT/"defined_law_parameter_mapping.csv",index=False)
 reg[["law_id","source_file","source_section_or_line_hint","current_forward_function","implementation_status","missing_mapping"]].to_csv(OUT/"source_law_traceability.csv",index=False)
 fixed=[];flux=[];prepared={};seed=None
 for mode in MODES:
  m=model(mode);hist={}
  for rate in (5,50):
   f,h=run_path(m,RampNoHold(rate,1500,start_C=950),initial(),300 if rate==5 else 60,300,mode);hist[rate]=h
   fixed.append({"mode":mode,"path":f"PDF_conditioned_{rate}C_min","final_rho":f.rho,"final_G_um":f.G_m*1e6,"closed_inventory_formed":max(h.closed_fraction*(1-h.rho)),"final_A_closed":f.A_closed,"PR_memory_final":f.PR_memory})
   flux.append({"mode":mode,"path":f"PDF_conditioned_{rate}C_min","Delta_rho_open":integ(h,"rho_dot_open_sinv"),"Delta_rho_closed":min(integ(h,"rho_dot_closed_sinv"),1-.66),"PR_topology_transfer":integ(h,"PR_coarsening_flux"),"closed_transition":integ(h,"closure_rate"),"state_density_gain":f.rho-.66})
  ratio=matched(hist[5],hist[50]).G_5_over_G_50.median()
  for r in fixed[-2:]:r["matched_density_G5_over_G50_median"]=ratio
  prepared[mode],row=prepare(mode)
  if seed is None:seed=row
 pd.DataFrame(fixed).to_csv(OUT/"fixed_path_mode_summary.csv",index=False);pd.DataFrame(flux).to_csv(OUT/"fixed_path_flux_integrals.csv",index=False)
 cand=pd.concat([scan(mode,injected_state(seed),"candidate_like_injected") for mode in MODES],ignore_index=True);nat=pd.concat([scan(mode,prepared[mode],"naturally_prepared_by_mode") for mode in MODES],ignore_index=True);cand.to_csv(OUT/"candidate_state_T2_scan_by_mode.csv",index=False);nat.to_csv(OUT/"natural_state_T2_scan_by_mode.csv",index=False)
 eligible=[]
 for mode,g in nat.groupby("mode"):
  low=(g[g.T2_C.le(950)].final_rho<TARGET_RHO).any();mid=((g[g.T2_C.between(975,1175)].final_rho>=TARGET_RHO-.005)&(g[g.T2_C.between(975,1175)].final_G_um<=TARGET_G*1.1)).any();high=(g[g.T2_C.ge(1200)].final_G_um>TARGET_G).any()
  if low and mid and high:eligible.append(mode)
 mm=mini_map(eligible);mm.to_csv(OUT/"mini_map_classification_points.csv",index=False)
 if len(mm):win=boundaries(mm)
 else:win=pd.DataFrame(columns=["mode","T1_C","switch_density","hold_h","lower_boundary_present","upper_boundary_present","success_count","success_min_T2_C","success_max_T2_C","finite_window","window_width_C","boundary_gap_C"])
 win.to_csv(OUT/"mini_map_window_boundaries.csv",index=False)
 if len(win):gap=win.groupby("mode").agg(groups=("mode","size"),strict_success_count=("success_count","sum"),finite_window_count=("finite_window","sum"),lower_boundary_count=("lower_boundary_present","sum"),upper_boundary_count=("upper_boundary_present","sum"),min_gap_C=("boundary_gap_C","min"),max_gap_C=("boundary_gap_C","max")).reset_index()
 else:gap=pd.DataFrame(columns=["mode","groups","strict_success_count","finite_window_count","lower_boundary_count","upper_boundary_count","min_gap_C","max_gap_C"])
 gap.to_csv(OUT/"mini_map_boundary_gap_summary.csv",index=False)
 comp=[]
 for mode,g in cand.groupby("mode"):
  s=g[g.strict_success];comp.append({"model":mode,"switch_density":.88,"switch_G_nm":117,"first_step_growth":117/50-1,"closed_fraction_at_switch":.649,"A_closed":.152,"PR_memory":1.,"closed_density_contribution":g.Delta_rho_closed.max(),"low_T2_failure":bool((g[g.T2_C.le(950)].final_rho<TARGET_RHO).any()),"success_band":f"{s.T2_C.min()}-{s.T2_C.max()}" if len(s) else "none","high_T2_growth_failure":bool((g[g.T2_C.ge(1200)].final_G_um>TARGET_G).any()),"finite_accommodation_role":"enabled and bounded","destructive_ablations":"not rerun; inherited Tier-B evidence","conditional_comparator":False,"validated":False})
 comp.append({"model":"candidate_693168","switch_density":.88,"switch_G_nm":117,"first_step_growth":.137,"closed_fraction_at_switch":.649,"A_closed":.152,"PR_memory":1.,"closed_density_contribution":.244,"low_T2_failure":True,"success_band":"conditional inherited interval","high_T2_growth_failure":True,"finite_accommodation_role":"required; infinite accommodation destructive","destructive_ablations":"PR damage, transition, shrinkage, finite accommodation required","conditional_comparator":True,"validated":False});pd.DataFrame(comp).to_csv(OUT/"candidate693168_defined_law_comparison.csv",index=False)
 status=pd.DataFrame([{"mapping":"renewal densification","status":"physical directly implemented","remaining":"eligibility/site multiplier globally calibrated"},{"mapping":"PR topology transfer","status":"conservative reduced phenomenological","remaining":"absolute excess-power-to-transfer coefficient"},{"mapping":"closed transition/accommodation","status":"reduced phenomenological plus mechanistic alternatives","remaining":"pressure, shrinkable fraction, capacity/recovery kinetics"},{"mapping":"Zener pinning","status":"geometry-derived form implemented","remaining":"k_Z transfer uncertainty"},{"mapping":"empirical closure","status":"diagnostic only","remaining":"not eligible as physical input"}]);status.to_csv(OUT/"phenomenological_to_mechanistic_mapping_status.csv",index=False)
 pm[pm.parameter_class.isin(["global_calibration","bounded_uncertainty","reduced_phenomenological","empirical_diagnostic","missing_physical_mapping"])].to_csv(OUT/"unresolved_parameter_mapping.csv",index=False)
 pd.DataFrame([{"decision":"defined-law transfer outcome","action":"retain trusted renewal/Onsager/open-pore laws; do not claim a Chen window"},{"decision":"closed mapping","action":"measure closed pressure, shrinkable fraction, transport length, and accommodation capacity/recovery before selecting a mechanistic port"},{"decision":"candidate 693168","action":"retain as conditional comparator only"},{"decision":"new mechanisms","action":"none proposed"}]).to_csv(OUT/"next_implementation_decision.csv",index=False)
 state={"branch":"codex/zro2-forward-port-defined-laws-parameter-mapping","source_commit":"83eab35f95a90e0d65db05d9cddee8dced8cf55b","barrier_sha256":hashlib.sha256(BARRIER.read_bytes()).hexdigest(),"discussion_doc_available":False,"discussion_rules_source":"continuation brief","comparison_modes":list(MODES),"mini_map_eligible_modes":eligible,"strict_success_count":int(mm.strict_success.sum()) if len(mm) else 0,"finite_window_count":int(win.finite_window.sum()) if len(win) else 0,"Q_closed_physical_input":False,"validation":False};(OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n")
 print(pd.DataFrame(fixed).to_string(index=False));print(state)
if __name__=="__main__":main()
