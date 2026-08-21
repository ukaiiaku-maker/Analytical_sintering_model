#!/usr/bin/env python3
"""Fixed-physics ZrO2 heating-rate and two-step processing-window predictions."""
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
import hashlib,json,math,os
import numpy as np,pandas as pd
from zro2_forward.conditioned_950c import make_pdf_conditioned_initial_state,run_path,BARRIER,matched
from zro2_forward.resolved_rules import ResolvedRuleModel,ResolvedRuleParameters,resolved_initial_state
from zro2_forward.schedules import RampNoHold,RampHold,Iso
from zro2_forward.pore_population import diagnostics

OUT=Path("results/zro2_forward_processing_window_prediction_figures");TARGETS=(.90,.95,.98)
RHO0=[.55,.60,.65,.66,.70,.75,.80];G0=[10.2,15,24.5,33,50,75,100,150,225,300];PORE=[10,15,24.5,33,40,60];WIDTH=[.45,.65,.85];T1=[1200,1250,1300,1350,1400,1450,1500,1550];RATE=[1,5,10,20,50,100];SWITCH=[.72,.75,.78,.80,.83,.85,.88,.90,.92,.95];T2=list(range(850,1351,25));HOLDS=[5,10,20,30,40,96]
def params():return replace(ResolvedRuleParameters(),site_density_multiplier=100.,surface_power_length2_m2=3e-20,mechanism_mode="defined_closed_laws_port",closed_mapping_mode="resolved_proxy_current",gb_mobility_mode="current_default")
def model():return ResolvedRuleModel(parameters=params())
def initial(r):return resolved_initial_state(make_pdf_conditioned_initial_state(rho=r.rho0,G_nm=r.G0_nm,T_C=950,pore_D50_nm=r.pore_D50_nm,pore_log_width=r.pore_ln_sigma,phi_iso_fraction=0,phi_closed_fraction=0))
class RampToHold:
 def __init__(self,rate,T,start=950,h=20):self.rate=rate/60;self.peak=T;self.start=start;self.t_end_s=(T-start)/self.rate+h*3600
 def temperature_K(self,t,rho):return min(self.peak,self.start+self.rate*t)+273.15
def run_switch(r,keep_history=False):
 m=model();s=initial(r);path=RampToHold(r.rate_T1_C_min,r.T1_C);rows=[];nextrec=0.;att=False
 while s.t_s<path.t_end_s and s.rho<r.rho_switch:
  T=path.temperature_K(s.t_s,s.rho);trial=min(600.,path.t_end_s-s.t_s);od,ii,cc,_,growth,_=m.rates(s,T);mx=max(float(np.max(np.abs(od))),float(np.max(np.abs(ii))),float(np.max(np.abs(cc))),1e-300);trial=max(min(trial,2e-3/mx,.01*s.G_m/max(growth["G_dot_m_s"],1e-300),path.t_end_s-s.t_s),1e-6);s,d=m.step(s,T,trial)
  if keep_history and (s.t_s>=nextrec or s.rho>=r.rho_switch):rows.append({"t_s":s.t_s,"T_C":T-273.15,"rho":s.rho,"G_nm":s.G_m*1e9,**d});nextrec=s.t_s+600
 att=s.rho>=r.rho_switch
 return s,rows,att
def screen_one(d):
 r=type("R",(),d);s,_,att=run_switch(r);q=diagnostics(s.pores);growth=s.G_m/(r.G0_nm*1e-9)-1;closed=float(s.pores.phi_closed.sum());return {**d,"switch_reached":att,"actual_rho_switch":s.rho,"G1_nm":s.G_m*1e9,"first_step_growth_fraction":growth,"phi_open_total":float(s.pores.phi_open.sum()),"phi_iso_total":float(s.pores.phi_iso.sum()),"phi_closed_total":closed,"A_closed":s.A_closed,"PR_memory":s.PR_memory,"pore_D50_nm_switch":q["pore_D50_m"]*1e9,"pore_D90_nm_switch":q["pore_D90_m"]*1e9,"remaining_open_fraction":float(s.pores.phi_open.sum()/max(s.pores.total,1e-300)),"already_at_final_target":s.rho>=.98,"promotable":bool(att and growth<=.20 and .75<=r.rho_switch<=.92 and (closed>0 or s.A_closed>0 or s.PR_memory>0) and s.rho<.98)}
def grid():
 rng=np.random.default_rng(224902);rows=[]
 for i in range(5040):
  vals=[x[rng.integers(len(x))] for x in (RHO0,G0,PORE,WIDTH,T1,RATE,SWITCH)];rows.append(dict(initial_state_id=f"S{i:05d}",rho0=vals[0],G0_nm=vals[1],pore_D50_nm=vals[2],pore_ln_sigma=vals[3],T1_C=vals[4],rate_T1_C_min=vals[5],rho_switch=vals[6]))
 return pd.DataFrame(rows).drop_duplicates(subset=["rho0","G0_nm","pore_D50_nm","pore_ln_sigma","T1_C","rate_T1_C_min","rho_switch"]).reset_index(drop=True)
def history_frame(h,meta):
 x=h.copy();x["physical_time_s"]=x.t_s;x["physical_time_h"]=x.t_s/3600;x["G_nm"]=x.G_m*1e9;x["G_um"]=x.G_m*1e6
 alias={"phi_open_total":"open_fraction","phi_iso_total":"isolated_fraction","phi_closed_total":"closed_fraction","pore_D50_nm":"pore_D50_m","pore_D90_nm":"pore_D90_m","connected_fine_pore_fraction":"fine_pore_fraction","large_pore_fraction":"large_pore_fraction","rho_dot_open":"rho_dot_open_sinv","rho_dot_closed":"rho_dot_closed_sinv","rho_dot_total":"rho_dot_total_sinv","PR_work":"cumulative_PR_work","A_closed_available":"A_closed","A_closed_used":"A_closed_used_total","A_closed_recovered":"A_closed_recovered_total","activity":"activity","Gstar_eV":"Gstar_eV","sigma_eff":"sigma_eff_Pa","tau_nuc":"tau_nuc_s","tau_exchange":"tau_exchange_s","tau_transport":"tau_transport_s","tau_cycle":"tau_cycle_s","Gamma_migration":"Gamma_migration","G_dot_intrinsic":"G_dot_intrinsic_m_s","G_dot_actual":"G_dot_actual_m_s","closed_law_mode":"closed_law_mode"}
 for dst,src in alias.items():x[dst]=x[src]*1e9 if dst.startswith("pore_D") and src in x else x[src] if src in x else np.nan
 x["cumulative_open_density_gain"]=np.cumsum(x.rho_dot_open.fillna(0)*x.physical_time_s.diff().fillna(0));x["cumulative_closed_density_gain"]=np.cumsum(x.rho_dot_closed.fillna(0)*x.physical_time_s.diff().fillna(0));x["mechanism_mode"]="defined_closed_laws_port";x["mobility_case"]="current_default"
 for k,v in meta.items():x[k]=v
 keep=["run_id","case_id","initial_state_id","schedule_id","path_type","physical_time_s","physical_time_h","T_C","rho","G_nm","G_um","phi_open_total","phi_iso_total","phi_closed_total","pore_D50_nm","pore_D90_nm","connected_fine_pore_fraction","large_pore_fraction","rho_dot_open","rho_dot_closed","rho_dot_total","cumulative_open_density_gain","cumulative_closed_density_gain","PR_memory","PR_work","A_closed_available","A_closed_used","A_closed_recovered","activity","Gstar_eV","sigma_eff","tau_nuc","tau_exchange","tau_transport","tau_cycle","Gamma_migration","G_dot_intrinsic","G_dot_actual","mechanism_mode","closed_law_mode","mobility_case"]
 return x[[c for c in keep if c in x]]
def classify(final,first,target,tol,att=True):
 if not att:return "UNATTAINABLE_FIRST_STEP"
 if first.rho>=target:return "INELIGIBLE_TARGET_ALREADY_REACHED"
 d=final.rho>=target;g=final.G_m/first.G_m-1<=tol
 return "SUCCESS" if d and g else "DENSIFICATION_EXHAUSTION_FAILURE" if not d and g else "GRAIN_GROWTH_FAILURE" if d else "MIXED_FAILURE"
def scan_job(a):
 rr,first,state_hash,target,hold,t2=a;f,h=run_path(model(),Iso(t2,hold),first,3600,7200,f"{rr['case_id']}_{t2}_{hold}");rows=[]
 for tier,tol in [("Tier_like_A",.05),("Tier_like_B",.20),("Diagnostic",.30)]:
  cl=classify(f,first,target,tol,True);rows.append(dict(case_id=rr["case_id"],initial_state_id=rr["initial_state_id"],switch_state_hash=state_hash,rho0=rr["rho0"],G0_nm=rr["G0_nm"],pore_D50_nm=rr["pore_D50_nm"],T1_C=rr["T1_C"],rate_T1_C_min=rr["rate_T1_C_min"],rho_switch=rr["rho_switch"],G1_nm=first.G_m*1e9,T2_C=t2,hold_h=hold,density_target=target,tolerance_tier=tier,grain_tolerance=tol,final_rho=f.rho,final_G_nm=f.G_m*1e9,second_step_growth_fraction=f.G_m/first.G_m-1,classification=cl,strict_success=cl=="SUCCESS" and tier!="Diagnostic",physical_start_time_s=first.t_s,physical_end_time_s=f.t_s))
 return rows
def second_scan(promoted):
 pts=[];states=[];reph=[];tasks=[];prepared=[]
 for rr in promoted.itertuples():
  first,fh,att=run_switch(rr,True);state_hash=hashlib.sha256(np.r_[first.rho,first.G_m,first.pores.phi_open,first.pores.phi_iso,first.pores.phi_closed,first.A_closed,first.PR_memory].tobytes()).hexdigest()[:16];states.append(dict(case_id=rr.case_id,initial_state_id=rr.initial_state_id,switch_state_hash=state_hash,rho_switch_requested=rr.rho_switch,rho_switch_actual=first.rho,G1_nm=first.G_m*1e9,first_step_end_time_s=first.t_s,first_step_growth_fraction=first.G_m/(rr.G0_nm*1e-9)-1,phi_closed_total=float(first.pores.phi_closed.sum()),A_closed=first.A_closed,PR_memory=first.PR_memory))
  rd=rr._asdict();prepared.append((rd,first,state_hash))
  for target in TARGETS:
   for hold in HOLDS:
    for t2 in T2:
     tasks.append((rd,first,state_hash,target,hold,t2))
  # exact first-step history retained once per promoted state
  if fh:reph.append(history_frame(pd.DataFrame(fh).rename(columns={"G_nm":"G_m"}).assign(G_m=lambda z:z.G_m*1e-9),dict(run_id=f"{rr.case_id}_first",case_id=rr.case_id,initial_state_id=rr.initial_state_id,schedule_id="first_step",path_type="first_step")))
 with ProcessPoolExecutor(max_workers=min(8,os.cpu_count() or 2)) as ex:
  for rows in ex.map(scan_job,tasks,chunksize=2):pts.extend(rows)
 for rd,first,state_hash in prepared[:3]:
  for t2 in (850,1050,1350):
   f,h=run_path(model(),Iso(t2,40),first,1800,3600,f"{rd['case_id']}_{t2}_40");h.t_s+=first.t_s;reph.append(history_frame(h,dict(run_id=f"{rd['case_id']}_{t2}",case_id=rd['case_id'],initial_state_id=rd['initial_state_id'],schedule_id=f"T2_{t2}_40h",path_type="two_step")))
 return pd.DataFrame(pts),pd.DataFrame(states),pd.concat(reph,ignore_index=True) if reph else pd.DataFrame()
def boundaries(pts):
 rows=[]
 for keys,g in pts.groupby(["case_id","density_target","tolerance_tier","hold_h"]):
  s=g[g.classification.eq("SUCCESS")].sort_values("T2_C");lo=g[g.classification.eq("DENSIFICATION_EXHAUSTION_FAILURE")];hi=g[g.classification.eq("GRAIN_GROWTH_FAILURE")];finite=len(s)>=2 and len(lo[lo.T2_C<s.T2_C.min()])>0 and len(hi[hi.T2_C>s.T2_C.max()])>0
  rows.append(dict(case_id=keys[0],density_target=keys[1],tolerance_tier=keys[2],hold_h=keys[3],G1_nm=g.G1_nm.iloc[0],lower_boundary_C=s.T2_C.min() if finite else np.nan,upper_boundary_C=s.T2_C.max() if finite else np.nan,lower_boundary_present=bool(finite),upper_boundary_present=bool(finite),success_count=len(s),finite_window=bool(finite),window_width_C=float(s.T2_C.max()-s.T2_C.min()) if finite else 0,boundary_gap_C=float(s.T2_C.max()-s.T2_C.min()) if finite else np.nan))
 return pd.DataFrame(rows)
def heating(cases):
 schedules=[];ends=[];hist=[];ratios=[]
 for rr in cases.itertuples():
  for peak in (1300,1350,1400,1450,1500,1550):
   for hold in (0,1,2,8,20):
    hh={}
    for rate in (.2,1,5,10,20,50,100):
     rec=3600 if hold==2 else 1e9;f,h=run_path(model(),RampHold(rate,peak,hold*60,start_C=950),initial(rr),300 if rate<=5 else 60,rec,f"H_{rr.case_id}_{rate}_{peak}_{hold}");meta=dict(run_id=f"H_{rr.case_id}_{rate}_{peak}_{hold}",case_id=rr.case_id,initial_state_id=rr.initial_state_id,schedule_id=f"r{rate}_p{peak}_h{hold}",path_type="heating_rate");x=history_frame(h,meta);hist.append(x);hh[rate]=x;schedules.append(dict(**meta,rate_C_min=rate,T_peak_C=peak,hold_h=hold,T_start_C=950));ends.append(dict(**meta,rate_C_min=rate,T_peak_C=peak,hold_h=hold,final_rho=f.rho,final_G_nm=f.G_m*1e9))
    ref=.2 if hh[.2].rho.max()>=.90 else 1.
    for fast in (5,10,20,50,100):
     a=hh[ref].sort_values("rho").drop_duplicates("rho");b=hh[fast].sort_values("rho").drop_duplicates("rho");lo=max(a.rho.min(),b.rho.min());hi=min(a.rho.max(),b.rho.max())
     for rho in np.linspace(lo,hi,50):
      gr=np.interp(rho,a.rho,a.G_nm);gf=np.interp(rho,b.rho,b.G_nm);ratio=gr/max(gf,1e-30);ratios.append(dict(case_id=rr.case_id,T_peak_C=peak,hold_h=hold,reference_rate_C_min=ref,fast_rate_C_min=fast,rho=rho,G_reference_nm=gr,G_fast_nm=gf,G_reference_over_G_fast=ratio,percent_reduction=100*(1-1/ratio),both_attained=hi>=.90))
 return pd.DataFrame(schedules),pd.DataFrame(ends),pd.concat(hist,ignore_index=True),pd.DataFrame(ratios)
def main():
 OUT.mkdir(parents=True,exist_ok=True);(OUT/"source_tables").mkdir(exist_ok=True)
 pd.DataFrame([dict(parameter="barrier",value=hashlib.sha256(BARRIER.read_bytes()).hexdigest(),units="sha256",status="fixed"),dict(parameter="D_GB0",value=.056,units="m2/s",status="fixed"),dict(parameter="Q_GB",value=380000,units="J/mol",status="fixed"),dict(parameter="D_s0",value=.10,units="m2/s",status="fixed"),dict(parameter="Q_s",value=380000,units="J/mol",status="fixed"),dict(parameter="M0_growth",value=.0058,units="m4/J/s",status="current default provisional"),dict(parameter="mechanism_mode",value="defined_closed_laws_port",units="",status="fixed"),dict(parameter="closed_mapping",value="resolved_proxy_current",units="",status="fixed")]).to_csv(OUT/"fixed_zro2_parameter_manifest.csv",index=False)
 prior=pd.DataFrame([dict(prior_candidate_id=693168,status="conditional Tier B; not validation",parameterized_as_ZrO2=False,state_values_copied=False,response_form_target_only=True,prior_rho_switch=.88,prior_G0_nm=103,prior_G1_nm=117,prior_success_band_C="925-1200",prior_mechanism="PR_prepared_closed_accommodation_memory",use_as_ZrO2_parameter=False)]);prior.to_csv(OUT/"prior_solution_manifest.csv",index=False);prior.to_csv(OUT/"prior_solution_reference_cases.csv",index=False)
 g=grid();g.to_csv(OUT/"initial_state_grid.csv",index=False);workers=min(8,os.cpu_count() or 2)
 if (OUT/"first_step_processing_state_screen.csv").exists():scr=pd.read_csv(OUT/"first_step_processing_state_screen.csv")
 else:
  with ProcessPoolExecutor(max_workers=workers) as ex:scr=pd.DataFrame(ex.map(screen_one,g.to_dict("records"),chunksize=8))
  scr.to_csv(OUT/"first_step_processing_state_screen.csv",index=False)
 eligible=scr[scr.promotable].sort_values(["G1_nm","first_step_growth_fraction"]);prom=eligible.drop_duplicates(["rho0","G0_nm","rho_switch","T1_C"]).head(24).copy();prom.insert(0,"case_id",[f"P{i:03d}" for i in range(len(prom))]);prom.to_csv(OUT/"promoted_first_step_states.csv",index=False)
 if (OUT/"twostep_second_step_classification_points.csv").exists():pts=pd.read_csv(OUT/"twostep_second_step_classification_points.csv");states=pd.read_csv(OUT/"twostep_first_step_state_table.csv");reph=pd.read_csv(OUT/"twostep_representative_histories.csv")
 else:pts,states,reph=second_scan(prom);pts.to_csv(OUT/"twostep_second_step_classification_points.csv",index=False);states.to_csv(OUT/"twostep_first_step_state_table.csv",index=False);reph.to_csv(OUT/"twostep_representative_histories.csv",index=False)
 pd.DataFrame([dict(schedule_id=f"T2_{t}_{h}h",T2_C=t,hold_h=h) for t in T2 for h in HOLDS]).to_csv(OUT/"twostep_schedule_definitions.csv",index=False);win=boundaries(pts);win.to_csv(OUT/"twostep_window_boundaries.csv",index=False);win.groupby(["density_target","tolerance_tier"]).agg(finite_window_count=("finite_window","sum"),max_window_width_C=("window_width_C","max"),minimum_G1_finite_nm=("G1_nm",lambda x:x[win.loc[x.index,"finite_window"]].min() if win.loc[x.index,"finite_window"].any() else np.nan)).reset_index().to_csv(OUT/"twostep_boundary_gap_summary.csv",index=False)
 heatcases=(prom.head(2) if len(prom) else scr.nsmallest(2,"first_step_growth_fraction")).copy()
 if (OUT/"heating_rate_histories.csv").exists():hs=pd.read_csv(OUT/"heating_rate_schedule_definitions.csv");he=pd.read_csv(OUT/"heating_rate_endpoint_summary.csv");hh=pd.read_csv(OUT/"heating_rate_histories.csv");hr=pd.read_csv(OUT/"heating_rate_matched_density_ratios.csv")
 else:hs,he,hh,hr=heating(heatcases);hs.to_csv(OUT/"heating_rate_schedule_definitions.csv",index=False);he.to_csv(OUT/"heating_rate_endpoint_summary.csv",index=False);hh.to_csv(OUT/"heating_rate_histories.csv",index=False);hr.to_csv(OUT/"heating_rate_matched_density_ratios.csv",index=False);hh.groupby("run_id").tail(1).to_csv(OUT/"heating_rate_state_diagnostics.csv",index=False)
 # Representative cases and comparator histories.
 finite=win[win.finite_window&win.tolerance_tier.eq("Tier_like_B")];sel=[]
 if len(hr):z=hr.loc[hr.G_reference_over_G_fast.idxmax()];sel.append(dict(case_id=z.case_id,reason_selected="best heating-rate response",T2_C=np.nan,hold_h=2,classification="heating_rate",finite_window_width_C=np.nan,density_target=z.rho,grain_tolerance=np.nan,notes=f"{z.reference_rate_C_min:g} vs {z.fast_rate_C_min:g} C/min"))
 if len(finite):
  z=finite.sort_values("window_width_C",ascending=False).iloc[0];sel.append(dict(case_id=z.case_id,reason_selected="best two-step window",T2_C=(z.lower_boundary_C+z.upper_boundary_C)/2,hold_h=z.hold_h,classification="SUCCESS",finite_window_width_C=z.window_width_C,density_target=z.density_target,grain_tolerance=.20,notes="fixed physics"))
  zn=finite.sort_values("G1_nm").iloc[0];sel.append(dict(case_id=zn.case_id,reason_selected="best nanoscale two-step case",T2_C=(zn.lower_boundary_C+zn.upper_boundary_C)/2,hold_h=zn.hold_h,classification="SUCCESS",finite_window_width_C=zn.window_width_C,density_target=zn.density_target,grain_tolerance=.20,notes="smallest Tier-like-B G1 with finite window"))
  for reason,temp,cl in [("low-T2 density-exhaustion path",zn.lower_boundary_C-25,"DENSIFICATION_EXHAUSTION_FAILURE"),("intermediate T2 success path",(zn.lower_boundary_C+zn.upper_boundary_C)/2,"SUCCESS"),("high-T2 growth-failure path",zn.upper_boundary_C+25,"GRAIN_GROWTH_FAILURE")]:sel.append(dict(case_id=zn.case_id,reason_selected=reason,T2_C=temp,hold_h=zn.hold_h,classification=cl,finite_window_width_C=zn.window_width_C,density_target=zn.density_target,grain_tolerance=.20,notes="boundary-order representative"))
 near=win[(~win.finite_window)&win.success_count.gt(0)].sort_values("G1_nm")
 if len(near):z=near.iloc[0];sel.append(dict(case_id=z.case_id,reason_selected="near-miss nanoscale case",T2_C=np.nan,hold_h=z.hold_h,classification="NEAR_MISS",finite_window_width_C=0,density_target=z.density_target,grain_tolerance=.20,notes="success points lack both finite boundaries"))
 reps=pd.DataFrame(sel);meta=prom[["case_id","rho0","G0_nm","pore_D50_nm","T1_C","rate_T1_C_min","rho_switch","G1_nm"]] if len(prom) else pd.DataFrame();reps=reps.merge(meta,on="case_id",how="left");reps.to_csv(OUT/"selected_representative_cases.csv",index=False)
 comp=[]
 if len(reps) and reps.reason_selected.str.contains("two-step").any():
  z=reps[reps.reason_selected.str.contains("two-step")].iloc[0];rr=prom[prom.case_id.eq(z.case_id)].iloc[0];first,_,_=run_switch(type("R",(),rr.to_dict()),False)
  for typ,temp,state0 in [("high_T_isothermal",rr.T1_C,initial(type("R",(),rr.to_dict()))),("low_T_isothermal",z.T2_C,initial(type("R",(),rr.to_dict()))),("two_step",z.T2_C,first)]:
   f,h=run_path(model(),Iso(temp,z.hold_h),state0,1800,1800,typ);comp.append(history_frame(h,dict(run_id=f"C_{typ}",case_id=z.case_id,initial_state_id=rr.initial_state_id,schedule_id=typ,path_type=typ)))
 comps=pd.concat(comp,ignore_index=True) if comp else pd.DataFrame()
 if not (OUT/"isothermal_comparator_histories.csv").exists():comps.to_csv(OUT/"isothermal_comparator_histories.csv",index=False)
 ii=pd.read_csv(OUT/"isothermal_comparator_histories.csv");mr=[]
 if len(ii) and {"two_step","high_T_isothermal","low_T_isothermal"}<=set(ii.path_type):
  curves={k:g.sort_values("rho").drop_duplicates("rho") for k,g in ii.groupby("path_type")};lo=max(g.rho.min() for g in curves.values());hi=min(g.rho.max() for g in curves.values())
  for rho in np.linspace(lo,hi,50):
   vals={k:np.interp(rho,g.rho,g.G_nm) for k,g in curves.items()};mr.append(dict(case_id=ii.case_id.iloc[0],rho=rho,G_two_step_nm=vals["two_step"],G_high_T_isothermal_nm=vals["high_T_isothermal"],G_low_T_isothermal_nm=vals["low_T_isothermal"],G_high_over_two=vals["high_T_isothermal"]/max(vals["two_step"],1e-30),G_low_over_two=vals["low_T_isothermal"]/max(vals["two_step"],1e-30)))
 pd.DataFrame(mr).to_csv(OUT/"twostep_vs_isothermal_matched_density_ratios.csv",index=False);reph.groupby("run_id").tail(1).to_csv(OUT/"twostep_state_diagnostics.csv",index=False)
 pd.DataFrame([dict(reference="candidate_693168",role="response-form only",zro2_result="independently searched",parameter_values_transferred=False,validation=False)]).to_csv(OUT/"prior_solution_reference_vs_zro2_forward_summary.csv",index=False)
 def span(th):
  spans=[]
  for _,g in hr.groupby(["case_id","T_peak_C","hold_h","reference_rate_C_min","fast_rate_C_min"]):
   q=g[g.G_reference_over_G_fast>=th];spans.append(float(q.rho.max()-q.rho.min()) if len(q) else 0.)
  return max(spans or [0.])
 maxratio=hr.G_reference_over_G_fast.max() if len(hr) else np.nan;nwin=int(win.finite_window.sum());besth=hr.loc[hr.G_reference_over_G_fast.idxmax()] if len(hr) else None;outcome=pd.DataFrame([dict(response="heating_rate",maximum_ratio=maxratio,median_ratio=float(hr.G_reference_over_G_fast.median()),density_span_ratio_ge_1p2=span(1.2),density_span_ratio_ge_1p5=span(1.5),density_span_ratio_ge_2p0=span(2.),best_initial_state=besth.case_id if besth is not None else "none",best_rates=f"{besth.reference_rate_C_min:g}/{besth.fast_rate_C_min:g}" if besth is not None else "none",finite_window_count=np.nan,maximum_window_width_C=np.nan,smallest_G1_nm=np.nan,lower_boundary_cause="not applicable",upper_boundary_cause="not applicable",scale_class="nanoscale",outcome_label="response_form_present" if maxratio>=1.2 else "response_form_absent",not_validated=True),dict(response="two_step",maximum_ratio=np.nan,median_ratio=np.nan,density_span_ratio_ge_1p2=np.nan,density_span_ratio_ge_1p5=np.nan,density_span_ratio_ge_2p0=np.nan,best_initial_state=finite.sort_values("window_width_C",ascending=False).case_id.iloc[0] if len(finite) else "none",best_rates="not applicable",finite_window_count=nwin,maximum_window_width_C=win.window_width_C.max() if len(win) else 0,smallest_G1_nm=finite.G1_nm.min() if len(finite) else np.nan,lower_boundary_cause="density exhaustion",upper_boundary_cause="grain growth",scale_class="nanoscale" if len(finite) and finite.G1_nm.min()<100 else "large-particle only",outcome_label="response_form_present" if nwin else "no_finite_window",not_validated=True)]);outcome.to_csv(OUT/"prediction_outcome_summary.csv",index=False)
 pd.DataFrame([dict(field="A_closed_used/A_closed_recovered",status="available only when closed-law diagnostics emit fields",handling="NaN, not fabricated")]).to_csv(OUT/"missing_history_fields.csv",index=False)
 state=dict(branch="codex/zro2-forward-processing-window-prediction-figures",source_branch="codex/zro2-forward-port-defined-closed-laws",source_commit="2249f02",screen_states=len(scr),promoted_states=len(prom),second_step_points=len(pts),finite_windows=nwin,failed_global_mobility_active=False,candidate_693168_parameter_values_used=False,model_physics_changed=False,validation=False);(OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n");print(state)
if __name__=="__main__":main()
