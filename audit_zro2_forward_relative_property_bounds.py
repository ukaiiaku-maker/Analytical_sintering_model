#!/usr/bin/env python3
"""Audit current ZrO2 properties against observed relative windows; not validation."""
from pathlib import Path
from dataclasses import replace
import hashlib,json
import numpy as np,pandas as pd
from zro2_forward.barrier_json import BarrierModel
from zro2_forward.conditioned_950c import BARRIER
from zro2_forward.integrator import ModelState
from zro2_forward.material_zro2 import MaterialParameters,R,EV_MOL
from zro2_forward.pore_population import PorePopulation
from zro2_forward.resolved_rules import ResolvedRuleModel
from run_zro2_forward_resolved_rules import params,initial

OUT=Path("results/zro2_forward_relative_property_bound_audit");TCS=np.arange(900,1501,100);NA=6.02214076e23;EVJ=1.602176634e-19
SCORE=Path("results/relative_material_property_window_attribution/source_tables/material_property_window_scorecard.csv")
BOUNDS={"Q_nuc_minus_Q_growth_kJ":(-50,96),"Q_nuc_minus_Q_PR_kJ":(51,251),"Q_nuc_minus_Q_transport_kJ":(15,202),"Q_closed_minus_Q_growth_kJ":(-252,-127),"Q_PR_minus_Q_closed_kJ":(-74,212),"log10_kclosed_over_kgrowth":(-1.52,1.50),"log10_kPR_over_kgrowth":(-1.48,1.48),"log10_knuc_over_ktransport":(-1.48,1.52)}
def make_states():
 s0=initial();states={"PDF_conditioned_nominal":s0}
 p=s0.pores.copy();scale=.12/p.total;p.phi_open*=scale;p.phi_iso[:]=0;p.phi_closed[:]=0
 states["T1_switch_like"]=ModelState(0,1150+273.15,.88,125e-9,p,0.,0.,0.)
 p=s0.pores.copy();p.phi_open*=.12/p.total;total=p.phi_open.sum();p.phi_closed=.65*total*(p.phi_open/total);p.phi_open*=.35;p.phi_iso[:]=0
 states["candidate_like_injected"]=ModelState(0,1100+273.15,.88,117e-9,p,.152178631,1.,0.)
 h=pd.read_csv("results/zro2_forward_open_closed_rate_handoff_audit/handoff_diagnostic_histories.csv").query("mode_id=='resolved_default'")
 r=h.loc[h.closed_fraction.idxmax()];p=PorePopulation(np.array(json.loads(r.pore_radii_m_json)),np.array(json.loads(r.phi_open_json)),np.array(json.loads(r.phi_iso_json)),np.array(json.loads(r.phi_closed_json)))
 states["resolved_best_prepared"]=ModelState(r.t_s,r.T_K,r.rho,r.G_m,p,r.A_closed,r.PR_memory,r.cumulative_PR_work)
 defs=[]
 for n,s in states.items():defs.append({"state_id":n,"rho":s.rho,"G_nm":s.G_m*1e9,"phi_open":s.pores.phi_open.sum(),"phi_precursor":s.pores.phi_iso.sum(),"phi_closed":s.pores.phi_closed.sum(),"closed_fraction_of_pores":s.pores.phi_closed.sum()/s.pores.total,"PR_memory":s.PR_memory,"A_closed":s.A_closed,"external_injection":n=="candidate_like_injected"})
 return states,pd.DataFrame(defs)
def effective_q(TK,k):
 x=1/np.asarray(TK);y=np.log(np.maximum(k,1e-300));return -R*np.gradient(y,x)/1000
def inside(v,b):return bool(np.isfinite(v) and b[0]<=v<=b[1])
def main():
 OUT.mkdir(parents=True,exist_ok=True);m=MaterialParameters();bar=BarrierModel.load(BARRIER);model=ResolvedRuleModel(parameters=params())
 wins=[]
 for k,(lo,hi) in BOUNDS.items():wins.append({"window":k,"lower":lo,"upper":hi,"window_type":"exact_promoted_both_pass","source":str(SCORE),"universal_constant":False})
 for k,lo,hi,note in [("OAT_Delta_Q_nuc",0,50,"-25 and +75 fail"),("OAT_Delta_Q_closed",-25,100,"-50 loses lower boundary"),("OAT_Delta_Q_growth",-100,100,"limit not found"),("OAT_PR_prefactor_factor",.3,np.inf,"0.1 fails"),("OAT_growth_prefactor_factor",.1,np.inf,"0.03 loses upper boundary")]:wins.append({"window":k,"lower":lo,"upper":hi,"window_type":"OAT_local","source":"docs/FINAL_MECHANISM_SYNTHESIS_CAPTIONS.md and source scorecard","universal_constant":False,"note":note})
 pd.DataFrame(wins).to_csv(OUT/"source_property_windows.csv",index=False)
 pd.DataFrame([("Q_transport",380,"kJ/mol","D_GB=0.056 exp(-380k/RT)"),("Q_surface",380,"kJ/mol","D_s=0.10 exp(-380k/RT)"),("Q_PR_separate",180,"kJ/mol","resolved Q_PR_J_mol"),("Q_growth",m.Q_M_J_mol/1000,"kJ/mol","4.2 eV high-temperature branch"),("Q_closed_physical",np.nan,"kJ/mol","proxy/unidentified; composite activity uses 0.35 Q_PR"),("D_GB0",m.D_GB0_m2_s,"m2/s","unchanged"),("D_s0",m.D_s0_m2_s,"m2/s","unchanged"),("M0",m.M0_m4_J_s,"m4/(J s)","unchanged")],columns=["property","value","unit","status_or_law"]).to_csv(OUT/"zro2_forward_property_values.csv",index=False)
 states,defs=make_states();defs.to_csv(OUT/"representative_state_definitions.csv",index=False)
 rows=[]
 for sid,s in states.items():
  for TC in TCS:
   TK=TC+273.15;_,_,_,_,_,d=model.rates(s,TK);GJ=bar.Gstar(d["sigma_eff_Pa"],TK);G0=bar.Gstar(0,TK)
   op=max(s.pores.phi_open.sum(),1e-300);cl=s.pores.phi_closed.sum();kg=max(d["G_dot_intrinsic_m_s"]/s.G_m,1e-300);kp=max(d["PR_coarsening_flux"]/op,1e-300);kc=d["rho_dot_closed_sinv"]/cl if cl>0 else np.nan
   rows.append({"state_id":sid,"source":"representative_state","T_C":TC,"rho":s.rho,"G_nm":s.G_m*1e9,"sigma_eff_Pa":d["sigma_eff_Pa"],"Gstar_eV":GJ/EVJ,"Gstar_kJ_mol":GJ*NA/1000,"Gstar_zero_stress_eV":G0/EVJ,"Gstar_zero_stress_kJ_mol":G0*NA/1000,"r_nuc_sinv":d["r_nuc_sinv"],"k_nuc_eff_sinv":d["r_nuc_sinv"],"k_transport_eff_sinv":1/d["tau_transport_s"],"k_PR_eff_sinv":kp,"k_closed_eff_sinv":kc,"k_growth_eff_sinv":kg,"closed_availability":d["closed_availability"],"A_closed":s.A_closed,"PR_memory":s.PR_memory,"barrier_temperature_extrapolated":not bar.temperature_in_fit_range(TK)})
 # Existing resolved histories, sampled without rerunning schedules.
 h=pd.read_csv("results/zro2_forward_open_closed_rate_handoff_audit/handoff_diagnostic_histories.csv").query("mode_id=='resolved_default' and path in ['5C','50C','low_T2','near_best','high_T2']")
 for path,g in h.groupby("path"):
  for _,r in g.iloc[::max(1,len(g)//30)].iterrows():
   GJ=bar.Gstar(r.sigma_eff_Pa,r.T_K);G0=bar.Gstar(0,r.T_K);cl=r.closed_fraction;op=max(r.open_fraction,1e-300)
   rows.append({"state_id":f"history_{path}","source":"existing_history","T_C":r.T_C,"rho":r.rho,"G_nm":r.G_m*1e9,"sigma_eff_Pa":r.sigma_eff_Pa,"Gstar_eV":GJ/EVJ,"Gstar_kJ_mol":GJ*NA/1000,"Gstar_zero_stress_eV":G0/EVJ,"Gstar_zero_stress_kJ_mol":G0*NA/1000,"r_nuc_sinv":r.r_nuc_sinv,"k_nuc_eff_sinv":r.r_nuc_sinv,"k_transport_eff_sinv":1/r.tau_transport_s,"k_PR_eff_sinv":max(r.PR_coarsening_flux/op,1e-300),"k_closed_eff_sinv":r.rho_dot_closed_sinv/cl if cl>0 else np.nan,"k_growth_eff_sinv":max(r.G_dot_intrinsic_m_s/r.G_m,1e-300),"closed_availability":r.closed_availability,"A_closed":r.A_closed,"PR_memory":r.PR_memory,"barrier_temperature_extrapolated":not bar.temperature_in_fit_range(r.T_K)})
 rates=pd.DataFrame(rows)
 # Effective closed Q is a composite proxy from finite-difference local rate.
 rates["Q_closed_eff_kJ_mol"]=np.nan
 for sid,g in rates[rates.source.eq("representative_state")].groupby("state_id"):
  if g.k_closed_eff_sinv.notna().all() and (g.k_closed_eff_sinv>0).all():rates.loc[g.index,"Q_closed_eff_kJ_mol"]=effective_q(g.T_C+273.15,g.k_closed_eff_sinv)
 rates["Q_closed_status"]=np.where(rates.Q_closed_eff_kJ_mol.notna(),"effective finite-difference proxy","proxy/unidentified")
 rates.to_csv(OUT/"effective_rate_values_vs_T.csv",index=False)
 rates[["state_id","source","T_C","rho","G_nm","sigma_eff_Pa","Gstar_eV","Gstar_kJ_mol","Gstar_zero_stress_eV","Gstar_zero_stress_kJ_mol","r_nuc_sinv","barrier_temperature_extrapolated"]].to_csv(OUT/"effective_barrier_values_vs_T.csv",index=False)
 qg=m.Q_M_J_mol/1000;qpr=m.Q_s_J_mol/1000;qpr_separate=params().Q_PR_J_mol/1000;qt=m.Q_GB_J_mol/1000
 q=rates.copy();q["Q_nuc_eff_minus_Q_growth_kJ"]=q.Gstar_kJ_mol-qg;q["Q_nuc_eff_minus_Q_PR_kJ"]=q.Gstar_kJ_mol-qpr;q["Q_nuc_eff_minus_Q_PR_separate_kJ"]=q.Gstar_kJ_mol-qpr_separate;q["Q_nuc_eff_minus_Q_transport_kJ"]=q.Gstar_kJ_mol-qt;q["Q_closed_eff_minus_Q_growth_kJ"]=q.Q_closed_eff_kJ_mol-qg;q["Q_PR_minus_Q_closed_eff_kJ"]=qpr-q.Q_closed_eff_kJ_mol;q["Q_PR_separate_minus_Q_closed_eff_kJ"]=qpr_separate-q.Q_closed_eff_kJ_mol
 mapping=[("inside_both_pass_nuc_growth","Q_nuc_eff_minus_Q_growth_kJ","Q_nuc_minus_Q_growth_kJ"),("inside_both_pass_nuc_PR","Q_nuc_eff_minus_Q_PR_kJ","Q_nuc_minus_Q_PR_kJ"),("inside_both_pass_nuc_transport","Q_nuc_eff_minus_Q_transport_kJ","Q_nuc_minus_Q_transport_kJ"),("inside_both_pass_closed_growth","Q_closed_eff_minus_Q_growth_kJ","Q_closed_minus_Q_growth_kJ"),("inside_both_pass_PR_closed","Q_PR_minus_Q_closed_eff_kJ","Q_PR_minus_Q_closed_kJ")]
 for out,col,b in mapping:q[out]=q[col].map(lambda v:inside(v,BOUNDS[b]))
 q["inside_OAT_Delta_Q_nuc"]=q.Q_nuc_eff_minus_Q_growth_kJ.between(0,50);q["inside_OAT_Delta_Q_closed"]=q.Q_closed_eff_minus_Q_growth_kJ.between(-25,100);q["inside_OAT_Delta_Q_growth"]=True;q["inside_PR_prefactor_threshold"]=True;q["inside_growth_prefactor_threshold"]=True
 q.to_csv(OUT/"relative_barrier_group_audit.csv",index=False)
 p=q[["state_id","source","T_C"]].copy();p["log10_kclosed_over_kgrowth"]=np.log10(np.maximum(q.k_closed_eff_sinv,1e-300)/np.maximum(q.k_growth_eff_sinv,1e-300));p["log10_kPR_over_kgrowth"]=np.log10(np.maximum(q.k_PR_eff_sinv,1e-300)/np.maximum(q.k_growth_eff_sinv,1e-300));p["log10_knuc_over_ktransport"]=np.log10(np.maximum(q.k_nuc_eff_sinv,1e-300)/np.maximum(q.k_transport_eff_sinv,1e-300))
 for c in ["log10_kclosed_over_kgrowth","log10_kPR_over_kgrowth","log10_knuc_over_ktransport"]:p[f"inside_{c}"]=p[c].map(lambda v:inside(v,BOUNDS[c]))
 p["inside_prefactor_envelope"]=p[[c for c in p if c.startswith("inside_")]].all(axis=1);p.to_csv(OUT/"relative_prefactor_group_audit.csv",index=False)
 flags=q[["state_id","source","T_C"]+[x[0] for x in mapping]+["inside_OAT_Delta_Q_nuc","inside_OAT_Delta_Q_closed","inside_OAT_Delta_Q_growth","inside_PR_prefactor_threshold","inside_growth_prefactor_threshold"]].copy();flags["inside_prefactor_envelope"]=p.inside_prefactor_envelope.to_numpy()
 flags["overall_classification"]=np.where(q.Q_closed_status.eq("proxy/unidentified"),"proxy_not_comparable",np.where(flags.filter(like="inside_").all(axis=1),"likely_inside",np.where(q.Q_closed_eff_minus_Q_growth_kJ<BOUNDS["Q_closed_minus_Q_growth_kJ"][0],"likely_outside_low","likely_outside_high")));flags.to_csv(OUT/"statewise_property_window_classification.csv",index=False)
 summ=[]
 for c in flags.columns[3:-1]:summ.append({"property_group":c,"inside_count":int(flags[c].sum()),"outside_count":int((~flags[c]).sum()),"outside_fraction":float((~flags[c]).mean()),"primary_unknown":c in ("inside_both_pass_closed_growth","inside_both_pass_PR_closed")})
 pd.DataFrame(summ).to_csv(OUT/"property_outside_bounds_summary.csv",index=False)
 pd.DataFrame([{"priority":1,"action":"identify/calibrate physical closed-channel activation and availability","reason":"Q_closed is composite proxy and closed/growth plus closed prefactor groups are outside or unidentifiable","change_model_now":False},{"priority":2,"action":"retain stress-resolved nucleation barrier and fixed transport/PR diffusion","reason":"evaluate statewise; do not replace Gstar with scalar","change_model_now":False},{"priority":3,"action":"do not tune GB mobility as rescue","reason":"separate mobility audit produced zero all-gate cases","change_model_now":False}]).to_csv(OUT/"recommended_next_action_from_property_audit.csv",index=False)
 state={"status":"complete_not_validated","barrier_sha256":hashlib.sha256(Path(BARRIER).read_bytes()).hexdigest(),"Q_transport_kJ_mol":qt,"Q_PR_surface_kJ_mol":m.Q_s_J_mol/1000,"Q_PR_separate_kJ_mol":qpr_separate,"Q_growth_kJ_mol":qg,"Q_closed_status":"proxy/unidentified; finite-difference effective values reported where closed inventory exists","effective_Gstar_eV_min":rates.Gstar_eV.min(),"effective_Gstar_eV_max":rates.Gstar_eV.max(),"model_physics_changed":False,"validation_claim":False};(OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n")
 print(json.dumps(state,indent=2));print(pd.DataFrame(summ).to_string(index=False))
if __name__=="__main__":main()
