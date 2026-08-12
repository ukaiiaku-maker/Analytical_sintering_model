#!/usr/bin/env python3
"""Reclassify nucleation candidates and dynamically test migration-only Chen closures."""
from __future__ import annotations
from dataclasses import asdict,replace
from pathlib import Path
import argparse,csv,math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import separated_fast_chen_model as model
import separated_mechanism_production_search as prior
import production_mechanism_assessment as protocols
import plot_style as ps

OLD=Path("results/separated_mechanism_production_search");OUT=Path("results/nucleation_limited_fast_firing_chen_production")

def write(path,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def reclassify(limit=3):
    OUT.mkdir(parents=True,exist_ok=True);near=pd.read_csv(OLD/"fast_firing_near_hit_material_sets.csv");rows=[]
    for _,r in near.iterrows():
        strict=bool(r.meaningful) and bool(r.attained) and not bool(r.promotion_blocked);nuc=not bool(r.no_nucleation_meaningful) if pd.notna(r.no_nucleation_meaningful) else False;pr=bool(r.no_PR_meaningful) if pd.notna(r.no_PR_meaningful) else False
        causal="nucleation_limited_PR_independent" if strict and nuc and pr else ("nucleation_limited_PR_assisted" if strict and nuc else ("stiff_unphysical" if r.promotion_blocked else "censored_or_unattained" if not r.attained else "not_strict"))
        plausible=400e3<=r.Q_disconnection_nucleation<=600e3 and 1.5<=r.max_ratio<=5 and not r.extreme_warning
        rows.append({**r.to_dict(),"causal_type":causal,"physical_plausibility":plausible,"selected":False,"reclassified_rejection_reason":"" if strict and nuc else "nucleation_ablation_or_strict_rule_failed"})
    candidates=sorted((r for r in rows if r["causal_type"].startswith("nucleation_limited") and r["physical_plausibility"]),key=lambda r:(-r["span_ge_1p5"],abs(r["max_ratio"]-2.2)))
    for r in candidates[:limit]:r["selected"]=True
    write(OUT/"nucleation_limited_reclassified_materials.csv",rows);write(OUT/"selected_nucleation_material_sets.csv",[r for r in rows if r["selected"]]);write(OUT/"reclassified_rejection_reasons.csv",[r for r in rows if r["reclassified_rejection_reason"]]);write(OUT/"physical_plausibility_filter.csv",rows);write(OUT/"selected_materials_for_dynamic_Chen.csv",[r for r in rows if r["selected"]]);write(OUT/"stiff_or_extreme_materials.csv",[r for r in rows if r["causal_type"]=="stiff_unphysical"])
    return [r for r in rows if r["selected"]]

def material(row,mode="full_material_model"):
    names=asdict(model.MaterialKinetics()).keys();return model.MaterialKinetics(**{k:row[k] for k in names if k in row and k not in ("ablation_mode","growth_activity_threshold")},ablation_mode=mode)

def timing(selected):
    summary=[];curves=[];hist=[]
    for r in selected:
      mid=r["material_id"]
      for mode in ("full_material_model","no_PR_redistribution","no_nucleation_limitation","transport_only","exchange_limited_variant","no_growth_before_activation"):
        p=material(r,mode);p=replace(p,growth_activity_threshold=1e-5) if mode=="no_growth_before_activation" else p;ref=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(1,r["peak_T_C"],r["hold_h"]));fast=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(r["fast_rate_C_min"],r["peak_T_C"],r["hold_h"]));m=prior.metrics(ref,fast);c=m.pop("curve",None);summary.append(dict(material_id=mid,ablation_mode=mode,**m))
        if c is not None:
          for x in c.to_dict("records"):curves.append(dict(material_id=mid,ablation_mode=mode,**x))
        for path,h in (("reference",ref),("fast",fast)):
          before=0.;during=0.;lastG=h["G"][0] if len(h["G"]) else p.G0
          for i in range(len(h["rho"])):
            dg=max(h["G"][i]-lastG,0);before+=dg if h["activity"][i]<1e-5 else 0;during+=dg if h["activity"][i]>=1e-5 else 0;lastG=h["G"][i]
            if i%max(1,len(h["rho"])//150)==0:hist.append(dict(material_id=mid,ablation_mode=mode,path=path,t_h=h["t"][i]/3600,T_C=h["T_C"][i],rho=h["rho"][i],G_nm=h["G"][i]*1e9,tau_nuc=h["tau_nuc"][i],tau_exchange=h["tau_exchange"][i],tau_transport=h["tau_transport"][i],activity=h["activity"][i],rho_dot=h["rho_dot"][i],G_dot=h["G_dot"][i],cumulative_growth_before_activation_nm=before*1e9,cumulative_growth_during_densification_nm=during*1e9,PR_exposure=h["PR_exposure"][i],pore_D50_nm=h["pore_D50"][i]*1e9,pore_D90_nm=h["pore_D90"][i]*1e9,connected_fine=h["connected_fine"][i]))
    write(OUT/"fast_firing_timing_audit.csv",summary);write(OUT/"fast_firing_ablation_ratio_curves.csv",curves);write(OUT/"fast_firing_timing_histories.csv",hist);return summary

def topology_registry():
    families=("topology_disabled","pore_drag_only","TJ_drag_only","persistent_junction_only","structural_TJ_multihit_q0","structural_TJ_multihit_q1","persistent_q0_multihit","persistent_q1_multihit","pore_pinned_structural","pore_relaxed_structural","mixed_relaxed_pinned")
    out=[]
    for i,f in enumerate(families):
        q=1 if "q1" in f else 0;disabled=f=="topology_disabled";tp=model.TopologyGrowthClosure(mode="disabled" if disabled else f,TJ_drag_strength=0 if disabled or f=="pore_drag_only" else (8 if "TJ_drag" in f else 16),pore_drag_strength=10 if "pore" in f and not disabled else 0,XJ_capacity=.5 if "persistent" in f else .25,lambda_ref=2,K_ref=2,q_TJ=q,pore_relax_fraction=1 if "relaxed" in f else (.5 if "mixed" in f else 0),pore_drag_fraction=1 if "pinned" in f else (.5 if "mixed" in f else 0))
        out.append((f"C{i:02d}",f,tp))
    return out

class RampToSwitch:
    def __init__(self,rate,T1):self.rate=rate/60;self.T1=T1;self.ramp=(T1-25)/self.rate;self.t_end=self.ramp+96*3600
    def T(self,t,rho):return min(self.T1,25+self.rate*t)

def classify(rho1,rho2,growth,tol,T2,T1,attained):
    if not attained:return "UNATTAINABLE_FIRST_STEP"
    if rho1>=.90-1e-12:return "INELIGIBLE_TARGET_ALREADY_REACHED"
    dense=rho2>=.90-1e-12;small=growth<=tol+1e-12
    if dense and small and T2<T1:return "SUCCESS"
    if dense and not small:return "GRAIN_GROWTH_FAILURE"
    if not dense and small:return "DENSIFICATION_EXHAUSTION_FAILURE"
    return "MIXED_FAILURE"

def dynamic_chen(selected):
    points=[];bounds=[];registry=[];failed=[]
    # Bounded exact reintegration: full topology families, representative G/T1/rate/switch axes.
    for tid,fam,tp in topology_registry():registry.append(dict(topology_id=tid,family=fam,**asdict(tp)))
    for r in selected:
      p=material(r)
      for tid,fam,tp in topology_registry():
       for G0 in (75,225,600):
        pm=replace(p,G0=G0*1e-9)
        for T1 in (1300,1500):
         for rate in (20,):
          for sw in (.80,.88):
            h1=model.run(pm,tp,RampToSwitch(rate,T1),stop_at_rho=sw);att=h1["rho_final"]>=sw-1e-10
            if not att:
                failed.append(dict(material_id=r["material_id"],topology_id=tid,G0_nm=G0,T1_C=T1,heating_rate_T1=rate,rho_switch=sw,rejection_reason="unattainable_first_step"));continue
            state=h1["final_state"];rho1=state["rho"];G1=state["G"]
            for T2 in range(800,1501,50):
                h2=model.run(pm,tp,protocols.FastSchedule(1,T2,96),initial=state);rho2=h2["rho_final"];growth=(h2["G_final"]-G1)/G1
                row=dict(material_id=r["material_id"],topology_id=tid,topology_family=fam,q_TJ=tp.q_TJ,G0_nm=G0,G1_nm=G1*1e9,T1_C=T1,heating_rate_T1=rate,rho_switch=sw,T2_C=T2,rho1=rho1,rho2=rho2,growth_fraction=growth,first_step_attained=att,numerical_censored=h2["numerical_censored"])
                for tol in (.05,.10):points.append({**row,"growth_tolerance":tol,"classification":classify(rho1,rho2,growth,tol,T2,T1,att)})
    df=pd.DataFrame(points)
    if len(df):
      keys=["material_id","topology_id","G0_nm","T1_C","heating_rate_T1","rho_switch","growth_tolerance"]
      for key,g in df.groupby(keys):
        practical=g[g.T2_C<g.T1_C];success=practical[practical.classification=="SUCCESS"];low=bool(len(success) and np.any((practical.T2_C<success.T2_C.min())&(practical.classification=="DENSIFICATION_EXHAUSTION_FAILURE")));up=bool(len(success) and np.any((practical.T2_C>success.T2_C.max())&(practical.classification=="GRAIN_GROWTH_FAILURE")));status="COMPLETE_WINDOW" if len(success) and low and up else "NO_COMPLETE_WINDOW"
        width=success.T2_C.max()-success.T2_C.min() if len(success) else math.nan
        if status=="COMPLETE_WINDOW" and width<=0:status="ISOLATED_SUCCESS_NOT_WINDOW"
        row=dict(zip(keys,key),boundary_status=status,lower_bracketed=low,upper_bracketed=up,T_first_success_C=success.T2_C.min() if len(success) else math.nan,T_last_success_C=success.T2_C.max() if len(success) else math.nan,window_width_C=width);bounds.append(row)
        if status!="COMPLETE_WINDOW":failed.append({**row,"rejection_reason":"missing_complete_practical_boundaries"})
    write(OUT/"dynamic_Chen_classification_points.csv",points);write(OUT/"dynamic_Chen_window_boundaries.csv",bounds);write(OUT/"dynamic_Chen_success_intervals.csv",[r for r in bounds if r["boundary_status"]=="COMPLETE_WINDOW"]);write(OUT/"dynamic_Chen_failed_or_censored_cases.csv",failed);write(OUT/"topology_parameter_registry.csv",registry)
    summary=[]
    for tid,fam,tp in topology_registry():summary.append(dict(topology_id=tid,family=fam,q_TJ=tp.q_TJ,complete_windows=sum(r["topology_id"]==tid and r["boundary_status"]=="COMPLETE_WINDOW" for r in bounds)))
    write(OUT/"topology_parameter_summary.csv",summary);return points,bounds,summary

def plots(selected,points,bounds):
    ps.apply_style();fd=OUT/"production_figures";fd.mkdir(exist_ok=True);inventory=[]
    if points:
      p=pd.DataFrame(points);b=pd.DataFrame(bounds);complete=b[b.boundary_status=="COMPLETE_WINDOW"]
      fig,ax=plt.subplots(figsize=(7,5));colors={"SUCCESS":"#009E73","DENSIFICATION_EXHAUSTION_FAILURE":"#0072B2","GRAIN_GROWTH_FAILURE":"#D55E00","MIXED_FAILURE":"#CC79A7"};q=p[(p.growth_tolerance==.10)&(p.topology_id==p.topology_id.iloc[0])]
      for c,col in colors.items():z=q[q.classification==c];ax.scatter(z.G1_nm,z.T2_C,c=col,s=8,label=c,alpha=.5)
      ax.set(xscale="log",xlabel=r"First-step grain size, $G_1$ [nm]",ylabel=r"Second-step temperature, $T_2$ [$^\circ$C]");ax.legend(fontsize=6);ps.clean(ax);ps.finish(fig,fd/"Figure6_dynamic_Chen_classification");inventory.append(ps.inventory_row("6","Figure6_dynamic_Chen_classification","Dynamic Chen classification","dynamic_Chen_classification_points.csv","Full exact-state classifications","Results"))
      if len(complete):
        fig,ax=plt.subplots(figsize=(7,5));q=complete[complete.growth_tolerance==.10]
        for tid,g in q.groupby("topology_id"):ax.vlines(g.G0_nm,g.T_first_success_C,g.T_last_success_C,alpha=.5,label=tid)
        ax.set(xscale="log",xlabel=r"Initial grain size, $G_0$ [nm]",ylabel=r"$T_2$ success interval [$^\circ$C]");ax.legend(fontsize=6,ncol=2);ps.clean(ax);ps.finish(fig,fd/"Figure7_Chen_filled_windows");inventory.append(ps.inventory_row("7","Figure7_Chen_filled_windows","Filled practical Chen windows","dynamic_Chen_window_boundaries.csv","Shows finite lower/upper intervals","Results"))
    ps.write_inventory(OUT/"production_figure_inventory.csv",inventory)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--selected",type=int,default=3);a=ap.parse_args();sel=reclassify(a.selected);timing(sel);points,bounds,top=dynamic_chen(sel);plots(sel,points,bounds);complete=[r for r in bounds if r["boundary_status"]=="COMPLETE_WINDOW"]
    overlap=[]
    for r in sel:
      tids=sorted({x["topology_id"] for x in complete if x["material_id"]==r["material_id"]})
      if tids:
       for tid in tids:overlap.append(dict(material_id=r["material_id"],topology_id=tid,classification="both_confirmed",material_parameters_frozen=True))
      else:overlap.append(dict(material_id=r["material_id"],topology_id="",classification="fast_only",material_parameters_frozen=True))
    write(OUT/"confirmed_overlap_scorecard.csv",overlap);write(OUT/"common_material_topology_sets.csv",[r for r in overlap if r["classification"]=="both_confirmed"]);write(OUT/"fast_only_cases.csv",[r for r in overlap if r["classification"]=="fast_only"]);write(OUT/"Chen_only_cases.csv",[]);write(OUT/"both_response_cases.csv",[r for r in overlap if r["classification"]=="both_confirmed"]);write(OUT/"conflict_cases.csv",[]);print("selected",[r["material_id"] for r in sel],"complete windows",len(complete),"overlap",sum(r["classification"]=="both_confirmed" for r in overlap))
if __name__=="__main__":main()
