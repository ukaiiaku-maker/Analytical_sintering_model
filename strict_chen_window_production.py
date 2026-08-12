#!/usr/bin/env python3
"""Strict tier recheck and bounded production confirmation of Chen windows."""
from __future__ import annotations
from pathlib import Path
import csv,math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import dynamic_chen_topology_search as search
import nucleation_fast_chen_production as nuc
import plot_style as ps

SRC=Path("results/dynamic_chen_topology_for_nucleation_materials");OUT=Path("results/strict_chen_window_production")

def write(path,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def tier(prep_growth,max_second,width,G1,practical=True,lower=True,upper=True):
    if not practical or not lower or not upper or width<=0:return "reject"
    if prep_growth<=.05 and max_second<=.05 and width>=25 and G1<=300:return "Tier_A"
    if prep_growth<=.10 and max_second<=.10 and width>=25 and G1<=450:return "Tier_B"
    return "Tier_C"

def lookup():
    reg=pd.read_csv(SRC/"dynamic_topology_parameter_registry.csv");design={x[0]:x for x in search.design(len(reg))};return design

def recheck():
    OUT.mkdir(parents=True,exist_ok=True);c=pd.read_csv(SRC/"confirmed_common_candidates.csv");mats=search.materials();tops=lookup();rows=[];reject=[]
    for _,r in c.iterrows():
      tid,fam,tp=tops[r.topology_id];mat=mats[r.material_id];p,h,att=search.first_state(mat,tp,r.G0_nm,r.T1_C,r.rho_switch)
      if not att:reject.append({**r.to_dict(),"rejection_reason":"first_step_unattainable"});continue
      state=h["final_state"];prep=(state["G"]-r.G0_nm*1e-9)/(r.G0_nm*1e-9);points=[search.point(r.material_id,tid,fam,p,tp,state,r.G0_nm,r.T1_C,r.rho_switch,T2) for T2 in range(800,int(r.T1_C),10)]
      for prep_tol in (.05,.10,.20):
       for sec_tol in (.05,.10):
        s=search.summarize(points,sec_tol,prep_tol);success=[x for x in points if x["prep_growth_fraction"]<=prep_tol and nuc.classify(x["rho1"],x["rho2"],x["growth_fraction"],sec_tol,x["T2_C"],x["T1_C"],True)=="SUCCESS"];max_second=max((x["growth_fraction"] for x in success),default=math.nan);G1=state["G"]*1e9;tr=tier(prep,max_second,s["window_width_C"],G1,True,s["lower_bracketed"],s["upper_bracketed"])
        row=dict(material_id=r.material_id,topology_id=tid,q_TJ=tp.q_TJ,G0_nm=r.G0_nm,G1_nm=G1,T1_C=r.T1_C,rho_switch=r.rho_switch,T2_first_success_C=s["T_first_success"],T2_last_success_C=s["T_last_success"],T_lower_density_C=s["T_first_success"]-10 if s["lower_bracketed"] else math.nan,T_upper_growth_C=s["T_last_success"]+10 if s["upper_bracketed"] else math.nan,window_width_C=s["window_width_C"],first_step_growth_fraction=prep,maximum_second_step_growth_fraction=max_second,prep_growth_tolerance=prep_tol,second_step_growth_tolerance=sec_tol,practical_T2_less_T1=True,tier=tr,rejection_reason="" if tr!="reject" else s["rejection_reason"]);rows.append(row)
        if tr=="reject":reject.append(row)
    write(OUT/"strict_window_recheck.csv",rows);write(OUT/"strict_window_rejected_cases.csv",reject)
    summary=[]
    unique=pd.DataFrame(rows).drop_duplicates(["material_id","topology_id","G0_nm","T1_C","rho_switch","tier"])
    for (mid,tr),g in unique.groupby(["material_id","tier"]):summary.append(dict(material_id=mid,tier=tr,n_unique_windows=len(g),max_width_C=g.window_width_C.max(),minimum_G1_nm=g.G1_nm.min()))
    write(OUT/"strict_window_tier_summary.csv",summary);return rows

def analysis(rows):
    q=[r for r in rows if r["tier"]!="reject"];boundary=[];sens=[];qcmp=[]
    for r in q:
        boundary.append(dict(material_id=r["material_id"],topology_id=r["topology_id"],tier=r["tier"],lower_controller="densification_exhaustion",upper_controller="migration_reactivation_Lambda_over_K_and_drag_release",preparation_growth_limiting=r["tier"]=="Tier_C",practical_truncation=False,T_lower_density_C=r["T_lower_density_C"],T_upper_growth_C=r["T_upper_growth_C"]))
    reg=pd.read_csv(SRC/"dynamic_topology_parameter_registry.csv");rdf=pd.DataFrame(q).merge(reg,on="topology_id") if q else pd.DataFrame()
    if len(rdf):
      for name in ("TJ_drag_strength","pore_drag_strength","XJ_capacity","lambda_ref","K_ref","Q_TJ_event"):
        for val,g in rdf.groupby(name):sens.append(dict(parameter=name,value=val,n_windows=len(g),median_width_C=g.window_width_C.median(),best_tier=sorted(g.tier)[0]))
      qcol="q_TJ_x" if "q_TJ_x" in rdf else "q_TJ"
      for (mid,qv),g in rdf.groupby(["material_id",qcol]):qcmp.append(dict(material_id=mid,q_TJ=qv,n_windows=len(g),max_width_C=g.window_width_C.max(),minimum_G1_nm=g.G1_nm.min(),best_tier=sorted(g.tier)[0]))
    write(OUT/"window_boundary_mechanism_analysis.csv",boundary);write(OUT/"topology_parameter_sensitivity.csv",sens);write(OUT/"q0_q1_dynamic_comparison.csv",qcmp);write(OUT/"confirmed_common_candidates.csv",[r for r in q if r["tier"] in ("Tier_A","Tier_B")]);write(OUT/"conditional_relaxed_candidates.csv",[r for r in q if r["tier"]=="Tier_C"]);return q

def plots(rows):
    ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);inv=[];d=pd.DataFrame(rows)
    fig,ax=plt.subplots(figsize=(7,4.5));colors={"Tier_A":"#009E73","Tier_B":"#0072B2","Tier_C":"#E69F00","reject":"#999999"}
    for tr,col in colors.items():g=d[d.tier==tr];ax.scatter(g.G1_nm,g.window_width_C,c=col,label=tr,alpha=.65,s=18)
    ax.axvline(300,color="#555",ls="--");ax.axvline(450,color="#777",ls=":");ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel=r"Finite $T_2$ window width [$^\circ$C]");ax.legend();ps.clean(ax);ps.finish(fig,fd/"strict_window_tier_map");inv.append(ps.inventory_row("1","strict_window_tier_map","Strict window tiers","strict_window_recheck.csv","Separates publishable and relaxed windows","Results"))
    best=d[d.tier!="reject"].sort_values(["tier","window_width_C"],ascending=[True,False]).head(1)
    if len(best):
      r=best.iloc[0];fig,ax=plt.subplots(figsize=(7,5));g=d[(d.material_id==r.material_id)&(d.topology_id==r.topology_id)&(d.prep_growth_tolerance==r.prep_growth_tolerance)&(d.second_step_growth_tolerance==r.second_step_growth_tolerance)];ax.vlines(g.G1_nm,g.T2_first_success_C,g.T2_last_success_C,color="#009E73",lw=6,alpha=.65);ax.scatter(g.G1_nm,g.T_lower_density_C,c="#0072B2",label="lower density boundary");ax.scatter(g.G1_nm,g.T_upper_growth_C,c="#D55E00",label="upper growth boundary");ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel=r"$T_2$ [$^\circ$C]");ax.legend();ps.clean(ax);ps.finish(fig,fd/"Chen_window_fill_map_best_candidate");inv.append(ps.inventory_row("2","Chen_window_fill_map_best_candidate","Best strict Chen window","strict_window_recheck.csv","Filled finite success interval","Results"))
    ps.write_inventory(OUT/"figure_inventory.csv",inv)

def main():
    rows=recheck();q=analysis(rows);plots(rows);print(pd.DataFrame(rows).groupby(["material_id","tier"]).size().to_string());print("strict common",sum(r["tier"] in ("Tier_A","Tier_B") for r in q),"relaxed",sum(r["tier"]=="Tier_C" for r in q))
if __name__=="__main__":main()
