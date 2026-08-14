#!/usr/bin/env python3
"""Adaptive migration-only topology search for frozen E0021/E0142 kinetics."""
from __future__ import annotations
from dataclasses import asdict,replace
from pathlib import Path
import argparse,csv,math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import separated_fast_chen_model as model
import nucleation_fast_chen_production as base
import production_mechanism_assessment as protocols
import plot_style as ps

OUT=Path("results/dynamic_chen_topology_for_nucleation_materials")

def write(path,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def materials():
    d=pd.read_csv("results/nucleation_limited_fast_firing_chen_production/selected_nucleation_material_sets.csv");return {r.material_id:base.material(r.to_dict()) for _,r in d.iterrows()}

def design(n=512):
    families=("persistent_junction_only","structural_q0","structural_q1","persistent_q0","persistent_q1","pore_pinned_q0","pore_pinned_q1","pore_relaxed_q0","pore_relaxed_q1","mixed_q0","mixed_q1");out=[]
    levels=((2,6,12,24,48,80,120),(1e4,5e4,1.5e5,4.5e5,1.35e6,3e6,1e7),(.25,.5,.75,1,1.5),(.25,.5,1,2,4,8),(0,.25,.5,1,2),(.1,.25,.5,1.5,4.5,13.5,30),(1,2,4,8,12),(280e3,320e3,340e3,360e3,380e3,420e3,460e3,500e3),(0,3,10,30,80),(0,.25,.5,.75,1),(0,.25,.5,.75,1));pr=(1,3,5,7,11,13,17,19,23,29,31)
    for i in range(n):
        v=[x[(i*p+i//len(x))%len(x)] for x,p in zip(levels,pr)];fam=families[i%len(families)];q=1 if "q1" in fam else 0;pinned="pinned" in fam or "mixed" in fam;relaxed="relaxed" in fam or "mixed" in fam;persistent="persistent" in fam
        tp=model.TopologyGrowthClosure(mode=fam,TJ_drag_strength=v[0] if persistent else 0,pore_drag_strength=v[8] if pinned else 0,XJ_capacity=v[2],XJ_relaxation=1/max(v[1],1),lambda_ref=v[5],K_ref=v[6],q_TJ=q,pore_relax_fraction=v[9] if relaxed else 0,pore_drag_fraction=v[10] if pinned else 0,XJ_prod_TJ=v[3],XJ_prod_sweep=v[4],Q_TJ_event=v[7]);out.append((f"D{i:04d}",fam,tp))
    return out

def first_state(mat,tp,G,T1,sw):
    p=replace(mat,G0=G*1e-9);h=base.model.run(p,tp,base.RampToSwitch(20,T1),stop_at_rho=sw);return p,h,h["rho_final"]>=sw-1e-10

def point(mid,tid,fam,p,tp,state,G,T1,sw,T2):
    h=model.run(p,tp,protocols.FastSchedule(1,T2,96),initial=state);rho1=state["rho"];G1=state["G"];rho2=h["rho_final"];growth=(h["G_final"]-G1)/G1;factor,td=model.topology_growth_factor({"G":G1,"X_J":state["X_J"],"connected_coverage":model.material_rates(rho1,G1,state["phi"],state["radii"],T2,p)["connected_fine"]},T2,tp)
    return dict(material_id=mid,topology_id=tid,family=fam,q_TJ=tp.q_TJ,G0_nm=G,G1_nm=G1*1e9,prep_growth_fraction=(G1-G*1e-9)/(G*1e-9),T1_C=T1,rho_switch=sw,T2_C=T2,rho1=rho1,rho2=rho2,growth_fraction=growth,X_J=h["final_state"]["X_J"],Lambda_over_K=h["Lambda_over_K"][-1],migration_factor=factor,numerical_censored=h["numerical_censored"])

def summarize(points,tol=.10,prep_tol=.20):
    q=sorted(points,key=lambda x:x["T2_C"])
    for r in q:r["classification"]="PREP_GROWTH_EXCEEDED" if r["prep_growth_fraction"]>prep_tol else base.classify(r["rho1"],r["rho2"],r["growth_fraction"],tol,r["T2_C"],r["T1_C"],True)
    success=[r for r in q if r["classification"]=="SUCCESS"];low=bool(success and any(r["T2_C"]<min(x["T2_C"] for x in success) and r["classification"]=="DENSIFICATION_EXHAUSTION_FAILURE" for r in q));up=bool(success and any(r["T2_C"]>max(x["T2_C"] for x in success) and r["classification"]=="GRAIN_GROWTH_FAILURE" for r in q));width=max((r["T2_C"] for r in success),default=math.nan)-min((r["T2_C"] for r in success),default=math.nan) if success else math.nan
    if q and all(r["prep_growth_fraction"]>prep_tol for r in q):reason="prep_growth_exceeded"
    elif not success:reason="no_success"
    elif len(success)==1:reason="isolated_success_not_window"
    elif not low:reason="missing_lower_boundary"
    elif not up:reason="missing_upper_boundary"
    elif width<=0:reason="isolated_success_not_window"
    else:reason=""
    return dict(n_success=len(success),lower_bracketed=low,upper_bracketed=up,T_first_success=min((r["T2_C"] for r in success),default=math.nan),T_last_success=max((r["T2_C"] for r in success),default=math.nan),window_width_C=width,outcome="finite_window" if not reason else "near_hit" if success and low and up else "rejected",rejection_reason=reason)

def c08_refinement(mat):
    tp=base.topology_registry()[8][2];rows=[]
    for T1 in (1300,1500):
      for sw in (.80,.88):
        p,h,att=first_state(mat,tp,600,T1,sw)
        if att:
          state=h["final_state"]
          for T2 in range(950,1251,10):rows.append(point("E0142","C08","pore_pinned_structural",p,tp,state,600,T1,sw,T2))
    for r in rows:r["classification_5pct"]=base.classify(r["rho1"],r["rho2"],r["growth_fraction"],.05,r["T2_C"],r["T1_C"],True);r["classification_10pct"]=base.classify(r["rho1"],r["rho2"],r["growth_fraction"],.10,r["T2_C"],r["T1_C"],True)
    write(OUT/"E0142_C08_near_hit_refinement.csv",rows);return rows

def screen(mats,n):
    registry=[];summary=[];near=[];reject=[];allpoints=[]
    states={}
    for mid,mat in mats.items():
      for G in (150,300,600):
       for T1 in (1300,1400,1500):
        for sw in (.80,.88):
         # First-step state depends on topology, so cache per closure below.
         states[(mid,G,T1,sw)]=(mat,)
    for j,(tid,fam,tp) in enumerate(design(n),1):
      registry.append(dict(topology_id=tid,family=fam,**asdict(tp)));best=None
      for mid,mat in mats.items():
       points=[];first_fail=0
       for G,T1,sw in ((300,1300,.80),(600,1300,.80),(300,1500,.88),(600,1500,.88)):
        p,h,att=first_state(mat,tp,G,T1,sw)
        if not att:first_fail+=1;continue
        state=h["final_state"]
        for T2 in range(900,min(1451,T1),50):points.append(point(mid,tid,fam,p,tp,state,G,T1,sw,T2))
       s=summarize(points);row=dict(material_id=mid,topology_id=tid,family=fam,q_TJ=tp.q_TJ,first_step_failures=first_fail,**s);summary.append(row)
       if s["outcome"]=="finite_window" or s["rejection_reason"]=="isolated_success_not_window":near.append(row)
       else:reject.append(row)
       if best is None or s["n_success"]>best[0]:best=(s["n_success"],points)
      if best and (best[0]>0 or j<=12):allpoints.extend(best[1])
      if j%25==0 or j==n:write(OUT/"dynamic_topology_parameter_registry.csv",registry);write(OUT/"dynamic_topology_screen_summary.csv",summary);write(OUT/"dynamic_topology_near_hits.csv",near);write(OUT/"dynamic_topology_rejected_cases.csv",reject)
    write(OUT/"dynamic_topology_classification_points.csv",allpoints);return summary,near,allpoints

def refine(mats,near):
    lookup={x[0]:x for x in design(max(int(r["topology_id"][1:]) for r in near)+1)} if near else {};rows=[];bounds=[]
    for r in sorted(near,key=lambda x:x["n_success"],reverse=True)[:12]:
      tid,fam,tp=lookup[r["topology_id"]];mat=mats[r["material_id"]]
      for G,T1,sw in ((300,1300,.80),(600,1300,.80),(300,1500,.88),(600,1500,.88)):
        p,h,att=first_state(mat,tp,G,T1,sw)
        if not att:continue
        q=[point(r["material_id"],tid,fam,p,tp,h["final_state"],G,T1,sw,T2) for T2 in range(900,min(1451,T1),10)];s=summarize(q);rows.extend(q);bounds.append(dict(material_id=r["material_id"],topology_id=tid,G0_nm=G,T1_C=T1,rho_switch=sw,prep_growth_tolerance=.20,second_step_growth_tolerance=.10,**s))
    write(OUT/"dynamic_topology_refined_boundaries.csv",bounds);write(OUT/"dynamic_topology_success_intervals.csv",[r for r in bounds if r["outcome"]=="finite_window"]);return rows,bounds

def plots(c08,summary,refined,bounds):
    ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);inv=[];q=pd.DataFrame(c08);fig,axs=plt.subplots(2,2,figsize=(9,7))
    for _,g in q.groupby(["T1_C","rho_switch"]):axs[0,0].plot(g.T2_C,g.rho2);axs[0,1].plot(g.T2_C,g.growth_fraction);axs[1,0].plot(g.T2_C,g.X_J);axs[1,1].plot(g.T2_C,g.Lambda_over_K)
    axs[0,0].axhline(.90,color="k",ls="--");axs[0,1].axhline(.05,color="k",ls="--");axs[0,0].set(ylabel=r"$\rho_2$");axs[0,1].set(ylabel="Growth fraction");axs[1,0].set(ylabel=r"$X_J$");axs[1,1].set(ylabel=r"$\Lambda/K$");[a.set(xlabel=r"$T_2$ [$^\circ$C]") for a in axs.flat];ps.panel_labels(axs);ps.finish(fig,fd/"E0142_C08_refined_T2_scan");inv.append(ps.inventory_row("1","E0142_C08_refined_T2_scan","E0142+C08 refinement","E0142_C08_near_hit_refinement.csv","Tests isolated point at 10 C","Results"))
    s=pd.DataFrame(summary);fig,ax=plt.subplots(figsize=(8,4));codes={"rejected":0,"near_hit":1,"finite_window":2};ax.scatter(range(len(s)),[codes[x] for x in s.outcome],c=[codes[x] for x in s.outcome],cmap="viridis",s=8);ax.set(xlabel="Topology-material evaluation",ylabel="Outcome code",yticks=(0,1,2),yticklabels=("rejected","near hit","finite window"));ps.clean(ax);ps.finish(fig,fd/"topology_search_outcome_map");inv.append(ps.inventory_row("2","topology_search_outcome_map","Topology search outcomes","dynamic_topology_screen_summary.csv","Shows broad search classifications","Results"));ps.write_inventory(OUT/"figure_inventory.csv",inv)

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--topologies",type=int,default=512);a=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True);mats=materials();c08=c08_refinement(mats["E0142"]);summary,near,points=screen(mats,a.topologies);refined,bounds=refine(mats,near);plots(c08,summary,refined,bounds);write(OUT/"E0021_topology_search_summary.csv",[r for r in summary if r["material_id"]=="E0021"]);write(OUT/"E0142_topology_search_summary.csv",[r for r in summary if r["material_id"]=="E0142"]);confirmed=[r for r in bounds if r["outcome"]=="finite_window"];write(OUT/"confirmed_common_candidates.csv",confirmed);write(OUT/"overlap_scorecard.csv",[{"material_id":mid,"classification":"both_confirmed" if any(r["material_id"]==mid for r in confirmed) else "fast_only"} for mid in mats]);print("screen",len(summary),"near",len(near),"refined finite",len(confirmed))
if __name__=="__main__":main()
