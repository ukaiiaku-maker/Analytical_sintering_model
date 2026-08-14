#!/usr/bin/env python3
"""Bounded pore-occupied versus structurally constrained TJ ablation."""
from __future__ import annotations
import csv,math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import adaptive_T2_boundary_search as adaptive
import agentic_mechanism_model as discovery
import joint_pr_desintering_search as joint
import pr_desintering_memory_model as memory
import preparation_window_search as preparation
import production_mechanism_assessment as production

MODES=discovery.TJ_CONSTRAINT_MODES


def write(path,rows):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def params(mid,mode):
    p=joint.frozen_base()[mid];p=replace(p,TJ_constraint_mode=mode);return memory.PRMemoryParams(p,"PR_plus_connected_fine_attrition",k_PR_ref_s=2e-4)


def task(item):
    mid,mode=item;p0=params(mid,mode);fast=[];history=[]
    for topo,frac in production.TOPOLOGIES.items():
      p=replace(p0,base=production.fast_params(p0.base,75,.70,frac))
      for peak in (1350.,1400.,1450.):
       for hold in (8.,20.):
        paths={r:memory.run(p,production.FastSchedule(r,peak,hold)) for r in (1.,20.,100.)}
        if topo=="baseline" and peak==1400 and hold==20:
            h=paths[20.];stride=max(1,len(h["rho"])//250)
            for i in list(range(0,len(h["rho"]),stride))+[len(h["rho"])-1]:history.append(dict(base_mechanism=mid,TJ_constraint_mode=mode,q_TJ=p0.base.q_TJ,t_s=float(h["t"][i]),T_C=float(h["T_C"][i]),rho=float(h["rho"][i]),G_nm=float(h["G"][i])*1e9,C_TJ_total=float(h["C_TJ_total"][i]),C_TJ_pore=float(h["C_TJ_pore"][i]),C_TJ_structural=float(h["C_TJ_structural"][i]),C_TJ_constraint=float(h["C_TJ_constraint"][i]),C_TJ_relaxed=float(h["C_TJ_relaxed"][i]),C_TJ_pinned=float(h["C_TJ_pinned"][i]),R_TJ_pore_drag=float(h["R_TJ_pore_drag"][i]),Lambda_TJ_structural=float(h["Lambda_TJ_structural"][i]),K_TJ_structural=float(h["K_TJ_structural"][i]),Lambda_over_K_TJ=float(h["Lambda_over_K_TJ"][i]),P_comp_TJ_structural=float(h["P_comp_TJ_structural"][i]),P_TJ_multihit=float(h["P_TJ_multihit"][i]),P_TJ_pore_drag=float(h["P_TJ_pore_drag"][i]),P_TJ_assisted_densification=float(h["P_TJ_assisted_densification"][i])))
        for target in (.85,.90):
         ids={r:np.flatnonzero(h["rho"]>=target-1e-12) for r,h in paths.items()}
         for rate in (20.,100.):
          if not len(ids[1.]) or not len(ids[rate]):fast.append(dict(initial_topology=topo,peak_T_C=peak,hold_h=hold,rho_target=target,rate=rate,attained=False));continue
          i=int(ids[rate][0]);j=int(ids[1.][0]);hr=100*(paths[1.]["G"][j]-paths[rate]["G"][i])/paths[1.]["G"][j]
          fast.append(dict(initial_topology=topo,peak_T_C=peak,hold_h=hold,rho_target=target,rate=rate,attained=True,HR_pct=hr,beneficial=hr>1,C_TJ_pore=float(paths[rate]["C_TJ_pore"][i]),C_TJ_constraint=float(paths[rate]["C_TJ_constraint"][i]),C_TJ_relaxed=float(paths[rate]["C_TJ_relaxed"][i]),C_TJ_pinned=float(paths[rate]["C_TJ_pinned"][i]),R_TJ_pore_drag=float(paths[rate]["R_TJ_pore_drag"][i]),Lambda_TJ=float(paths[rate]["Lambda_TJ"][i]),K_TJ=float(paths[rate]["K_TJ"][i]),Lambda_over_K_TJ=float(paths[rate]["Lambda_over_K_TJ"][i]),P_comp_TJ=float(paths[rate]["P_comp_TJ"][i]),P_TJ_multihit=float(paths[rate]["P_TJ_multihit"][i]),P_TJ_pore_drag=float(paths[rate]["P_TJ_pore_drag"][i]),P_TJ_assisted_densification=float(paths[rate]["P_TJ_assisted_densification"][i])))
    b=replace(p0.base.action.location.base,G0=225e-9);pc=replace(p0,base=replace(p0.base,action=replace(p0.base.action,location=replace(p0.base.action.location,base=b))));h1=memory.run(pc,preparation.FixedBudgetRamp(20,1400),stop_at_rho=.84);chen_status="UNATTAINABLE_FIRST_STEP";complete=False
    if h1["rho"][-1]>=.84-1e-12:
        st=memory.final_state(h1,pc);pts=joint.adaptive_points(f"{mid}_{mode}",pc,225,1400,.84,st);s=adaptive.status(pts,.10,1400,True);chen_status=s["boundary_status"];complete=chen_status=="COMPLETE_WINDOW" and (st.base.pore.G*1e9-225)/225<=.10
    valid=[r for r in fast if r.get("attained")];benef=[r for r in valid if r.get("beneficial")];return dict(base_mechanism=mid,TJ_constraint_mode=mode,q_TJ=p0.base.q_TJ,chen_boundary_status=chen_status,complete_practical_chen=complete,n_fast_attained=len(valid),n_fast_beneficial=len(benef),n_fast_harmful=sum(r["HR_pct"]<-1 for r in valid),n_fast_neutral=sum(abs(r["HR_pct"])<=1 for r in valid),HR_pct_median=np.median([r["HR_pct"] for r in benef]) if benef else math.nan,C_TJ_pore_median=np.median([r["C_TJ_pore"] for r in benef]) if benef else math.nan,C_TJ_constraint_median=np.median([r["C_TJ_constraint"] for r in benef]) if benef else math.nan,P_TJ_multihit_median=np.median([r["P_TJ_multihit"] for r in benef]) if benef else math.nan,P_TJ_pore_drag_median=np.median([r["P_TJ_pore_drag"] for r in benef]) if benef else math.nan,joint_positive=bool(complete and benef and len(benef)<len(valid))),fast,history


def plots(out,summary,history):
    fd=out/"tj_constraint_figures";fd.mkdir(exist_ok=True);modes=list(MODES);joint=[sum(r["joint_positive"] for r in summary if r["TJ_constraint_mode"]==m) for m in modes];chen=[sum(r["complete_practical_chen"] for r in summary if r["TJ_constraint_mode"]==m) for m in modes];benef=[sum(r["n_fast_beneficial"] for r in summary if r["TJ_constraint_mode"]==m) for m in modes]
    for vals,name,ylabel in ((joint,"figure1_joint_scorecard.png","joint-positive bases"),(chen,"figure2_chen_window_count.png","complete canonical Chen bases"),(benef,"figure3_beneficial_fast_count.png","beneficial fast cases")):
        fig,ax=plt.subplots(figsize=(9,4));ax.bar(modes,vals);ax.set_ylabel(ylabel);ax.tick_params(axis="x",rotation=25);fig.tight_layout();fig.savefig(fd/name,dpi=180);plt.close(fig)
    sample=[r for r in history if r["base_mechanism"]=="mech_009"]
    fig,ax=plt.subplots();
    for mode in modes:
        q=[r for r in sample if r["TJ_constraint_mode"]==mode];ax.plot([r["rho"] for r in q],[r["C_TJ_pore"] for r in q],ls="--",label=f"{mode}: pore");ax.plot([r["rho"] for r in q],[r["C_TJ_constraint"] for r in q],label=f"{mode}: constraint")
    ax.set(xlabel="density",ylabel="TJ coverage");ax.legend(fontsize=5,ncol=2);fig.tight_layout();fig.savefig(fd/"figure4_pore_vs_constraint_histories.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots();
    for mode in modes:
        q=[r for r in sample if r["TJ_constraint_mode"]==mode];ax.plot([r["rho"] for r in q],[r["P_TJ_multihit"] for r in q],label=f"{mode}: multihit");ax.plot([r["rho"] for r in q],[r["P_TJ_pore_drag"] for r in q],ls="--",label=f"{mode}: pore drag")
    ax.set(xlabel="density",ylabel="migration dissipation");ax.legend(fontsize=5,ncol=2);fig.tight_layout();fig.savefig(fd/"figure5_multihit_vs_pore_drag.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4));x=np.arange(len(modes));q0=[sum(r["n_fast_beneficial"] for r in summary if r["TJ_constraint_mode"]==m and r["q_TJ"]==0) for m in modes];q1=[sum(r["n_fast_beneficial"] for r in summary if r["TJ_constraint_mode"]==m and r["q_TJ"]==1) for m in modes];ax.bar(x-.2,q0,.4,label="q0");ax.bar(x+.2,q1,.4,label="q1");ax.set_xticks(x,modes,rotation=25);ax.set_ylabel("beneficial fast cases");ax.legend();fig.tight_layout();fig.savefig(fd/"figure6_q0_vs_q1.png",dpi=180);plt.close(fig)


def main():
    items=[(mid,mode) for mid in joint.frozen_base() for mode in MODES];summary=[];points=[];history=[]
    with ProcessPoolExecutor(max_workers=8) as pool:
        for s,p,h in pool.map(task,items,chunksize=1):summary.append(s);points.extend([{**r,"base_mechanism":s["base_mechanism"],"TJ_constraint_mode":s["TJ_constraint_mode"]} for r in p]);history.extend(h)
    out=Path("results/production_pr_desintering_assessment");write(out/"TJ_constraint_ablation.csv",summary);write(out/"TJ_constraint_fast_points.csv",points);write(out/"tj_constraint_mode_summary.csv",summary);write(out/"tj_constraint_mode_fast_firing.csv",points);write(out/"tj_constraint_mode_chen_windows.csv",[{k:r[k] for k in ("base_mechanism","TJ_constraint_mode","q_TJ","chen_boundary_status","complete_practical_chen","joint_positive")} for r in summary]);write(out/"tj_constraint_mode_diagnostics.csv",history);write(out/"tj_constraint_mode_rejected_cases.csv",[r for r in summary if not r["joint_positive"]]);plots(out,summary,history);print("DONE TJ ablation",flush=True)


if __name__=="__main__":main()
