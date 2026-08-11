#!/usr/bin/env python3
"""Bounded pore-occupied versus structurally constrained TJ ablation."""
from __future__ import annotations
import csv,math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
import numpy as np

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
    mid,mode=item;p0=params(mid,mode);fast=[]
    for topo,frac in production.TOPOLOGIES.items():
      p=replace(p0,base=production.fast_params(p0.base,75,.70,frac))
      for peak in (1350.,1400.,1450.):
       for hold in (8.,20.):
        paths={r:memory.run(p,production.FastSchedule(r,peak,hold)) for r in (1.,20.,100.)}
        for target in (.85,.90):
         ids={r:np.flatnonzero(h["rho"]>=target-1e-12) for r,h in paths.items()}
         for rate in (20.,100.):
          if not len(ids[1.]) or not len(ids[rate]):fast.append(dict(initial_topology=topo,peak_T_C=peak,hold_h=hold,rho_target=target,rate=rate,attained=False));continue
          i=int(ids[rate][0]);j=int(ids[1.][0]);hr=100*(paths[1.]["G"][j]-paths[rate]["G"][i])/paths[1.]["G"][j]
          fast.append(dict(initial_topology=topo,peak_T_C=peak,hold_h=hold,rho_target=target,rate=rate,attained=True,HR_pct=hr,beneficial=hr>1,C_TJ_pore=float(paths[rate]["C_TJ_pore"][i]),C_TJ_constraint=float(paths[rate]["C_TJ_constraint"][i]),C_TJ_relaxed=float(paths[rate]["C_TJ_relaxed"][i]),C_TJ_pinned=float(paths[rate]["C_TJ_pinned"][i]),R_TJ_pore_drag=float(paths[rate]["R_TJ_pore_drag"][i]),Lambda_TJ=float(paths[rate]["Lambda_TJ"][i]),K_TJ=float(paths[rate]["K_TJ"][i]),Lambda_over_K_TJ=float(paths[rate]["Lambda_over_K_TJ"][i]),P_comp_TJ=float(paths[rate]["P_comp_TJ"][i]),P_TJ_multihit=float(paths[rate]["P_TJ_multihit"][i]),P_TJ_pore_drag=float(paths[rate]["P_TJ_pore_drag"][i]),P_TJ_assisted_densification=float(paths[rate]["P_TJ_assisted_densification"][i])))
    b=replace(p0.base.action.location.base,G0=225e-9);pc=replace(p0,base=replace(p0.base,action=replace(p0.base.action,location=replace(p0.base.action.location,base=b))));h1=memory.run(pc,preparation.FixedBudgetRamp(20,1400),stop_at_rho=.84);chen_status="UNATTAINABLE_FIRST_STEP";complete=False
    if h1["rho"][-1]>=.84-1e-12:
        st=memory.final_state(h1,pc);pts=joint.adaptive_points(f"{mid}_{mode}",pc,225,1400,.84,st);s=adaptive.status(pts,.10,1400,True);chen_status=s["boundary_status"];complete=chen_status=="COMPLETE_WINDOW" and (st.base.pore.G*1e9-225)/225<=.10
    valid=[r for r in fast if r.get("attained")];benef=[r for r in valid if r.get("beneficial")];return dict(base_mechanism=mid,TJ_constraint_mode=mode,q_TJ=p0.base.q_TJ,chen_boundary_status=chen_status,complete_practical_chen=complete,n_fast_attained=len(valid),n_fast_beneficial=len(benef),n_fast_harmful=sum(r["HR_pct"]<-1 for r in valid),n_fast_neutral=sum(abs(r["HR_pct"])<=1 for r in valid),HR_pct_median=np.median([r["HR_pct"] for r in benef]) if benef else math.nan,C_TJ_pore_median=np.median([r["C_TJ_pore"] for r in benef]) if benef else math.nan,C_TJ_constraint_median=np.median([r["C_TJ_constraint"] for r in benef]) if benef else math.nan,P_TJ_multihit_median=np.median([r["P_TJ_multihit"] for r in benef]) if benef else math.nan,P_TJ_pore_drag_median=np.median([r["P_TJ_pore_drag"] for r in benef]) if benef else math.nan,joint_positive=bool(complete and benef and len(benef)<len(valid))),fast


def main():
    items=[(mid,mode) for mid in joint.frozen_base() for mode in MODES];summary=[];points=[]
    with ProcessPoolExecutor(max_workers=8) as pool:
        for s,p in pool.map(task,items,chunksize=1):summary.append(s);points.extend([{**r,"base_mechanism":s["base_mechanism"],"TJ_constraint_mode":s["TJ_constraint_mode"]} for r in p])
    out=Path("results/production_pr_desintering_assessment");write(out/"TJ_constraint_ablation.csv",summary);write(out/"TJ_constraint_fast_points.csv",points);print("DONE TJ ablation",flush=True)


if __name__=="__main__":main()
