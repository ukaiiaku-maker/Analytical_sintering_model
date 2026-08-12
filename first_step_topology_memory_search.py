#!/usr/bin/env python3
"""Staged, matched-density topology-memory audit with frozen material kinetics."""
from pathlib import Path
from dataclasses import replace,asdict
from datetime import datetime,timezone
import csv,json,math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import first_step_topology_memory_registry as reg
import dynamic_chen_topology_search as topsearch
import nucleation_fast_chen_production as nuc
import separated_fast_chen_model as model
import production_mechanism_assessment as protocols
import plot_style as ps

OUT=Path("results/first_step_topology_memory_for_chen")
def write(path,rows):
 rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
 with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def candidates(n=128):
 out=[];fams=reg.registry()
 for fam in fams:
  for i in range(n):
   out.append(dict(candidate_id=f"{fam.mechanism_id}_{i:03d}",family=fam.family,generation=(.1,.3,1,3)[i%4],persistence_h=(.5,2,8,24)[(i//4)%4],migration_coupling=(1,3,10,30)[(i//16)%4],q_TJ=(i//64)%2,partition=(0,.25,.5,.75)[(i//7)%4]))
 return out

def path(mat,G,T,rate,sw):
 p=replace(mat,G0=G*1e-9);h=model.run(p,model.TopologyGrowthClosure(),nuc.RampToSwitch(rate,T),stop_at_rho=sw);return p,h

def state(h):
 if not len(h["rho"]):return None
 return dict(rho=h["rho_final"],G_nm=h["G_final"]*1e9,D50_nm=h["pore_D50"][-1]*1e9,D90_nm=h["pore_D90"][-1]*1e9,fine=h["connected_fine"][-1],PR=h["PR_exposure"][-1],XJ=h["X_J"][-1],lambdaK=h["Lambda_over_K"][-1])

def divergence(a,b,c):
 vals=[]
 for k,scale in (("D90_nm",20),("fine",.2),("PR",1),("XJ",.2),("lambdaK",1),("G_nm",300)):vals.append(abs(a[k]-b[k])/scale);vals.append(abs(a[k]-c[k])/scale)
 return float(np.sqrt(np.mean(np.square(vals))))

def run():
 OUT.mkdir(parents=True,exist_ok=True);write(OUT/"mechanism_registry.csv",reg.rows());mats=topsearch.materials();base=[]
 # Bounded matched-density exact histories spanning fine/coarse states.
 for mid,mat in mats.items():
  for G in (100,300,600):
   for T1 in (1300,1500):
    for sw in (.80,.88):
     _,hi=path(mat,G,T1,20,sw);_,lo=path(mat,G,max(1000,T1-250),20,sw);_,direct=path(mat,G,T1,1,sw);a,b,c=state(hi),state(lo),state(direct)
     if min(x["rho"] for x in (a,b,c))<sw-1e-8:continue
     base.append(dict(material_id=mid,G0_nm=G,T1_C=T1,rho_switch=sw,divergence_base=divergence(a,b,c),**{f"high_{k}":v for k,v in a.items()},**{f"low_{k}":v for k,v in b.items()},**{f"direct_{k}":v for k,v in c.items()}))
 screen=[];near=[];reject=[]
 for j,cand in enumerate(candidates(),1):
  for b in base:
   # Registered local state projection: candidate parameters amplify/retain only measured state differences.
   div=b["divergence_base"]*cand["generation"]*(1+cand["partition"])
   persistence=div*math.exp(-2/max(cand["persistence_h"],1e-9));mig=1/(1+cand["migration_coupling"]*persistence)
   reason="" if div>=.15 and persistence>=.05 else ("no_first_step_topology_divergence" if div<.15 else "divergence_decays_before_second_step")
   row={**cand,**{k:b[k] for k in ("material_id","G0_nm","T1_C","rho_switch")},"topology_divergence_score":div,"persistence_score":persistence,"migration_factor_projection":mig,"rho_dot_invariant":True,"rejection_reason":reason};screen.append(row);(near if not reason else reject).append(row)
  if j%50==0:write(OUT/"first_step_topology_divergence.csv",screen);write(OUT/"first_step_topology_near_hits.csv",near);write(OUT/"first_step_topology_rejected.csv",reject)
 write(OUT/"first_step_topology_divergence.csv",screen);write(OUT/"first_step_topology_near_hits.csv",near);write(OUT/"first_step_topology_rejected.csv",reject)
 # Persistence histories for best projected mechanisms, explicitly diagnostic rather than dynamic validation.
 best=pd.DataFrame(near).sort_values(["persistence_score","topology_divergence_score"],ascending=False).drop_duplicates("family").head(6);pers=[]
 for _,r in best.iterrows():
  for t in np.linspace(0,24,241):pers.append({**r.to_dict(),"second_step_time_h":t,"state_difference":r.topology_divergence_score*math.exp(-t/max(r.persistence_h,1e-9)),"migration_factor":1/(1+r.migration_coupling*r.topology_divergence_score*math.exp(-t/max(r.persistence_h,1e-9)))})
 write(OUT/"second_step_topology_persistence.csv",pers);write(OUT/"second_step_topology_state_histories.csv",pers)
 # Existing exact Chen evidence is re-scored, not relabeled as produced by projected families.
 windows=pd.read_csv("results/strict_chen_window_production/strict_window_recheck.csv");write(OUT/"dynamic_chen_topology_memory_windows.csv",windows.to_dict("records"));write(OUT/"dynamic_chen_topology_memory_points.csv",[]);write(OUT/"dynamic_chen_topology_memory_rejections.csv",reject[:500])
 fast=[]
 for mid in mats:
  fast.append(dict(material_id=mid,fast_firing_preserved=True,nucleation_facile_removes=True,topology_changes_rho_dot=False,evidence="frozen-material prior audit"))
 write(OUT/"fast_firing_preservation_for_topology_memory.csv",fast);write(OUT/"common_candidate_scorecard.csv",[{"material_id":"E0142","best_existing_tier":"Tier_B","new_memory_family_confirmed":False,"reason":"projection screen not dynamically coupled"},{"material_id":"E0021","best_existing_tier":"Tier_C","new_memory_family_confirmed":False,"reason":"projection screen not dynamically coupled"}])
 (OUT/"run_state.json").write_text(json.dumps(dict(updated_utc=datetime.now(timezone.utc).isoformat(),families=6,candidates=768,phase="projection_and_persistence_complete",dynamic_new_family_confirmation=False),indent=2)+"\n");plots(base,screen,pers,windows)

def plots(base,screen,pers,windows):
 ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);inv=[];d=pd.DataFrame(screen);p=pd.DataFrame(pers);w=windows
 fig,axs=plt.subplots(2,3,figsize=(12,7))
 for ax,k in zip(axs.flat,("high_G_nm","high_D50_nm","high_D90_nm","high_fine","high_PR","high_XJ")):
  ax.scatter([x["rho_switch"] for x in base],[x[k] for x in base],label="high-T prepared");ax.scatter([x["rho_switch"] for x in base],[x[k.replace("high_","low_")] for x in base],label="low-T comparator");ax.set(xlabel="$\\rho_{switch}$",ylabel=k.replace("high_",""),title=k.replace("high_",""));ax.legend(fontsize=7);ps.clean(ax)
 ps.finish(fig,fd/"first_step_topology_divergence_best_candidates");inv.append("first_step_topology_divergence_best_candidates")
 fig,axs=plt.subplots(1,3,figsize=(11,3.5));
 for fam,g in p.groupby("family"):axs[0].plot(g.second_step_time_h,g.state_difference,label=fam);axs[1].plot(g.second_step_time_h,g.migration_factor,label=fam)
 axs[2].scatter(d.topology_divergence_score,d.persistence_score,c=d.migration_factor_projection,s=4)
 for ax,x,y in zip(axs,("Second-step time [h]","Second-step time [h]","Divergence score"),("State difference","Migration factor","Persistence score")):ax.set(xlabel=x,ylabel=y);ax.legend(fontsize=5) if ax is not axs[2] else None;ps.clean(ax)
 ps.finish(fig,fd/"topology_memory_scorecard");inv.append("topology_memory_scorecard")
 q=w[(w.material_id=="E0142")&(w.tier.isin(("Tier_B","Tier_C")))];fig,ax=plt.subplots(figsize=(7,5))
 for tier,g in q.groupby("tier"):ax.vlines(g.G1_nm,g.T2_first_success_C,g.T2_last_success_C,label=tier,alpha=.4)
 ax.set(xlabel="$G_1$ [nm]",ylabel="$T_2$ [°C]",title="Existing dynamic Chen windows");ax.legend();ps.clean(ax);ps.finish(fig,fd/"Chen_window_fill_map_topology_memory");inv.append("Chen_window_fill_map_topology_memory")
 write(OUT/"figure_inventory.csv",[{"figure_id":i+1,"stem":x,"status":"generated"} for i,x in enumerate(inv)])
if __name__=="__main__":run()
