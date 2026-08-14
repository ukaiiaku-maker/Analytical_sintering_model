#!/usr/bin/env python3
"""Production confirmation of the frozen moderate PR/de-sintering candidate."""
from __future__ import annotations
import argparse,csv,math,time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import adaptive_T2_boundary_search as adaptive
import joint_pr_desintering_search as joint
import pr_desintering_memory_model as memory
import preparation_window_search as preparation
import production_mechanism_assessment as prior
import topology_constrained_sintering as aggregate

G0S=(75.,100.,125.,150.,175.,200.,225.,250.,300.);T1S=(1325.,1350.,1375.,1400.,1425.,1450.,1475.,1500.);RATES=(1.,5.,20.,50.,100.);SWITCHES=(.76,.78,.80,.82,.84,.86,.88,.90);PREP_TOLS=(.05,.10,.20);SECOND_TOLS=(.05,.10);TARGET=.90;BUDGET=96*3600


def write(path,rows,empty=("candidate_id","status")):
    fields=[]
    for row in rows:
        for key in row:
            if key not in fields:fields.append(key)
    if not fields:fields=list(empty)
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def candidates():
    return {mid:memory.PRMemoryParams(p,"PR_plus_connected_fine_attrition",k_PR_ref_s=2e-4) for mid,p in joint.frozen_base().items()}


def _first_task(task):
    mid,p0,G,T1,rate=task;b=replace(p0.base.action.location.base,G0=G*1e-9);p=replace(p0,base=replace(p0.base,action=replace(p0.base.action,location=replace(p0.base.action.location,base=b))));h=memory.run(p,preparation.FixedBudgetRamp(rate,T1));rows=[]
    for sw in SWITCHES:
        ids=np.flatnonzero(h["rho"]>=sw-1e-12)
        if not len(ids):rows.append((dict(candidate_id=mid,G0_nm=G,T1_C=T1,heating_rate_T1_C_min=rate,rho_switch=sw,first_step_attained=False,reason="switch_density_unattainable"),None,p));continue
        i=int(ids[0]);state=memory.final_state(h,p,i);G1=state.base.pore.G*1e9;phi=[float(np.sum(h[k][i])) for k in ("phi_GBseg","phi_TJ","phi_iso")];z=max(sum(phi),1e-300);meta=dict(candidate_id=mid,G0_nm=G,T1_C=T1,heating_rate_T1_C_min=rate,rho_switch=sw,first_step_attained=True,rho1=state.base.pore.rho,G1_nm=G1,first_step_growth_fraction=(G1-G)/G,T_at_switch_C=float(h["T_C"][i]),f_GBseg=phi[0]/z,f_TJ=phi[1]/z,f_iso_location=phi[2]/z,X_J=state.base.X_J,cumulative_PR_desintering_work=state.cumulative_PR_desintering_work)
        reason="target_already_reached_first_step" if meta["rho1"]>=TARGET-1e-12 else ("first_step_growth_above_20pct" if meta["first_step_growth_fraction"]>.20 else "")
        rows.append(({**meta,"reason":reason},state if not reason else None,p))
    return rows


def first_states(workers):
    tasks=[(mid,p,G,T1,r) for mid,p in candidates().items() for G in G0S for T1 in T1S for r in RATES];allrows=[];states=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for n,result in enumerate(pool.map(_first_task,tasks,chunksize=1),1):
            for meta,state,p in result:
                allrows.append(meta)
                if state is not None:states.append((meta,state,p))
            if n%50==0:print("first groups",n,len(tasks),flush=True)
    return allrows,states


def state_key(item):
    meta,s,_=item;p=s.base.pore
    return (meta["candidate_id"],round(p.rho,13),round(p.G,18),round(s.base.X_J,13),round(s.cumulative_PR_desintering_work,4),*(a.round(15).tobytes() for a in (p.phi_GBseg,p.phi_TJ,p.phi_iso)))


def _adaptive_task(group):
    members=group;meta,state,p=members[0];points=joint.adaptive_points(meta["candidate_id"],p,meta["G0_nm"],meta["T1_C"],meta["rho_switch"],state);rows=[]
    for route,_,_ in members:
        for second_tol in SECOND_TOLS:
            for practical,map_type in ((False,"kinetic"),(True,"practical")):
                status=adaptive.status(points,second_tol,route["T1_C"],practical)
                for prep_tol in PREP_TOLS:
                    complete=(status["boundary_status"]=="COMPLETE_WINDOW" and route["first_step_growth_fraction"]<=prep_tol+1e-12)
                    rows.append({**route,"map_type":map_type,"prep_growth_tolerance":prep_tol,"second_step_growth_tolerance":second_tol,"complete_window":complete,**status})
    return rows


def chen_map(states,workers):
    groups={}
    for item in states:groups.setdefault(state_key(item),[]).append(item)
    rows=[];values=list(groups.values());print("unique PR states",len(values),"routes",len(states),flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for n,result in enumerate(pool.map(_adaptive_task,values,chunksize=1),1):
            rows.extend(result)
            if n%50==0:print("adaptive states",n,len(values),flush=True)
    return rows


def production_fast_rows():
    path=Path("results/pr_desintering_fast_firing_memory/raw_fast_firing_response_map_full.csv")
    if not path.exists():raise FileNotFoundError("full fast-firing point table is required")
    return pd.read_csv(path)


def sensitivity_design():
    base=dict(early_memory_mode="PR_plus_connected_fine_attrition",k_PR_ref_s=2e-4,Q_PR_J_mol=180e3,renewal_gate_mid=.35,renewal_power=2.,smoothing_share=.65,GB_to_TJ_share=.25,TJ_to_iso_share=.10)
    levels=[]
    for name,vals in (("k_PR_ref_s",(1e-4,2e-4,4e-4)),("Q_PR_J_mol",(150e3,180e3,220e3)),("renewal_gate_mid",(.25,.35,.45)),("renewal_power",(1.,2.,3.))):
        for v in vals:levels.append((f"{name}_{v:g}",{**base,name:v},name,v))
    shares=(("baseline",.65,.25,.10),("more_smoothing",.80,.15,.05),("more_GB_to_TJ",.50,.45,.05),("more_isolation",.50,.20,.30))
    for label,a,b,c in shares:levels.append((f"partition_{label}",{**base,"smoothing_share":a,"GB_to_TJ_share":b,"TJ_to_iso_share":c},"partition",label))
    unique={label:(kw,group,value) for label,kw,group,value in levels};return [(label,*x) for label,x in unique.items()]


def _sensitivity_task(task):
    mid,label,kw,group,value=task;p0=joint.frozen_base()[mid];p=memory.PRMemoryParams(p0,**kw)
    fast_cases={"mech_009":(50,.70,"baseline",1450,8,.90,5),"mech_009_q0":(50,.70,"GBseg_rich",1450,8,.92,100),"mech_019":(50,.75,"TJ_rich",1400,20,.88,100),"mech_019_q0":(50,.70,"mixed_GBseg_TJ",1450,8,.90,50)};Gf,rho0,topo,peak,hold,target,rate=fast_cases[mid];basep=prior.fast_params(p.base,Gf,rho0,prior.TOPOLOGIES[topo]);pf=replace(p,base=basep);paths={r:memory.run(pf,prior.FastSchedule(r,peak,hold)) for r in (1.,rate)};ids={r:np.flatnonzero(h["rho"]>=target-1e-12) for r,h in paths.items()};fast_ok=all(len(x) for x in ids.values());hr=math.nan
    if fast_ok:hr=100*(paths[1.]["G"][ids[1.][0]]-paths[rate]["G"][ids[rate][0]])/paths[1.]["G"][ids[1.][0]]
    chen_cases={"mech_009":(150,1325,1,.76),"mech_019":(100,1325,5,.82),"mech_009_q0":(175,1325,1,.82),"mech_019_q0":(100,1325,5,.82)};Gc,T1,prep_rate,sw=chen_cases[mid];b=replace(p.base.action.location.base,G0=Gc*1e-9);pc=replace(p,base=replace(p.base,action=replace(p.base.action,location=replace(p.base.action.location,base=b))));h1=memory.run(pc,preparation.FixedBudgetRamp(prep_rate,T1),stop_at_rho=sw);chen=False;status="UNATTAINABLE_FIRST_STEP"
    if h1["rho"][-1]>=sw-1e-12:
        st=memory.final_state(h1,pc);pts=joint.adaptive_points(label,pc,Gc,T1,sw,st);s=adaptive.status(pts,.10,T1,True);status=s["boundary_status"];chen=status=="COMPLETE_WINDOW" and (st.base.pore.G*1e9-Gc)/Gc<=.10
    return dict(base_mechanism=mid,sensitivity_id=label,parameter_group=group,parameter_value=value,fast_both_attained=fast_ok,HR_pct=hr,beneficial_fast=fast_ok and hr>1,chen_boundary_status=status,complete_practical_chen=chen,joint_positive=bool(chen and fast_ok and hr>1))


def parameter_sensitivity(workers):
    tasks=[(mid,*level) for mid in candidates() for level in sensitivity_design()]
    with ProcessPoolExecutor(max_workers=workers) as pool:return list(pool.map(_sensitivity_task,tasks,chunksize=1))


def histories():
    rows_fast=[];rows_two=[];p0=candidates()["mech_009_q0"]
    pf=replace(p0,base=prior.fast_params(p0.base,75,.70,prior.TOPOLOGIES["baseline"]))
    for label,rate in (("slow",1.),("fast",20.)):
        h=memory.run(pf,prior.FastSchedule(rate,1400,20));rows_fast+=history_rows(label,"fast_firing",h)
    b=replace(p0.base.action.location.base,G0=225e-9);pt=replace(p0,base=replace(p0.base,action=replace(p0.base.action,location=replace(p0.base.action.location,base=b))));h1=memory.run(pt,preparation.FixedBudgetRamp(20,1400),stop_at_rho=.84);rows_two+=history_rows("two_step","first_step",h1);state=memory.final_state(h1,pt);h2=memory.run(pt,aggregate.Iso(1250,BUDGET),initial=state);rows_two+=history_rows("two_step","second_step",h2);return rows_fast,rows_two


def history_rows(label,phase,h):
    rows=[];stride=max(1,len(h["rho"])//400);r=h["pore_radii"]
    for i in list(range(0,len(h["rho"]),stride))+[len(h["rho"])-1]:
        phi=[float(np.sum(h[k][i])) for k in ("phi_GBseg","phi_TJ","phi_iso")];z=max(sum(phi),1e-300);conn=h["phi_GBseg"][i]+h["phi_TJ"][i];zc=max(float(np.sum(conn)),1e-300)
        rows.append(dict(path=label,phase=phase,t_s=float(h["t"][i]),T_C=float(h["T_C"][i]),rho=float(h["rho"][i]),G_nm=float(h["G"][i])*1e9,cumulative_PR_desintering_work=float(h.get("cumulative_PR_desintering_work",np.zeros(len(h["rho"])))[i]),cumulative_densifying_work=float(h.get("cumulative_densifying_work",np.zeros(len(h["rho"])))[i]),cumulative_non_densifying_work=float(h.get("cumulative_non_densifying_work",np.zeros(len(h["rho"])))[i]),connected_fine_pore_fraction=float(h.get("connected_fine_pore_fraction",np.zeros(len(h["rho"])))[i]),connected_mean_radius_nm=float(np.sum(conn*r)/zc)*1e9,large_pore_fraction=float(h.get("large_pore_fraction",np.zeros(len(h["rho"])))[i]),f_GBseg=phi[0]/z,f_TJ=phi[1]/z,f_iso_location=phi[2]/z,X_J=float(h["X_J"][i]),Lambda_over_K_TJ=float(h["Lambda_over_K_TJ"][i]),P_comp_TJ=float(h["P_comp_TJ"][i]),C_GBseg=float(h["C_GBseg"][i]),C_TJ=float(h["C_TJ"][i]),f_clean_GB=float(h["f_clean_GB"][i]),f_iso=float(h["f_iso"][i]),P_GBseg_dens=float(h["P_GBseg_dens"][i]),P_TJ_dens=float(h["P_TJ_dens"][i]),P_clean_GB=float(h["P_clean_GB"][i]),P_persistent_junction_drag=float(h["P_persistent_junction_drag"][i]),P_TJ_multihit=float(h["P_TJ_multihit"][i]),sigma_base=float(h["sigma_base"][i]),sigma_GBseg=float(h["sigma_GBseg_pore"][i]),sigma_TJ=float(h["sigma_TJ_pore"][i]),sigma_total=float(h["sigma_act_total"][i])))
    return rows


def figures(out,chen,fast,hfast,htwo,score):
    fd=out/"figures";fd.mkdir(exist_ok=True);complete=[r for r in chen if r["map_type"]=="practical" and r["complete_window"] and r["prep_growth_tolerance"]==.05 and r["second_step_growth_tolerance"]==.05]
    fig,ax=plt.subplots();ax.scatter([r["G1_nm"] for r in complete],[r["T_first_success_C"] for r in complete],c=[r["T1_C"] for r in complete],s=8,cmap="viridis");ax.set(xlabel="G1 [nm]",ylabel="first success T2 [C]");fig.tight_layout();fig.savefig(fd/"figure1_practical_chen_PR.png",dpi=180);plt.close(fig)
    valid=fast[fast.comparison_attained==True];fig,ax=plt.subplots();sc=ax.scatter(valid.heating_rate_C_min,valid.peak_T_C,c=valid.HR_pct,cmap="coolwarm",vmin=-15,vmax=15,s=5,alpha=.25);ax.set(xscale="log",xlabel="rate [C/min]",ylabel="peak T [C]");fig.colorbar(sc,ax=ax,label="HR_pct");fig.tight_layout();fig.savefig(fd/"figure2_fast_firing_HR.png",dpi=180);plt.close(fig)
    for field,name,ylabel in (("connected_fine_pore_fraction","figure3_connected_fine_memory.png","connected fine fraction"),("cumulative_PR_desintering_work","figure4_PR_exposure.png","cumulative PR work")):
        fig,ax=plt.subplots();
        for label in ("slow","fast"):
            q=[r for r in hfast if r["path"]==label];ax.plot([r["rho"] for r in q],[r[field] for r in q],label=label)
        ax.set(xlabel="density",ylabel=ylabel);ax.legend();fig.tight_layout();fig.savefig(fd/name,dpi=180);plt.close(fig)
    b=valid[valid.response_class=="beneficial"];fig,ax=plt.subplots();ax.scatter(b.PR_exposure_difference,b.HR_pct,s=5,alpha=.25);ax.set(xlabel="slow-reference minus fast PR work",ylabel="HR_pct");fig.tight_layout();fig.savefig(fd/"figure5_HR_vs_PR_difference.png",dpi=180);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(14,4));
    for field,ax in zip(("X_J","Lambda_over_K_TJ","P_comp_TJ"),axs):ax.plot([r["rho"] for r in htwo],[r[field] for r in htwo]);ax.set(xlabel="density",ylabel=field)
    fig.tight_layout();fig.savefig(fd/"figure6_two_step_junction_histories.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,2.5));ax.axis("off");labels=("PR/de-sintering\nearly memory","persistent junction\ndrag","TJ multihit\nreactivation","joint processing\nresponse");
    for i,label in enumerate(labels):ax.text(.12+i*.25,.5,label,ha="center",va="center",bbox=dict(boxstyle="round",fc="white"));
    for i in range(3):ax.annotate("",xy=(.23+i*.25,.5),xytext=(.17+i*.25,.5),arrowprops=dict(arrowstyle="->"))
    fig.tight_layout();fig.savefig(fd/"figure7_mechanism_ingredients.png",dpi=180);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--outdir",default="results/production_pr_desintering_assessment");ap.add_argument("--workers",type=int,default=8);ap.add_argument("--resume-first",action="store_true");ap.add_argument("--resume-chen",action="store_true");args=ap.parse_args();out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    if args.resume_chen:
        chen=pd.read_csv(out/"raw_chen_boundaries.csv").to_dict("records")
        for r in chen:r["complete_window"]=bool(r["complete_window"])
    else:
        first,states=first_states(args.workers);write(out/"raw_first_states.csv",first);chen=chen_map(states,args.workers);write(out/"raw_chen_boundaries.csv",chen)
    fast=production_fast_rows();sens=parameter_sensitivity(args.workers);hfast,htwo=histories();summary=[];success=[];failed=[]
    for mid in candidates():
        q=[r for r in chen if r["candidate_id"]==mid];success.extend([r for r in q if r["map_type"]=="practical" and r["complete_window"]]);failed.extend([r for r in q if not r["complete_window"] or "CENSORED" in r["boundary_status"]]);summary.append(dict(candidate_id=mid,complete_practical_5_5=sum(r["map_type"]=="practical" and r["complete_window"] and r["prep_growth_tolerance"]==.05 and r["second_step_growth_tolerance"]==.05 for r in q),complete_practical_all_tolerances=sum(r["map_type"]=="practical" and r["complete_window"] for r in q)))
    score=[]
    for row in summary:
        fq=fast[fast.candidate_id==f"{row['candidate_id']}_PR_attrition_moderate"];benef=fq[fq.response_class=="beneficial"];valid=fq[fq.comparison_attained==True];score.append({**row,"beneficial_fast_count":len(benef),"harmful_count":sum(valid.response_class=="harmful"),"neutral_count":sum(valid.response_class=="neutral"),"unattainable_count":sum(fq.response_class=="unattainable"),"universal_fast":len(benef)==len(valid),"joint_positive":bool(row["complete_practical_all_tolerances"] and len(benef) and len(benef)<len(valid))})
    score.append(dict(candidate_id="disabled_control",complete_practical_5_5="production_negative_control",complete_practical_all_tolerances="production_negative_control",beneficial_fast_count=0,harmful_count="see prior production",neutral_count="see prior production",unattainable_count="see prior production",universal_fast=False,joint_positive=False))
    fast_summary=fast.groupby(["candidate_id","rho_target","initial_topology","heating_rate_C_min","peak_T_C","response_class"],dropna=False).agg(n_cases=("G0_nm","size"),HR_pct_median=("HR_pct","median"),PR_difference_median=("PR_exposure_difference","median"),fine_difference_median=("connected_fine_difference","median")).reset_index();write(out/"production_joint_scorecard.csv",score);write(out/"production_chen_window_summary.csv",summary);fast_summary.to_csv(out/"production_fast_firing_summary.csv",index=False,lineterminator="\n");write(out/"successful_practical_windows.csv",success);fast[fast.response_class=="beneficial"].to_csv(out/"successful_fast_firing_cases.csv",index=False,lineterminator="\n");write(out/"failed_or_censored_cases.csv",failed);write(out/"representative_slow_fast_histories.csv",hfast);write(out/"representative_two_step_histories.csv",htwo);write(out/"PR_parameter_sensitivity.csv",sens);write(out/"mechanism_ingredient_table.csv",[{"ingredient":"PR/de-sintering memory","role":"early slow-ramp pore-removability loss"},{"ingredient":"persistent junction drag","role":"low-T2 migration suppression"},{"ingredient":"TJ multihit reactivation","role":"upper grain-growth boundary"},{"ingredient":"censor-aware adaptive map","role":"honest finite-window classification"}]);figures(out,chen,fast,hfast,htwo,score);write(out/"runtime_summary.csv",[{"wall_s":time.perf_counter()-start,"rho_target":TARGET,"budget_h":96,"frozen_k_PR_ref_s":2e-4,"n_states":len(chen)}]);print("DONE production PR",flush=True)


if __name__=="__main__":main()
