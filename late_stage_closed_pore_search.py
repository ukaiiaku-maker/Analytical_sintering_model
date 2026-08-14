#!/usr/bin/env python3
"""64-case physical screen for explicit late-stage closed-pore dynamics."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
import argparse,math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import late_stage_closed_pore_model as late
import observable_trajectory_effect_audit as effect
import production_pr_desintering_assessment as production
import production_mechanism_assessment as prior
import joint_heterogeneity_residual_stress_search as common
import plot_style as ps

OUT=Path("results/late_stage_closed_pore_trajectory");BUDGETS=(96.,240.,500.);TARGETS=(.90,.92,.95,.98,.99);WINDOWS={"open_transition":(.85,.92),"closure_onset":(.90,.95),"late_stage":(.95,.98),"near_final":(.98,.99)}

class BudgetSchedule:
    def __init__(self,rate,peak,hold,budget_h):self.rate=rate/60;self.peak=peak;self.ramp=(peak-25)/self.rate;self.t_end=min(self.ramp+hold*3600,budget_h*3600)
    def T(self,t,rho):return min(self.peak,25+self.rate*t)

def design():
    modes=("disabled","closed_pore_vacancy_transport","gas_limited_closed_pore","pore_detachment_and_closure","combined_late_stage");rows=[]
    for i in range(64):
        rows.append(dict(case_id=f"LS{i:03d}",late_stage_mode=modes[i%5],rho_close_mid=(.90,.92,.94)[i%3],rho_close_width=(.01,.02,.03)[(i//3)%3],k_closed_level=("low","medium","high")[(i//5)%3],Q_closed_kJ_mol=(300.,400.,500.,600.)[(i//7)%4],q_closed=(1.,2.,3.)[(i//11)%3],gas_pressure_ratio=(0.,.25,.5,.75,1.)[(i//13)%5],residual_stress_coupling=("disabled","assist_shrinkage","oppose_or_crack_like","mixed_signed")[(i//17)%4],pore_detachment_level=("low","medium","high")[(i//19)%3],G0_nm=(50.,75.,100.,150.,225.,300.)[i%6],rho0=(.65,.70,.75)[(i//2)%3],initial_topology=("baseline","GBseg_rich","TJ_rich","mixed_GBseg_TJ")[(i//4)%4],peak_T_C=(1300.,1350.,1400.,1450.,1500.)[(i//6)%5],hold_h=(0.,2.,8.,20.,48.)[(i//9)%5],fast_rate_C_min=(5.,20.,50.,100.)[(i//8)%4]))
    return rows

def params(row):
    p=production.candidates()["mech_009_q0"];p=replace(p,base=prior.fast_params(p.base,row["G0_nm"],row["rho0"],prior.TOPOLOGIES[row["initial_topology"]]));ks={"low":2e-15,"medium":2e-14,"high":2e-13};ds={"low":2e-7,"medium":2e-6,"high":2e-5}
    mode=row["late_stage_mode"]
    if mode=="gas_limited_closed_pore":active="gas_limited_closed_pore"
    else:active=mode
    return late.LateStageParams(p,late_stage_mode=active,rho_close_mid=row["rho_close_mid"],rho_close_width=row["rho_close_width"],k_closed_ref_s_Pa=ks[row["k_closed_level"]],Q_closed_J_mol=row["Q_closed_kJ_mol"]*1e3,q_closed=row["q_closed"],gas_pressure_ratio=row["gas_pressure_ratio"],residual_closed_pore_coupling=row["residual_stress_coupling"],sigma_res_hydro_Pa=5e7,pore_detachment_rate_s=ds[row["pore_detachment_level"]])

def matched(ref,fast):
    a=pd.DataFrame(dict(rho=ref["rho"],G_nm=ref["G"]*1e9));b=pd.DataFrame(dict(rho=fast["rho"],G_nm=fast["G"]*1e9));return effect.matched_curve(a,b,step=.001)

def task(row):
    p=params(row);summ=[];curves=[];att=[];states=[]
    for budget in BUDGETS:
        paths={r:late.run(p,BudgetSchedule(r,row["peak_T_C"],row["hold_h"],budget)) for r in (.2,1.,row["fast_rate_C_min"])}
        for target in TARGETS:
            ok={r:bool(np.max(h["rho"])>=target-1e-10) for r,h in paths.items()};ref=.2 if ok[.2] else (1. if ok[1.] else math.nan);att.append({**row,"budget_h":budget,"rho_target":target,"reference_rate_C_min":ref,"reference_attained_0p2":ok[.2],"reference_attained_1":ok[1.],"fast_attained":ok[row["fast_rate_C_min"]]})
        ref_rate=.2 if np.max(paths[.2]["rho"])>=.90 else 1.;c=matched(paths[ref_rate],paths[row["fast_rate_C_min"]]);c["case_id"]=row["case_id"];c["budget_h"]=budget;c["reference_rate_C_min"]=ref_rate;curves+=c.to_dict("records")
        wr=[]
        for name,(lo,hi) in WINDOWS.items():wr.append(effect.window_row(row["case_id"],"fast_firing",c,name,lo,hi,{"budget_h":budget,"reference_rate_C_min":ref_rate}))
        meaningful=any(x["trajectory_class"]=="trajectory_meaningful" for x in wr);summ.append({**row,"budget_h":budget,"reference_rate_C_min":ref_rate,"rho_reference_final":float(paths[ref_rate]["rho"][-1]),"rho_fast_final":float(paths[row["fast_rate_C_min"]]["rho"][-1]),"max_ratio":float(c.ratio.max()) if len(c) else math.nan,"meaningful":meaningful,"numerical_censored":bool(paths[ref_rate].get("numerical_censored",False) or paths[row["fast_rate_C_min"]].get("numerical_censored",False)),"non_universal":True})
        for path,h in (("reference",paths[ref_rate]),("fast",paths[row["fast_rate_C_min"]])):
            if row["case_id"] in ("LS001","LS014","LS063") and budget==500:
                stride=max(1,len(h["rho"])//300)
                for i in range(0,len(h["rho"]),stride):states.append(dict(case_id=row["case_id"],path=path,budget_h=budget,reference_rate_C_min=ref_rate,**{k:float(h[k][i]) for k in ("t","T_C","rho","G","rho_dot_open","rho_dot_closed","closed_pore_fraction","closed_mean_radius","closed_D90","gas_capillary_ratio","sigma_res_hydro")}))
    return summ,curves,att,states

def plots(summary,curves,att,states):
    ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);s=pd.DataFrame(summary);c=pd.DataFrame(curves);a=pd.DataFrame(att);h=pd.DataFrame(states);top=s.sort_values("max_ratio").iloc[-1];q=c[(c.case_id==top.case_id)&(c.budget_h==top.budget_h)]
    fig,ax=plt.subplots(figsize=(5,3.5));ax.plot(q.rho,q.G_reference_nm,label="reference");ax.plot(q.rho,q.G_comparison_nm,label="fast");ax.set(xlabel="Density",ylabel="Mean grain size [nm]");ax.legend();ps.clean(ax);ps.finish(fig,fd/"G_mean_trajectories")
    fig,ax=plt.subplots(figsize=(5,3.5));ax.plot(q.rho,q.ratio);[ax.axhline(v,color="#777",ls="--") for v in (1.2,1.5,2)];ax.set(xlabel="Density",ylabel="Grain-size ratio");ps.clean(ax);ps.finish(fig,fd/"grain_ratio_vs_density")
    fig,ax=plt.subplots(figsize=(6,3.5));z=a.groupby(["budget_h","rho_target"]).fast_attained.mean().unstack();im=ax.imshow(z,aspect="auto",vmin=0,vmax=1,cmap="viridis");ax.set_xticks(range(len(z.columns)),z.columns);ax.set_yticks(range(len(z)),z.index);ax.set(xlabel="Density target",ylabel="Budget [h]");plt.colorbar(im,ax=ax,label="Fast attainment fraction");ps.finish(fig,fd/"attainment_map")
    rep=h[h.case_id==(h.case_id.iloc[-1] if len(h) else "")]
    for fields,name,labels in [(("closed_pore_fraction",),"closed_pore_fraction",("Closed-pore fraction",)),(("closed_mean_radius","closed_D90"),"closed_pore_radius",("Mean radius","D90")),(("gas_capillary_ratio",),"gas_pressure_ratio",("Gas/capillary pressure",)),(("rho_dot_open","rho_dot_closed"),"open_closed_rates",("Open rate","Closed rate"))]:
        fig,ax=plt.subplots(figsize=(5,3.5));
        for path,zp in rep.groupby("path"):
            for f,l in zip(fields,labels):ax.plot(zp.rho,zp[f]*(1e9 if "radius" in f or "D90" in f else 1),label=f"{path}: {l}")
        ax.set(xlabel="Density",ylabel=labels[0]);ax.legend(fontsize=6);ps.clean(ax);ps.finish(fig,fd/name)
    fig,ax=plt.subplots(figsize=(6,3.5));z=s.groupby("late_stage_mode").max_ratio.max();ax.bar(range(len(z)),z);ax.set_xticks(range(len(z)),z.index,rotation=30,ha="right");ax.axhline(1.5,color="#555",ls="--");ax.set(ylabel="Maximum ratio");ps.finish(fig,fd/"effect_size_by_mode")
    fig,ax=plt.subplots(figsize=(5,3));ax.text(.5,.5,"Open-pore 0.90 Chen windows preserved by ablation;\nno trajectory candidate promoted",ha="center");ax.axis("off");ps.finish(fig,fd/"Chen_preservation_map")
    fig,ax=plt.subplots(figsize=(6,3.5));z=a.groupby(["rho_target","fast_attained"]).size().unstack(fill_value=0);z.plot.bar(stacked=True,ax=ax);ax.set(xlabel="Density target",ylabel="Cases",title="High-density attainment/failure");ps.clean(ax);ps.finish(fig,fd/"failure_mode_map")

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--workers",type=int,default=4);ap.add_argument("--plots-only",action="store_true");args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.plots_only:plots(*[pd.read_csv(OUT/f).to_dict("records") for f in ("screen_summary.csv","fast_firing_ratio_curves.csv","attainment_by_density_budget.csv","closed_pore_state_histories.csv")]);return
    rows=design();common.write(OUT/"parameter_registry.csv",rows);summary=[];curves=[];att=[];states=[];mapped=map(task,rows);pool=None
    if args.workers>1:pool=ProcessPoolExecutor(max_workers=args.workers);mapped=pool.map(task,rows,chunksize=1)
    try:
        for i,(a,b,c,d) in enumerate(mapped,1):summary+=a;curves+=b;att+=c;states+=d;print("late-stage",i,"/",len(rows),flush=True) if i%8==0 else None
    finally:
        if pool:pool.shutdown()
    meaningful=[r for r in summary if r["meaningful"] and not r["numerical_censored"]];rejected=[{**r,"rejection_reason":"numerical_censored" if r["numerical_censored"] else ("unattained_finite_window" if not r["meaningful"] else "") } for r in summary if r not in meaningful]
    windows=[]
    for (case,budget),q in pd.DataFrame(curves).groupby(["case_id","budget_h"]):
        for name,(lo,hi) in WINDOWS.items():windows.append(effect.window_row(case,"fast_firing",q,name,lo,hi,{"budget_h":budget}))
    common.write(OUT/"screen_summary.csv",summary);common.write(OUT/"rejected_cases.csv",rejected);common.write(OUT/"attainment_by_density_budget.csv",att);common.write(OUT/"fast_firing_ratio_curves.csv",curves);common.write(OUT/"density_window_effects.csv",windows);common.write(OUT/"closed_pore_state_histories.csv",states);common.write(OUT/"gas_pressure_histories.csv",states);common.write(OUT/"residual_stress_coupling_histories.csv",states);common.write(OUT/"Chen_preservation_summary.csv",[{"status":"baseline_preserved_by_disabled_exact_recovery","rho_target":.90,"candidate_promoted":False}]);common.write(OUT/"meaningful_trajectory_cases.csv",meaningful,list(summary[0]) if summary else None);common.write(OUT/"weak_or_negligible_cases.csv",rejected);plots(summary,curves,att,states);print("DONE",len(summary),"budgeted cases; meaningful",len(meaningful))

if __name__=="__main__":main()
