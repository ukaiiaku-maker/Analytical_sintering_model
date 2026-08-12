#!/usr/bin/env python3
"""96-case fractional screen for local connected-sink defect memory."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import local_connected_sink_defect_model as local
import joint_heterogeneity_residual_stress_search as parent
import plot_style as ps

OUT=Path("results/local_connected_sink_defect_memory")


def design():
    modes=local.MODES;rows=[]
    for i in range(96):
        rows.append(dict(case_id=f"L{i:03d}",mode=modes[i%len(modes)],rho0=(.65,.70,.75)[i%3],G0_nm=(75.,150.,225.)[(i//3)%3],defect_weight=(.02,.05,.10,.20)[i%4],defect_large_pore_factor=(4.,8.,12.)[(i//4)%3],matrix_connected_fraction=(.5,.7,.9)[(i//5)%3],stress_scale=(0.,.5,1.,2.)[(i//6)%4],tau_defect_over_matrix=(1.,3.,10.,30.)[(i//7)%4],peak_T_C=(1350.,1400.,1450.,1500.)[(i//8)%4],hold_h=(0.,2.,8.,20.)[(i//12)%4]))
    return rows


def params(row):
    base=parent.base_params("mech_019_q0",row["G0_nm"],row["rho0"],"GBseg_rich")
    return local.LocalMixtureParams(base,**{k:row[k] for k in ("mode","G0_nm","rho0","defect_weight","defect_large_pore_factor","matrix_connected_fraction","stress_scale","tau_defect_over_matrix")})


def at_density(h,rho):return parent.at_density(h,rho)


def task(row):
    p=params(row);ref,refl=local.run(p,parent.prior.FastSchedule(1,row["peak_T_C"],row["hold_h"]));summ=[];curves=[];grains=[];pores=[];stress=[];hist=[]
    for rate in (20.,100.):
        fast,fastl=local.run(p,parent.prior.FastSchedule(rate,row["peak_T_C"],row["hold_h"]));cs,score=parent.classify(ref,fast);censored=bool(ref.get("numerical_censored",False) or fast.get("numerical_censored",False))
        if censored:score={**score,"both_paths_attained":False,"meaningful":False,"meaningful_span":0.,"classification":"numerical_censored"}
        summ.append({**row,"fast_rate_C_min":rate,**score,"numerical_censored":censored,"trajectory_gate_passed":score["meaningful"]})
        for metric,c in cs.items():curves += [{"case_id":row["case_id"],"mode":row["mode"],"fast_rate_C_min":rate,"metric":metric,**x} for x in c.to_dict("records")]
        for rho in (.85,.88,.90,.92):
            a=at_density(ref,rho);b=at_density(fast,rho)
            if a and b:
                grains.append(dict(case_id=row["case_id"],fast_rate_C_min=rate,rho=rho,**{f"{k}_{side}_nm":v[k]*1e9 for side,v in (("reference",a),("fast",b)) for k in ("G_mean","G50","G90")}))
                pores.append(dict(case_id=row["case_id"],fast_rate_C_min=rate,rho=rho,pore_D90_reference_nm=a["pore_D90"]*1e9,pore_D90_fast_nm=b["pore_D90"]*1e9,large_pore_reference=a["large_pore_fraction"],large_pore_fast=b["large_pore_fraction"],fine_reference=a["connected_fine_pore_fraction"],fine_fast=b["connected_fine_pore_fraction"]))
        if row["case_id"] in ("L000","L004","L094") and rate==100:
            for path,agg,lh in (("reference",ref,refl),("fast",fast,fastl)):
                for name,z in lh.items():
                    stride=max(1,len(z["t"])//250)
                    for j in range(0,len(z["t"]),stride):
                        rec=dict(case_id=row["case_id"],mode=row["mode"],path=path,subpopulation=name,t_s=z["t"][j],rho_i=z["rho"][j],G_i_nm=z["G"][j]*1e9,pore_D50_i_nm=z["pore_mean_radius"][j]*1e9,pore_D90_i_nm=z["pore_mean_radius"][j]*4,connected_fine_i=z["connected_fine_pore_fraction"][j],large_pore_i=z["large_pore_fraction"][j],rho_dot_i=z["rho_dot"][j],G_dot_i_nm_s=z["G_dot"][j]*1e9,PR_work_i=z["cumulative_PR_desintering_work"][j]);
                        for k in ("sigma_res_GBseg","sigma_res_TJ","sigma_res_large_pore","f_defect_large_pore","stored_PR_work"):rec[k]=z[k][j] if k in z else 0
                        hist.append(rec);stress.append({k:rec[k] for k in ("case_id","mode","path","subpopulation","t_s","rho_i","sigma_res_GBseg","sigma_res_TJ","sigma_res_large_pore")})
    return summ,curves,grains,pores,stress,hist


def plots(summary,curves,hist,grains=None,stress=None):
    ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);s=pd.DataFrame(summary);c=pd.DataFrame(curves);h=pd.DataFrame(hist);g=pd.DataFrame(grains or []);st=pd.DataFrame(stress or []);top=s.sort_values("max_ratio").iloc[-1];q=c[(c.case_id==top.case_id)&(c.fast_rate_C_min==top.fast_rate_C_min)&(c.metric=="G_mean")]
    fig,ax=plt.subplots(figsize=(5,3.5));ax.plot(q.rho,q.ratio);[ax.axhline(v,color="#777",ls="--") for v in (1.2,1.5,2)];ax.set(xlabel="Density",ylabel="Mean-grain ratio",title="Strongest rejected local mixture");ps.clean(ax);ps.finish(fig,fd/"grain_ratio_vs_density")
    rep=h[h.case_id==("L094" if "L094" in set(h.case_id) else h.case_id.iloc[0])]
    fig,axs=plt.subplots(1,3,figsize=(7.2,3),constrained_layout=True)
    for path,z in rep.groupby("path"):
        for ax,k,y in zip(axs,("G_i_nm","rho_i","large_pore_i"),("Local grain size [nm]","Local density","Large-pore fraction")):
            for name,qh in z.groupby("subpopulation"):ax.plot(qh.rho_i,qh[k],label=f"{path}: {name}")
    for ax,y in zip(axs,("Local grain size [nm]","Local density","Large-pore fraction")):ax.set(xlabel="Local density",ylabel=y);ps.clean(ax)
    axs[0].legend(fontsize=5);ps.panel_labels(axs);ps.finish(fig,fd/"local_subpopulation_paths")
    fig,axs=plt.subplots(1,2,figsize=(7.2,3),constrained_layout=True)
    for path,z in rep.groupby("path"):
        for name,qh in z.groupby("subpopulation"):axs[0].plot(qh.rho_i,qh.sigma_res_large_pore/1e6,label=f"{path}: {name}");axs[1].plot(qh.rho_i,qh.PR_work_i,label=f"{path}: {name}")
    axs[0].set(xlabel="Local density",ylabel="Large-pore stress [MPa]");axs[1].set(xlabel="Local density",ylabel="PR work");axs[0].legend(fontsize=5);ps.panel_labels(axs);ps.finish(fig,fd/"stress_PR_by_subpopulation")
    fig,ax=plt.subplots(figsize=(6,3.5));z=s.groupby("mode").agg(max_ratio=("max_ratio","max"),attained=("both_paths_attained","sum"));ax.scatter(z.max_ratio,z.attained);[ax.text(r.max_ratio,r.attained,i,fontsize=6) for i,r in z.iterrows()];ax.set(xlabel="Maximum ratio",ylabel="Jointly attained comparisons");ps.clean(ax);ps.finish(fig,fd/"failure_mode_summary")
    fig,ax=plt.subplots(figsize=(5,3));ax.text(.5,.5,"No trajectory-qualified candidates\nChen solver not run",ha="center");ax.axis("off");ps.finish(fig,fd/"Chen_preservation_map")
    fig,ax=plt.subplots(figsize=(5,3.5));
    if len(g):
        case=g.case_id.iloc[0];zg=g[g.case_id==case]
        for metric in ("G_mean","G50","G90"):ax.plot(zg.rho,zg[f"{metric}_reference_nm"],"--",label=f"reference {metric}");ax.plot(zg.rho,zg[f"{metric}_fast_nm"],label=f"fast {metric}")
    ax.set(xlabel="Density",ylabel="Grain metric [nm]");ax.legend(fontsize=5);ps.clean(ax);ps.finish(fig,fd/"grain_metrics_vs_density")
    fig,axs=plt.subplots(1,2,figsize=(7.2,3),constrained_layout=True)
    for path,z in rep.groupby("path"):
        for name,qh in z.groupby("subpopulation"):axs[0].plot(qh.rho_i,qh.pore_D90_i_nm,label=f"{path}: {name}");axs[1].plot(qh.rho_i,qh.large_pore_i,label=f"{path}: {name}")
    axs[0].set(xlabel="Local density",ylabel="Local pore D90 [nm]");axs[1].set(xlabel="Local density",ylabel="Large-pore fraction");axs[0].legend(fontsize=5);ps.panel_labels(axs);ps.finish(fig,fd/"defect_pore_metrics")
    fig,ax=plt.subplots(figsize=(5,3.5));
    if len(st):
        zs=st[st.case_id==rep.case_id.iloc[0]]
        for (path,name),qh in zs.groupby(["path","subpopulation"]):ax.plot(qh.rho_i,qh.sigma_res_large_pore/1e6,label=f"{path}: {name}")
    ax.set(xlabel="Local density",ylabel="Large-pore residual stress [MPa]");ax.legend(fontsize=5);ps.clean(ax);ps.finish(fig,fd/"residual_stress_by_subpopulation")
    fig,ax=plt.subplots(figsize=(5,3.5));
    for (path,name),qh in rep.groupby(["path","subpopulation"]):ax.plot(qh.rho_i,qh.PR_work_i,label=f"{path}: {name}")
    ax.set(xlabel="Local density",ylabel="Cumulative PR work");ax.legend(fontsize=5);ps.clean(ax);ps.finish(fig,fd/"PR_work_by_subpopulation")


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--workers",type=int,default=4);ap.add_argument("--plots-only",action="store_true");args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.plots_only:
        plots(pd.read_csv(OUT/"reduced_screen_summary.csv").to_dict("records"),pd.read_csv(OUT/"matched_density_ratio_curves.csv").to_dict("records"),pd.read_csv(OUT/"local_subpopulation_histories.csv").to_dict("records"),pd.read_csv(OUT/"grain_metric_summary.csv").to_dict("records"),pd.read_csv(OUT/"residual_stress_summary.csv").to_dict("records"));return
    rows=design();parent.write(OUT/"parameter_registry.csv",rows);summary=[];curves=[];grains=[];pores=[];stress=[];hist=[];mapped=map(task,rows);pool=None
    if args.workers>1:pool=ProcessPoolExecutor(max_workers=args.workers);mapped=pool.map(task,rows,chunksize=1)
    try:
        for i,x in enumerate(mapped,1):a,b,c,d,e,f=x;summary+=a;curves+=b;grains+=c;pores+=d;stress+=e;hist+=f;print("local mixture",i,"/",len(rows),flush=True) if i%12==0 else None
    finally:
        if pool:pool.shutdown()
    meaningful=[r for r in summary if r["meaningful"]];rejected=[{**r,"rejection_reason":"numerical_censored" if r.get("numerical_censored") else ("reference_nonattainment" if not r["both_paths_attained"] else ("short_span_below_0.03" if r["max_ratio"]>=1.5 else "ratio_below_1.5"))} for r in summary if not r["meaningful"]]
    parent.write(OUT/"reduced_screen_summary.csv",summary);parent.write(OUT/"rejected_cases.csv",rejected);parent.write(OUT/"meaningful_trajectory_cases.csv",meaningful,list(summary[0]) if summary else None);parent.write(OUT/"weak_or_transient_cases.csv",rejected);parent.write(OUT/"local_subpopulation_histories.csv",hist);parent.write(OUT/"matched_density_ratio_curves.csv",curves);parent.write(OUT/"grain_metric_summary.csv",grains);parent.write(OUT/"pore_metric_summary.csv",pores);parent.write(OUT/"residual_stress_summary.csv",stress);parent.write(OUT/"Chen_preservation_summary.csv",[],["case_id","trajectory_gate_passed","complete_practical_window","reason"]);plots(summary,curves,hist,grains,stress);print("DONE",len(summary),"comparisons meaningful",len(meaningful))

if __name__=="__main__":main()
