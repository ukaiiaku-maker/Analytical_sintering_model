#!/usr/bin/env python3
"""Reduced, fixed-parameter screen of persistent defect memory."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import heterogeneous_initial_state_model as hetero
import residual_stress_memory_model as residual
import persistent_defect_topology_stress_model as persistent
import joint_heterogeneity_residual_stress_search as search
import plot_style as ps

OUT=Path("results/persistent_defect_topology_stress_memory")


def design():
    rows=[]
    for mode in persistent.MODES:
        for rho0 in (.65,.70,.75):
            for G in (75.,150.,225.):
                for peak in (1350.,1400.,1450.,1500.):
                    rows.append(dict(case_id=f"{mode[:3]}_r{rho0:.2f}_G{G:g}_T{peak:g}",mode=mode,rho0=rho0,G0_nm=G,peak_T_C=peak))
    return rows


def run_path(row,rate):
    # C091 is the physical seed, but rho0/G0/T are independently crossed.
    seed=search.combined_design()[91];seed={**seed,"rho0":row["rho0"],"G0_mean_nm":row["G0_nm"],"peak_T_C":row["peak_T_C"]}
    hp,rp=search.params_from(seed);pp=persistent.PersistentParams(mode=row["mode"]);items=[]
    for spec,cp in hetero.cohort_params(hp):items.append((spec,persistent.run(cp,search.prior.FastSchedule(rate,row["peak_T_C"],8),rp,pp,spec.stress_sign)))
    return hetero.aggregate_histories(items)


def task(row):
    ref=run_path(row,1);out=[];curves=[];hist=[]
    for rate in (20.,100.):
        fast=run_path(row,rate);cs,score=search.classify(ref,fast);out.append({**row,"fast_rate_C_min":rate,**score})
        for metric,c in cs.items():curves += [{"case_id":row["case_id"],"mode":row["mode"],"fast_rate_C_min":rate,"metric":metric,**x} for x in c.to_dict("records")]
        if row["rho0"]==.65 and row["G0_nm"]==150 and row["peak_T_C"]==1500 and rate==100:
            for label,h in (("reference",ref),("fast",fast)):
                stride=max(1,len(h["t"])//300)
                keys=("rho","G_mean","G50","G90","pore_D90","large_pore_fraction","connected_fine_pore_fraction","cumulative_PR_desintering_work","f_defect_large_pore","f_crack_like_pore","defect_D90","defect_connectedness","stored_PR_work","stored_shear_coupled_stress","persistent_eligibility","persistent_growth_factor")
                hist += [{"case_id":row["case_id"],"mode":row["mode"],"path":label,"t_s":h["t"][i],**{k:h[k][i] for k in keys if k in h}} for i in range(0,len(h["t"]),stride)]
    return out,curves,hist


def plots(summary,curves,hist):
    ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);s=pd.DataFrame(summary);c=pd.DataFrame(curves);h=pd.DataFrame(hist)
    top=s.sort_values("max_ratio").iloc[-1];q=c[(c.case_id==top.case_id)&(c.fast_rate_C_min==top.fast_rate_C_min)&(c.metric=="G_mean")]
    fig,ax=plt.subplots(figsize=(5,3.5));ax.plot(q.rho,q.ratio);[ax.axhline(v,color="#777",ls="--") for v in (1.2,1.5,2)];ax.set(xlabel="Density",ylabel="$G_{reference}/G_{fast}$",title=f"{top.case_id}: unattained finite window");ps.clean(ax);ps.finish(fig,fd/"strongest_rejected_ratio")
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.4),constrained_layout=True);ph=h[h["mode"]=="persistent_defect_memory"]
    for path,z in ph.groupby("path"):
        for ax,k,y in zip(axs.flat,("f_defect_large_pore","stored_PR_work","persistent_eligibility","persistent_growth_factor"),("Large-pore defect fraction","Stored PR work","Densification eligibility","Growth factor")):ax.plot(z.rho,z[k],label=path);ax.set(xlabel="Density",ylabel=y);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);ps.finish(fig,fd/"persistent_state_histories")
    fig,ax=plt.subplots(figsize=(5,3.5));
    for mode,z in s.groupby("mode"):ax.scatter(z.max_ratio,z.meaningful_span,label=mode,alpha=.6,s=15)
    ax.axhline(.03,color="#555",ls="--");ax.set(xlabel="Maximum mean-grain ratio",ylabel="Qualified density span");ax.legend(fontsize=7);ps.clean(ax);ps.finish(fig,fd/"effect_vs_span")
    fig,ax=plt.subplots(figsize=(5.5,3.5));z=s.groupby(["mode","both_paths_attained"]).size().unstack(fill_value=0);z.plot.bar(stacked=True,ax=ax);ax.set(xlabel="Mode",ylabel="Comparisons",title="Joint interval attainment");ps.clean(ax);ps.finish(fig,fd/"attainment_by_mode")


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--workers",type=int,default=4);ap.add_argument("--plots-only",action="store_true");args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.plots_only:
        plots(pd.read_csv(OUT/"screen_summary.csv").to_dict("records"),pd.read_csv(OUT/"ratio_curves.csv").to_dict("records"),pd.read_csv(OUT/"representative_histories.csv").to_dict("records"));return
    rows=design();summary=[];curves=[];hist=[]
    mapped=map(task,rows);pool=None
    if args.workers>1:pool=ProcessPoolExecutor(max_workers=args.workers);mapped=pool.map(task,rows,chunksize=1)
    try:
        for i,(a,b,c) in enumerate(mapped,1):summary+=a;curves+=b;hist+=c;print("persistent",i,"/",len(rows),flush=True) if i%12==0 else None
    finally:
        if pool:pool.shutdown()
    rejected=[{**r,"rejection_reason":"unattained_interval" if not r["both_paths_attained"] else "ratio_below_1.5_over_Drho_0.03"} for r in summary if not r["meaningful"]];meaningful=[r for r in summary if r["meaningful"]]
    search.write(OUT/"screen_summary.csv",summary);search.write(OUT/"ratio_curves.csv",curves);search.write(OUT/"representative_histories.csv",hist);search.write(OUT/"rejected_cases.csv",rejected);search.write(OUT/"meaningful_cases.csv",meaningful,list(summary[0]) if summary else None);plots(summary,curves,hist)
    print("DONE",len(summary),"comparisons; meaningful",len(meaningful))

if __name__=="__main__":main()
