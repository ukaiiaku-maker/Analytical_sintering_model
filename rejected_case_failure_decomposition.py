#!/usr/bin/env python3
"""Decompose transient high-ratio rejected cases without changing physics."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import joint_heterogeneity_residual_stress_search as search
import plot_style as ps

OUT=Path("results/rejected_case_failure_decomposition")


def history(case):
    row={**search.combined_design()[int(case[1:])],"persist_history":True}
    hp,rp=search.params_from(row);paths={}
    for label,rate in (("reference",1),("fast",row["rate"])):
        h=search.run_path(hp,rp,rate,row["peak_T_C"],row["hold_h"]);t=h["t"]
        h["rho_dot"]=np.gradient(h["rho"],t,edge_order=1);h["G_mean_dot_nm_s"]=np.gradient(h["G_mean"]*1e9,t,edge_order=1)
        paths[label]=h
    return row,paths


def sampled_rows(case,row,paths):
    rows=[]
    for label,h in paths.items():
        stride=max(1,len(h["t"])//500)
        for i in list(range(0,len(h["t"]),stride))+[len(h["t"])-1]:
            rows.append(dict(case_id=case,path=label,t_s=h["t"][i],T_C=search.prior.FastSchedule(1 if label=="reference" else row["rate"],row["peak_T_C"],row["hold_h"]).T(h["t"][i],h["rho"][i]),**{k:h[k][i] for k in ("rho","G_mean","G50","G90","pore_D50","pore_D90","large_pore_fraction","connected_fine_pore_fraction","cumulative_PR_desintering_work","sigma_res_GBseg","sigma_res_TJ","sigma_res_large_pore","sigma_res_crack_like","rho_dot","G_mean_dot_nm_s") if k in h}))
    return rows


def plots(case,row,paths,ratio):
    ps.apply_style();OUT.mkdir(parents=True,exist_ok=True);ref,fast=paths["reference"],paths["fast"]
    fig,axs=plt.subplots(3,2,figsize=(7.2,8),constrained_layout=True)
    for label,h,col in (("reference",ref,ps.COLORS["slow"]),("fast",fast,ps.COLORS["fast"])):
        th=h["t"]/3600;T=[search.prior.FastSchedule(1 if label=="reference" else row["rate"],row["peak_T_C"],row["hold_h"]).T(t,r) for t,r in zip(h["t"],h["rho"])]
        for ax,x,y in ((axs[0,0],th,h["rho"]),(axs[0,1],th,T),(axs[1,0],th,h["G_mean"]*1e9),(axs[1,1],h["rho"],h["G50"]*1e9),(axs[2,0],h["rho"],h["G90"]*1e9),(axs[2,1],h["rho"],h["G_mean"]*1e9)):ax.plot(x,y,label=label,color=col)
    labels=(("Time [h]","Density"),("Time [h]","Temperature [°C]"),("Time [h]","$G_{mean}$ [nm]"),("Density","$G_{50}$ [nm]"),("Density","$G_{90}$ [nm]"),("Density","$G_{mean}$ [nm]"))
    for ax,(x,y) in zip(axs.flat,labels):ax.set(xlabel=x,ylabel=y);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);ps.finish(fig,OUT/f"{case}_trajectory_decomposition")
    fig,axs=plt.subplots(3,2,figsize=(7.2,8),constrained_layout=True)
    fields=(("pore_D50",1e9,"Pore D50 [nm]"),("pore_D90",1e9,"Pore D90 [nm]"),("large_pore_fraction",1,"Large-pore fraction"),("connected_fine_pore_fraction",1,"Connected fine-pore fraction"),("cumulative_PR_desintering_work",1,"Cumulative PR work"),("rho_dot",1,"Densification rate [s$^{-1}$]"))
    for label,h,col in (("reference",ref,ps.COLORS["slow"]),("fast",fast,ps.COLORS["fast"])):
        for ax,(field,scale,ylabel) in zip(axs.flat,fields):ax.plot(h["rho"],h[field]*scale,label=label,color=col);ax.set(xlabel="Density",ylabel=ylabel);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);ps.finish(fig,OUT/f"{case}_pore_rate_decomposition")
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.4),constrained_layout=True)
    for label,h,col in (("reference",ref,ps.COLORS["slow"]),("fast",fast,ps.COLORS["fast"])):
        axs[0,0].plot(h["rho"],h["sigma_res_GBseg"]/1e6,label=label,color=col);axs[0,1].plot(h["rho"],h["sigma_res_TJ"]/1e6,label=label,color=col);axs[1,0].plot(h["rho"],h["sigma_res_large_pore"]/1e6,label=label,color=col);axs[1,1].plot(h["rho"],h["G_mean_dot_nm_s"],label=label,color=col)
    for ax,y in zip(axs.flat,("GB-segment stress [MPa]","TJ stress [MPa]","Large-pore stress [MPa]","Grain-growth rate [nm/s]")):ax.set(xlabel="Density",ylabel=y);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);ps.finish(fig,OUT/f"{case}_stress_growth_decomposition")
    fig,ax=plt.subplots(figsize=(5,3.5));ax.plot(ratio.rho,ratio.ratio,color=ps.COLORS["fast"]);[ax.axhline(v,color="#777",ls="--") for v in (1.2,1.5,2)];ax.axvspan(.85,min(.88,ratio.rho.max()),color="#E69F00",alpha=.15);ax.set(xlabel="Density",ylabel="$G_{reference}/G_{fast}$",title="High ratio terminates at reference nonattainment");ps.clean(ax);ps.finish(fig,OUT/f"{case}_ratio_termination")


def main():
    OUT.mkdir(parents=True,exist_ok=True);screen=pd.read_csv("results/heterogeneity_residual_stress_search/combined_screen.csv");cases=screen[screen.max_ratio>=1.5].case_id.tolist();summ=[];allhist=[];allratio=[]
    for case in cases:
        row,paths=history(case);ratio=search.ratio_curve(paths["reference"],paths["fast"],"G_mean");ratio["case_id"]=case;allratio+=ratio.to_dict("records");allhist+=sampled_rows(case,row,paths);r=ratio[ratio.ratio>=1.5];span=search.effect.longest_span(ratio.rho.to_numpy(),ratio.ratio.to_numpy(),1.5)
        ref,fast=paths["reference"],paths["fast"];summ.append(dict(case_id=case,max_ratio=ratio.ratio.max(),ratio_ge_1p5_span=span,rho_reference_final=ref["rho"][-1],rho_fast_final=fast["rho"][-1],reference_attains_0p92=ref["rho"].max()>=.92,fast_attains_0p92=fast["rho"].max()>=.92,reference_final_rho_dot=ref["rho_dot"][-1],fast_final_rho_dot=fast["rho_dot"][-1],termination_reason="reference_target_nonattainment_at_fixed_schedule_end",stress_persistence_fraction=ref["sigma_res_large_pore"][-1]/max(np.max(ref["sigma_res_large_pore"]),1e-300),late_stage_required=False,recommended_next_physics="persistent_defect_topology_stress_memory"));plots(case,row,paths,ratio)
    search.write(OUT/"failure_summary.csv",summ);search.write(OUT/"path_histories.csv",allhist);search.write(OUT/"ratio_curves.csv",allratio);print("decomposed",len(cases),"transient high-ratio cases")

if __name__=="__main__":main()
