#!/usr/bin/env python3
"""Generate exhaustive supplementary figures from compact production tables."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plot_style as ps

ROOT=Path("results/production_pr_desintering_assessment");OUT=Path("results/paper_figures/supplement")

def save_inventory(rows):ps.write_inventory(OUT/"figure_inventory.csv",rows)

def main():
    ps.apply_style();OUT.mkdir(parents=True,exist_ok=True);inv=[];ff=pd.read_csv(ROOT/"production_fast_firing_summary.csv");tj=pd.read_csv(ROOT/"tj_constraint_mode_summary.csv");sens=pd.read_csv(ROOT/"PR_parameter_sensitivity.csv");diag=pd.read_csv(ROOT/"tj_constraint_mode_diagnostics.csv");sf=pd.read_csv(ROOT/"representative_slow_fast_histories.csv");ts=pd.read_csv(ROOT/"representative_two_step_histories.csv");fail=pd.read_csv(ROOT/"failed_or_censored_cases.csv")
    # S1: each frozen base.
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.4),sharex=True,sharey=True,constrained_layout=True)
    for ax,(mid,q) in zip(axs.flat,ff.groupby("candidate_id")):
        z=q[(q.rho_target==.90)&(q.response_class!="unattainable")].pivot_table(index="peak_T_C",columns="heating_rate_C_min",values="HR_pct_median",aggfunc="median");im=ax.imshow(z.values,origin="lower",aspect="auto",cmap="RdBu_r",vmin=-8,vmax=8,extent=(np.log10(z.columns.min()),np.log10(z.columns.max()),z.index.min()-25,z.index.max()+25));ax.set_title(mid.replace("_PR_attrition_moderate",""));ax.set_xticks(np.log10(z.columns),[f"{x:g}" for x in z.columns]);ax.set(xlabel="Heating rate [°C min$^{-1}$]",ylabel="Peak $T$ [°C]")
    fig.colorbar(im,ax=axs,label=r"$HR$ improvement [%]",shrink=.8);ps.panel_labels(axs);stem="S1_frozen_base_fast_maps";ps.finish(fig,OUT/stem);inv.append(ps.inventory_row("S1",stem,"Frozen-base fast maps","production_fast_firing_summary.csv","Show response for each base","supplement"))
    # S2 q0/q1 by mode.
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),constrained_layout=True)
    for ax,(family,bases) in zip(axs,(("mech_009",["mech_009","mech_009_q0"]),("mech_019",["mech_019","mech_019_q0"]))):
        q=tj[tj.base_mechanism.isin(bases)];x=np.arange(4);modes=list(ps.TJ_COLORS);q1=[q[(q.TJ_constraint_mode==m)&(q.q_TJ==1)].n_fast_beneficial.iloc[0] for m in modes];q0=[q[(q.TJ_constraint_mode==m)&(q.q_TJ==0)].n_fast_beneficial.iloc[0] for m in modes];ax.bar(x-.2,q0,.4,color=ps.COLORS["q0"],label="$q=0$");ax.bar(x+.2,q1,.4,color=ps.COLORS["q1"],label="$q=1$");ax.set_xticks(x,[m.replace("_","\n") for m in modes],rotation=20);ax.set(title=family,ylabel="Beneficial fast cases")
    axs[0].legend();ps.panel_labels(axs);stem="S2_q0_q1_TJ_modes";ps.finish(fig,OUT/stem);inv.append(ps.inventory_row("S2",stem,"q0 versus q1","tj_constraint_mode_summary.csv","Resolve size-exponent sensitivity by TJ mode","supplement"))
    # S3 PR OAT sensitivity.
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.4),sharey=True,constrained_layout=True)
    for ax,(mid,q) in zip(axs.flat,sens.groupby("base_mechanism")):
        x=np.arange(len(q));ax.scatter(x,q.HR_pct,c=q.joint_positive.map({True:ps.COLORS["positive"],False:ps.COLORS["negative"]}),s=24);ax.axhline(1,color="#444",ls="--");ax.set(title=mid,ylabel=r"$HR$ improvement [%]",xlabel="OAT variant index")
    ps.panel_labels(axs);stem="S3_PR_parameter_sensitivity";ps.finish(fig,OUT/stem);inv.append(ps.inventory_row("S3",stem,"PR parameter sensitivity","PR_parameter_sensitivity.csv","Show bounded OAT robustness","supplement"))
    # S4 full TJ population histories.
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.4),constrained_layout=True)
    fields=(("C_TJ_pore",r"$C_{\mathrm{TJ}}^{\mathrm{pore}}$"),("C_TJ_constraint",r"$C_{\mathrm{TJ}}^{\mathrm{constraint}}$"),("C_TJ_relaxed",r"$C_{\mathrm{TJ}}^{\mathrm{relaxed}}$"),("C_TJ_pinned",r"$C_{\mathrm{TJ}}^{\mathrm{pinned}}$"))
    q=diag[diag.base_mechanism=="mech_009"]
    for ax,(field,label) in zip(axs.flat,fields):
        for mode,col in ps.TJ_COLORS.items():z=q[q.TJ_constraint_mode==mode];ax.plot(z.rho,z[field],label=mode,color=col)
        ax.set(xlabel=ps.LABELS["rho"],ylabel=label);ps.clean(ax)
    axs[0,0].legend(fontsize=6);ps.panel_labels(axs);stem="S4_TJ_population_histories";ps.finish(fig,OUT/stem);inv.append(ps.inventory_row("S4",stem,"TJ population histories","tj_constraint_mode_diagnostics.csv","Show full population decomposition","supplement"))
    # S5 detailed power and stress histories.
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.4),constrained_layout=True)
    for path,col in (("slow",ps.COLORS["slow"]),("fast",ps.COLORS["fast"])):
        q=sf[sf.path==path];axs[0,0].plot(q.rho,q.P_GBseg_dens+q.P_TJ_dens,color=col,label=path);axs[0,1].plot(q.rho,q.sigma_total,color=col,label=path)
    axs[1,0].plot(ts.rho,ts.P_persistent_junction_drag,color=ps.COLORS["two_step"]);axs[1,1].plot(ts.rho,ts.P_TJ_multihit,color=ps.COLORS["two_step"])
    for ax,l in zip(axs.flat,("Densification power [model units]","Activation stress [Pa]","Persistent-junction dissipation","TJ multihit dissipation")):ax.set(xlabel=ps.LABELS["rho"],ylabel=l);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);stem="S5_power_stress_histories";ps.finish(fig,OUT/stem);inv.append(ps.inventory_row("S5",stem,"Power and stress histories","representative histories","Expose detailed diagnostic channels","supplement"))
    # S6 failure/censoring.
    f=fail.groupby("boundary_status").n_cases.sum().sort_values();fig,ax=plt.subplots(figsize=(6.4,3.8));ax.barh([x.replace("_"," ") for x in f.index],f.values,color="#7A7A7A");ax.set(xlabel="Number of classified routes",ylabel="Boundary status");ps.clean(ax,False);stem="S6_failure_censoring";ps.finish(fig,OUT/stem);inv.append(ps.inventory_row("S6",stem,"Failure and censoring","failed_or_censored_cases.csv","Document honest boundary failures","supplement"))
    # S7 attainment completeness.
    a=ff.groupby(["peak_T_C","heating_rate_C_min","response_class"]).n_cases.sum().reset_index();tot=a.groupby(["peak_T_C","heating_rate_C_min"]).n_cases.transform("sum");a["fraction"]=a.n_cases/tot;u=a[a.response_class=="unattainable"].pivot(index="peak_T_C",columns="heating_rate_C_min",values="fraction").fillna(0);fig,ax=plt.subplots(figsize=(5.5,3.8));im=ax.imshow(u.values,origin="lower",aspect="auto",cmap="magma",vmin=0,vmax=1);ax.set_xticks(range(len(u.columns)),[f"{x:g}" for x in u.columns]);ax.set_yticks(range(len(u.index)),[f"{x:g}" for x in u.index]);ax.set(xlabel="Heating rate [°C min$^{-1}$]",ylabel="Peak $T$ [°C]");fig.colorbar(im,ax=ax,label="Unattainable fraction [-]");stem="S7_attainment_map";ps.finish(fig,OUT/stem);inv.append(ps.inventory_row("S7",stem,"Attainment map","production_fast_firing_summary.csv","Show parameter-space completeness","supplement"))
    save_inventory(inv);print("generated",len(inv),"supplement figures")

if __name__=="__main__":main()
