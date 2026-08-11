#!/usr/bin/env python3
"""Generate deterministic main-text figures from committed production results."""
from pathlib import Path
from dataclasses import replace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

import plot_style as ps
import production_pr_desintering_assessment as prod
import pr_desintering_memory_model as memory
import topology_constrained_sintering as aggregate

ROOT=Path("results/production_pr_desintering_assessment");OUT=Path("results/paper_figures/main")

def data():
    names=("representative_slow_fast_histories","representative_two_step_histories","production_fast_firing_summary","successful_practical_windows","production_joint_scorecard","tj_constraint_mode_summary","tj_constraint_mode_diagnostics","tj_constraint_mode_fast_firing")
    return {n:pd.read_csv(ROOT/f"{n}.csv") for n in names}

def highT_history():
    p0=prod.candidates()["mech_009_q0"];b=replace(p0.base.action.location.base,G0=225e-9);p=replace(p0,base=replace(p0.base,action=replace(p0.base.action,location=replace(p0.base.action.location,base=b))));h=memory.run(p,aggregate.Iso(1400,prod.BUDGET));return pd.DataFrame(prod.history_rows("highT","highT_isothermal",h))

def add_box(ax,xy,w,h,text,color):
    patch=FancyBboxPatch(xy,w,h,boxstyle="round,pad=0.02,rounding_size=0.025",facecolor=color,edgecolor="#333333",linewidth=1.1);ax.add_patch(patch);ax.text(xy[0]+w/2,xy[1]+h/2,text,ha="center",va="center",fontsize=9)

def M1(out):
    fig,ax=plt.subplots(figsize=(7.2,4.5));ax.set(xlim=(0,1),ylim=(0,1));ax.axis("off");add_box(ax,(.04,.66),.22,.18,"Connected pore state\n$f_{\\mathrm{fine}}^{c}$, $\\bar r_p^{c}$","#DCECF6");add_box(ax,(.39,.70),.22,.15,"Renewal densification\n$\\dot\\rho$","#D8F0E5");add_box(ax,(.73,.70),.22,.15,"TJ-assisted\ndensification","#D8F0E5");add_box(ax,(.04,.24),.22,.18,"PR/de-sintering\n$W_{\\mathrm{PR}}$","#FBE3D5");add_box(ax,(.39,.24),.22,.18,"Persistent junction\n$X_J$ + multihit","#E8E0F3");add_box(ax,(.73,.24),.22,.18,"Migration / growth\n$\\dot G$","#F7E6EE")
    for a,b in (((.26,.75),(.39,.77)),((.61,.77),(.73,.77)),((.15,.66),(.15,.42)),((.26,.33),(.39,.33)),((.61,.33),(.73,.33)),((.50,.70),(.50,.42))):ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=12,color="#444"))
    ax.text(.50,.08,"Pore-filled TJs partition into relaxed accommodation and Class-A pore drag;\nClass-B multihit applies to the structurally constrained population.",ha="center",fontsize=9);ps.finish(fig,out)

def M2(out):
    items=[(r"$\rho$","solid fraction"),(r"$G$","mean grain size"),(r"$f_{\mathrm{fine}}^{c}$","removable connected pores"),(r"$\bar r_p^{\,c}$","connected mean pore radius"),(r"$W_{\mathrm{PR}}$","non-densifying work"),(r"$X_J$","persistent junction state"),(r"$\Lambda_{\mathrm{TJ}}/K_{\mathrm{TJ}}$","multihit activity ratio"),(r"$C_{\mathrm{TJ}}^{\mathrm{pore,constraint}}$","distinct TJ populations")];fig,axs=plt.subplots(2,4,figsize=(7.2,3.5),constrained_layout=True)
    for ax,(symbol,text) in zip(axs.flat,items):ax.axis("off");ax.text(.5,.62,symbol,ha="center",fontsize=17,color="#2F4B7C");ax.text(.5,.30,text,ha="center",fontsize=8,wrap=True);ax.add_patch(FancyBboxPatch((.04,.08),.92,.82,transform=ax.transAxes,boxstyle="round,pad=.02",fill=False,edgecolor="#BBBBBB"))
    ps.panel_labels(axs);ps.finish(fig,out)

def M3(out):
    modes=[("current_all_TJ_multihit",(.0,.0,1.0)),("pore_relaxed_TJ",(.6,.0,.4)),("pore_pinned_drag",(.0,.8,.2)),("mixed_relaxed_pinned",(.6,.4,.0))];fig,axs=plt.subplots(1,4,figsize=(7.2,2.6),sharey=True,constrained_layout=True)
    for ax,(mode,(rel,pin,rem)) in zip(axs,modes):ax.bar([0],[rel],color="#72B7B2",label="relaxed");ax.bar([0],[pin],bottom=rel,color="#F58518",label="pinned");ax.bar([0],[rem],bottom=rel+pin,color="#4C78A8",label="constraint");ax.set(ylim=(0,1),xticks=[],title=mode.replace("_","\n"));ps.clean(ax)
    axs[0].set_ylabel(r"Partition of $C_{\mathrm{TJ}}^{\mathrm{pore}}$ [-]");axs[-1].legend(loc="center left",bbox_to_anchor=(1.02,.5));ps.panel_labels(axs);ps.finish(fig,out)

def P1(d,out):
    sf=d["representative_slow_fast_histories"];ts=d["representative_two_step_histories"];fig,axs=plt.subplots(1,2,figsize=(7.2,3),constrained_layout=True)
    for path,col in (("slow",ps.COLORS["slow"]),("fast",ps.COLORS["fast"])):
        q=sf[sf.path==path];axs[0].plot(q.t_s/3600,q.T_C,color=col,label=f"{path} heating")
    axs[1].plot([0,96],[1400,1400],color=ps.COLORS["highT"],label="high-$T$ isothermal");axs[1].plot(ts.t_s/3600,ts.T_C,color=ps.COLORS["two_step"],label="two-step")
    for ax in axs:ax.set(xlabel="Time [h]",ylabel="Temperature [°C]");ax.legend();ps.clean(ax)
    ps.panel_labels(axs);ps.finish(fig,out)

def protocol_frames(d):
    sf=d["representative_slow_fast_histories"].copy();ts=d["representative_two_step_histories"].copy();ht=highT_history();return {"slow heating":sf[sf.path=="slow"],"fast heating":sf[sf.path=="fast"],"high-$T$ isothermal":ht,"two-step":ts}

def P2(frames,out):
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.6),constrained_layout=True);fields=(("rho",ps.LABELS["rho"]),("G_nm",ps.LABELS["G"]),("connected_mean_radius_nm",ps.LABELS["radius"]),("connected_fine_pore_fraction",ps.LABELS["fine"]));cols=[ps.COLORS[k] for k in ("slow","fast","highT","two_step")]
    for ax,(field,label) in zip(axs.flat,fields):
        for (name,q),col in zip(frames.items(),cols):ax.plot(q.t_s/3600,q[field],label=name,color=col)
        ax.set(xlabel="Time [h]",ylabel=label);ps.clean(ax)
    axs[0,0].legend(fontsize=8);ps.panel_labels(axs);ps.finish(fig,out)

def P3(frames,out):
    fig,ax=plt.subplots(figsize=(4.3,3.6));
    for (name,q),col in zip(frames.items(),[ps.COLORS[k] for k in ("slow","fast","highT","two_step")]):ax.plot(q.G_nm,q.rho,label=name,color=col)
    ax.set(xlabel=ps.LABELS["G"],ylabel=ps.LABELS["rho"]);ax.legend();ps.clean(ax);ps.finish(fig,out)

def P4(frames,out):
    fig,axs=plt.subplots(1,2,figsize=(7.2,3),constrained_layout=True)
    for (name,q),col in zip(frames.items(),[ps.COLORS[k] for k in ("slow","fast","highT","two_step")]):axs[0].plot(q.rho,q.connected_mean_radius_nm,label=name,color=col);axs[1].plot(q.rho,q.large_pore_fraction,label=name,color=col)
    axs[0].set(xlabel=ps.LABELS["rho"],ylabel=ps.LABELS["radius"]);axs[1].set(xlabel=ps.LABELS["rho"],ylabel="Large-pore fraction [-]");axs[0].legend(fontsize=8)
    for ax in axs:ps.clean(ax)
    ps.panel_labels(axs);ps.finish(fig,out)

def P5(d,out):
    q=d["production_fast_firing_summary"];q=q[(q.candidate_id=="mech_009_PR_attrition_moderate")&(q.rho_target==.90)];tops=sorted(q.initial_topology.unique());fig,axs=plt.subplots(2,2,figsize=(7.2,5.3),sharex=True,sharey=True,constrained_layout=True)
    for ax,top in zip(axs.flat,tops):z=q[q.initial_topology==top].pivot_table(index="peak_T_C",columns="heating_rate_C_min",values="HR_pct_median",aggfunc="median");im=ax.imshow(z.values,origin="lower",aspect="auto",cmap="RdBu_r",vmin=-8,vmax=8,extent=(np.log10(z.columns.min()),np.log10(z.columns.max()),z.index.min()-25,z.index.max()+25));ax.set_title(top.replace("_"," "));ax.set_xticks(np.log10(z.columns),[f"{x:g}" for x in z.columns]);ax.set(xlabel="Heating rate [°C min$^{-1}$]",ylabel="Peak temperature [°C]")
    fig.colorbar(im,ax=axs,label=r"$HR$ improvement [%]",shrink=.8);ps.panel_labels(axs);ps.finish(fig,out)

def P6(d,out):
    q=d["successful_practical_windows"];q=q[(q.candidate_id=="mech_009")&(q.prep_growth_tolerance==.05)&(q.second_step_growth_tolerance==.05)];q=q.sort_values("G1_nm");fig,ax=plt.subplots(figsize=(5.2,4));ax.vlines(q.G1_nm,q.T_lower_density_C,q.T_upper_growth_C,color="#B9D9CF",alpha=.12);sc=ax.scatter(q.G1_nm,q.T_first_success_C,c=q.window_width_C,cmap="viridis",s=12);ax.scatter(q.G1_nm,q.T_lower_density_C,s=4,color="#0072B2",label="lower density boundary");ax.scatter(q.G1_nm,q.T_upper_growth_C,s=4,color="#D55E00",label="upper growth boundary");ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel="Second-step temperature, $T_2$ [°C]");fig.colorbar(sc,ax=ax,label="Window width [°C]");ax.legend(fontsize=7);ps.clean(ax);ps.finish(fig,out)

def P7(d,out):
    q=d["production_joint_scorecard"].iloc[:4].copy();vals=np.column_stack([pd.to_numeric(q.complete_practical_5_5)>0,pd.to_numeric(q.beneficial_fast_count)>0,q.joint_positive.astype(bool)]);fig,ax=plt.subplots(figsize=(5.5,2.8));ax.imshow(vals,cmap=plt.matplotlib.colors.ListedColormap(["#E6E6E6",ps.COLORS["positive"]]),vmin=0,vmax=1,aspect="auto");ax.set_xticks(range(3),["Practical Chen","Fast-firing benefit","Joint positive"]);ax.set_yticks(range(4),q.candidate_id);[ax.text(j,i,"✓" if vals[i,j] else "—",ha="center",va="center",color="white" if vals[i,j] else "#555",fontweight="bold") for i in range(4) for j in range(3)];ps.finish(fig,out)

def P8(d,out):
    q=d["representative_slow_fast_histories"];fig,axs=plt.subplots(2,2,figsize=(7.2,5.5),constrained_layout=True);fields=(("cumulative_PR_desintering_work",ps.LABELS["WPR"]),("connected_fine_pore_fraction",ps.LABELS["fine"]),("connected_mean_radius_nm",ps.LABELS["radius"]),("large_pore_fraction","Large-pore fraction [-]"))
    for ax,(f,l) in zip(axs.flat,fields):
        for path,col in (("slow",ps.COLORS["slow"]),("fast",ps.COLORS["fast"])):z=q[q.path==path];ax.plot(z.rho,z[f],color=col,label=f"{path} heating")
        ax.set(xlabel=ps.LABELS["rho"],ylabel=l);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);ps.finish(fig,out)

def P9(d,out):
    q=d["representative_two_step_histories"];fig,axs=plt.subplots(2,2,figsize=(7.2,5.4),constrained_layout=True);fields=(("X_J",ps.LABELS["XJ"]),("Lambda_over_K_TJ",ps.LABELS["lambdaK"]),("P_comp_TJ",r"Multihit completion, $P_{\mathrm{comp}}^{\mathrm{TJ}}$"),("P_TJ_multihit",r"TJ multihit dissipation [model units]"))
    for ax,(f,l) in zip(axs.flat,fields):
        for phase,ls in (("first_step","-"),("second_step","--")):z=q[q.phase==phase];ax.plot(z.rho,z[f],ls=ls,color=ps.COLORS["two_step"],label=phase.replace("_"," "))
        ax.set(xlabel=ps.LABELS["rho"],ylabel=l);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);ps.finish(fig,out)

def P10(d,out):
    q=d["tj_constraint_mode_summary"];g=q.groupby("TJ_constraint_mode").agg(chen=("complete_practical_chen","sum"),fast=("n_fast_beneficial","sum"),joint=("joint_positive","sum")).reindex(ps.TJ_COLORS);fig,ax=plt.subplots(figsize=(7.2,3.5));x=np.arange(len(g));ax.bar(x-.25,g.chen,.25,label="Chen-positive bases");ax.bar(x,g.joint,.25,label="joint-positive bases");ax.bar(x+.25,g.fast/20,.25,label="beneficial fast cases / 20");ax.set_xticks(x,[m.replace("_","\n") for m in g.index]);ax.set_ylabel("Count [scaled where noted]");ax.legend();ps.clean(ax);ps.finish(fig,out)

def P11(d,out):
    ff=d["production_fast_firing_summary"];tj=d["tj_constraint_mode_fast_firing"];ch=d["successful_practical_windows"];fig,axs=plt.subplots(2,2,figsize=(7.2,5.5),constrained_layout=True);v=ff[ff.response_class!="unattainable"];axs[0,0].scatter(v.PR_difference_median,v.HR_pct_median,s=7,alpha=.25);axs[0,1].scatter(v.fine_difference_median,v.HR_pct_median,s=7,alpha=.25);a=tj[tj.attained==True];axs[1,0].scatter(a.Lambda_over_K_TJ,a.HR_pct,c=a.beneficial.map({True:ps.COLORS["positive"],False:ps.COLORS["negative"]}),s=8,alpha=.35);c=ch[(ch.prep_growth_tolerance==.05)&(ch.second_step_growth_tolerance==.05)];axs[1,1].scatter(c.G1_nm,c.window_width_C,s=7,alpha=.25)
    labels=((r"$\Delta W_{\mathrm{PR}}$ [model units]",r"$HR$ improvement [%]"),(r"Retained $\Delta f_{\mathrm{fine}}^c$",r"$HR$ improvement [%]"),(ps.LABELS["lambdaK"],r"$HR$ improvement [%]"),(r"Prepared grain size, $G_1$ [nm]","Window width [°C]"))
    for ax,(x,y) in zip(axs.flat,labels):ax.set(xlabel=x,ylabel=y);ps.clean(ax)
    ps.panel_labels(axs);ps.finish(fig,out)

def main():
    ps.apply_style();OUT.mkdir(parents=True,exist_ok=True);d=data();frames=protocol_frames(d);spec=[("M1","M1_mechanism_overview","Mechanism overview","generate_paper_figures.py","Separate densification, PR memory, migration, and TJ channels",M1), ("M2","M2_variable_definitions","Variable definitions","generate_paper_figures.py","Define manuscript observables and notation",M2),("M3","M3_TJ_mode_schematic","TJ population modes","tj_constraint_mode_summary.csv","Explain relaxed, pinned, and constrained TJ partitions",M3)]
    inventory=[]
    for fid,stem,title,source,purpose,fn in spec:fn(OUT/stem);inventory.append(ps.inventory_row(fid,stem,title,source,purpose,"main text: model formulation"))
    funcs=[("P1","P1_thermal_histories","Thermal histories","representative histories","Define protocol families",lambda o:P1(d,o)),("P2","P2_evolution_histories","Evolution histories","representative histories","Compare density, grain, and pore evolution",lambda o:P2(frames,o)),("P3","P3_density_grain_trajectories","Density–grain trajectories","representative histories","Show sintering efficiency directly",lambda o:P3(frames,o)),("P4","P4_density_pore_trajectories","Density–pore trajectories","representative histories","Show pore-distribution pathway memory",lambda o:P4(frames,o)),("P5","P5_fast_firing_map","Fast-firing response map","production_fast_firing_summary.csv","Locate beneficial, neutral, and harmful response",lambda o:P5(d,o)),("P6","P6_chen_window_map","Practical Chen window","successful_practical_windows.csv","Show finite lower and upper boundaries",lambda o:P6(d,o)),("P7","P7_joint_scorecard","Joint scorecard","production_joint_scorecard.csv","Summarize joint-positive frozen bases",lambda o:P7(d,o)),("P8","P8_PR_memory","PR-memory comparison","representative_slow_fast_histories.csv","Connect fast firing to observable pore memory",lambda o:P8(d,o)),("P9","P9_two_step_migration","Two-step migration suppression","representative_two_step_histories.csv","Explain persistent-junction and multihit histories",lambda o:P9(d,o)),("P10","P10_TJ_ablation","TJ-mode ablation","tj_constraint_mode_summary.csv","Test pore relaxation, pinning, and constraint",lambda o:P10(d,o)),("P11","P11_success_correlations","Success correlations","compact production tables","Relate response to PR, topology, and TJ activity",lambda o:P11(d,o))]
    for fid,stem,title,source,purpose,fn in funcs:fn(OUT/stem);inventory.append(ps.inventory_row(fid,stem,title,source,purpose,"main text: results"))
    ps.write_inventory(OUT/"figure_inventory.csv",inventory);print("generated",len(inventory),"main figures")

if __name__=="__main__":main()
