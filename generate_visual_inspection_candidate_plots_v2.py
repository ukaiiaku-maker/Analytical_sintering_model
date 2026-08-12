#!/usr/bin/env python3
"""Dense, auditable visual-inspection rebuild; does not alter model physics."""
from pathlib import Path
from dataclasses import replace
import csv,json,math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import plot_style as ps
import dynamic_chen_topology_search as chen
import nucleation_fast_chen_production as nuc
import production_mechanism_assessment as protocols
import separated_fast_chen_model as model
import separated_mechanism_production_search as prod

OUT=Path("results/visual_inspection_candidate_plots_v2");FIG=OUT/"figures";TAB=OUT/"tables";HIST=OUT/"histories"
COL={1:"#0072B2",20:"#E69F00",50:"#009E73",100:"#D55E00"};STY={"lower_failure":"--","success":"-","upper_failure":"-."}

def write(path,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def select():
    d=pd.read_csv("results/strict_chen_window_production/strict_window_recheck.csv");eb=d[(d.material_id=="E0142")&(d.tier=="Tier_B")].sort_values(["window_width_C","G1_nm"],ascending=[False,True]).drop_duplicates(["topology_id","G0_nm","T1_C","rho_switch"]);ec=d[(d.material_id=="E0021")&(d.tier=="Tier_C")].sort_values(["window_width_C","G1_nm"],ascending=[False,True]).head(1);q=pd.concat([eb,ec],ignore_index=True);q["heating_rate_T1_C_min"]=20.;write(TAB/"selected_candidates_for_visualization.csv",q.to_dict("records"));return q

def dense(h,n=1200,offset=0,stage="single",extra=None):
    if not len(h["t"]):return []
    idx=np.unique(np.r_[np.arange(len(h["t"])),0,len(h["t"])-1]);t=h["t"][idx];grid=np.unique(np.r_[t,np.linspace(t.min(),t.max(),max(n,len(t)))])
    keys=("T_C","rho","G","tau_nuc","tau_exchange","tau_transport","activity","rho_dot","G_dot","PR_exposure","pore_D50","pore_D90","connected_fine","large_pore_fraction","X_J","Lambda_over_K","pore_drag")
    rows=[]
    for x in grid:
      r=dict(physical_time_s=offset+x,physical_time_h=(offset+x)/3600,local_time_s=x,stage=stage)
      for k in keys:
       if k in h and len(h[k]):r[k]=float(np.interp(x,h["t"],h[k]))
      if "G" in r:r["G_nm"]=r.pop("G")*1e9
      if "G_dot" in r:r["G_dot_nm_s"]=r.pop("G_dot")*1e9
      if "pore_D50" in r:r["pore_D50_nm"]=r.pop("pore_D50")*1e9
      if "pore_D90" in r:r["pore_D90_nm"]=r.pop("pore_D90")*1e9
      r.update(extra or {});rows.append(r)
    return rows

def fast_histories(mats,sel):
    rows=[];labels=[];rob=[]
    schedules={"E0021":(1550,20),"E0142":(1550,20)}
    for mid,mat in mats.items():
      peak,hold=schedules[mid]
      for rate in (1,20,50,100):
        h=model.run(mat,model.TopologyGrowthClosure(),protocols.FastSchedule(rate,peak,hold));rows+=dense(h,extra=dict(material_id=mid,rate_C_min=rate,protocol_label=f"{mid} | full kinetics | topology disabled | {rate} C/min",ablation="full"))
      for mode in ("no_PR_redistribution","no_nucleation_limitation","no_growth_before_activation","transport_only","exchange_limited_variant"):
        p=replace(mat,ablation_mode=mode,growth_activity_threshold=1e-5);h=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(100,peak,hold));rows+=dense(h,extra=dict(material_id=mid,rate_C_min=100,protocol_label=f"{mid} | {mode} | topology disabled | 100 C/min",ablation=mode));labels.append(dict(material_id=mid,ablation=mode,PR_redistribution_enabled=mode!="no_PR_redistribution",nucleation_limitation_enabled=mode!="no_nucleation_limitation",topology_enabled=False,growth_before_activation_enabled=mode!="no_growth_before_activation",peak_T_C=peak,hold_time_h=hold,reference_rate_C_min=1,fast_rate_C_min=100))
      for rho0 in (.60,.65,.70,.75,.80):
       for G0 in (50,75,100,150,225,300):
        p=replace(mat,rho0=rho0,G0=G0*1e-9);ref=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(1,peak,hold))
        for rate in (20,50,100):
          fast=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(rate,peak,hold));m=prod.metrics(ref,fast);m.pop("curve",None);rob.append(dict(material_id=mid,rho0=rho0,G0_nm=G0,peak_T_C=peak,hold_time_h=hold,fast_rate_C_min=rate,reference_rate_C_min=1,both_paths_attained=m.get("attained",False),max_ratio=m.get("max_ratio",math.nan),span_ratio_ge_1p2=m.get("span_ge_1p2",0),span_ratio_ge_1p5=m.get("span_ge_1p5",0),span_ratio_ge_2p0=m.get("span_ge_2p0",0),density_interval_start=m.get("effect_rho_min",math.nan),density_interval_end=m.get("effect_rho_max",math.nan),trajectory_class="meaningful" if m.get("meaningful") else "rejected",rejection_reason=m.get("rejection_reason","")))
    write(HIST/"dense_fast_histories.csv",rows);write(HIST/"dense_ablation_histories.csv",[r for r in rows if r["ablation"]!="full"]);write(TAB/"dense_fast_histories_index.csv",pd.DataFrame(rows).groupby(["material_id","rate_C_min","ablation"]).size().reset_index(name="n_points").to_dict("records"));write(TAB/"fast_firing_ablation_labels.csv",labels);write(TAB/"fast_firing_initial_condition_robustness.csv",rob);return pd.DataFrame(rows),pd.DataFrame(rob)

def two_step(mats,sel):
    r=sel[sel.material_id=="E0142"].sort_values(["window_width_C","G1_nm"],ascending=[False,True]).iloc[0];tops={x[0]:x for x in chen.design(128)};tid,fam,tp=tops[r.topology_id];p,h1,att=chen.first_state(mats["E0142"],tp,r.G0_nm,r.T1_C,r.rho_switch);state=h1["final_state"];end=float(h1["t"][-1]);Ts=(int(r.T_lower_density_C),int((r.T2_first_success_C+r.T2_last_success_C)/2),int(r.T_upper_growth_C));roles=("lower_failure","success","upper_failure");rows=[];labels=[]
    first=dense(h1,1200,0,"first_step",dict(material_id="E0142",topology_id=tid,T1_C=r.T1_C,rho_switch=r.rho_switch))
    for role,T2 in zip(roles,Ts):
      h2=model.run(p,tp,protocols.FastSchedule(1,T2,96),initial=state);combined=[{**x,"role":role,"T2_C":T2} for x in first]+dense(h2,1200,end,"second_step",dict(material_id="E0142",topology_id=tid,T1_C=r.T1_C,T2_C=T2,rho_switch=r.rho_switch,role=role));rows+=combined;labels.append(dict(material_id="E0142",topology_id=tid,q_TJ=tp.q_TJ,topology_mode=fam,T1_C=r.T1_C,T2_C=T2,rho_switch=r.rho_switch,role=role))
    write(HIST/"dense_two_step_histories.csv",rows);write(HIST/"dense_chen_triplet_histories.csv",rows);write(TAB/"dense_two_step_histories_index.csv",pd.DataFrame(rows).groupby(["role","stage"]).size().reset_index(name="n_points").to_dict("records"));write(TAB/"two_step_ablation_labels.csv",labels);write(TAB/"E0142_TierB_triplet_selection.csv",labels);return pd.DataFrame(rows),r

def save(fig,stem,inventory,category,candidate,sources,purpose):
    for ax in fig.axes:
      if ax.get_visible() and ax.has_data():
       if not ax.get_xlabel():ax.set_xlabel("Independent variable")
       if not ax.get_ylabel():ax.set_ylabel("Response")
    ps.finish(fig,FIG/stem);inventory.append(dict(figure_id=f"F{len(inventory)+1:02d}",filename_pdf=f"{stem}.pdf",filename_png=f"{stem}.png",category=category,candidate_id=candidate,source_tables=sources,rerun_required=True,purpose=purpose,status="generated"))

def fast_plots(fast,rob,inv):
    for mid in ("E0021","E0142"):
      q=fast[(fast.material_id==mid)&(fast.ablation=="full")];fig,axs=plt.subplots(2,3,figsize=(12,7))
      for rate,g in q.groupby("rate_C_min"):
       lab=f"{rate:g} °C/min";axs[0,0].plot(g.physical_time_h,g.T_C,label=lab,color=COL[rate]);axs[0,1].plot(g.physical_time_h,g.rho,label=lab,color=COL[rate]);axs[0,2].plot(g.physical_time_h,g.G_nm,label=lab,color=COL[rate]);axs[1,0].plot(g.rho,g.G_nm,label=lab,color=COL[rate]);axs[1,2].plot(g.rho,g.activity,label=lab,color=COL[rate])
      ref=q[q.rate_C_min==1]
      for rate in (20,50,100):
       g=q[q.rate_C_min==rate];lo=max(ref.rho.min(),g.rho.min());hi=min(ref.rho.max(),g.rho.max());x=np.linspace(lo,hi,300);axs[1,1].plot(x,np.interp(x,ref.rho,ref.G_nm)/np.interp(x,g.rho,g.G_nm),label=f"1/{rate}",color=COL[rate])
      for y in (1.2,1.5,2):axs[1,1].axhline(y,color="#777",ls="--",lw=.8)
      labs=(("Time, $t$ [h]","Temperature, $T$ [°C]"),("Time, $t$ [h]","Relative density, $\\rho$"),("Time, $t$ [h]","Grain size, $G$ [nm]"),("Relative density, $\\rho$","Grain size, $G$ [nm]"),("Relative density, $\\rho$",r"$G_{ref}/G_{fast}$"),("Relative density, $\\rho$","Activity, $a$"))
      for ax,(x,y) in zip(axs.flat,labs):ax.set(xlabel=x,ylabel=y,title=y);ax.legend(fontsize=7);ps.clean(ax)
      ps.panel_labels(axs);save(fig,f"fast_firing/{mid}_time_histories",inv,"fast_firing",mid,"dense_fast_histories.csv","Dense fast-firing histories")
      fig,axs=plt.subplots(2,3,figsize=(12,7));
      for rate,g in q.groupby("rate_C_min"):
       for ax,k in zip(axs.flat,("tau_nuc","tau_exchange","tau_transport","activity","rho_dot","G_dot_nm_s")):ax.plot(g.T_C,g[k],label=f"{rate:g} °C/min",color=COL[rate])
      yl=(r"$\tau_{nuc}$ [s]",r"$\tau_{ex}$ [s]",r"$\tau_{tr}$ [s]","Activity, $a$",r"$\dot\rho$ [s$^{-1}$]",r"$\dot G$ [nm s$^{-1}$]")
      for i,(ax,y) in enumerate(zip(axs.flat,yl)):ax.set(xlabel="Temperature, $T$ [°C]",ylabel=y,title=y);ax.legend(fontsize=6);ps.clean(ax);ax.set_yscale("log" if i<3 else "linear")
      save(fig,f"fast_firing/{mid}_kinetic_times",inv,"fast_firing",mid,"dense_fast_histories.csv","Serial kinetics")
      a=fast[(fast.material_id==mid)&(fast.rate_C_min==100)];fig,axs=plt.subplots(2,2,figsize=(9,7))
      for mode,g in a.groupby("ablation"):
       lab=f"{mid} | {mode} | 100 °C/min";axs[0,0].plot(g.rho,g.G_nm,label=lab);axs[1,0].plot(g.rho,g.PR_exposure,label=lab);axs[1,1].plot(g.rho,g.activity,label=lab)
      ref=q[q.rate_C_min==1];full=a[a.ablation=="full"];x=np.linspace(max(ref.rho.min(),full.rho.min()),min(ref.rho.max(),full.rho.max()),300);axs[0,1].plot(x,np.interp(x,ref.rho,ref.G_nm)/np.interp(x,full.rho,full.G_nm),label="full: 1/100");[axs[0,1].axhline(y,color="#777",ls="--") for y in (1.2,1.5,2)]
      for ax,xl,yl in zip(axs.flat,("Relative density, $\\rho$",)*4,("Grain size, $G$ [nm]",r"$G_{ref}/G_{fast}$","PR exposure","Activity, $a$")):ax.set(xlabel=xl,ylabel=yl,title=yl);ax.legend(fontsize=5);ps.clean(ax)
      save(fig,f"fast_firing/{mid}_ablation_comparison",inv,"ablation",mid,"dense_ablation_histories.csv","Fully labeled fast ablations")
    # Robustness heat maps.
    fig,axs=plt.subplots(2,3,figsize=(12,7))
    for ax,((mid,rate),g) in zip(axs.flat,rob.groupby(["material_id","fast_rate_C_min"])):
      p=g.pivot(index="rho0",columns="G0_nm",values="span_ratio_ge_1p5");im=ax.imshow(p,origin="lower",aspect="auto",vmin=0,vmax=.17,cmap="viridis");ax.set(xticks=range(len(p.columns)),xticklabels=p.columns,yticks=range(len(p.index)),yticklabels=p.index,xlabel="$G_0$ [nm]",ylabel="$\\rho_0$",title=f"{mid}, {rate:g} °C/min")
    fig.colorbar(im,ax=axs.ravel().tolist(),label=r"Span with ratio $\geq1.5$");save(fig,"robustness/fast_firing_ratio_heatmap_by_rho0_G0",inv,"robustness","both","fast_firing_initial_condition_robustness.csv","Initial-state robustness")
    fig,axs=plt.subplots(2,3,figsize=(12,7))
    for ax,((mid,rate),g) in zip(axs.flat,rob.groupby(["material_id","fast_rate_C_min"])):
      p=g.pivot(index="rho0",columns="G0_nm",values="max_ratio");im=ax.imshow(p,origin="lower",aspect="auto",vmin=1,vmax=min(5,np.nanmax(p)),cmap="magma");ax.set(xticks=range(len(p.columns)),xticklabels=p.columns,yticks=range(len(p.index)),yticklabels=p.index,xlabel="$G_0$ [nm]",ylabel="$\\rho_0$",title=f"{mid}, {rate:g} °C/min")
    fig.colorbar(im,ax=axs.ravel().tolist(),label="Maximum grain ratio");save(fig,"robustness/fast_firing_max_ratio_heatmap_by_rho0_G0",inv,"robustness","both","fast_firing_initial_condition_robustness.csv","Maximum ratio robustness")

def two_step_plots(ts,r,inv):
    fig,axs=plt.subplots(2,3,figsize=(12,7))
    for role,g in ts.groupby("role"):
      g=g.sort_values("physical_time_h");ls=STY[role];lab=role.replace("_"," ");axs[0,0].plot(g.physical_time_h,g.T_C,ls=ls,label=lab);axs[0,1].plot(g.physical_time_h,g.rho,ls=ls,label=lab);axs[0,2].plot(g.physical_time_h,g.G_nm,ls=ls,label=lab);axs[1,0].plot(g.rho,g.G_nm,ls=ls,label=lab);axs[1,1].plot(g.physical_time_h,g.rho_dot,ls=ls,label=lab);axs[1,2].plot(g.physical_time_h,g.G_dot_nm_s,ls=ls,label=lab)
      sw=g[g.stage=="second_step"].physical_time_h.min()
      for ax in axs.flat:ax.axvline(sw,color="#777",lw=.7,alpha=.5)
    labs=(("Time, $t$ [h]","Temperature, $T$ [°C]"),("Time, $t$ [h]","Relative density, $\\rho$"),("Time, $t$ [h]","Grain size, $G$ [nm]"),("Relative density, $\\rho$","Grain size, $G$ [nm]"),("Time, $t$ [h]",r"$\dot\rho$ [s$^{-1}$]"),("Time, $t$ [h]",r"$\dot G$ [nm s$^{-1}$]"))
    for ax,(x,y) in zip(axs.flat,labs):ax.set(xlabel=x,ylabel=y,title=y);ax.legend(fontsize=7);ps.clean(ax)
    ps.panel_labels(axs);save(fig,"two_step/E0142_TierB_two_step_time_histories",inv,"two_step","E0142","dense_two_step_histories.csv","Continuous physical time triplet")
    fig,axs=plt.subplots(2,2,figsize=(9,7))
    for role,g in ts.groupby("role"):
      for ax,k in zip(axs.flat,("X_J","Lambda_over_K","pore_drag","connected_fine")):ax.plot(g.physical_time_h,g[k],ls=STY[role],label=role)
    for ax,y in zip(axs.flat,(r"$X_J$",r"$\Lambda_{TJ}/K_{TJ}$","Pore drag",r"$f_{fine}^c$")):ax.set(xlabel="Time, $t$ [h]",ylabel=y,title=y);ax.legend(fontsize=7);ps.clean(ax)
    save(fig,"two_step/E0142_TierB_two_step_topology_histories",inv,"two_step","E0142","dense_two_step_histories.csv","Available topology histories")
    fig,axs=plt.subplots(1,3,figsize=(11,3.5))
    for role,g in ts.groupby("role"):
      for ax,k in zip(axs,("pore_D50_nm","pore_D90_nm","connected_fine")):ax.plot(g.rho,g[k],ls=STY[role],label=role)
    for ax,y in zip(axs,(r"$D_{50}$ [nm]",r"$D_{90}$ [nm]",r"$f_{fine}^c$")):ax.set(xlabel="Relative density, $\\rho$",ylabel=y,title=y);ax.legend(fontsize=7);ps.clean(ax)
    save(fig,"two_step/E0142_TierB_two_step_microstructure_histories",inv,"two_step","E0142","dense_two_step_histories.csv","Available pore histories")

def chen_plots(sel,inv):
    raw=pd.read_csv("results/dynamic_chen_topology_for_nucleation_materials/dynamic_topology_classification_points.csv");b=sel[(sel.material_id=="E0142")&(sel.tier=="Tier_B")].drop_duplicates(["topology_id","G0_nm","T1_C","rho_switch"]);ids=set(b.topology_id);q=raw[(raw.material_id=="E0142")&raw.topology_id.isin(ids)];write(TAB/"E0142_TierB_Chen_window_boundaries_v2.csv",b.to_dict("records"));write(TAB/"E0142_TierB_Chen_classification_points_v2.csv",q.to_dict("records"));write(TAB/"E0021_TierC_failure_summary_v2.csv",sel[sel.material_id=="E0021"].to_dict("records"))
    colors={"SUCCESS":"#009E73","DENSIFICATION_EXHAUSTION_FAILURE":"#0072B2","GRAIN_GROWTH_FAILURE":"#D55E00","MIXED_FAILURE":"#CC79A7"};fig,ax=plt.subplots(figsize=(7,5))
    for c,col in colors.items():g=q[q.classification==c];ax.scatter(g.G1_nm,g.T2_C,c=col,label=c,s=22,alpha=.7)
    ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel=r"Second-step temperature, $T_2$ [°C]",title="E0142 Tier B dynamic classifications");ax.legend(fontsize=7);ps.clean(ax);save(fig,"chen_maps/E0142_TierB_Chen_classification_map_v2",inv,"chen_map","E0142","E0142_TierB_Chen_classification_points_v2.csv","Distinct failure classes")
    fig,ax=plt.subplots(figsize=(7,5))
    for _,x in b.iterrows():ax.vlines(x.G1_nm,x.T2_first_success_C,x.T2_last_success_C,color="#009E73",lw=9,alpha=.5);ax.scatter(x.G1_nm,x.T_lower_density_C,c="#0072B2",marker="v");ax.scatter(x.G1_nm,x.T_upper_growth_C,c="#D55E00",marker="^");ax.annotate(f"{x.window_width_C:.0f}°",(x.G1_nm,x.T2_last_success_C),fontsize=7)
    ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel=r"Second-step temperature, $T_2$ [°C]",title="E0142 Tier B finite success bands");ax.legend(handles=[mlines.Line2D([],[],color="#009E73",lw=7,label="Tier B success band"),mlines.Line2D([],[],color="#0072B2",marker="v",ls="",label="density boundary"),mlines.Line2D([],[],color="#D55E00",marker="^",ls="",label="growth boundary")]);ps.clean(ax);save(fig,"chen_maps/E0142_TierB_Chen_filled_window_map_v2",inv,"chen_map","E0142","E0142_TierB_Chen_window_boundaries_v2.csv","Filled finite windows only")
    fig,ax=plt.subplots(figsize=(6,4));
    for qv,g in b.groupby("q_TJ"):ax.scatter(g.G1_nm,g.window_width_C,marker="o" if qv==0 else "s",label=f"qTJ={qv}",s=55)
    best=b.sort_values("window_width_C",ascending=False).iloc[0];ax.annotate("best",(best.G1_nm,best.window_width_C));ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel="Window width [°C]",title="Tier B window width");ax.legend();ps.clean(ax);save(fig,"chen_maps/E0142_TierB_window_width_map_v2",inv,"chen_map","E0142","E0142_TierB_Chen_window_boundaries_v2.csv","q0/q1 widths")
    fig,axs=plt.subplots(1,4,figsize=(13,3.5));
    for qv,g in b.groupby("q_TJ"):axs[0].bar(str(qv),len(g));axs[1].scatter([qv]*len(g),g.window_width_C,label=f"q={qv}");axs[2].scatter([qv]*len(g),g.G1_nm,label=f"q={qv}");axs[3].scatter([qv]*len(g),g.T2_first_success_C,label=f"q={qv}")
    for ax,y in zip(axs,("Count","Width [°C]",r"$G_1$ [nm]",r"First-success $T_2$ [°C]")):
        ax.set(xlabel=r"$q_{TJ}$",ylabel=y,title=y)
        if ax is not axs[0]:ax.legend(fontsize=7)
        ps.clean(ax)
    save(fig,"chen_maps/E0142_TierB_q0_q1_comparison_v2",inv,"chen_map","E0142","E0142_TierB_Chen_window_boundaries_v2.csv","q0 and q1 both Tier B")

def main():
    for p in (OUT,TAB,HIST,FIG/"fast_firing",FIG/"two_step",FIG/"chen_maps",FIG/"robustness",FIG/"ablations",FIG/"comparison"):p.mkdir(parents=True,exist_ok=True)
    ps.apply_style();sel=select();mats=chen.materials();fast,rob=fast_histories(mats,sel);ts,best=two_step(mats,sel);inv=[];fast_plots(fast,rob,inv);two_step_plots(ts,best,inv);chen_plots(sel,inv)
    # Comparison and ablation matrices reuse actual plotted dense histories.
    fig,axs=plt.subplots(1,2,figsize=(9,4));
    for mid,g in fast[(fast.ablation=="full")&(fast.rate_C_min.isin((1,100)))].groupby("material_id"):
      for rate,z in g.groupby("rate_C_min"):axs[0].plot(z.rho,z.G_nm,label=f"{mid}, {rate:g}°/min");axs[1].plot(z.rho,z.activity,label=f"{mid}, {rate:g}°/min")
    for ax,y in zip(axs,("Grain size, $G$ [nm]","Activity, $a$")):ax.set(xlabel="Relative density, $\\rho$",ylabel=y,title=y);ax.legend(fontsize=7);ps.clean(ax)
    save(fig,"comparison/E0021_vs_E0142_fast_firing_comparison",inv,"comparison","both","dense_fast_histories.csv","Material comparison")
    fig,axs=plt.subplots(2,6,figsize=(18,6),sharex=False)
    modes=("full","no_PR_redistribution","no_nucleation_limitation","no_growth_before_activation","transport_only","exchange_limited_variant")
    for row,mid in enumerate(("E0021","E0142")):
      for ax,mode in zip(axs[row],modes):
       g=fast[(fast.material_id==mid)&(fast.ablation==mode)];ax.plot(g.rho,g.G_nm,label=f"{mid} | {mode} | 100°C/min");ax.set(xlabel="$\\rho$",ylabel="$G$ [nm]",title=mode.replace("_"," "));ax.legend(fontsize=5);ps.clean(ax)
    save(fig,"ablations/fast_firing_ablation_matrix_v2",inv,"ablation","both","dense_ablation_histories.csv","Separated fast ablations")
    fig,axs=plt.subplots(1,3,figsize=(11,3.5));
    for role,g in ts.groupby("role"):axs[0].plot(g.rho,g.G_nm,ls=STY[role],label=role);axs[1].plot(g.physical_time_h,g.rho,ls=STY[role],label=role);axs[2].plot(g.physical_time_h,g.G_nm,ls=STY[role],label=role)
    for ax,x,y in zip(axs,("$\\rho$","Time, $t$ [h]","Time, $t$ [h]"),("$G$ [nm]","$\\rho$","$G$ [nm]")):ax.set(xlabel=x,ylabel=y,title=y);ax.legend(fontsize=7);ps.clean(ax)
    save(fig,"ablations/two_step_topology_ablation_matrix_v2",inv,"ablation","E0142","dense_two_step_histories.csv","Representative topology-controlled triplet")
    write(TAB/"candidate_visual_summary_v2.csv",sel.to_dict("records"));write(TAB/"figure_inventory.csv",inv);(OUT/"run_metadata.json").write_text(json.dumps(dict(n_figures=len(inv),model_files_changed=False),indent=2)+"\n")
if __name__=="__main__":main()
