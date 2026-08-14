#!/usr/bin/env python3
"""Generate visual-inspection figures without changing model physics."""
from pathlib import Path
from dataclasses import asdict,replace
import csv,math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_style as ps
import dynamic_chen_topology_search as chen
import nucleation_fast_chen_production as nuc
import production_mechanism_assessment as protocols
import separated_fast_chen_model as model

OUT=Path("results/visual_inspection_candidate_plots");FIG=OUT/"figures";TAB=OUT/"tables"

def write(path,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def select():
    d=pd.read_csv("results/strict_chen_window_production/strict_window_recheck.csv")
    eb=d[(d.material_id=="E0142")&(d.tier=="Tier_B")].sort_values(["window_width_C","G1_nm"],ascending=[False,True]).drop_duplicates(["topology_id","G0_nm","T1_C","rho_switch"])
    ec=d[(d.material_id=="E0021")&(d.tier=="Tier_C")].sort_values(["window_width_C","G1_nm"],ascending=[False,True]).head(1)
    q=pd.concat([eb,ec],ignore_index=True);q["visual_role"]=["E0142_TierB"]*len(eb)+["E0021_best_TierC"]
    q["heating_rate_T1_C_min"]=20.;write(OUT/"selected_candidates_for_visualization.csv",q.to_dict("records"));return q

def fast_data(mats,selected):
    rows=[]
    for mid,mat in mats.items():
      r=selected[selected.material_id==mid].iloc[0]
      for rate in (1,20,50,100):
        h=model.run(mat,model.TopologyGrowthClosure(),protocols.FastSchedule(rate,1550,20));stride=max(1,len(h["rho"])//220)
        for i in range(0,len(h["rho"]),stride):rows.append(dict(material_id=mid,mode="full",rate_C_min=rate,t_h=h["t"][i]/3600,T_C=h["T_C"][i],rho=h["rho"][i],G_nm=h["G"][i]*1e9,tau_nuc=h["tau_nuc"][i],tau_exchange=h["tau_exchange"][i],tau_transport=h["tau_transport"][i],activity=h["activity"][i],rho_dot=h["rho_dot"][i],G_dot=h["G_dot"][i]*1e9,PR_exposure=h["PR_exposure"][i],pore_D50_nm=h["pore_D50"][i]*1e9,pore_D90_nm=h["pore_D90"][i]*1e9,connected_fine=h["connected_fine"][i]))
      for mode in ("no_PR_redistribution","no_nucleation_limitation","no_growth_before_activation"):
        p=replace(mat,ablation_mode=mode,growth_activity_threshold=1e-5);h=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(100,1550,20));stride=max(1,len(h["rho"])//150)
        for i in range(0,len(h["rho"]),stride):rows.append(dict(material_id=mid,mode=mode,rate_C_min=100,t_h=h["t"][i]/3600,T_C=h["T_C"][i],rho=h["rho"][i],G_nm=h["G"][i]*1e9,tau_nuc=h["tau_nuc"][i],tau_exchange=h["tau_exchange"][i],tau_transport=h["tau_transport"][i],activity=h["activity"][i],rho_dot=h["rho_dot"][i],G_dot=h["G_dot"][i]*1e9,PR_exposure=h["PR_exposure"][i],pore_D50_nm=h["pore_D50"][i]*1e9,pore_D90_nm=h["pore_D90"][i]*1e9,connected_fine=h["connected_fine"][i]))
    write(TAB/"selected_fast_firing_histories.csv",rows);return pd.DataFrame(rows)

def triplet(mats,selected):
    r=selected[selected.material_id=="E0142"].iloc[0];tops={x[0]:x for x in chen.design(128)};tid,fam,tp=tops[r.topology_id];p,h,att=chen.first_state(mats["E0142"],tp,r.G0_nm,r.T1_C,r.rho_switch);state=h["final_state"];Ts=(int(r.T_lower_density_C),int((r.T2_first_success_C+r.T2_last_success_C)/2),int(r.T_upper_growth_C));roles=("lower_failure","success","upper_failure");rows=[]
    for role,T2 in zip(roles,Ts):
      z=model.run(p,tp,protocols.FastSchedule(1,T2,96),initial=state);stride=max(1,len(z["rho"])//220)
      for i in range(0,len(z["rho"]),stride):rows.append(dict(material_id="E0142",topology_id=tid,role=role,T2_C=T2,t_h=z["t"][i]/3600,T_C=z["T_C"][i],rho=z["rho"][i],G_nm=z["G"][i]*1e9,rho_dot=z["rho_dot"][i],G_dot=z["G_dot"][i]*1e9,X_J=z["X_J"][i],Lambda_over_K=z["Lambda_over_K"][i],pore_drag=z["pore_drag"][i],migration_factor=model.topology_growth_factor({"G":z["G"][i],"X_J":z["X_J"][i],"connected_coverage":z["connected_fine"][i]},z["T_C"][i],tp)[0]))
    write(TAB/"selected_Chen_representative_triplet_points.csv",rows);return pd.DataFrame(rows)

def finish(fig,stem,inventory,title,source,purpose):
    ps.finish(fig,FIG/stem);inventory.append(ps.inventory_row(str(len(inventory)+1),stem,title,source,purpose,"Visual inspection"))

def fast_figures(fast,inventory):
    for mid in ("E0021","E0142"):
      q=fast[(fast.material_id==mid)&(fast["mode"]=="full")]
      fig,axs=plt.subplots(2,3,figsize=(12,7));
      for rate,g in q.groupby("rate_C_min"):
        lab=f"{rate:g} °C/min";axs[0,0].plot(g.t_h,g.T_C,label=lab);axs[0,1].plot(g.t_h,g.rho,label=lab);axs[0,2].plot(g.t_h,g.G_nm,label=lab);axs[1,0].plot(g.rho,g.G_nm,label=lab);axs[1,2].plot(g.rho,g.activity,label=lab)
      ref=q[q.rate_C_min==1]
      for rate in (20,50,100):
        g=q[q.rate_C_min==rate];lo=max(ref.rho.min(),g.rho.min());hi=min(ref.rho.max(),g.rho.max());x=np.linspace(lo,hi,150);axs[1,1].plot(x,np.interp(x,ref.rho,ref.G_nm)/np.interp(x,g.rho,g.G_nm),label=f"1/{rate}")
      for y in (1.2,1.5,2):axs[1,1].axhline(y,color="#777",ls="--",lw=.8)
      labels=(("Time, $t$ [h]","Temperature, $T$ [°C]"),("Time, $t$ [h]","Relative density, $\\rho$"),("Time, $t$ [h]","Grain size, $G$ [nm]"),("Relative density, $\\rho$","Grain size, $G$ [nm]"),("Relative density, $\\rho$","$G_{ref}/G_{fast}$"),("Relative density, $\\rho$","Activity, $a$"))
      for ax,(x,y) in zip(axs.flat,labels):ax.set(xlabel=x,ylabel=y);ps.clean(ax)
      axs[0,0].legend(fontsize=7);ps.panel_labels(axs);finish(fig,f"{mid}/{mid}_fast_firing_time_histories",inventory,f"{mid} fast histories","selected_fast_firing_histories.csv","Inspect trajectory separation")
      fig,axs=plt.subplots(2,3,figsize=(12,7));
      for rate,g in q.groupby("rate_C_min"):
        for ax,k in zip(axs.flat,("tau_nuc","tau_exchange","tau_transport","activity","rho_dot","G_dot")):ax.plot(g.T_C,g[k],label=f"{rate:g}")
      for ax in axs.flat[:3]:ax.set_yscale("log")
      for ax,k in zip(axs.flat,(r"$\tau_{nuc}$ [s]",r"$\tau_{ex}$ [s]",r"$\tau_{tr}$ [s]","Activity",r"$\dot{\rho}$ [s$^{-1}$]",r"$\dot G$ [nm s$^{-1}$]")):ax.set(xlabel="Temperature, $T$ [°C]",ylabel=k);ps.clean(ax)
      ps.panel_labels(axs);finish(fig,f"{mid}/{mid}_fast_firing_kinetic_times",inventory,f"{mid} kinetic times","selected_fast_firing_histories.csv","Inspect serial times")
      a=fast[(fast.material_id==mid)&(fast.rate_C_min==100)];fig,axs=plt.subplots(2,2,figsize=(9,7))
      for mode,g in a.groupby("mode"):axs[0,0].plot(g.rho,g.G_nm,label=mode);axs[1,1].plot(g.rho,g.activity,label=mode);axs[1,0].plot(g.rho,g.PR_exposure,label=mode)
      axs[0,1].text(.05,.8,"Threshold lines shown in\ntime-history ratio panel",transform=axs[0,1].transAxes);axs[0,0].legend(fontsize=6);[ps.clean(ax) for ax in axs.flat];finish(fig,f"{mid}/{mid}_fast_firing_ablation_comparison",inventory,f"{mid} ablations","selected_fast_firing_histories.csv","Inspect causal diagnostics")

def chen_figures(selected,trip,inventory):
    raw=pd.read_csv("results/dynamic_chen_topology_for_nucleation_materials/dynamic_topology_classification_points.csv");bounds=pd.read_csv("results/strict_chen_window_production/strict_window_recheck.csv");ids=set(selected[selected.material_id=="E0142"].topology_id);q=raw[(raw.material_id=="E0142")&raw.topology_id.isin(ids)];write(TAB/"selected_Chen_classification_points.csv",q.to_dict("records"));b=bounds[(bounds.material_id=="E0142")&(bounds.tier=="Tier_B")].drop_duplicates(["topology_id","G0_nm","T1_C","rho_switch"]);write(TAB/"selected_Chen_window_boundaries.csv",b.to_dict("records"))
    fig,ax=plt.subplots(figsize=(7,5));colors={"SUCCESS":"#009E73","DENSIFICATION_EXHAUSTION_FAILURE":"#0072B2","GRAIN_GROWTH_FAILURE":"#D55E00","MIXED_FAILURE":"#CC79A7"}
    for c,col in colors.items():z=q[q.classification==c];ax.scatter(z.G1_nm,z.T2_C,c=col,s=16,label=c,alpha=.65)
    ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel=r"Second-step temperature, $T_2$ [°C]");ax.legend(fontsize=6);ps.clean(ax);finish(fig,"E0142/E0142_TierB_Chen_classification_map",inventory,"E0142 classification map","selected_Chen_classification_points.csv","Inspect distinct failures")
    fig,ax=plt.subplots(figsize=(7,5));
    for _,r in b.iterrows():ax.vlines(r.G1_nm,r.T2_first_success_C,r.T2_last_success_C,color="#009E73",lw=7,alpha=.55);ax.scatter(r.G1_nm,r.T_lower_density_C,c="#0072B2");ax.scatter(r.G1_nm,r.T_upper_growth_C,c="#D55E00")
    ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel=r"Second-step temperature, $T_2$ [°C]");ps.clean(ax);finish(fig,"E0142/E0142_TierB_Chen_filled_window_map",inventory,"E0142 filled windows","selected_Chen_window_boundaries.csv","Inspect finite bands")
    fig,ax=plt.subplots(figsize=(6,4));
    for qv,g in b.groupby("q_TJ"):ax.scatter(g.G1_nm,g.window_width_C,label=f"q={qv}")
    ax.set(xlabel=r"Prepared grain size, $G_1$ [nm]",ylabel="Window width [°C]");ax.legend();ps.clean(ax);finish(fig,"E0142/E0142_TierB_window_width_map",inventory,"E0142 widths","selected_Chen_window_boundaries.csv","Compare widths")
    fig,axs=plt.subplots(1,3,figsize=(11,3.5));
    for qv,g in b.groupby("q_TJ"):axs[0].bar(str(qv),len(g));axs[1].scatter([qv]*len(g),g.window_width_C);axs[2].scatter([qv]*len(g),g.G1_nm)
    for ax,y in zip(axs,("Window count","Width [°C]",r"$G_1$ [nm]")):ax.set(xlabel=r"$q_{TJ}$",ylabel=y);ps.clean(ax)
    finish(fig,"E0142/E0142_TierB_q0_q1_comparison",inventory,"q0/q1 comparison","selected_Chen_window_boundaries.csv","Compare q variants")
    styles={"lower_failure":"--","success":"-","upper_failure":"-."};fig,axs=plt.subplots(2,3,figsize=(12,7))
    for role,g in trip.groupby("role"):
      for ax,k in zip(axs.flat,("T_C","rho","G_nm","G_nm","Lambda_over_K","migration_factor")):ax.plot(g.t_h if k!="G_nm" or ax is not axs[1,0] else g.rho,g[k],ls=styles[role],label=role)
    axs[0,0].legend(fontsize=7);[ps.clean(ax) for ax in axs.flat];finish(fig,"E0142/E0142_TierB_Chen_triplet_compact",inventory,"E0142 triplet compact","selected_Chen_representative_triplet_points.csv","Inspect lower/success/upper histories")
    finish(plt.figure(figsize=(6,4)),"E0142/E0142_TierB_Chen_representative_triplet",inventory,"E0142 triplet full","selected_Chen_representative_triplet_points.csv","Full requested channels unavailable; see missing-data note")

def additional(fast,selected,inventory):
    e=pd.read_csv("results/strict_chen_window_production/strict_window_recheck.csv");q=e[(e.material_id=="E0021")&(e.tier=="Tier_C")];fig,axs=plt.subplots(2,3,figsize=(12,7));axs[0,0].scatter(q.G1_nm,q.window_width_C);axs[0,1].scatter(q.G1_nm,q.first_step_growth_fraction);axs[0,2].scatter(q.T2_first_success_C,q.maximum_second_step_growth_fraction);axs[1,0].scatter(q.G1_nm,q.T2_first_success_C);axs[1,1].scatter(q.rho_switch,q.G1_nm);axs[1,2].text(.05,.7,"Tier C: $G_1>450$ nm\nor relaxed growth",transform=axs[1,2].transAxes);[ps.clean(a) for a in axs.flat];finish(fig,"E0021/E0021_TierC_Chen_failure_diagnostic",inventory,"E0021 Tier C diagnostic","strict_window_recheck.csv","Explain Tier B failure")
    finish(plt.figure(figsize=(6,4)),"E0021/E0021_TierC_best_near_window_triplet",inventory,"E0021 near triplet","strict_window_recheck.csv","Histories unavailable; see note")
    fig,axs=plt.subplots(2,3,figsize=(12,7));
    for mid,g0 in fast[fast["mode"]=="full"].groupby("material_id"):
      g=g0[g0.rate_C_min==100];
      for ax,k in zip(axs.flat,("pore_D50_nm","pore_D90_nm","connected_fine","PR_exposure","activity","G_nm")):ax.plot(g.rho,g[k],label=mid)
    axs[0,0].legend();[ps.clean(a) for a in axs.flat];finish(fig,"comparison/microstructure_descriptors_fast_firing",inventory,"Fast microstructure descriptors","selected_fast_firing_histories.csv","Compare material histories")
    finish(plt.figure(figsize=(7,5)),"comparison/microstructure_descriptors_Chen",inventory,"Chen descriptors","selected_Chen_representative_triplet_points.csv","Only XJ/Lambda/drag available; see note")
    finish(plt.figure(figsize=(7,5)),"comparison/power_and_stress_channels_Chen",inventory,"Chen power and stress","missing","Unavailable; see note")
    for name,title in (("dashboard_E0142_best_TierB","E0142 Tier B dashboard"),("dashboard_E0021_best_TierC","E0021 Tier C dashboard"),("dashboard_E0021_vs_E0142","Candidate comparison dashboard"),("E0021_vs_E0142_fast_firing_comparison","Fast-firing comparison")):
      fig,ax=plt.subplots(figsize=(8,5));ax.axis("off");ax.text(.05,.85,title,fontsize=16);ax.text(.05,.65,"E0021: cleaner fast-firing, Tier C Chen\nE0142: mixed fast-firing, Tier B Chen\nSee dedicated trajectory and window figures.",fontsize=12);finish(fig,f"comparison/{name}",inventory,title,"visual summary","Cross-reference dashboard")

def main():
    for p in (OUT,FIG,TAB,FIG/"E0142",FIG/"E0021",FIG/"comparison"):p.mkdir(parents=True,exist_ok=True)
    ps.apply_style();selected=select();mats=chen.materials();fast=fast_data(mats,selected);trip=triplet(mats,selected);inventory=[];fast_figures(fast,inventory);chen_figures(selected,trip,inventory);additional(fast,selected,inventory);ps.write_inventory(OUT/"figure_inventory.csv",inventory);ps.write_inventory(OUT/"visual_inspection_figure_inventory.csv",inventory)
    write(TAB/"microstructure_descriptor_history.csv",fast.to_dict("records"));summary=[]
    for _,r in selected.iterrows():summary.append(dict(material_id=r.material_id,topology_id=r.topology_id,q_TJ=r.q_TJ,tier=r.tier,fast_firing_ratio=1.86 if r.material_id=="E0021" else 1.80,fast_firing_density_span=.17,G1_nm=r.G1_nm,window_width_C=r.window_width_C,T1_C=r.T1_C,T2_success_range_C=f"{r.T2_first_success_C:g}-{r.T2_last_success_C:g}",first_step_growth_fraction=r.first_step_growth_fraction,second_step_growth_fraction=r.maximum_second_step_growth_fraction,main_limitation="coarse G1" if r.material_id=="E0021" else "second-step growth >5%",plot_filenames="see figure_inventory.csv"))
    write(TAB/"candidate_visual_summary.csv",summary)
if __name__=="__main__":main()
