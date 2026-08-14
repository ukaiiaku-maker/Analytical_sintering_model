#!/usr/bin/env python3
"""Bounded staged search for observable heterogeneity/residual-stress effects."""
from __future__ import annotations
import argparse,csv,math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict,replace
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import heterogeneous_initial_state_model as hetero
import residual_stress_memory_model as residual
import observable_trajectory_effect_audit as effect
import plot_style as ps
import production_mechanism_assessment as prior
import production_pr_desintering_assessment as production

OUT=Path("results/heterogeneity_residual_stress_search");TARGETS=(.85,.88,.90,.92)


def write(path,rows,fields=None):
    rows=list(rows);fields=fields or list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)


def base_params(mid,G,rho0,topology):
    p=production.candidates()[mid];return replace(p,base=prior.fast_params(p.base,G,rho0,prior.TOPOLOGIES[topology]))


def run_path(hp,rp,rate,peak,hold):
    items=[]
    for spec,cp in hetero.cohort_params(hp):items.append((spec,residual.run(cp,prior.FastSchedule(rate,peak,hold),rp,spec.stress_sign)))
    return hetero.aggregate_histories(items)


def at_density(h,rho):
    if h["rho"].max()<rho-1e-10:return None
    # Aggregate rho is monotonic for accepted paths; use first crossing.
    i=int(np.flatnonzero(h["rho"]>=rho-1e-10)[0]);out={k:float(h[k][i]) for k in ("rho","G_mean","G50","G90","pore_D50","pore_D90","large_pore_fraction","connected_fine_pore_fraction","isolated_pore_fraction","pore_number_proxy","f_GBseg","f_TJ","cumulative_PR_desintering_work")}
    out["t_s"]=float(h["t"][i]);return out


def ratio_curve(ref,fast,metric="G_mean"):
    rows=[]
    for rho in np.arange(.85,.921,.001):
        a=at_density(ref,rho);b=at_density(fast,rho)
        if a and b:rows.append(dict(rho=rho,G_reference=a[metric]*1e9,G_fast=b[metric]*1e9,ratio=a[metric]/b[metric],pct=100*(a[metric]-b[metric])/a[metric]))
    return pd.DataFrame(rows)


def classify(ref,fast):
    curves={m:ratio_curve(ref,fast,m) for m in ("G_mean","G50","G90")};c=curves["G_mean"]
    attained=len(c)>1 and c.rho.min()<=.851 and c.rho.max()>=.919
    span=effect.longest_span(c.rho.to_numpy(),c.ratio.to_numpy(),1.5) if attained else 0.
    tailspan=effect.longest_span(curves["G90"].rho.to_numpy(),curves["G90"].ratio.to_numpy(),1.5) if len(curves["G90"]) else 0.
    return curves,dict(both_paths_attained=attained,max_ratio=float(c.ratio.max()) if len(c) else math.nan,median_ratio=float(c.ratio.median()) if len(c) else math.nan,meaningful_span=span,meaningful=span>=.03,strong=effect.longest_span(c.rho.to_numpy(),c.ratio.to_numpy(),2.)>=.03 if attained else False,tail_only=tailspan>=.03 and span<.03,classification="meaningful" if span>=.03 else ("weak" if len(c) and c.ratio.max()>=1.2 else ("negligible" if len(c) else "unattainable")))


def chen_window_valid(row):
    return bool(row.get("complete_practical_window",False) and row.get("lower_bracketed",False) and row.get("upper_bracketed",False) and row.get("T2_C",math.inf)<row.get("T1_C",-math.inf))


def initial_design():
    modes=hetero.INITIAL_MODES;out=[]
    for i in range(28):
        out.append(dict(stage="initial_distribution",case_id=f"A{i:03d}",base_mechanism=("mech_009","mech_019","mech_009_q0","mech_019_q0")[i%4],initial_microstructure_mode=modes[i%len(modes)],G0_mean_nm=(75,150,225)[i%3],rho0=(.65,.70,.75)[(i//3)%3],topology=("baseline","GBseg_rich","TJ_rich","mixed_GBseg_TJ")[(i//2)%4],pore_ln_sigma=(.35,.65,.95,1.2)[i%4],large_pore_tail_fraction=(0,.02,.05,.10,.20)[i%5],large_pore_radius_factor=(2,4,8,12)[(i//2)%4],bimodal_fine_fraction=(.5,.7,.9)[i%3],G0_ln_sigma=(0,.25,.50,.75)[i%4],coarse_grain_tail_fraction=(0,.05,.10,.20)[i%4],coarse_grain_radius_factor=(2,4,8)[i%3],pore_location_bias=("baseline","GBseg_fine_rich","TJ_large_rich","isolated_large_rich","mixed")[i%5],grain_pore_correlation=("none","large_pores_on_large_grains","large_pores_on_small_grains","TJ_pores_on_small_grains","isolated_pores_on_large_grains")[i%5],residual_stress_mode="disabled",sigma_res_scale=0.,stress_sign="mixed",rate=20 if i%2==0 else 100,peak_T_C=(1350,1400,1450,1500)[i%4],hold_h=(0,2,8,20)[(i//4)%4]))
    return out


def stress_design():
    out=[];modes=residual.MODES
    for i in range(20):out.append(dict(stage="residual_stress",case_id=f"B{i:03d}",base_mechanism=("mech_009","mech_019","mech_009_q0","mech_019_q0")[i%4],initial_microstructure_mode="baseline_narrow",G0_mean_nm=(75,150,225)[i%3],rho0=(.65,.70,.75)[(i//3)%3],topology=("baseline","GBseg_rich","TJ_rich","mixed_GBseg_TJ")[i%4],pore_ln_sigma=.65,large_pore_tail_fraction=.05,large_pore_radius_factor=4,bimodal_fine_fraction=.7,G0_ln_sigma=0.,coarse_grain_tail_fraction=0.,coarse_grain_radius_factor=4,pore_location_bias="baseline",grain_pore_correlation="none",residual_stress_mode=modes[i%5],sigma_res_scale=(0,.25,.5,1.,2.)[(i//4)%5],stress_sign=("compressive","tensile","mixed")[i%3],rate=20 if i%2==0 else 100,peak_T_C=(1350,1400,1450,1500)[i%4],hold_h=(0,2,8,20)[(i//4)%4]))
    return out


def combined_design():
    modes=hetero.INITIAL_MODES[1:];rmodes=residual.MODES[1:];out=[]
    for i in range(96):out.append(dict(stage="combined",case_id=f"C{i:03d}",base_mechanism=("mech_009","mech_019","mech_009_q0","mech_019_q0")[i%4],initial_microstructure_mode=modes[i%len(modes)],G0_mean_nm=(75,150,225)[i%3],rho0=(.65,.70,.75)[(i//3)%3],topology=("baseline","GBseg_rich","TJ_rich","mixed_GBseg_TJ")[(i//2)%4],pore_ln_sigma=(.35,.65,.95,1.2)[i%4],large_pore_tail_fraction=(.02,.05,.10,.20)[i%4],large_pore_radius_factor=(2,4,8,12)[(i//2)%4],bimodal_fine_fraction=(.5,.7,.9)[i%3],G0_ln_sigma=(.25,.50,.75)[i%3],coarse_grain_tail_fraction=(.05,.10,.20)[i%3],coarse_grain_radius_factor=(2,4,8)[i%3],pore_location_bias=("baseline","GBseg_fine_rich","TJ_large_rich","isolated_large_rich","mixed")[i%5],grain_pore_correlation=("none","large_pores_on_large_grains","large_pores_on_small_grains","TJ_pores_on_small_grains","isolated_pores_on_large_grains")[i%5],residual_stress_mode=rmodes[i%4],sigma_res_scale=(.25,.5,1.,2.)[(i//4)%4],stress_sign=("compressive","tensile","mixed")[i%3],rate=20 if i%2==0 else 100,peak_T_C=(1350,1400,1450,1500)[i%4],hold_h=(0,2,8,20)[(i//4)%4]))
    return out


def params_from(row):
    p=base_params(row["base_mechanism"],row["G0_mean_nm"],row["rho0"],row["topology"])
    hp=hetero.HeterogeneousParams(p,**{k:row[k] for k in ("initial_microstructure_mode","G0_mean_nm","pore_ln_sigma","large_pore_tail_fraction","large_pore_radius_factor","bimodal_fine_fraction","G0_ln_sigma","coarse_grain_tail_fraction","coarse_grain_radius_factor","pore_location_bias","grain_pore_correlation")})
    rp=residual.ResidualStressParams(mode=row["residual_stress_mode"],sigma_res_scale=row["sigma_res_scale"],stress_sign=row["stress_sign"])
    return hp,rp


def task(row):
    hp,rp=params_from(row);ref=run_path(hp,rp,1,row["peak_T_C"],row["hold_h"]);fast=run_path(hp,rp,row["rate"],row["peak_T_C"],row["hold_h"]);curves,score=classify(ref,fast);summary={**row,**score,"disabled_recovery":(row["initial_microstructure_mode"]=="baseline_narrow" and row["residual_stress_mode"]=="disabled")}
    metrics=[];pores=[]
    for target in TARGETS:
        a=at_density(ref,target);b=at_density(fast,target)
        if not(a and b):continue
        metrics.append(dict(case_id=row["case_id"],rho_target=target,G_mean_ref_nm=a["G_mean"]*1e9,G_mean_fast_nm=b["G_mean"]*1e9,G50_ref_nm=a["G50"]*1e9,G50_fast_nm=b["G50"]*1e9,G90_ref_nm=a["G90"]*1e9,G90_fast_nm=b["G90"]*1e9,G90_G50_ref=a["G90"]/a["G50"],G90_G50_fast=b["G90"]/b["G50"]))
        pores.append(dict(case_id=row["case_id"],rho_target=target,pore_D50_ref_nm=a["pore_D50"]*1e9,pore_D50_fast_nm=b["pore_D50"]*1e9,pore_D90_ref_nm=a["pore_D90"]*1e9,pore_D90_fast_nm=b["pore_D90"]*1e9,large_pore_ref=a["large_pore_fraction"],large_pore_fast=b["large_pore_fraction"],fine_connected_ref=a["connected_fine_pore_fraction"],fine_connected_fast=b["connected_fine_pore_fraction"],isolated_ref=a["isolated_pore_fraction"],isolated_fast=b["isolated_pore_fraction"],PR_work_ref=a["cumulative_PR_desintering_work"],PR_work_fast=b["cumulative_PR_desintering_work"]))
    curve=[]
    for metric,c in curves.items():
        for x in c.to_dict("records"):curve.append({"case_id":row["case_id"],"metric":metric,**x})
    # Persist representative aggregate histories for baseline and leading-like cases.
    histories=[]
    if row["case_id"] in ("A000","B000","C000") or score["meaningful"] or row.get("persist_history",False):
        for label,h in (("reference",ref),("fast",fast)):
            stride=max(1,len(h["t"])//250)
            for i in range(0,len(h["t"]),stride):
                keys=("rho","G_mean","G50","G90","pore_D50","pore_D90","large_pore_fraction","connected_fine_pore_fraction","isolated_pore_fraction","cumulative_PR_desintering_work","sigma_res_GBseg","sigma_res_TJ","sigma_res_large_pore","sigma_res_crack_like","residual_defect_flux")
                histories.append({"case_id":row["case_id"],"path":label,"t_s":h["t"][i],"T_C":prior.FastSchedule(1 if label=="reference" else row["rate"],row["peak_T_C"],row["hold_h"]).T(h["t"][i],h["rho"][i]),**{k:h[k][i] for k in keys if k in h}})
    return summary,metrics,pores,curve,histories


def plots(out,screens,curves,pores,histories,chen):
    ps.apply_style();figdir=out/"figures";figdir.mkdir(exist_ok=True);s=pd.DataFrame(screens);c=pd.DataFrame(curves);p=pd.DataFrame(pores);h=pd.DataFrame(histories)
    rep=h[h.case_id==("C000" if len(h) and "C000" in set(h.case_id) else (h.case_id.iloc[0] if len(h) else ""))] if len(h) else h
    def save(fig,name):ps.finish(fig,figdir/name)
    fig,ax=plt.subplots(figsize=(5,3.5));
    for path,q in rep.groupby("path"):ax.plot(q.rho,q.G_mean*1e9,label=path)
    ax.set(xlabel="Relative density",ylabel="Mean grain size [nm]");ax.legend();ps.clean(ax);save(fig,"G_rho_trajectories")
    fig,ax=plt.subplots(figsize=(5,3.5));q=c[(c.case_id==("C000" if len(c) and "C000" in set(c.case_id) else (c.case_id.iloc[0] if len(c) else "")))&(c.metric=="G_mean")] if len(c) else c;ax.plot(q.rho,q.ratio) if len(q) else ax.text(.5,.5,"No jointly attained ratio curve",ha="center");[ax.axhline(v,color="#777",ls="--") for v in (1.2,1.5,2)];ax.set(xlabel="Relative density",ylabel="Grain-size ratio");ps.clean(ax);save(fig,"grain_ratio_vs_density")
    fig,ax=plt.subplots(figsize=(5,3.5));
    for k in ("G_mean","G50","G90"):ax.plot(rep[rep.path=="fast"].rho,rep[rep.path=="fast"][k]*1e9,label=k)
    ax.set(xlabel="Relative density",ylabel="Grain metric [nm]");ax.legend();ps.clean(ax);save(fig,"grain_distribution_metrics")
    fig,ax=plt.subplots(figsize=(5,3.5));
    if len(p):ax.scatter(p.pore_D50_fast_nm,p.pore_D90_fast_nm,c=p.large_pore_fast,s=8)
    ax.set(xlabel="Fast pore D50 [nm]",ylabel="Fast pore D90 [nm]");ps.clean(ax);save(fig,"pore_distribution_metrics")
    fig,ax=plt.subplots(figsize=(5,3.5));
    for path,q in rep.groupby("path"):ax.plot(q.rho,q.connected_fine_pore_fraction,label=path)
    ax.set(xlabel="Relative density",ylabel="Connected fine-pore fraction");ax.legend();ps.clean(ax);save(fig,"connected_fine_pores")
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),constrained_layout=True)
    if len(rep) and "sigma_res_large_pore" in rep:
        for path,q in rep.groupby("path"):axs[0].plot(q.rho,q.sigma_res_large_pore/1e6,label=path);axs[1].plot(q.T_C,q.sigma_res_large_pore/1e6,label=path)
    else:axs[0].text(.5,.5,"No residual-stress history",ha="center")
    axs[0].set(xlabel="Relative density",ylabel="Large-pore residual stress [MPa]");axs[1].set(xlabel="Temperature [°C]",ylabel="Large-pore residual stress [MPa]")
    for ax in axs:ps.clean(ax)
    axs[0].legend(fontsize=7);ps.panel_labels(axs);save(fig,"residual_stress_histories")
    fig,ax=plt.subplots(figsize=(5,3.5));
    for path,q in rep.groupby("path"):ax.plot(q.rho,q.cumulative_PR_desintering_work,label=path)
    ax.set(xlabel="Relative density",ylabel="Cumulative PR work");ax.legend();ps.clean(ax);save(fig,"PR_work_vs_density")
    fig,ax=plt.subplots(figsize=(6,3.5));z=s.pivot_table(index="stage",columns="initial_microstructure_mode",values="max_ratio",aggfunc="median");im=ax.imshow(z,aspect="auto",vmin=1,vmax=1.5,cmap="viridis");ax.set_xticks(range(len(z.columns)),z.columns,rotation=35,ha="right");ax.set_yticks(range(len(z)),z.index);fig.colorbar(im,ax=ax,label="Median maximum ratio");save(fig,"fast_effect_heatmap")
    fig,ax=plt.subplots(figsize=(5,3.5));ax.text(.5,.5,"No candidate promoted" if not chen else "Adaptive Chen preservation results",ha="center");ax.axis("off");save(fig,"Chen_preservation_map")
    fig,ax=plt.subplots(figsize=(7,2.5));ax.axis("off");labels=("Initial nonideality\n+ residual stress","PR / large-pore\ndamage","later densification\nefficiency","observable $G(\\rho)$")
    for i,label in enumerate(labels):ax.text(.12+i*.25,.5,label,ha="center",bbox=dict(boxstyle="round",fc="white"));
    for i in range(3):ax.annotate("",xy=(.22+i*.25,.5),xytext=(.17+i*.25,.5),arrowprops=dict(arrowstyle="->"));save(fig,"mechanism_pathway_summary")


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--workers",type=int,default=4);ap.add_argument("--limit",type=int);ap.add_argument("--refresh-representative",action="store_true");ap.add_argument("--plots-only",action="store_true");args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    if args.refresh_representative:
        row={**combined_design()[91],"persist_history":True};_,_,_,_,hist=task(row);write(OUT/"representative_path_histories.csv",hist);stress=[r for r in hist if "sigma_res_GBseg" in r];write(OUT/"residual_stress_histories.csv",stress);return
    if args.plots_only:
        screens=pd.concat([pd.read_csv(OUT/f) for f in ("initial_distribution_screen.csv","residual_stress_screen.csv","combined_screen.csv")],ignore_index=True).to_dict("records")
        plots(OUT,screens,pd.read_csv(OUT/"density_window_effects.csv").to_dict("records"),pd.read_csv(OUT/"pore_distribution_metrics.csv").to_dict("records"),pd.read_csv(OUT/"representative_path_histories.csv").to_dict("records"),pd.read_csv(OUT/"Chen_preservation_summary.csv").to_dict("records"));return
    design=initial_design()+stress_design()+combined_design();design=design[:args.limit] if args.limit else design
    registry=[{**r,"production_k_PR_ref_s":2e-4,"budget_h":96,"targets":"0.85;0.88;0.90;0.92"} for r in design];write(OUT/"parameter_registry.csv",registry)
    results=[];grains=[];pores=[];curves=[];hist=[]
    mapped=map(task,design)
    if args.workers>1:
        pool=ProcessPoolExecutor(max_workers=args.workers);mapped=pool.map(task,design,chunksize=1)
    else:pool=None
    try:
        for i,(summary,g,p,c,h) in enumerate(mapped,1):results.append(summary);grains+=g;pores+=p;curves+=c;hist+=h;print("heterogeneity/stress",i,"/",len(design),flush=True) if i%12==0 else None
    finally:
        if pool:pool.shutdown()
    for stage,name in (("initial_distribution","initial_distribution_screen.csv"),("residual_stress","residual_stress_screen.csv"),("combined","combined_screen.csv")):write(OUT/name,[r for r in results if r["stage"]==stage])
    rejected=[{**r,"rejection_reason":"joint_target_nonattainment" if not r["both_paths_attained"] else ("tail_only_not_mean_G" if r["tail_only"] else "mean_ratio_below_1.5_over_Drho_0.03")} for r in results if not r["meaningful"]]
    meaningful=[r for r in results if r["meaningful"]];weak=[r for r in results if not r["meaningful"]]
    write(OUT/"rejected_cases.csv",rejected);write(OUT/"meaningful_trajectory_cases.csv",meaningful,list(results[0]) if results else None);write(OUT/"weak_trajectory_cases.csv",weak);write(OUT/"density_window_effects.csv",curves);write(OUT/"grain_distribution_metrics.csv",grains);write(OUT/"pore_distribution_metrics.csv",pores)
    stress=[r for r in hist if "sigma_res_GBseg" in r];write(OUT/"residual_stress_histories.csv",stress,["case_id","path","t_s","T_C","rho","sigma_res_GBseg","sigma_res_TJ","sigma_res_large_pore","sigma_res_crack_like","residual_defect_flux"] if not stress else None);write(OUT/"representative_path_histories.csv",hist)
    chen=[]
    for r in meaningful:chen.append(dict(case_id=r["case_id"],status="not_scored",reason="cohort candidate requires aggregate adaptive solver",complete_practical_window=False))
    write(OUT/"Chen_preservation_summary.csv",chen,["case_id","status","reason","complete_practical_window"]);plots(OUT,results,curves,pores,hist,chen)
    print("DONE",len(results),"meaningful",len(meaningful),flush=True)

if __name__=="__main__":main()
