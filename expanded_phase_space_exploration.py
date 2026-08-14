#!/usr/bin/env python3
"""Expanded fixed-model phase-space exploration with calculated refinements."""
from __future__ import annotations

import argparse,csv,gzip,math,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import fields,replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from density_window_processing_map import sample_at_density,write_csv
from two_step_window_map import classify

G0_GRID=(25.,35.,50.,75.,100.,150.,225.,300.,450.,600.,800.,1000.,1500.,2000.)
STYLES=("density","connectivity","hybrid_topology")
T2_COARSE=tuple(float(x) for x in range(900,1351,25))
TARGETS=(.88,.90,.92,.94,.95,.96,.98)
TOLERANCES=(.025,.05,.075,.10,.15,.20)
RHO0_GRID=(.60,.65,.70,.75,.80)
RATES=(.1,.2,.5,1.,2.,5.,10.,20.,50.,100.)
PEAKS=(1300.,1400.,1450.)
STEP_BUDGET_S=96*3600
VARYING={"rho0","G0","smoothing_gate_mode"}


def base_params()->model.Params:
    return model.Params(memory_model="pore_bin_redistribution",rho0=.70,G0=150e-9,pore_radius0=25e-9,pore_ln_sigma=.65)


def assert_fixed(params_list:list[model.Params],base:model.Params)->None:
    fixed={f.name:getattr(base,f.name) for f in fields(base) if f.name not in VARYING}
    for p in params_list:
        if {f.name:getattr(p,f.name) for f in fields(p) if f.name not in VARYING}!=fixed:raise AssertionError("fixed material parameter changed")


def first_metrics(result:dict)->dict:
    i=-1
    return {"rho1":float(result["rho"][i]),"G1_nm":float(result["G"][i])*1e9,"first_time_h":float(result["t"][i])/3600,
            "first_mean_radius_nm":float(result["pore_mean_radius"][i])*1e9,"first_large_fraction":float(result["large_pore_fraction"][i]),
            "first_fine_fraction":float(result["removable_fine_pore_fraction"][i]),"first_f_pore":float(result["f_pore"][i]),
            "first_connected_coverage":float(result["f_pore"][i]*result["connectivity"][i]),"first_isolation":float(result["isolated_pore_fraction"][i]),
            "first_activity":float(result["activity"][i]),"first_E_G":float(result["E_G"][i])}


def second_metrics(result:dict)->dict:
    i=-1
    return {"rho2":float(result["rho"][i]),"G2_nm":float(result["G"][i])*1e9,"second_time_h":float(result["t"][i])/3600,
            "second_initial_rho_dot":float(result["rho_dot"][0]),"second_final_rho_dot":float(result["rho_dot"][i]),
            "second_initial_dGdt_nm_s":float(result["dGdt"][0])*1e9,"second_final_dGdt_nm_s":float(result["dGdt"][i])*1e9,
            "second_median_E_G":float(np.median(result["E_G"])),"second_mean_radius_nm":float(result["pore_mean_radius"][i])*1e9,
            "second_large_fraction":float(result["large_pore_fraction"][i]),"second_fine_fraction":float(result["removable_fine_pore_fraction"][i]),
            "second_f_pore":float(result["f_pore"][i]),"second_connected_coverage":float(result["f_pore"][i]*result["connectivity"][i]),
            "second_isolation":float(result["isolated_pore_fraction"][i]),"second_activity":float(result["activity"][i])}


def empty_first_metrics()->dict:
    return {key:np.nan for key in ("rho1","G1_nm","first_time_h","first_mean_radius_nm","first_large_fraction","first_fine_fraction","first_f_pore","first_connected_coverage","first_isolation","first_activity","first_E_G")}


def empty_second_metrics()->dict:
    return {key:np.nan for key in ("rho2","G2_nm","second_time_h","second_initial_rho_dot","second_final_rho_dot","second_initial_dGdt_nm_s","second_final_dGdt_nm_s","second_median_E_G","second_mean_radius_nm","second_large_fraction","second_fine_fraction","second_f_pore","second_connected_coverage","second_isolation","second_activity")}


def two_step_group_task(args:tuple)->list[dict]:
    design,style,g0,rho0,T1,switch,T2_values=args;p=replace(base_params(),smoothing_gate_mode=style,G0=g0*1e-9,rho0=rho0)
    first=model.run(p,model.Iso(T1,STEP_BUDGET_S),stop_at_rho=switch);attained=float(np.max(first["rho"]))>=switch-1e-12;fm=first_metrics(first)
    state=model.state_from_result(first,p) if attained else None;rows=[]
    for T2 in T2_values:
        if attained:
            second=model.run(p,model.Iso(T2,STEP_BUDGET_S),initial=state);sm=second_metrics(second);growth=(sm["G2_nm"]-fm["G1_nm"])/fm["G1_nm"]
            density_gain=sm["rho2"]-fm["rho1"]
        else:sm=empty_second_metrics();growth=density_gain=np.nan
        rows.append({"design":design,"model_style":style,"G0_nm":g0,"rho0":rho0,"T1_C":T1,"rho_switch":switch,"T2_C":T2,
                     "first_step_attained":attained,"first_step_budget_h":STEP_BUDGET_S/3600,"second_step_budget_h":STEP_BUDGET_S/3600,
                     "density_gain":density_gain,"growth_fraction":growth,**fm,**sm})
    return rows


def group_design()->list[tuple]:
    groups=set()
    for style in STYLES:
        for g0 in G0_GRID:
            for T1 in (1200.,1300.,1400.):
                for switch in (.75,.825,.875,.925):groups.add(("core_size",style,g0,.70,T1,switch,T2_COARSE))
    for style in ("density","connectivity"):
        for g0 in (75.,450.,1000.,2000.):
            for T1 in tuple(float(x) for x in range(1150,1451,50)):
                for switch in (.70,.75,.80,.825,.85,.875,.90,.925):groups.add(("T1_switch_extension",style,g0,.70,T1,switch,T2_COARSE))
        for g0 in (75.,450.,1000.,2000.):
            for rho0 in RHO0_GRID:groups.add(("rho0_probe",style,g0,rho0,1300.,.85,T2_COARSE))
    return sorted(groups,key=lambda x:x[:-1])


def run_groups(groups:list[tuple],workers:int,label:str)->list[dict]:
    rows=[];started=time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(two_step_group_task,g) for g in groups]
        for index,future in enumerate(as_completed(futures),1):
            rows.extend(future.result())
            if index%25==0 or index==len(groups):print(f"[{label} {index}/{len(groups)}] trajectories={len(rows)}",flush=True)
    print(f"{label} wall_s={time.perf_counter()-started:.2f}",flush=True);return rows


def failure_detail(row:dict,target:float,tolerance:float)->str:
    category=expanded_classify(row,target,tolerance)
    if category not in ("DENSIFICATION_EXHAUSTION_FAILURE","MIXED_FAILURE"):return category
    if row.get("second_connected_coverage",1)<.05 or row.get("second_isolation",0)>.8:return "TOPOLOGY_EXHAUSTION"
    remaining=target-row.get("rho2",0);projected=max(row.get("second_final_rho_dot",0),0)*STEP_BUDGET_S
    return "KINETIC_EXHAUSTION" if projected<remaining else "TIME_BUDGET_EXHAUSTION"


def expanded_classify(row:dict,target:float,tolerance:float)->str:
    if row["first_step_attained"] and target<=row["rho1"]+1e-12:return "INELIGIBLE_TARGET_ALREADY_REACHED"
    return classify(row["first_step_attained"],row.get("rho2",np.nan),row["growth_fraction"],target,tolerance)


def write_classifications(path:Path,trajectories:list[dict])->None:
    fields_out=("design","model_style","G0_nm","rho0","T1_C","rho_switch","T2_C","rho_target","growth_tolerance","eligible_second_step_target","classification","dominant_failure")
    with gzip.open(path,"wt",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields_out,lineterminator="\n");writer.writeheader()
        for row in trajectories:
            for target in TARGETS:
                for tolerance in TOLERANCES:
                    category=expanded_classify(row,target,tolerance)
                    writer.writerow({**{k:row[k] for k in fields_out[:7]},"rho_target":target,"growth_tolerance":tolerance,
                                     "eligible_second_step_target":target>row["rho1"]+1e-12,
                                     "classification":category,"dominant_failure":failure_detail(row,target,tolerance)})


def boundary_rows(trajectories:list[dict])->list[dict]:
    output=[];keys={(r["design"],r["model_style"],r["G0_nm"],r["rho0"],r["T1_C"],r["rho_switch"]) for r in trajectories}
    for key in sorted(keys):
        group=[r for r in trajectories if (r["design"],r["model_style"],r["G0_nm"],r["rho0"],r["T1_C"],r["rho_switch"])==key]
        for target in TARGETS:
            eligible=target>group[0]["rho1"]+1e-12
            dense=[r for r in group if eligible and r["first_step_attained"] and r.get("rho2",0)>=target-1e-12]
            min_required=min((r["growth_fraction"] for r in dense),default=np.nan)
            for tolerance in TOLERANCES:
                nogrowth=[r for r in group if eligible and r["first_step_attained"] and r["growth_fraction"]<=tolerance+1e-12]
                success=[r for r in dense if r["growth_fraction"]<=tolerance+1e-12]
                output.append({"design":key[0],"model_style":key[1],"G0_nm":key[2],"rho0":key[3],"T1_C":key[4],"rho_switch":key[5],
                    "rho_target":target,"growth_tolerance":tolerance,"first_step_attained":any(r["first_step_attained"] for r in group),
                    "eligible_second_step_target":eligible,
                    "rho1":group[0]["rho1"],"G1_nm":group[0]["G1_nm"],"first_connected_coverage":group[0]["first_connected_coverage"],
                    "T_lower_density_C":min((r["T2_C"] for r in dense),default=np.nan),"T_upper_no_growth_C":max((r["T2_C"] for r in nogrowth),default=np.nan),
                    "T_success_lower_C":min((r["T2_C"] for r in success),default=np.nan),"T_success_upper_C":max((r["T2_C"] for r in success),default=np.nan),
                    "window_width_C":max((r["T2_C"] for r in success),default=np.nan)-min((r["T2_C"] for r in success),default=np.nan) if success else np.nan,
                    "optimum_T2_C":min(success,key=lambda r:r["growth_fraction"])["T2_C"] if success else np.nan,
                    "minimum_required_growth_tolerance":min_required})
    return output


def refinement_groups(boundaries:list[dict])->list[tuple]:
    groups={}
    for row in boundaries:
        if row["rho_target"]!=.90 or row["growth_tolerance"] not in (.05,.10):continue
        candidates=[]
        for value in (row["T_lower_density_C"],row["T_upper_no_growth_C"]):
            if np.isfinite(value):candidates.extend(value+d for d in (-20,-10,10,20))
        refined=tuple(sorted({float(x) for x in candidates if 900<=x<=1350 and x not in T2_COARSE}))
        if refined:
            key=(row["design"],row["model_style"],row["G0_nm"],row["rho0"],row["T1_C"],row["rho_switch"])
            groups[key]=tuple(sorted(set(groups.get(key,()))|set(refined)))
    return [(*key,values) for key,values in groups.items()]


def fast_task(args:tuple)->dict:
    style,g0,rho0,rate,peak=args;p=replace(base_params(),smoothing_gate_mode=style,G0=g0*1e-9,rho0=rho0);protocol=model.RampHold(rate,target_C=peak)
    result=model.run(p,protocol);sample=sample_at_density(result,.90,rho0);i=-1;peak_i=int(np.argmax(result["T_C"]))
    return {"model_style":style,"G0_nm":g0,"rho0":rho0,"heating_rate_C_min":rate,"peak_target_C":peak,"actual_peak_C":float(result["T_C"][peak_i]),
            "time_budget_h":min(protocol.t_end,p.t_max_s)/3600,"reached_0p90":sample["reached_target"],"G_at_0p90_nm":sample["G_at_target_nm"],
            "time_to_0p90_h":sample["time_to_target_h"],"final_rho":float(result["rho"][i]),"final_G_nm":float(result["G"][i])*1e9,
            "final_mean_radius_nm":float(result["pore_mean_radius"][i])*1e9,"final_fine_fraction":float(result["removable_fine_pore_fraction"][i]),
            "final_connected_coverage":float(result["f_pore"][i]*result["connectivity"][i]),"peak_rho":float(result["rho"][peak_i]),
            "peak_G_nm":float(result["G"][peak_i])*1e9,"peak_connected_coverage":float(result["f_pore"][peak_i]*result["connectivity"][peak_i])}


def run_fast(workers:int)->list[dict]:
    tasks=[(s,g,rho,rate,peak) for s in STYLES for g in G0_GRID for rho in (.60,.70,.80) for rate in RATES for peak in PEAKS]
    rows=[];started=time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(fast_task,t) for t in tasks]
        for i,f in enumerate(as_completed(futures),1):
            rows.append(f.result())
            if i%250==0 or i==len(tasks):print(f"[fast {i}/{len(tasks)}]",flush=True)
    lookup={(r["model_style"],r["G0_nm"],r["rho0"],r["peak_target_C"],r["heating_rate_C_min"]):r for r in rows}
    for row in rows:
        ref=lookup[(row["model_style"],row["G0_nm"],row["rho0"],row["peak_target_C"],.1)]
        ref_0p2=lookup[(row["model_style"],row["G0_nm"],row["rho0"],row["peak_target_C"],.2)]
        row["HR_pct_vs_0p1"]=model.percent_gain(ref["G_at_0p90_nm"],row["G_at_0p90_nm"]) if ref["reached_0p90"] and row["reached_0p90"] else np.nan
        row["HR_pct_vs_0p2"]=model.percent_gain(ref_0p2["G_at_0p90_nm"],row["G_at_0p90_nm"]) if ref_0p2["reached_0p90"] and row["reached_0p90"] else np.nan
    print(f"fast wall_s={time.perf_counter()-started:.2f}",flush=True);return sorted(rows,key=lambda r:(r["model_style"],r["G0_nm"],r["rho0"],r["peak_target_C"],r["heating_rate_C_min"]))


def matched_history()->list[dict]:
    rows=[]
    for style in STYLES:
        for g0 in (35.,150.,600.,1000.,2000.):
            for target in (.80,.85):
                for rate in (.1,100.):
                    p=replace(base_params(),smoothing_gate_mode=style,G0=g0*1e-9);history=model.run(p,model.RampHold(rate),stop_at_rho=target);reached=float(np.max(history["rho"]))>=target-1e-12
                    if reached:
                        state=model.state_from_result(history,p);follow=model.run(p,model.Iso(1200.,STEP_BUDGET_S),initial=state);fm=first_metrics(history);sm=second_metrics(follow)
                    else:fm=empty_first_metrics();sm=empty_second_metrics()
                    rows.append({"model_style":style,"G0_nm":g0,"matched_density":target,"history_rate_C_min":rate,"history_reached":reached,**fm,**sm})
    return rows


def combined_history()->list[dict]:
    groups=[]
    for style in ("density","connectivity"):
        for g0 in (75.,800.,2000.):
            for rate in (.2,5.,20.):
                p=replace(base_params(),smoothing_gate_mode=style,G0=g0*1e-9);protocol=model.RampHold(rate,target_C=1300.,hold_s=STEP_BUDGET_S)
                first=model.run(p,protocol,stop_at_rho=.85);attained=float(np.max(first["rho"]))>=.85-1e-12;fm=first_metrics(first);state=model.state_from_result(first,p) if attained else None
                for T2 in tuple(float(x) for x in range(1000,1301,25)):
                    if attained:
                        second=model.run(p,model.Iso(T2,STEP_BUDGET_S),initial=state);sm=second_metrics(second);growth=(sm["G2_nm"]-fm["G1_nm"])/fm["G1_nm"]
                    else:sm=empty_second_metrics();growth=np.nan
                    groups.append({"model_style":style,"G0_nm":g0,"T1_heating_rate_C_min":rate,"rho_switch":.85,"T2_C":T2,"first_step_attained":attained,
                                   "classification_5pct":classify(attained,sm.get("rho2",np.nan),growth,.90,.05),"growth_fraction":growth,**fm,**sm})
    return groups


def make_plots(outdir:Path,trajectories:list[dict],boundaries:list[dict],fast:list[dict])->None:
    b=[r for r in boundaries if r["design"]=="core_size" and r["rho_target"]==.90 and r["growth_tolerance"]==.05]
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharey=True)
    for axis,style in zip(axes,STYLES):
        g=[r for r in b if r["model_style"]==style]
        axis.scatter([r["G0_nm"] for r in g],[r["window_width_C"] for r in g],alpha=.6);axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("5% window width [C]");fig.tight_layout();fig.savefig(outdir/"window_width_vs_G0.png",dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharey=True)
    for axis,style in zip(axes,STYLES):
        g=[r for r in b if r["model_style"]==style]
        axis.scatter([r["G0_nm"] for r in g],[100*r["minimum_required_growth_tolerance"] for r in g],c=[r["T1_C"] for r in g],cmap="viridis",alpha=.6)
        axis.axhline(5,color="black",linestyle="--");axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("Minimum growth tolerance [%]");fig.tight_layout();fig.savefig(outdir/"minimum_growth_tolerance_vs_G0.png",dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharey=True)
    for axis,style in zip(axes,STYLES):
        g=[r for r in b if r["model_style"]==style]
        axis.scatter([r["G0_nm"] for r in g],[r["T_lower_density_C"] for r in g],label="lower density")
        axis.scatter([r["G0_nm"] for r in g],[r["T_upper_no_growth_C"] for r in g],label="upper 5% growth")
        axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 boundary [C]");axes[-1].legend(fontsize=8);fig.tight_layout();fig.savefig(outdir/"T2_boundaries_vs_G0.png",dpi=150);plt.close(fig)
    subset=[r for r in trajectories if r["design"]=="core_size" and r["T1_C"]==1300 and r["rho_switch"]==.825]
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharey=True)
    for axis,style in zip(axes,STYLES):
        g=[r for r in subset if r["model_style"]==style];c=[100*r["growth_fraction"] if r.get("rho2",0)>=.90 else -1 for r in g]
        im=axis.scatter([r["G0_nm"] for r in g],[r["T2_C"] for r in g],c=c,cmap="coolwarm",vmin=-1,vmax=20);axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 [C]");fig.colorbar(im,ax=axes,label="growth % if rho>=0.90; -1 otherwise");fig.subplots_adjust(right=.9,wspace=.15);fig.savefig(outdir/"G0_T2_success_topology.png",dpi=150);plt.close(fig)
    subset=[r for r in fast if r["rho0"]==.70 and r["peak_target_C"]==1400 and np.isfinite(r["HR_pct_vs_0p2"])]
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharey=True)
    for axis,style in zip(axes,STYLES):
        g=[r for r in subset if r["model_style"]==style];im=axis.scatter([r["G0_nm"] for r in g],[r["heating_rate_C_min"] for r in g],c=[r["HR_pct_vs_0p2"] for r in g],cmap="coolwarm",vmin=-10,vmax=10)
        axis.set_xscale("log");axis.set_yscale("log");axis.set(xlabel="G0 [nm]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("Heating rate [C/min]");fig.colorbar(im,ax=axes,label="HR_pct vs 0.2 C/min");fig.subplots_adjust(right=.9,wspace=.15);fig.savefig(outdir/"fast_firing_G0_rate_surface.png",dpi=150);plt.close(fig)


def main()->None:
    parser=argparse.ArgumentParser();parser.add_argument("--outdir",default="results/expanded_phase_space");parser.add_argument("--workers",type=int,default=4);args=parser.parse_args()
    outdir=Path(args.outdir);outdir.mkdir(parents=True,exist_ok=True);manifest=group_design();started=time.perf_counter()
    trajectories=run_groups(manifest,args.workers,"coarse");write_csv(outdir/"coarse_checkpoint.csv",trajectories)
    boundaries=boundary_rows(trajectories);ref_groups=refinement_groups(boundaries)
    refinements=run_groups(ref_groups,args.workers,"refine") if ref_groups else [];write_csv(outdir/"refinement_checkpoint.csv",refinements) if refinements else None
    all_trajectories=trajectories+refinements;all_boundaries=boundary_rows(all_trajectories)
    fast=run_fast(args.workers);write_csv(outdir/"fast_checkpoint.csv",fast)
    matched=matched_history();write_csv(outdir/"matched_checkpoint.csv",matched)
    combined=combined_history();write_csv(outdir/"combined_checkpoint.csv",combined)
    write_csv(outdir/"two_step_trajectories.csv",all_trajectories);write_classifications(outdir/"two_step_classifications.csv.gz",all_trajectories)
    write_csv(outdir/"two_step_boundaries.csv",all_boundaries);write_csv(outdir/"refined_two_step_trajectories.csv",refinements)
    write_csv(outdir/"fast_firing_surface.csv",fast);write_csv(outdir/"matched_history_expanded.csv",matched);write_csv(outdir/"combined_history.csv",combined)
    write_csv(outdir/"sweep_manifest.csv",[{"design":g[0],"model_style":g[1],"G0_nm":g[2],"rho0":g[3],"T1_C":g[4],"rho_switch":g[5],"n_T2":len(g[6])} for g in manifest])
    make_plots(outdir,all_trajectories,all_boundaries,fast)
    write_csv(outdir/"runtime_summary.csv",[{"coarse_groups":len(manifest),"coarse_trajectories":len(trajectories),"refined_groups":len(ref_groups),"refined_trajectories":len(refinements),"fast_trajectories":len(fast),"matched_histories":len(matched),"combined_trajectories":len(combined),"total_wall_s":time.perf_counter()-started}])
    print(f"DONE coarse={len(trajectories)} refined={len(refinements)} fast={len(fast)} matched={len(matched)} combined={len(combined)} wall_s={time.perf_counter()-started:.2f}")


if __name__=="__main__":main()
