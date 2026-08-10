#!/usr/bin/env python3
"""Bounded Chen-style ablation of nanoscale GB-migration mobility closures."""
from __future__ import annotations
import argparse,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from dataclasses import fields,replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from density_window_processing_map import write_csv
from two_step_window_map import classify

GROWTH_MODES=("baseline","junction_limited","threshold_mobility")
G0_VALUES=(25.,35.,50.,75.,100.,150.,225.,300.,450.,600.)
T1_VALUES=(1250.,1300.,1350.)
SWITCHES=(.75,.80,.85)
T2_VALUES=tuple(float(x) for x in range(900,1301,25))
TOLERANCES=(.05,.10)
TARGET=.90
RHO0=.70
STEP_BUDGET_S=96*3600
VARYING={"growth_mode","G0"}


def base_params()->model.Params:
    return model.Params(memory_model="pore_bin_redistribution",smoothing_gate_mode="density",growth_mode="baseline",
                        rho0=RHO0,G0=150e-9,pore_radius0=25e-9,pore_ln_sigma=.65)


def assert_only_growth_mode_and_size_vary(params_list:list[model.Params],base:model.Params)->None:
    fixed={f.name:getattr(base,f.name) for f in fields(base) if f.name not in VARYING}
    for p in params_list:
        if {f.name:getattr(p,f.name) for f in fields(p) if f.name not in VARYING}!=fixed:raise AssertionError("non-growth material parameter changed")


def group_task(args:tuple)->list[dict]:
    mode,g0,T1,switch=args;p=replace(base_params(),growth_mode=mode,G0=g0*1e-9)
    first=model.run(p,model.Iso(T1,STEP_BUDGET_S),stop_at_rho=switch);attained=float(np.max(first["rho"]))>=switch-1e-12;i=-1
    common={"growth_mode":mode,"G0_nm":g0,"T1_C":T1,"rho_switch":switch,"rho0":RHO0,"first_step_attained":attained,
            "rho1":float(first["rho"][i]),"G1_nm":float(first["G"][i])*1e9,"first_connected_coverage":float(first["f_pore"][i]*first["connectivity"][i]),
            "first_fine_fraction":float(first["removable_fine_pore_fraction"][i]),"first_isolation":float(first["isolated_pore_fraction"][i]),
            "first_growth_mobility_factor":float(first["growth_mobility_factor"][i]),"first_step_budget_h":STEP_BUDGET_S/3600,"second_step_budget_h":STEP_BUDGET_S/3600}
    state=model.state_from_result(first,p) if attained else None;rows=[]
    for T2 in T2_VALUES:
        if attained:
            second=model.run(p,model.Iso(T2,STEP_BUDGET_S),initial=state);j=-1;rho2=float(second["rho"][j]);G2=float(second["G"][j])*1e9;growth=(G2-common["G1_nm"])/common["G1_nm"]
            diagnostics={"rho2":rho2,"G2_nm":G2,"density_gain":rho2-common["rho1"],"growth_fraction":growth,
                "second_initial_rho_dot":float(second["rho_dot"][0]),"second_initial_dGdt_nm_s":float(second["dGdt"][0])*1e9,
                "second_final_rho_dot":float(second["rho_dot"][j]),"second_final_dGdt_nm_s":float(second["dGdt"][j])*1e9,
                "median_E_G":float(np.median(second["E_G"])),"initial_growth_mobility_factor":float(second["growth_mobility_factor"][0]),
                "final_growth_mobility_factor":float(second["growth_mobility_factor"][j]),"initial_junction_time_ratio":float(second["junction_time_ratio"][0]),
                "initial_migration_drive_ratio":float(second["migration_drive_ratio"][0]),"final_connected_coverage":float(second["f_pore"][j]*second["connectivity"][j]),
                "final_mean_radius_nm":float(second["pore_mean_radius"][j])*1e9,"final_fine_fraction":float(second["removable_fine_pore_fraction"][j]),
                "final_isolation":float(second["isolated_pore_fraction"][j])}
        else:
            diagnostics={key:np.nan for key in ("rho2","G2_nm","density_gain","growth_fraction","second_initial_rho_dot","second_initial_dGdt_nm_s","second_final_rho_dot","second_final_dGdt_nm_s","median_E_G","initial_growth_mobility_factor","final_growth_mobility_factor","initial_junction_time_ratio","initial_migration_drive_ratio","final_connected_coverage","final_mean_radius_nm","final_fine_fraction","final_isolation")}
        rows.append({**common,"T2_C":T2,**diagnostics})
    return rows


def run_map(workers:int)->tuple[list[dict],float]:
    settings=[(mode,g,T1,switch) for mode in GROWTH_MODES for g in G0_VALUES for T1 in T1_VALUES for switch in SWITCHES]
    params=[replace(base_params(),growth_mode=m,G0=g*1e-9) for m,g,_,_ in settings];assert_only_growth_mode_and_size_vary(params,base_params())
    rows=[];started=time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(group_task,s) for s in settings]
        for index,future in enumerate(as_completed(futures),1):
            rows.extend(future.result())
            if index%20==0 or index==len(settings):print(f"[{index}/{len(settings)}] trajectories={len(rows)}",flush=True)
    return rows,time.perf_counter()-started


def classification_rows(trajectories:list[dict])->list[dict]:
    rows=[]
    for r in trajectories:
        for tolerance in TOLERANCES:
            rows.append({"growth_mode":r["growth_mode"],"G0_nm":r["G0_nm"],"T1_C":r["T1_C"],"rho_switch":r["rho_switch"],"T2_C":r["T2_C"],
                         "rho_target":TARGET,"growth_tolerance":tolerance,"first_step_attained":r["first_step_attained"],
                         "classification":classify(r["first_step_attained"],r["rho2"],r["growth_fraction"],TARGET,tolerance)})
    return rows


def boundary_rows(trajectories:list[dict])->list[dict]:
    output=[]
    for mode in GROWTH_MODES:
        for g0 in G0_VALUES:
            for T1 in T1_VALUES:
                for switch in SWITCHES:
                    group=[r for r in trajectories if r["growth_mode"]==mode and r["G0_nm"]==g0 and r["T1_C"]==T1 and r["rho_switch"]==switch]
                    dense=[r for r in group if r["first_step_attained"] and r["rho2"]>=TARGET-1e-12];minimum=min((r["growth_fraction"] for r in dense),default=np.nan)
                    for tolerance in TOLERANCES:
                        no_growth=[r for r in group if r["first_step_attained"] and r["growth_fraction"]<=tolerance+1e-12];success=[r for r in dense if r["growth_fraction"]<=tolerance+1e-12]
                        output.append({"growth_mode":mode,"G0_nm":g0,"T1_C":T1,"rho_switch":switch,"rho_target":TARGET,"growth_tolerance":tolerance,
                            "first_step_attained":any(r["first_step_attained"] for r in group),"rho1":group[0]["rho1"],"G1_nm":group[0]["G1_nm"],
                            "T_lower_density_C":min((r["T2_C"] for r in dense),default=np.nan),"T_upper_no_growth_C":max((r["T2_C"] for r in no_growth),default=np.nan),
                            "T_success_lower_C":min((r["T2_C"] for r in success),default=np.nan),"T_success_upper_C":max((r["T2_C"] for r in success),default=np.nan),
                            "window_width_C":max((r["T2_C"] for r in success),default=np.nan)-min((r["T2_C"] for r in success),default=np.nan) if success else np.nan,
                            "n_success_T2":len(success),"minimum_required_growth_tolerance":minimum})
    return output


def baseline_recovery(trajectories:list[dict])->list[dict]:
    from expanded_phase_space_analysis import load
    previous=load(Path("results/expanded_phase_space/two_step_trajectories.csv"));rows=[]
    for g0 in G0_VALUES:
        old=[r for r in previous if r["design"]=="core_size" and r["model_style"]=="density" and r["G0_nm"]==g0 and r["T1_C"]==1300 and r["rho_switch"]==.75 and r["T2_C"] in T2_VALUES]
        new=[r for r in trajectories if r["growth_mode"]=="baseline" and r["G0_nm"]==g0 and r["T1_C"]==1300 and r["rho_switch"]==.75]
        lookup={r["T2_C"]:r for r in old}
        for r in new:
            o=lookup[r["T2_C"]];rows.append({"G0_nm":g0,"T2_C":r["T2_C"],"rho2_abs_difference":abs(r["rho2"]-o["rho2"]),"G2_nm_abs_difference":abs(r["G2_nm"]-o["G2_nm"])})
    return rows


def make_plots(out:Path,trajectories:list[dict],classes:list[dict],boundaries:list[dict])->None:
    colors={"SUCCESS":"tab:green","GRAIN_GROWTH_FAILURE":"tab:red","DENSIFICATION_EXHAUSTION_FAILURE":"tab:blue","MIXED_FAILURE":"tab:purple","UNATTAINABLE_FIRST_STEP":"0.5"}
    markers={"SUCCESS":"o","GRAIN_GROWTH_FAILURE":"^","DENSIFICATION_EXHAUSTION_FAILURE":"v","MIXED_FAILURE":"X","UNATTAINABLE_FIRST_STEP":"x"}
    fig,axes=plt.subplots(1,3,figsize=(17,5),sharey=True)
    selected=[r for r in classes if r["growth_tolerance"]==.05]
    for axis,mode in zip(axes,GROWTH_MODES):
        group=[r for r in selected if r["growth_mode"]==mode]
        for category in colors:
            g=[r for r in group if r["classification"]==category];axis.scatter([r["G0_nm"] for r in g],[r["T2_C"] for r in g],c=colors[category],marker=markers[category],s=24,alpha=.55,label=category)
        axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=mode);axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 [C]");axes[-1].legend(fontsize=7);fig.tight_layout();fig.savefig(out/"chen_style_T2_vs_G_growth_modes_5pct.png",dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(17,5),sharey=True)
    for axis,mode in zip(axes,GROWTH_MODES):
        group=[r for r in boundaries if r["growth_mode"]==mode and r["growth_tolerance"]==.05]
        axis.scatter([r["G0_nm"] for r in group],[r["window_width_C"] for r in group],c=[r["T1_C"] for r in group],cmap="viridis",alpha=.65);axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=mode);axis.grid(alpha=.2)
    axes[0].set_ylabel("5% window width [C]");fig.tight_layout();fig.savefig(out/"window_width_vs_G_growth_modes.png",dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(17,5),sharey=True)
    for axis,mode in zip(axes,GROWTH_MODES):
        group=[r for r in boundaries if r["growth_mode"]==mode and r["growth_tolerance"]==.05]
        axis.scatter([r["G0_nm"] for r in group],[r["T_lower_density_C"] for r in group],label="lower density")
        axis.scatter([r["G0_nm"] for r in group],[r["T_upper_no_growth_C"] for r in group],label="upper growth")
        axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=mode);axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 boundary [C]");axes[-1].legend(fontsize=8);fig.tight_layout();fig.savefig(out/"T2_boundaries_growth_modes.png",dpi=150);plt.close(fig)
    fig,axis=plt.subplots(figsize=(8,5.5))
    for mode in GROWTH_MODES:
        g=[r for r in trajectories if r["growth_mode"]==mode and r["G0_nm"]==100 and r["T1_C"]==1300 and r["rho_switch"]==.80]
        axis.plot([r["T2_C"] for r in g],[r["initial_growth_mobility_factor"] for r in g],"o-",label=mode)
    axis.set(xlabel="T2 [C]",ylabel="Initial GB migration mobility factor",yscale="log");axis.grid(alpha=.2);axis.legend();fig.tight_layout();fig.savefig(out/"mobility_activation_vs_T2.png",dpi=150);plt.close(fig)


def main()->None:
    parser=argparse.ArgumentParser();parser.add_argument("--outdir",default="results/growth_mechanism_sensitivity");parser.add_argument("--workers",type=int,default=4);args=parser.parse_args()
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);trajectories,wall=run_map(args.workers);classes=classification_rows(trajectories);boundaries=boundary_rows(trajectories);recovery=baseline_recovery(trajectories)
    write_csv(out/"growth_mode_trajectories.csv",trajectories);write_csv(out/"growth_mode_classifications.csv",classes);write_csv(out/"growth_mode_boundaries.csv",boundaries);write_csv(out/"baseline_recovery.csv",recovery)
    write_csv(out/"runtime_summary.csv",[{"first_step_groups":270,"second_step_trajectories":len(trajectories),"classifications":len(classes),"wall_s":wall}]);make_plots(out,trajectories,classes,boundaries)
    print(f"DONE trajectories={len(trajectories)} classifications={len(classes)} wall_s={wall:.2f}")


if __name__=="__main__":main()
