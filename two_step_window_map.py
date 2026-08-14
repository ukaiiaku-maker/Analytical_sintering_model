#!/usr/bin/env python3
"""Bounded Chen/Wang-style second-step process-window mechanism audit."""
from __future__ import annotations

import argparse
from dataclasses import fields, replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from density_window_processing_map import write_csv

G0_NM=(50.,75.,100.,150.,225.,300.,450.)
GATE_MODES=("density","connectivity")
T1_VALUES=(1250.,1300.,1350.)
RHO_SWITCHES=(.75,.80,.85,.88)
T2_VALUES=(1000.,1050.,1100.,1150.,1200.,1250.,1300.)
RHO_TARGETS=(.90,.95,.98)
GROWTH_TOLERANCES=(.02,.05,.10)
RHO0=.70
FIRST_STEP_BUDGET_S=96*3600
SECOND_STEP_BUDGET_S=96*3600
CLASSIFICATIONS=("SUCCESS","GRAIN_GROWTH_FAILURE","DENSIFICATION_EXHAUSTION_FAILURE","MIXED_FAILURE","UNATTAINABLE_FIRST_STEP")
VARYING={"G0","smoothing_gate_mode"}


def base_params()->model.Params:
    return model.Params(memory_model="pore_bin_redistribution",rho0=RHO0,G0=150e-9,
                        pore_radius0=25e-9,pore_ln_sigma=.65)


def classify(first_attained:bool,rho_final:float,growth_fraction:float,target:float,tolerance:float)->str:
    if not first_attained:return "UNATTAINABLE_FIRST_STEP"
    dense=rho_final>=target-1e-12; bounded=growth_fraction<=tolerance+1e-12
    if dense and bounded:return "SUCCESS"
    if dense:return "GRAIN_GROWTH_FAILURE"
    if bounded:return "DENSIFICATION_EXHAUSTION_FAILURE"
    return "MIXED_FAILURE"


def assert_fixed_parameters(parameter_sets:list[model.Params],base:model.Params)->None:
    reference={f.name:getattr(base,f.name) for f in fields(base) if f.name not in VARYING}
    for params in parameter_sets:
        current={f.name:getattr(params,f.name) for f in fields(params) if f.name not in VARYING}
        if current!=reference:raise AssertionError("material parameter drift")


def last_diagnostics(result:dict)->dict:
    i=-1
    return {
        "rho_final":float(result["rho"][i]),"G_final_nm":float(result["G"][i])*1e9,
        "second_step_initial_rho_dot":float(result["rho_dot"][0]),
        "second_step_initial_dGdt_nm_s":float(result["dGdt"][0])*1e9,
        "second_step_median_rho_dot":float(np.median(result["rho_dot"])),
        "second_step_median_dGdt_nm_s":float(np.median(result["dGdt"]))*1e9,
        "second_step_median_E_G":float(np.median(result["E_G"])),
        "pore_mean_radius_nm":float(result["pore_mean_radius"][i])*1e9,
        "large_pore_fraction":float(result["large_pore_fraction"][i]),
        "removable_fine_pore_fraction":float(result["removable_fine_pore_fraction"][i]),
        "f_pore":float(result["f_pore"][i]),
        "connected_coverage":float(result["f_pore"][i]*result["connectivity"][i]),
        "isolated_pore_fraction":float(result["isolated_pore_fraction"][i]),
        "renewal_activity":float(result["activity"][i]),
    }


def first_diagnostics(result:dict)->dict:
    i=-1
    return {
        "rho1":float(result["rho"][i]),"G1_nm":float(result["G"][i])*1e9,
        "first_step_time_h":float(result["t"][i])/3600,
        "first_pore_mean_radius_nm":float(result["pore_mean_radius"][i])*1e9,
        "first_large_pore_fraction":float(result["large_pore_fraction"][i]),
        "first_removable_fine_pore_fraction":float(result["removable_fine_pore_fraction"][i]),
        "first_f_pore":float(result["f_pore"][i]),
        "first_connectivity":float(result["connectivity"][i]),
        "first_connected_coverage":float(result["f_pore"][i]*result["connectivity"][i]),
        "first_isolated_pore_fraction":float(result["isolated_pore_fraction"][i]),
        "first_activity":float(result["activity"][i]),
        "first_E_G":float(result["E_G"][i]),
    }


def run_map()->list[dict]:
    base=base_params(); parameter_sets=[replace(base,G0=g*1e-9,smoothing_gate_mode=mode) for mode in GATE_MODES for g in G0_NM]
    assert_fixed_parameters(parameter_sets,base)
    rows=[]; total=len(parameter_sets)*len(T1_VALUES)*len(RHO_SWITCHES); done=0
    for params in parameter_sets:
        for T1 in T1_VALUES:
            for switch in RHO_SWITCHES:
                first_protocol=model.Iso(T1,FIRST_STEP_BUDGET_S)
                first=model.run(params,first_protocol,stop_at_rho=switch)
                attained=float(np.max(first["rho"]))>=switch-1e-12
                fd=first_diagnostics(first)
                for T2 in T2_VALUES:
                    if attained:
                        state=model.state_from_result(first,params)
                        second=model.run(params,model.Iso(T2,SECOND_STEP_BUDGET_S),initial=state)
                        sd=last_diagnostics(second);growth=(sd["G_final_nm"]-fd["G1_nm"])/fd["G1_nm"]
                    else:
                        sd={key:np.nan for key in ("rho_final","G_final_nm","second_step_initial_rho_dot","second_step_initial_dGdt_nm_s","second_step_median_rho_dot","second_step_median_dGdt_nm_s","second_step_median_E_G","pore_mean_radius_nm","large_pore_fraction","removable_fine_pore_fraction","f_pore","connected_coverage","isolated_pore_fraction","renewal_activity")};growth=np.nan
                    for target in RHO_TARGETS:
                        for tolerance in GROWTH_TOLERANCES:
                            category=classify(attained,sd["rho_final"],growth,target,tolerance)
                            rows.append({
                                "smoothing_gate_mode":params.smoothing_gate_mode,"G0_nm":params.G0*1e9,
                                "rho0":params.rho0,"T1_C":T1,"rho_switch":switch,"T2_C":T2,
                                "rho_target":target,"growth_tolerance":tolerance,
                                "first_step_attained":attained,"classification":category,
                                "growth_fraction":growth,"first_step_budget_h":FIRST_STEP_BUDGET_S/3600,
                                "second_step_budget_h":SECOND_STEP_BUDGET_S/3600,**fd,**sd,
                            })
                done+=1;print(f"[{done}/{total}] {params.smoothing_gate_mode} G0={params.G0*1e9:g} T1={T1:g} switch={switch:.2f}",flush=True)
    return rows


def boundary_rows(points:list[dict])->list[dict]:
    rows=[]
    keys={(r["smoothing_gate_mode"],r["G0_nm"],r["T1_C"],r["rho_switch"],r["rho_target"],r["growth_tolerance"]) for r in points}
    for key in sorted(keys):
        mode,g0,T1,switch,target,tol=key;group=[r for r in points if (r["smoothing_gate_mode"],r["G0_nm"],r["T1_C"],r["rho_switch"],r["rho_target"],r["growth_tolerance"])==key]
        evaluated=[r for r in group if r["first_step_attained"]]
        density_t=[r["T2_C"] for r in evaluated if r["rho_final"]>=target-1e-12]
        growth_t=[r["T2_C"] for r in evaluated if r["growth_fraction"]<=tol+1e-12]
        success_t=[r["T2_C"] for r in evaluated if r["classification"]=="SUCCESS"]
        rows.append({
            "smoothing_gate_mode":mode,"G0_nm":g0,"T1_C":T1,"rho_switch":switch,"rho_target":target,"growth_tolerance":tol,
            "first_step_attained":bool(evaluated),"G1_nm":group[0]["G1_nm"],"rho1":group[0]["rho1"],
            "first_connected_coverage":group[0]["first_connected_coverage"],
            "T_lower_density_C":min(density_t) if density_t else np.nan,
            "T_upper_no_growth_C":max(growth_t) if growth_t else np.nan,
            "T_success_min_C":min(success_t) if success_t else np.nan,
            "T_success_max_C":max(success_t) if success_t else np.nan,
            "window_width_C":max(success_t)-min(success_t) if success_t else np.nan,
            "n_success":len(success_t),
        })
    return rows


def size_summary(points:list[dict],boundaries:list[dict])->list[dict]:
    rows=[]
    for mode in GATE_MODES:
        for g0 in G0_NM:
            b=[r for r in boundaries if r["smoothing_gate_mode"]==mode and np.isclose(r["G0_nm"],g0) and r["rho_target"]==.90 and r["growth_tolerance"]==.05]
            successes=[r for r in b if r["n_success"]>0]; widths=[r["window_width_C"] for r in successes]
            first_ok=any(r["first_step_attained"] for r in b)
            p=[r for r in points if r["smoothing_gate_mode"]==mode and np.isclose(r["G0_nm"],g0) and r["rho_target"]==.90 and r["growth_tolerance"]==.05]
            if not first_ok: label="FIRST_STEP_UNATTAINABLE"
            elif successes and sum(r["n_success"]>=3 for r in successes)>=2 and max(widths)>=100:label="ROBUST_TWO_STEP_WINDOW"
            elif successes:label="NARROW_TWO_STEP_WINDOW"
            else:
                counts={c:sum(r["classification"]==c for r in p) for c in CLASSIFICATIONS}
                if counts["GRAIN_GROWTH_FAILURE"]+counts["MIXED_FAILURE"]>=max(1,counts["DENSIFICATION_EXHAUSTION_FAILURE"]):label="NO_TWO_STEP_WINDOW_GRAIN_GROWTH"
                elif np.nanmedian([r["first_connected_coverage"] for r in p])<.1:label="NO_TWO_STEP_WINDOW_TOPOLOGY"
                else:label="NO_TWO_STEP_WINDOW_DENSIFICATION_EXHAUSTION"
            rows.append({"smoothing_gate_mode":mode,"G0_nm":g0,"size_window_classification":label,
                         "n_first_step_states":sum(r["first_step_attained"] for r in b),"n_windows":len(successes),
                         "max_window_width_C":max(widths) if widths else np.nan,
                         "median_window_center_C":np.median([(r["T_success_min_C"]+r["T_success_max_C"])/2 for r in successes]) if successes else np.nan})
    return rows


def make_plots(outdir:Path,points:list[dict],boundaries:list[dict])->None:
    selected=[r for r in points if r["rho_target"]==.90 and r["growth_tolerance"]==.05]
    colors={"SUCCESS":"tab:green","GRAIN_GROWTH_FAILURE":"tab:red","DENSIFICATION_EXHAUSTION_FAILURE":"tab:blue","MIXED_FAILURE":"tab:purple","UNATTAINABLE_FIRST_STEP":"0.6"}
    markers={"SUCCESS":"o","GRAIN_GROWTH_FAILURE":"^","DENSIFICATION_EXHAUSTION_FAILURE":"v","MIXED_FAILURE":"X","UNATTAINABLE_FIRST_STEP":"x"}
    for xkey,xlabel,filename in (("G1_nm","G1 after first step [nm]","chen_style_T2_vs_G1_classification.png"),("rho1","rho1","T2_vs_rho1_classification.png")):
        fig,axes=plt.subplots(1,2,figsize=(13,5),sharey=True)
        for axis,mode in zip(axes,GATE_MODES):
            group=[r for r in selected if r["smoothing_gate_mode"]==mode]
            for category in CLASSIFICATIONS:
                g=[r for r in group if r["classification"]==category]
                axis.scatter([r[xkey] for r in g],[r["T2_C"] for r in g],c=colors[category],marker=markers[category],s=28,alpha=.65,label=category)
            axis.set(xlabel=xlabel,title=mode);axis.grid(alpha=.2)
        axes[0].set_ylabel("T2 [C]");axes[1].legend(fontsize=7,loc="best");fig.tight_layout();fig.savefig(outdir/filename,dpi=150);plt.close(fig)
    b=[r for r in boundaries if r["rho_target"]==.90 and r["growth_tolerance"]==.05 and r["first_step_attained"]]
    fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
    for axis,mode in zip(axes,GATE_MODES):
        g=[r for r in b if r["smoothing_gate_mode"]==mode]
        axis.scatter([r["G1_nm"] for r in g],[r["T_lower_density_C"] for r in g],label="lowest T2 reaching density")
        axis.scatter([r["G1_nm"] for r in g],[r["T_upper_no_growth_C"] for r in g],label="highest T2 within 5% growth")
        axis.set(xlabel="G1 [nm]",title=mode);axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 [C]");axes[1].legend(fontsize=8);fig.tight_layout();fig.savefig(outdir/"window_boundaries_vs_G1.png",dpi=150);plt.close(fig)
    relaxed=[r for r in points if r["rho_target"]==.90 and r["growth_tolerance"]==.10]
    fig,axes=plt.subplots(1,2,figsize=(13,5),sharey=True)
    for axis,mode in zip(axes,GATE_MODES):
        group=[r for r in relaxed if r["smoothing_gate_mode"]==mode]
        for category in CLASSIFICATIONS:
            g=[r for r in group if r["classification"]==category]
            axis.scatter([r["G1_nm"] for r in g],[r["T2_C"] for r in g],c=colors[category],marker=markers[category],s=28,alpha=.65,label=category)
        axis.set(xlabel="G1 after first step [nm]",title=f"{mode}, 10% tolerance");axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 [C]");axes[1].legend(fontsize=7);fig.tight_layout();fig.savefig(outdir/"chen_style_T2_vs_G1_classification_growth10.png",dpi=150);plt.close(fig)
    for zkey,label,filename,cmap in (("growth_fraction","Second-step grain growth fraction","grain_growth_fraction_map.png","magma"),("rho_final","Final density","density_attainment_map.png","viridis")):
        fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
        for axis,mode in zip(axes,GATE_MODES):
            g=[r for r in selected if r["smoothing_gate_mode"]==mode and r["first_step_attained"]]
            im=axis.scatter([r["G1_nm"] for r in g],[r["T2_C"] for r in g],c=[r[zkey] for r in g],cmap=cmap,s=32)
            axis.set(xlabel="G1 [nm]",title=mode);axis.grid(alpha=.2)
        axes[0].set_ylabel("T2 [C]");fig.colorbar(im,ax=axes,label=label);fig.subplots_adjust(right=.9,wspace=.15);fig.savefig(outdir/filename,dpi=150);plt.close(fig)


def main()->None:
    parser=argparse.ArgumentParser();parser.add_argument("--outdir",default="results/two_step_window_map");args=parser.parse_args()
    outdir=Path(args.outdir);outdir.mkdir(parents=True,exist_ok=True)
    points=run_map();boundaries=boundary_rows(points);summary=size_summary(points,boundaries)
    failures=[r for r in points if r["classification"]!="SUCCESS"]
    write_csv(outdir/"two_step_window_points.csv",points);write_csv(outdir/"window_boundaries.csv",boundaries);write_csv(outdir/"failure_modes.csv",failures);write_csv(outdir/"size_window_summary.csv",summary)
    make_plots(outdir,points,boundaries)
    print(f"points={len(points)} boundaries={len(boundaries)} failures={len(failures)}")


if __name__=="__main__":main()
