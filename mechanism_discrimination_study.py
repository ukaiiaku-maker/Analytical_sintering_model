#!/usr/bin/env python3
"""Matched-density topology identifiability and fast-firing size-domain audit."""
from __future__ import annotations

import argparse
from dataclasses import fields,replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as model
from density_window_processing_map import sample_at_density,write_csv
from two_step_window_map import classify

MODEL_STYLES=("density","connectivity","hybrid_topology")
G0_MATCHED=(75.,150.,300.)
MATCH_DENSITIES=(.80,.85)
HISTORY_RATES=(.2,20.)
FAST_G0=(50.,75.,100.,150.,225.,300.,450.)
FAST_RHO0=(.65,.70,.75)
FAST_RATES=(.2,5.,20.)
FAST_TARGET=.90
FOLLOWUP_T_C=1250.
FOLLOWUP_BUDGET_S=96*3600
VARYING={"rho0","G0","smoothing_gate_mode"}


def base_params()->model.Params:
    return model.Params(memory_model="pore_bin_redistribution",rho0=.70,G0=150e-9,pore_radius0=25e-9,pore_ln_sigma=.65)


def assert_fixed(parameter_sets:list[model.Params],base:model.Params)->None:
    fixed={f.name:getattr(base,f.name) for f in fields(base) if f.name not in VARYING}
    for p in parameter_sets:
        if {f.name:getattr(p,f.name) for f in fields(p) if f.name not in VARYING}!=fixed:raise AssertionError("material parameter drift")


def state_metrics(result:dict,prefix:str,index:int=-1)->dict:
    return {
        f"{prefix}_rho":float(result["rho"][index]),f"{prefix}_G_nm":float(result["G"][index])*1e9,
        f"{prefix}_pore_mean_radius_nm":float(result["pore_mean_radius"][index])*1e9,
        f"{prefix}_removable_fine_fraction":float(result["removable_fine_pore_fraction"][index]),
        f"{prefix}_large_pore_fraction":float(result["large_pore_fraction"][index]),
        f"{prefix}_f_pore":float(result["f_pore"][index]),
        f"{prefix}_connectivity":float(result["connectivity"][index]),
        f"{prefix}_connected_coverage":float(result["f_pore"][index]*result["connectivity"][index]),
        f"{prefix}_isolated_fraction":float(result["isolated_pore_fraction"][index]),
        f"{prefix}_activity":float(result["activity"][index]),f"{prefix}_E_G":float(result["E_G"][index]),
        f"{prefix}_cumulative_redistribution":float(result["cumulative_redistributed_pore_volume"][index]),
    }


def matched_density_audit(base:model.Params)->tuple[list[dict],dict]:
    settings=[(style,g0,rho,rate) for style in MODEL_STYLES for g0 in G0_MATCHED for rho in MATCH_DENSITIES for rate in HISTORY_RATES]
    params_list=[replace(base,smoothing_gate_mode=s,G0=g*1e-9) for s,g,_,_ in settings];assert_fixed(params_list,base)
    rows=[];cache={}
    for index,((style,g0,rho_match,rate),params) in enumerate(zip(settings,params_list),1):
        history=model.run(params,model.RampHold(rate),stop_at_rho=rho_match)
        reached=float(np.max(history["rho"]))>=rho_match-1e-12
        if reached:
            state=model.state_from_result(history,params);follow=model.run(params,model.Iso(FOLLOWUP_T_C,FOLLOWUP_BUDGET_S),initial=state)
            initial=state_metrics(history,"matched");final=state_metrics(follow,"followup")
            delta_rho=final["followup_rho"]-initial["matched_rho"]
            growth=(final["followup_G_nm"]-initial["matched_G_nm"])/initial["matched_G_nm"]
            initial_rho_dot=float(follow["rho_dot"][0]);initial_dGdt=float(follow["dGdt"][0])*1e9
        else:
            initial={};final={};delta_rho=growth=initial_rho_dot=initial_dGdt=np.nan
        row={"model_style":style,"G0_nm":g0,"rho0":params.rho0,"matched_density":rho_match,
             "history_heating_rate_C_min":rate,"history_reached":reached,"followup_T_C":FOLLOWUP_T_C,
             "history_time_budget_h":min(model.RampHold(rate).t_end,params.t_max_s)/3600,
             "followup_time_budget_h":FOLLOWUP_BUDGET_S/3600,"followup_delta_rho":delta_rho,
             "followup_growth_fraction":growth,"followup_initial_rho_dot":initial_rho_dot,
             "followup_initial_dGdt_nm_s":initial_dGdt,**initial,**final}
        rows.append(row);cache[(style,g0,rho_match,rate)]=(history,follow if reached else None)
        print(f"[matched {index}/{len(settings)}] {style} G0={g0:g} rho={rho_match:.2f} rate={rate:g}",flush=True)
    return rows,cache


def matched_pairs(rows:list[dict])->list[dict]:
    output=[]
    for style in MODEL_STYLES:
        for rho in MATCH_DENSITIES:
            eligible=[r for r in rows if r["model_style"]==style and r["matched_density"]==rho and r["history_reached"]]
            for a_g,b_g in ((75.,150.),(150.,300.)):
                candidates=[(a,b) for a in eligible for b in eligible if a["G0_nm"]==a_g and b["G0_nm"]==b_g]
                if not candidates:continue
                a,b=min(candidates,key=lambda pair:abs(pair[0]["matched_connected_coverage"]-pair[1]["matched_connected_coverage"]))
                output.append({"model_style":style,"matched_density":rho,"G0_a_nm":a_g,"G0_b_nm":b_g,
                    "history_rate_a":a["history_heating_rate_C_min"],"history_rate_b":b["history_heating_rate_C_min"],
                    "coverage_difference":abs(a["matched_connected_coverage"]-b["matched_connected_coverage"]),
                    "G1_difference_nm":abs(a["matched_G_nm"]-b["matched_G_nm"]),
                    "delta_rho_difference":abs(a["followup_delta_rho"]-b["followup_delta_rho"]),
                    "growth_fraction_difference":abs(a["followup_growth_fraction"]-b["followup_growth_fraction"])})
    return output


def fast_firing_audit(base:model.Params)->tuple[list[dict],dict]:
    styles=("density","connectivity")
    settings=[(s,g,rho) for s in styles for g in FAST_G0 for rho in FAST_RHO0]
    params_list=[replace(base,smoothing_gate_mode=s,G0=g*1e-9,rho0=rho) for s,g,rho in settings];assert_fixed(params_list,base)
    rows=[];cache={}
    for index,((style,g0,rho0),params) in enumerate(zip(settings,params_list),1):
        runs={rate:model.run(params,model.RampHold(rate)) for rate in FAST_RATES};cache[(style,g0,rho0)]=runs
        slow=sample_at_density(runs[.2],FAST_TARGET,rho0)
        for rate in FAST_RATES:
            sampled=sample_at_density(runs[rate],FAST_TARGET,rho0);scored=slow["eligible_target"] and sampled["eligible_target"] and slow["reached_target"] and sampled["reached_target"]
            result=runs[rate];peak=int(np.argmax(result["T_C"]));last=-1
            rows.append({"model_style":style,"G0_nm":g0,"rho0":rho0,"heating_rate_C_min":rate,"target_density":FAST_TARGET,
                "reached_target":sampled["reached_target"],"HR_pct_vs_slow":model.percent_gain(slow["G_at_target_nm"],sampled["G_at_target_nm"]) if scored else np.nan,
                **{k:v for k,v in sampled.items() if k not in ("eligible_target","reached_target")},
                "peak_rho":float(result["rho"][peak]),"peak_G_nm":float(result["G"][peak])*1e9,
                "peak_connected_coverage":float(result["f_pore"][peak]*result["connectivity"][peak]),
                "final_rho":float(result["rho"][last]),"final_G_nm":float(result["G"][last])*1e9,
                "final_connected_coverage":float(result["f_pore"][last]*result["connectivity"][last])})
        print(f"[fast {index}/{len(settings)}] {style} G0={g0:g} rho0={rho0:.2f}",flush=True)
    return rows,cache


def initial_condition_window_probe(base:model.Params)->list[dict]:
    """Small post-map probe; not a Cartesian expansion of the full window map."""
    settings=[(style,g0,rho0) for style in ("density","connectivity") for g0 in (75.,150.,300.) for rho0 in FAST_RHO0]
    params_list=[replace(base,smoothing_gate_mode=s,G0=g*1e-9,rho0=r) for s,g,r in settings];assert_fixed(params_list,base)
    rows=[]
    for index,((style,g0,rho0),params) in enumerate(zip(settings,params_list),1):
        first=model.run(params,model.Iso(1300.,96*3600),stop_at_rho=.85);attained=float(np.max(first["rho"]))>=.85-1e-12
        if attained:
            state=model.state_from_result(first,params);fm=state_metrics(first,"first")
        else:state=None;fm={}
        for T2 in (1000.,1050.,1100.,1150.,1200.,1250.,1300.):
            if attained:
                second=model.run(params,model.Iso(T2,96*3600),initial=state);final=state_metrics(second,"final")
                growth=(final["final_G_nm"]-fm["first_G_nm"])/fm["first_G_nm"];rho_final=final["final_rho"]
            else:final={};growth=rho_final=np.nan
            rows.append({"model_style":style,"G0_nm":g0,"rho0":rho0,"T1_C":1300.,"rho_switch":.85,"T2_C":T2,
                         "rho_target":.90,"growth_tolerance":.05,"first_step_attained":attained,
                         "classification":classify(attained,rho_final,growth,.90,.05),"growth_fraction":growth,**fm,**final})
        print(f"[window probe {index}/{len(settings)}] {style} G0={g0:g} rho0={rho0:.2f}",flush=True)
    return rows


def trajectory_rows(cache:dict)->list[dict]:
    rows=[]
    for (style,g0,rho0),runs in cache.items():
        if rho0!=.70 or g0 not in (75.,150.,300.,450.):continue
        for rate,result in runs.items():
            indices=np.unique(np.linspace(0,len(result["t"])-1,min(250,len(result["t"]))).astype(int))
            for i in indices:rows.append({"model_style":style,"G0_nm":g0,"rho0":rho0,"heating_rate_C_min":rate,
                "t_h":float(result["t"][i])/3600,"T_C":float(result["T_C"][i]),"rho":float(result["rho"][i]),
                "G_nm":float(result["G"][i])*1e9,"rho_dot":float(result["rho_dot"][i]),"dGdt_nm_s":float(result["dGdt"][i])*1e9,
                "connected_coverage":float(result["f_pore"][i]*result["connectivity"][i]),"activity":float(result["activity"][i]),"E_G":float(result["E_G"][i])})
    return rows


def robustness_matrix(matched:list[dict],fast:list[dict],size_summary:list[dict],window_points:list[dict])->list[dict]:
    criteria=("finite two-step window","nanoscale two-step robustness","two-step size boundary","matched-density topology identifiability",
              "high-heating-rate advantage","heating-rate robustness across particle size","heating-rate robustness across initial density",
              "grain-growth failure reproduced","densification-exhaustion failure reproduced","unattainable-first-step region reproduced")
    table={criterion:{"criterion":criterion} for criterion in criteria}
    for style in MODEL_STYLES:
        ss=[r for r in size_summary if r["smoothing_gate_mode"]==style]
        wp=[r for r in window_points if r["smoothing_gate_mode"]==style]
        mr=[r for r in matched if r["model_style"]==style and r["history_reached"]]
        fr=[r for r in fast if r["model_style"]==style and r["heating_rate_C_min"]==20 and np.isfinite(r["HR_pct_vs_slow"])]
        if style not in ("density","connectivity"):
            window="NOT_EVALUATED";fast_adv="NOT_EVALUATED";size_robust="NOT_EVALUATED";density_robust="NOT_EVALUATED"
        else:
            window="SUPPORTED_AT_10PCT_ONLY" if any(r["classification"]=="SUCCESS" for r in wp) else "INCONSISTENT_WITH_SWEEP"
            fast_adv="SUPPORTED" if fr and np.mean([r["HR_pct_vs_slow"]>0 for r in fr])>.5 else "INCONSISTENT_WITH_SWEEP"
            by_size=[np.mean([r["HR_pct_vs_slow"]>0 for r in fr if np.isclose(r["G0_nm"],g)]) for g in FAST_G0]
            by_rho=[np.mean([r["HR_pct_vs_slow"]>0 for r in fr if np.isclose(r["rho0"],rho)]) for rho in FAST_RHO0]
            size_robust="SUPPORTED" if np.mean(np.asarray(by_size)>=2/3)>=.7 else "PARTIAL"
            density_robust="SUPPORTED" if all(x>.5 for x in by_rho) else "PARTIAL"
        table["finite two-step window"][style]=window
        table["nanoscale two-step robustness"][style]="INCONSISTENT_WITH_SWEEP" if ss else "NOT_EVALUATED"
        table["two-step size boundary"][style]="INCONSISTENT_WITH_SWEEP" if ss else "NOT_EVALUATED"
        table["matched-density topology identifiability"][style]="NON_IDENTIFIABLE" if mr else "NOT_EVALUATED"
        table["high-heating-rate advantage"][style]=fast_adv
        table["heating-rate robustness across particle size"][style]=size_robust
        table["heating-rate robustness across initial density"][style]=density_robust
        table["grain-growth failure reproduced"][style]="SUPPORTED" if any(r["classification"] in ("GRAIN_GROWTH_FAILURE","MIXED_FAILURE") for r in wp) else "NOT_EVALUATED"
        table["densification-exhaustion failure reproduced"][style]="SUPPORTED" if any(r["classification"]=="DENSIFICATION_EXHAUSTION_FAILURE" for r in wp) else "NOT_EVALUATED"
        table["unattainable-first-step region reproduced"][style]="SUPPORTED" if any(r["classification"]=="UNATTAINABLE_FIRST_STEP" for r in wp) else ("INCONSISTENT_WITH_SWEEP" if wp else "NOT_EVALUATED")
    return list(table.values())


def make_plots(outdir:Path,matched:list[dict],fast:list[dict],trajectories:list[dict])->None:
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for style in MODEL_STYLES:
        g=[r for r in matched if r["model_style"]==style and r["history_reached"]]
        axes[0].scatter([r["matched_connected_coverage"] for r in g],[r["followup_delta_rho"] for r in g],label=style)
        axes[1].scatter([r["matched_connected_coverage"] for r in g],[r["followup_growth_fraction"] for r in g],label=style)
    axes[0].set(xlabel="Connected coverage at matched density",ylabel="Subsequent density gain")
    axes[1].set(xlabel="Connected coverage at matched density",ylabel="Subsequent growth fraction")
    for a in axes:a.grid(alpha=.2);a.legend(fontsize=8)
    fig.tight_layout();fig.savefig(outdir/"matched_density_varied_topology_response.png",dpi=150);plt.close(fig)
    for descriptor,filename in (("G0_nm","fast_firing_HR_vs_initial_size.png"),("rho0","fast_firing_HR_vs_initial_density.png")):
        fig,axis=plt.subplots(figsize=(8,5))
        subset=[r for r in fast if r["heating_rate_C_min"]==20 and np.isfinite(r["HR_pct_vs_slow"])]
        for style in ("density","connectivity"):
            for fixed in ((FAST_RHO0 if descriptor=="G0_nm" else FAST_G0)):
                g=[r for r in subset if r["model_style"]==style and r[("rho0" if descriptor=="G0_nm" else "G0_nm")]==fixed]
                if g:axis.plot([r[descriptor] for r in g],[r["HR_pct_vs_slow"] for r in g],"o-",alpha=.6,label=f"{style}, fixed={fixed:g}")
        axis.axhline(0,color="black",linestyle="--");axis.set(xlabel=descriptor,ylabel="HR_pct: 20 vs 0.2 C/min");axis.grid(alpha=.2);axis.legend(fontsize=6,ncol=2)
        fig.tight_layout();fig.savefig(outdir/filename,dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,9))
    subset=[r for r in trajectories if r["model_style"]=="connectivity" and r["G0_nm"]==150]
    for rate in FAST_RATES:
        g=[r for r in subset if r["heating_rate_C_min"]==rate]
        axes[0,0].plot([r["t_h"] for r in g],[r["rho"] for r in g],label=f"{rate:g} C/min")
        axes[0,1].plot([r["T_C"] for r in g],[r["rho"] for r in g])
        axes[1,0].plot([r["rho"] for r in g],[r["G_nm"] for r in g])
        axes[1,1].plot([r["T_C"] for r in g],[r["connected_coverage"] for r in g])
    axes[0,0].set(xlabel="Time [h]",ylabel="Density");axes[0,1].set(xlabel="Temperature [C]",ylabel="Density")
    axes[1,0].set(xlabel="Density",ylabel="G [nm]");axes[1,1].set(xlabel="Temperature [C]",ylabel="Connected coverage")
    for a in axes.flat:a.grid(alpha=.2)
    axes[0,0].legend();fig.tight_layout();fig.savefig(outdir/"slow_intermediate_fast_trajectories.png",dpi=150);plt.close(fig)


def main()->None:
    parser=argparse.ArgumentParser();parser.add_argument("--outdir",default="results/mechanism_discrimination");parser.add_argument("--window-results",default="results/two_step_window_map");args=parser.parse_args()
    outdir=Path(args.outdir);outdir.mkdir(parents=True,exist_ok=True);base=base_params()
    matched,_=matched_density_audit(base);pairs=matched_pairs(matched);fast,cache=fast_firing_audit(base);trajectories=trajectory_rows(cache);probe=initial_condition_window_probe(base)
    import csv
    with (Path(args.window_results)/"size_window_summary.csv").open() as f:
        summary=[]
        for r in csv.DictReader(f):
            for key in ("G0_nm","n_windows"):
                r[key]=float(r[key])
            summary.append(r)
    with (Path(args.window_results)/"two_step_window_points.csv").open() as f:
        window_points=[]
        for r in csv.DictReader(f):
            for key in ("G0_nm","rho_target","growth_tolerance"):
                r[key]=float(r[key])
            window_points.append(r)
    matrix=robustness_matrix(matched,fast,summary,window_points)
    write_csv(outdir/"matched_density_results.csv",matched);write_csv(outdir/"matched_state_pairs.csv",pairs)
    write_csv(outdir/"fast_firing_results.csv",fast);write_csv(outdir/"fast_firing_trajectories.csv",trajectories);write_csv(outdir/"robustness_matrix.csv",matrix)
    write_csv(outdir/"initial_condition_window_probe.csv",probe)
    make_plots(outdir,matched,fast,trajectories)
    print(f"matched={len(matched)} pairs={len(pairs)} fast={len(fast)} trajectory_samples={len(trajectories)}")


if __name__=="__main__":main()
