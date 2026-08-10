#!/usr/bin/env python3
"""Bounded Chen-style ablation of observable pore/junction GB pinning."""
from __future__ import annotations
import argparse
from pathlib import Path

import growth_mechanism_sensitivity as study
from density_window_processing_map import write_csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODES=("baseline","junction_limited","pore_junction_pinning")


def diagnostic_plots(out:Path,trajectories:list[dict],classes:list[dict])->None:
    colors={"SUCCESS":"tab:green","GRAIN_GROWTH_FAILURE":"tab:red","DENSIFICATION_EXHAUSTION_FAILURE":"tab:blue","MIXED_FAILURE":"tab:purple","UNATTAINABLE_FIRST_STEP":"0.5"}
    markers={"SUCCESS":"o","GRAIN_GROWTH_FAILURE":"^","DENSIFICATION_EXHAUSTION_FAILURE":"v","MIXED_FAILURE":"X","UNATTAINABLE_FIRST_STEP":"x"}
    fig,axes=plt.subplots(1,3,figsize=(17,5),sharey=True)
    for axis,mode in zip(axes,MODES):
        group=[r for r in classes if r["growth_mode"]==mode and r["growth_tolerance"]==.10]
        for category in colors:
            rows=[r for r in group if r["classification"]==category]
            axis.scatter([r["G0_nm"] for r in rows],[r["T2_C"] for r in rows],c=colors[category],marker=markers[category],s=24,alpha=.55,label=category)
        axis.set(xscale="log",xlabel="G0 [nm]",title=mode);axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 [C]");axes[-1].legend(fontsize=7);fig.tight_layout();fig.savefig(out/"chen_style_T2_vs_G_pinning_10pct.png",dpi=150);plt.close(fig)
    rows=[r for r in trajectories if r["growth_mode"]=="pore_junction_pinning"]
    fig,axis=plt.subplots(figsize=(8,5.5))
    axis.scatter([r["initial_connected_pinning_coverage"] for r in rows],[r["initial_growth_mobility_factor"] for r in rows],s=12,alpha=.25,label="second-step start")
    axis.scatter([r["final_connected_coverage"] for r in rows],[r["final_growth_mobility_factor"] for r in rows],s=12,alpha=.25,label="second-step end")
    axis.set(xlabel="Connected pore-boundary coverage",ylabel="GB migration mobility factor",yscale="log")
    axis.grid(alpha=.2);axis.legend();fig.tight_layout();fig.savefig(out/"pinning_release_vs_connected_coverage.png",dpi=150);plt.close(fig)


def main()->None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",default="results/pore_junction_pinning")
    parser.add_argument("--workers",type=int,default=4)
    args=parser.parse_args()
    study.GROWTH_MODES=MODES
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    trajectories,wall=study.run_map(args.workers)
    classes=study.classification_rows(trajectories)
    boundaries=study.boundary_rows(trajectories)
    recovery=study.baseline_recovery(trajectories)
    write_csv(out/"pinning_trajectories.csv",trajectories)
    write_csv(out/"pinning_classifications.csv",classes)
    write_csv(out/"pinning_boundaries.csv",boundaries)
    write_csv(out/"baseline_recovery.csv",recovery)
    write_csv(out/"runtime_summary.csv",[{"first_step_groups":len(MODES)*len(study.G0_VALUES)*len(study.T1_VALUES)*len(study.SWITCHES),
        "second_step_trajectories":len(trajectories),"classifications":len(classes),"wall_s":wall}])
    study.make_plots(out,trajectories,classes,boundaries)
    diagnostic_plots(out,trajectories,classes)
    for old,new in {
        "chen_style_T2_vs_G_growth_modes_5pct.png":"chen_style_T2_vs_G_pinning_5pct.png",
        "window_width_vs_G_growth_modes.png":"pinning_window_width_vs_G.png",
        "T2_boundaries_growth_modes.png":"pinning_T2_boundaries.png",
        "mobility_activation_vs_T2.png":"pinning_mobility_activation_vs_T2.png",
    }.items():(out/old).rename(out/new)
    print(f"DONE trajectories={len(trajectories)} classifications={len(classes)} wall_s={wall:.2f}")


if __name__=="__main__":main()
