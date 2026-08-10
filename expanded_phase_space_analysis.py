#!/usr/bin/env python3
"""Boundary-focused analysis products for the expanded fixed-model campaign."""
from __future__ import annotations
import csv,math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from density_window_processing_map import write_csv

STYLES=("density","connectivity","hybrid_topology")


def load(path:Path)->list[dict]:
    rows=[]
    with path.open() as stream:
        for row in csv.DictReader(stream):
            for key,value in list(row.items()):
                if value=="True":row[key]=True
                elif value=="False":row[key]=False
                elif key not in ("design","model_style","classification_5pct"):
                    try:row[key]=float(value)
                    except ValueError:pass
            rows.append(row)
    return rows


def finite(value)->bool:return isinstance(value,(float,int)) and np.isfinite(value)


def matched_pairing(rows:list[dict])->list[dict]:
    output=[]
    for style in STYLES:
        for density in (.80,.85):
            group=[r for r in rows if r["model_style"]==style and r["matched_density"]==density and r["history_reached"]]
            pairs=[(a,b) for i,a in enumerate(group) for b in group[i+1:] if a["G0_nm"]!=b["G0_nm"]]
            scores={
                "max_topology_contrast":lambda p:-abs(p[0]["first_connected_coverage"]-p[1]["first_connected_coverage"]),
                "nearest_rho_connectivity":lambda p:abs(p[0]["first_connected_coverage"]-p[1]["first_connected_coverage"]),
                "nearest_rho_connectivity_G1":lambda p:abs(p[0]["first_connected_coverage"]-p[1]["first_connected_coverage"])/.1+abs(p[0]["G1_nm"]-p[1]["G1_nm"])/100,
            }
            for tier,score in scores.items():
                if not pairs:continue
                a,b=min(pairs,key=score)
                output.append({"model_style":style,"matched_density":density,"matching_tier":tier,
                    "G0_a_nm":a["G0_nm"],"G0_b_nm":b["G0_nm"],"rate_a":a["history_rate_C_min"],"rate_b":b["history_rate_C_min"],
                    "coverage_difference":abs(a["first_connected_coverage"]-b["first_connected_coverage"]),"G1_difference_nm":abs(a["G1_nm"]-b["G1_nm"]),
                    "rho2_difference":abs(a["rho2"]-b["rho2"]),"G2_difference_nm":abs(a["G2_nm"]-b["G2_nm"])})
    return output


def size_summary(boundaries:list[dict])->list[dict]:
    output=[]
    for style in STYLES:
        for g0 in sorted({r["G0_nm"] for r in boundaries if r["model_style"]==style}):
            group=[r for r in boundaries if r["model_style"]==style and r["G0_nm"]==g0 and r["rho_target"]==.90 and r["growth_tolerance"]==.05 and r.get("eligible_second_step_target",True)]
            wins=[r for r in group if finite(r["window_width_C"])];required=[r["minimum_required_growth_tolerance"] for r in group if finite(r["minimum_required_growth_tolerance"])]
            output.append({"model_style":style,"G0_nm":g0,"n_first_states":len(group),"n_5pct_windows":len(wins),
                "max_5pct_window_width_C":max((r["window_width_C"] for r in wins),default=np.nan),
                "minimum_required_growth_tolerance":min(required,default=np.nan),
                "minimum_T_lower_C":min((r["T_lower_density_C"] for r in group if finite(r["T_lower_density_C"])),default=np.nan),
                "maximum_T_upper_5pct_C":max((r["T_upper_no_growth_C"] for r in group if finite(r["T_upper_no_growth_C"])),default=np.nan)})
    return output


def fast_boundaries(fast:list[dict])->list[dict]:
    output=[]
    for style in STYLES:
        for rho0 in (.60,.70,.80):
            for peak in (1300.,1400.,1450.):
                for rate in (.5,1.,2.,5.,10.,20.,50.,100.):
                    group=[r for r in fast if r["model_style"]==style and r["rho0"]==rho0 and r["peak_target_C"]==peak and r["heating_rate_C_min"]==rate and finite(r["HR_pct_vs_0p2"])]
                    positive=sorted(r["G0_nm"] for r in group if r["HR_pct_vs_0p2"]>0)
                    output.append({"model_style":style,"rho0":rho0,"peak_C":peak,"heating_rate_C_min":rate,"n_scored":len(group),"n_positive":len(positive),
                        "lowest_positive_G0_nm":min(positive,default=np.nan),"highest_positive_G0_nm":max(positive,default=np.nan),
                        "max_HR_pct":max((r["HR_pct_vs_0p2"] for r in group),default=np.nan)})
    return output


def mechanism_assessment()->list[dict]:
    return [
        {"mechanism":"fractional grain-growth scaling controls nanoscale penalty","assessment":"SUPPORTED"},
        {"mechanism":"topology independently creates large-size onset","assessment":"INCONSISTENT_WITH_SWEEP"},
        {"mechanism":"topology moves T2 boundary surfaces","assessment":"SUPPORTED"},
        {"mechanism":"initial density controls fundamental size threshold","assessment":"INCONSISTENT_WITH_SWEEP"},
        {"mechanism":"T1/switch history modulates window width/location","assessment":"SUPPORTED"},
        {"mechanism":"closed-pore/late-stage physics above rho=0.92","assessment":"REQUIRED"},
        {"mechanism":"topology has independent matched-state predictive content","assessment":"NON_IDENTIFIABLE"},
        {"mechanism":"fast firing is monotonic across size","assessment":"INCONSISTENT_WITH_SWEEP"},
        {"mechanism":"fast-firing benefit saturates by high rate","assessment":"SUPPORTED"},
        {"mechanism":"rapid T1 heating universally creates two-step eligibility","assessment":"INCONSISTENT_WITH_SWEEP"},
    ]


def plots(out:Path,boundaries:list[dict],fast:list[dict],combined:list[dict])->None:
    extension=[r for r in boundaries if r["design"]=="T1_switch_extension" and r["rho_target"]==.90 and r["growth_tolerance"]==.05 and r.get("eligible_second_step_target",True)]
    fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
    for axis,style in zip(axes,("density","connectivity")):
        g=[r for r in extension if r["model_style"]==style];im=axis.scatter([r["G0_nm"] for r in g],[r["rho_switch"] for r in g],c=[r["window_width_C"] if finite(r["window_width_C"]) else -10 for r in g],cmap="viridis",vmin=-10,vmax=50)
        axis.set_xscale("log");axis.set(xlabel="G0 [nm]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("rho switch");fig.colorbar(im,ax=axes,label="5% window width [C]; -10=no window");fig.subplots_adjust(right=.9,wspace=.15);fig.savefig(out/"G0_rho_switch_success_region.png",dpi=150);plt.close(fig)
    fig,axis=plt.subplots(figsize=(8,5))
    for style in STYLES:
        xs=[];reach=[];success=[]
        for target in (.88,.90,.92,.94,.95,.96,.98):
            g=[r for r in boundaries if r["design"]=="core_size" and r["model_style"]==style and r["rho_target"]==target and r["growth_tolerance"]==.05 and r.get("eligible_second_step_target",True)]
            xs.append(target);reach.append(np.mean([finite(r["T_lower_density_C"]) for r in g]));success.append(np.mean([finite(r["window_width_C"]) for r in g]))
        axis.plot(xs,reach,"o-",label=f"{style}: density reached");axis.plot(xs,success,"s--",label=f"{style}: 5% success")
    axis.set(xlabel="rho target",ylabel="Fraction of eligible first states");axis.grid(alpha=.2);axis.legend(fontsize=7,ncol=2);fig.tight_layout();fig.savefig(out/"rho_target_accessibility.png",dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for style in STYLES:
        g=[r for r in boundaries if r["design"]=="core_size" and r["model_style"]==style and r["T1_C"]==1300 and r["rho_switch"]==.825 and r["rho_target"]==.90 and r["growth_tolerance"]==.05]
        axes[0].plot([r["G0_nm"] for r in g],[r["T_lower_density_C"] for r in g],"o-",label=style)
        axes[1].plot([r["G0_nm"] for r in g],[r["T_upper_no_growth_C"] for r in g],"o-",label=style)
    for axis,title in zip(axes,("Topology displacement: lower density boundary","Topology displacement: upper 5% boundary")):
        axis.set_xscale("log");axis.set(xlabel="G0 [nm]",ylabel="T2 [C]",title=title);axis.grid(alpha=.2);axis.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/"topology_gate_boundary_displacement.png",dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharey=True)
    for axis,style in zip(axes,STYLES):
        g=[r for r in fast if r["model_style"]==style and r["rho0"]==.70 and r["G0_nm"]==150]
        im=axis.scatter([r["heating_rate_C_min"] for r in g],[r["peak_target_C"] for r in g],c=[r["final_rho"] for r in g],cmap="viridis");axis.set_xscale("log");axis.set(xlabel="Heating rate [C/min]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("Peak target [C]");fig.colorbar(im,ax=axes,label="Final density");fig.subplots_adjust(right=.9,wspace=.15);fig.savefig(out/"heating_rate_peak_density_surface.png",dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
    for axis,style in zip(axes,("density","connectivity")):
        for g0 in (75.,800.,2000.):
            g=[r for r in combined if r["model_style"]==style and r["G0_nm"]==g0]
            axis.scatter([r["T1_heating_rate_C_min"] for r in g],[r["T2_C"] for r in g],c=[1 if r["classification_5pct"]=="SUCCESS" else 0 for r in g],cmap="RdYlGn",vmin=0,vmax=1,label=f"G0={g0:g}",alpha=.7)
        axis.set_xscale("log");axis.set(xlabel="T1 heating rate [C/min]",title=style);axis.grid(alpha=.2)
    axes[0].set_ylabel("T2 [C]");axes[1].legend(fontsize=7);fig.tight_layout();fig.savefig(out/"rapid_T1_two_step_eligibility.png",dpi=150);plt.close(fig)


def main()->None:
    out=Path("results/expanded_phase_space")
    boundaries=load(out/"two_step_boundaries.csv")+load(out/"upper_size_extension_boundaries.csv")+load(out/"size_onset_refinement_boundaries.csv")
    fast=load(out/"fast_firing_surface.csv");matched=load(out/"matched_history_expanded.csv");combined=load(out/"combined_history.csv")
    write_csv(out/"expanded_size_summary.csv",size_summary(boundaries));write_csv(out/"matched_history_pairing.csv",matched_pairing(matched))
    write_csv(out/"fast_sign_boundaries.csv",fast_boundaries(fast));write_csv(out/"mechanism_assessment.csv",mechanism_assessment())
    plots(out,boundaries,fast,combined);print("expanded analysis complete")


if __name__=="__main__":main()
