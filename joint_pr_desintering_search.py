#!/usr/bin/env python3
"""Bounded joint screen for local PR/de-sintering early-stage memory."""
from __future__ import annotations
import argparse, csv, inspect, math, time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import adaptive_T2_boundary_search as adaptive
import agentic_mechanism_search as classification
import preparation_window_search as preparation
import production_mechanism_assessment as production
import pr_desintering_memory_model as memory
import topology_constrained_sintering as aggregate

TARGET = .90
BUDGET = 96 * 3600
TOLS = (.05, .10)
G_CHEN = (150., 225., 300.)
T1S = (1350., 1400., 1450.)
SWITCHES = (.80, .84, .88)
G_FAST_REDUCED = (75., 150., 300.)
RHO0_REDUCED = (.65, .70, .75)
RATES_REDUCED = (1., 20., 100.)
PEAKS_REDUCED = (1350., 1400., 1450.)
HOLDS_REDUCED = (2., 8., 20.)
TARGETS = (.85, .88, .90, .92)


def write(path, records, empty=("candidate_id", "status")):
    fields = []
    for row in records:
        for key in row:
            if key not in fields:
                fields.append(key)
    if not fields:
        fields = list(empty)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(records)


def frozen_base():
    result = {}
    for mid, p in preparation.mechanisms().items():
        b = replace(p.action.location.base, dt_max_s=900.)
        result[mid] = replace(p, action=replace(p.action, location=replace(p.action.location, base=b)))
    return result


def design():
    variants = (
        ("disabled", dict(early_memory_mode="disabled")),
        ("PR_comp", dict(early_memory_mode="PR_desintering_competition", k_PR_ref_s=2e-4)),
        ("PR_attrition_moderate", dict(early_memory_mode="PR_plus_connected_fine_attrition", k_PR_ref_s=2e-4)),
        ("PR_attrition_strong", dict(early_memory_mode="PR_plus_connected_fine_attrition", k_PR_ref_s=5e-4)),
    )
    return [(f"{mid}_{tag}", memory.PRMemoryParams(p, **kw), mid, tag)
            for mid, p in frozen_base().items() for tag, kw in variants]


def fractions(h, i):
    values = [float(np.sum(h[k][i])) for k in ("phi_GBseg", "phi_TJ", "phi_iso")]
    z = max(sum(values), 1e-300)
    return tuple(v / z for v in values)


def state_metrics(h, i):
    radii = h["pore_radii"]
    connected = h["phi_GBseg"][i] + h["phi_TJ"][i]
    z = max(float(np.sum(connected)), 1e-300)
    rref = radii[len(radii)//2]
    fine = float(np.sum(connected[radii <= rref])) / z
    return dict(connected_fine_pore_fraction=float(h.get("connected_fine_pore_fraction", [fine])[i]) if "connected_fine_pore_fraction" in h else fine,
                pore_mean_radius_nm=float(np.sum(connected * radii) / z) * 1e9,
                large_pore_fraction=float(np.sum(connected[radii > rref])) / z,
                cumulative_PR_desintering_work=float(h.get("cumulative_PR_desintering_work", [0])[i]) if "cumulative_PR_desintering_work" in h else 0.)


def adaptive_points(cid, p, G, T1, sw, state):
    cache = {}
    def point(T2, stage):
        if T2 in cache:
            return cache[T2]
        h = memory.run(p, aggregate.Iso(T2, BUDGET), initial=state)
        rho2 = float(h["rho"][-1]); G1 = state.base.pore.G * 1e9; G2 = float(h["G"][-1]) * 1e9
        cache[T2] = dict(candidate_id=cid, G0_nm=G, T1_C=T1, rho_switch=sw, T2_C=T2,
                         scan_stage=stage, rho1=state.base.pore.rho, G1_nm=G1, rho2=rho2,
                         G2_nm=G2, growth_fraction=(G2-G1)/G1, **state_metrics(h, -1))
        return cache[T2]
    points = [point(T, "coarse") for T in adaptive.COARSE]
    if points[0]["rho2"] >= TARGET - 1e-12:
        for T in reversed(adaptive.DOWN):
            q = point(T, "downward_extension"); points.append(q)
            if q["rho2"] < TARGET - 1e-12: break
    for T in adaptive.UP:
        if not any(adaptive.status(points, tol, T1).get("n_success", 0) and
                   not adaptive.status(points, tol, T1).get("upper_bracketed", False) for tol in TOLS):
            break
        points.append(point(T, "upward_extension"))
    intervals = set()
    for tol in TOLS: intervals.update(adaptive.intervals_for_refinement(points, tol))
    for lo, hi in intervals:
        for T in np.arange(lo + 10, hi, 10): point(float(T), "local_refinement")
    return list(cache.values())


def chen_screen(candidates):
    boundaries = []
    for n, (cid, p, mid, tag) in enumerate(candidates, 1):
        print("chen", n, len(candidates), cid, flush=True)
        for G in G_CHEN:
            basep = p.base
            b = replace(basep.action.location.base, G0=G*1e-9)
            pp = replace(p, base=replace(basep, action=replace(basep.action, location=replace(basep.action.location, base=b))))
            for T1 in T1S:
                for sw in SWITCHES:
                    h1 = memory.run(pp, aggregate.Iso(T1, BUDGET), stop_at_rho=sw)
                    if float(h1["rho"][-1]) < sw - 1e-12:
                        for tol in TOLS:
                            boundaries.append(dict(candidate_id=cid, base_mechanism=mid, variant=tag,
                                G0_nm=G, T1_C=T1, rho_switch=sw, growth_tolerance=tol,
                                boundary_status="UNATTAINABLE_FIRST_STEP", complete_practical=False))
                        continue
                    state = memory.final_state(h1, pp)
                    growth1 = (state.base.pore.G*1e9-G)/G
                    pts = adaptive_points(cid, pp, G, T1, sw, state)
                    for tol in TOLS:
                        status = adaptive.status(pts, tol, T1, practical=True)
                        complete = status["boundary_status"] == "COMPLETE_WINDOW" and growth1 <= tol + 1e-12
                        boundaries.append(dict(candidate_id=cid, base_mechanism=mid, variant=tag,
                            G0_nm=G, T1_C=T1, rho_switch=sw, growth_tolerance=tol,
                            first_step_growth_fraction=growth1, complete_practical=complete, **status))
    return boundaries


def _fast_task(task):
    cid,p0,mid,tag,G,rho0,topo,frac,peak,hold,rates=task
    basep=production.fast_params(p0.base,G,rho0,frac);p=replace(p0,base=basep);paths={r:memory.run(p,production.FastSchedule(r,peak,hold)) for r in rates};rows=[]
    for target in TARGETS:
        indices={r:np.flatnonzero(h["rho"]>=target-1e-12) for r,h in paths.items()};ref_rate=.2 if .2 in paths and len(indices[.2]) else 1.;ref_ok=len(indices[ref_rate])>0
        for rate,h in paths.items():
            ok=len(indices[rate])>0;common=dict(candidate_id=cid,base_mechanism=mid,variant=tag,G0_nm=G,rho0=rho0,initial_topology=topo,heating_rate_C_min=rate,reference_rate_C_min=ref_rate,peak_T_C=peak,hold_time_h=hold,rho_target=target,comparison_attained=bool(ok and ref_ok))
            if not ok or not ref_ok:rows.append({**common,"response_class":"unattainable"});continue
            i=int(indices[rate][0]);j=int(indices[ref_rate][0]);Gx=float(h["G"][i])*1e9;Gr=float(paths[ref_rate]["G"][j])*1e9;hr=100*(Gr-Gx)/Gr;cls="beneficial" if hr>1 else ("harmful" if hr<-1 else "neutral");mx=state_metrics(h,i);mr=state_metrics(paths[ref_rate],j)
            rows.append({**common,"response_class":cls,"HR_pct":hr,"G_at_target_nm":Gx,"G_reference_nm":Gr,"time_to_target_h":float(h["t"][i])/3600,**mx,"reference_PR_work":mr["cumulative_PR_desintering_work"],"PR_exposure_difference":mr["cumulative_PR_desintering_work"]-mx["cumulative_PR_desintering_work"],"connected_fine_difference":mx["connected_fine_pore_fraction"]-mr["connected_fine_pore_fraction"],"mean_radius_difference_nm":mx["pore_mean_radius_nm"]-mr["pore_mean_radius_nm"]})
    return rows


def fast_screen(candidates, full=False, workers=1):
    Gs = production.G_FAST if full else G_FAST_REDUCED
    rho0s = production.RHO0S if full else RHO0_REDUCED
    rates = production.FAST_RATES if full else RATES_REDUCED
    peaks = production.PEAKS if full else PEAKS_REDUCED
    holds = production.HOLDS if full else HOLDS_REDUCED
    tasks=[(cid,p0,mid,tag,G,rho0,topo,frac,peak,hold,rates) for cid,p0,mid,tag in candidates for G in Gs for rho0 in rho0s for topo,frac in production.TOPOLOGIES.items() for peak in peaks for hold in holds]
    rows=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for n,result in enumerate(pool.map(_fast_task,tasks,chunksize=2),1):
            rows.extend(result)
            if n%100==0:print("fast groups",n,len(tasks),flush=True)
    return rows


def histories(best):
    cid,p,mid,tag=best; rows=[]
    basep=production.fast_params(p.base,75,.70,production.TOPOLOGIES["baseline"]);p=replace(p,base=basep)
    for label,rate in (("slow",1.),("fast",20.)):
        h=memory.run(p,production.FastSchedule(rate,1400,20));stride=max(1,len(h["rho"])//400)
        for i in list(range(0,len(h["rho"]),stride))+[len(h["rho"])-1]:
            rows.append(dict(candidate_id=cid,path=label,t_s=float(h["t"][i]),T_C=float(h["T_C"][i]),rho=float(h["rho"][i]),G_nm=float(h["G"][i])*1e9,**state_metrics(h,i)))
    return rows


def plots(out,fast,chen,hist,score):
    figdir=out/"figures";figdir.mkdir(exist_ok=True)
    for field,name,ylabel in (("connected_fine_pore_fraction","connected_fine_vs_density.png","connected fine-pore fraction"),("pore_mean_radius_nm","pore_size_memory_vs_density.png","connected mean radius [nm]"),("cumulative_PR_desintering_work","PR_exposure_vs_density.png","cumulative PR work")):
        fig,ax=plt.subplots()
        for label,q in _groups(hist,"path"):ax.plot([r["rho"] for r in q],[r[field] for r in q],label=label)
        ax.set(xlabel="density",ylabel=ylabel);ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(figdir/name,dpi=180);plt.close(fig)
    valid=[r for r in fast if r.get("comparison_attained") and math.isfinite(r.get("HR_pct",math.nan))]
    fig,ax=plt.subplots();ax.scatter([r["PR_exposure_difference"] for r in valid],[r["HR_pct"] for r in valid],s=8,alpha=.3);ax.axhline(1,color="k",ls="--");ax.set(xlabel="reference - fast PR work",ylabel="HR_pct");fig.tight_layout();fig.savefig(figdir/"HR_vs_PR_exposure_difference.png",dpi=180);plt.close(fig)
    complete=[r for r in chen if r["complete_practical"]];fig,ax=plt.subplots();ax.scatter([r["G0_nm"] for r in complete],[r.get("T_first_success_C",math.nan) for r in complete],c=[r["growth_tolerance"] for r in complete],cmap="viridis",s=20);ax.set(xlabel="G0 [nm]",ylabel="first success T2 [C]");fig.tight_layout();fig.savefig(figdir/"chen_window_preservation_map.png",dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,5));labels=[r["candidate_id"] for r in score];x=np.arange(len(labels));ax.bar(x-.2,[r["complete_chen_count"] for r in score],.4,label="Chen");ax.bar(x+.2,[r["beneficial_fast_count"] for r in score],.4,label="fast benefit");ax.set_xticks(x,labels,rotation=70);ax.legend();fig.tight_layout();fig.savefig(figdir/"joint_scorecard.png",dpi=180);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,5));
    for label,q in _groups(hist,"path"):axs[0].plot([r["rho"] for r in q],[r["G_nm"] for r in q],label=label);axs[1].plot([r["T_C"] for r in q],[r["connected_fine_pore_fraction"] for r in q],label=label)
    axs[0].set(xlabel="density",ylabel="G [nm]");axs[1].set(xlabel="temperature [C]",ylabel="connected fine fraction");axs[0].legend();fig.tight_layout();fig.savefig(figdir/"representative_path_diagnostics.png",dpi=180);plt.close(fig)


def _groups(rows,key):
    for value in sorted(set(r[key] for r in rows)):yield value,[r for r in rows if r[key]==value]


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--outdir",default="results/pr_desintering_fast_firing_memory");ap.add_argument("--full",action="store_true");ap.add_argument("--workers",type=int,default=4);ap.add_argument("--resume-chen",action="store_true");ap.add_argument("--selected-variant");args=ap.parse_args()
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();all_candidates=design();candidates=[c for c in all_candidates if args.selected_variant is None or c[3]==args.selected_variant]
    registry=[]
    for cid,p,mid,tag in all_candidates:
        row=asdict(p);row.pop("base");registry.append(dict(candidate_id=cid,base_mechanism=mid,variant=tag,**row))
    if args.resume_chen:
        with (out/"chen_window_preservation.csv").open(newline="") as stream:chen=list(csv.DictReader(stream))
        for row in chen:
            row["complete_practical"]=row["complete_practical"]=="True"
    else:
        chen=chen_screen(candidates);write(out/"chen_window_preservation.csv",chen)
    fast=fast_screen(candidates,args.full,args.workers)
    score=[]
    for cid,p,mid,tag in candidates:
        cq=[r for r in chen if r["candidate_id"]==cid and r["complete_practical"]]
        fq=[r for r in fast if r["candidate_id"]==cid and r.get("response_class")=="beneficial"]
        universal=bool(fq) and len(fq)==sum(r.get("comparison_attained",False) for r in fast if r["candidate_id"]==cid)
        joint=bool(cq and fq and not universal)
        score.append(dict(candidate_id=cid,base_mechanism=mid,variant=tag,complete_chen_count=len(cq),beneficial_fast_count=len(fq),universal_fast_benefit=universal,joint_positive=joint,rejection_reason="" if joint else ("no_complete_chen_window" if not cq else "no_beneficial_fast_case")))
    winners=[c for c in candidates if next(r for r in score if r["candidate_id"]==c[0])["joint_positive"]]
    best=max(winners,key=lambda c:next(r for r in score if r["candidate_id"]==c[0])["beneficial_fast_count"],default=candidates[0]);hist=histories(best)
    suffix="_full" if args.full else ""
    write(out/"parameter_registry.csv",registry);write(out/("full_joint_scorecard.csv" if args.full else "reduced_joint_screen.csv"),score);write(out/"chen_window_preservation.csv",chen);write(out/f"fast_firing_response_map{suffix}.csv",fast);write(out/("successful_joint_candidates_full.csv" if args.full else "successful_joint_candidates.csv"),[r for r in score if r["joint_positive"]]);write(out/"failed_joint_candidates.csv",[r for r in score if not r["joint_positive"]]);write(out/"rejected_parameter_sets.csv",[r for r in score if not r["joint_positive"]]);write(out/"representative_slow_fast_histories.csv",hist);write(out/f"PR_desintering_exposure{suffix}.csv",[{k:r.get(k,"") for k in ("candidate_id","G0_nm","rho0","initial_topology","heating_rate_C_min","peak_T_C","hold_time_h","rho_target","comparison_attained","HR_pct","cumulative_PR_desintering_work","reference_PR_work","PR_exposure_difference")} for r in fast]);write(out/f"connected_fine_pore_memory{suffix}.csv",[{k:r.get(k,"") for k in ("candidate_id","G0_nm","rho0","initial_topology","heating_rate_C_min","peak_T_C","hold_time_h","rho_target","comparison_attained","connected_fine_pore_fraction","connected_fine_difference","pore_mean_radius_nm","large_pore_fraction")} for r in fast]);plots(out,fast,chen,hist,score);write(out/f"runtime_summary{suffix}.csv",[{"wall_s":time.perf_counter()-start,"full_grid":args.full,"n_candidates":len(candidates),"local_signature":str(inspect.signature(memory.local_competition))}]);print("DONE",len(winners),"joint candidates",flush=True)


if __name__ == "__main__": main()
