#!/usr/bin/env python3
"""Audit observable matched-density trajectory separation without changing physics."""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_style as ps
import production_mechanism_assessment as prior
import production_pr_desintering_assessment as production
import pr_desintering_memory_model as memory
import topology_constrained_sintering as aggregate

SOURCE = Path("results/production_pr_desintering_assessment")
FAST_SOURCE = Path("results/pr_desintering_fast_firing_memory/raw_fast_firing_response_map_full.csv")
OUT = Path("results/observable_trajectory_effect_audit")
WINDOWS = {"early": (.75, .85), "intermediate": (.85, .92),
           "late_open_pore": (.92, .95), "near_final": (.95, .99)}
THRESHOLDS = (1.2, 1.5, 2.0)


def write_csv(path: Path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields or ["status"], lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def monotonic_xy(frame, x="rho", y="G_nm"):
    q = frame[[x, y]].dropna().sort_values(x).groupby(x, as_index=False).first()
    return q[x].to_numpy(float), q[y].to_numpy(float)


def matched_curve(reference, comparison, lo=None, hi=None, step=.001):
    """Return matched-density interpolation only over jointly attained density."""
    xr, yr = monotonic_xy(reference); xc, yc = monotonic_xy(comparison)
    if not len(xr) or not len(xc):
        return pd.DataFrame(columns=["rho", "G_reference_nm", "G_comparison_nm", "ratio", "pct", "delta_ln_G"])
    lower = max(xr.min(), xc.min(), -np.inf if lo is None else lo)
    upper = min(xr.max(), xc.max(), np.inf if hi is None else hi)
    if upper < lower:
        return pd.DataFrame(columns=["rho", "G_reference_nm", "G_comparison_nm", "ratio", "pct", "delta_ln_G"])
    grid = np.arange(math.ceil(lower / step) * step, upper + step * .25, step)
    gr = np.interp(grid, xr, yr); gc = np.interp(grid, xc, yc)
    ratio = gr / gc
    return pd.DataFrame(dict(rho=grid, G_reference_nm=gr, G_comparison_nm=gc,
                             ratio=ratio, pct=100 * (gr - gc) / gr,
                             delta_ln_G=np.log(ratio)))


def longest_span(rho, ratio, threshold):
    good = np.asarray(ratio) >= threshold
    best = 0.; start = None
    for i, ok in enumerate(good):
        if ok and start is None: start = i
        if start is not None and (not ok or i == len(good)-1):
            end = i if ok and i == len(good)-1 else i-1
            best = max(best, float(rho[end] - rho[start])); start = None
    return best


def tier(max_ratio):
    if not np.isfinite(max_ratio) or max_ratio < 1.2: return "negligible"
    if max_ratio < 1.5: return "weak"
    if max_ratio < 2.: return "meaningful"
    return "strong"


def classify(curve, supported_hi=.95):
    if curve.empty: return "unattainable"
    meaningful = curve[curve.ratio >= 1.5]
    span = longest_span(curve.rho.to_numpy(), curve.ratio.to_numpy(), 1.5)
    if span >= .03:
        if meaningful.rho.min() >= .95: return "unsupported_high_density" if supported_hi <= .95 else "late_only"
        return "trajectory_meaningful"
    weak_span = longest_span(curve.rho.to_numpy(), curve.ratio.to_numpy(), 1.2)
    return "trajectory_weak" if weak_span > 0 else "trajectory_negligible"


def window_row(cid, kind, curve, name, lo, hi, meta=None):
    q = curve[(curve.rho >= lo-1e-12) & (curve.rho <= hi+1e-12)]
    attained = len(q) >= 2 and q.rho.min() <= lo+.002 and q.rho.max() >= hi-.002
    row = dict(comparison_id=cid, comparison_type=kind, density_window=name,
               rho_lo=lo, rho_hi=hi, both_paths_attained=attained)
    if not attained:
        return {**row, "uncertainty_class": "unattainable", "trajectory_class": "unattainable"}
    ratios = q.ratio.to_numpy(); rho = q.rho.to_numpy()
    area = float(np.trapezoid(np.maximum(0, np.log(ratios)), rho))
    row.update(max_ratio=float(ratios.max()), median_ratio=float(np.median(ratios)),
               integrated_separation=area, uncertainty_class=tier(float(ratios.max())),
               trajectory_class=classify(q))
    for threshold in THRESHOLDS: row[f"span_ratio_ge_{threshold:g}"] = longest_span(rho, ratios, threshold)
    if meta: row.update(meta)
    return row


def fast_curves():
    raw = pd.read_csv(FAST_SOURCE)
    raw = raw[(raw.comparison_attained == True) & raw.G_at_target_nm.notna() & raw.G_reference_nm.notna()].copy()
    keys = ["candidate_id","G0_nm","rho0","initial_topology","heating_rate_C_min",
            "reference_rate_C_min","peak_T_C","hold_time_h"]
    curves=[]; windows=[]; summaries=[]
    for key, q in raw.groupby(keys, dropna=False):
        q=q.sort_values("rho_target"); meta=dict(zip(keys,key)); cid="|".join(map(str,key))
        c=pd.DataFrame(dict(rho=q.rho_target, G_reference_nm=q.G_reference_nm,
                            G_comparison_nm=q.G_at_target_nm))
        c["ratio"]=c.G_reference_nm/c.G_comparison_nm
        c["pct"]=100*(c.G_reference_nm-c.G_comparison_nm)/c.G_reference_nm
        c["delta_ln_G"]=np.log(c.ratio); c["comparison_id"]=cid
        for k,v in meta.items(): c[k]=v
        curves.append(c)
        wr=[window_row(cid,"fast_firing",c,n,*bounds,{"reference_fallback":meta["reference_rate_C_min"]==1.0}) for n,bounds in WINDOWS.items()]
        windows.extend(wr)
        valid=[r for r in wr if r["both_paths_attained"]]
        summaries.append({**meta,"comparison_id":cid,"comparison_type":"fast_firing",
                          "reference_fallback":meta["reference_rate_C_min"]==1.0,
                          "max_ratio":max((r["max_ratio"] for r in valid),default=math.nan),
                          "classification":"trajectory_meaningful" if any(r["trajectory_class"]=="trajectory_meaningful" for r in valid) else ("trajectory_weak" if any(r["trajectory_class"]=="trajectory_weak" for r in valid) else "trajectory_negligible")})
    return pd.concat(curves,ignore_index=True), windows, summaries


def representative_frames():
    sf=pd.read_csv(SOURCE/"representative_slow_fast_histories.csv")
    ts=pd.read_csv(SOURCE/"representative_two_step_histories.csv")
    p0=production.candidates()["mech_009_q0"]
    b=replace(p0.base.action.location.base,G0=225e-9)
    p=replace(p0,base=replace(p0.base,action=replace(p0.base.action,location=replace(p0.base.action.location,base=b))))
    ht=pd.DataFrame(production.history_rows("highT","highT_isothermal",memory.run(p,aggregate.Iso(1400,production.BUDGET))))
    return sf[sf.path=="slow"],sf[sf.path=="fast"],ht,ts


def representative_curves():
    slow,fast,high,two=representative_frames()
    hr=matched_curve(slow,fast); ts=matched_curve(high,two)
    hr["comparison_id"]="representative_fast_1_vs_20";hr["reference_rate_C_min"]=1.;hr["reference_fallback"]=True
    ts["comparison_id"]="representative_highT_vs_two_step"
    windows=[]
    for kind,c in (("fast_firing",hr),("two_step",ts)):
        for name,bounds in WINDOWS.items():windows.append(window_row(c.comparison_id.iloc[0],kind,c,name,*bounds))
    return hr,ts,windows


def rescue_design():
    base=dict(k_PR_ref_s=2e-4,Q_PR_J_mol=180e3,renewal_gate_mid=.35,renewal_power=2.,
              smoothing_share=.65,GB_to_TJ_share=.25,TJ_to_iso_share=.10,topology_power=1.)
    design=[("baseline",base)]
    for name,values in (("k_PR_ref_s",(2e-4,5e-4,1e-3,2e-3)),
                        ("Q_PR_J_mol",(140e3,180e3,220e3,260e3)),
                        ("renewal_gate_mid",(.25,.35,.45,.55)),
                        ("renewal_power",(1.,2.,3.,4.)),
                        ("topology_power",(1.,2.,3.))):
        for value in values: design.append((f"{name}_{value:g}",{**base,name:value}))
    for label,shares in (("more_smoothing",(.80,.15,.05)),("more_GB_to_TJ",(.45,.50,.05)),
                         ("more_TJ_to_iso",(.45,.15,.40)),("balanced",(1/3,1/3,1/3))):
        design.append((f"partition_{label}",{**base,"smoothing_share":shares[0],"GB_to_TJ_share":shares[1],"TJ_to_iso_share":shares[2]}))
    unique={label:kw for label,kw in design};return list(unique.items())


def rescue_screen(run=True):
    if not run:return [],[]
    pbase=production.candidates()["mech_009_q0"]
    pbase=replace(pbase,base=prior.fast_params(pbase.base,75,.70,prior.TOPOLOGIES["baseline"]))
    rows=[];rejected=[]
    for label,kw in rescue_design():
        p=replace(pbase,**kw); histories={r:memory.run(p,prior.FastSchedule(r,1400,20)) for r in (1.,20.)}
        c=matched_curve(pd.DataFrame(production.history_rows("reference","fast",histories[1.])),
                        pd.DataFrame(production.history_rows("fast","fast",histories[20.])),.85,.95)
        mid=c[(c.rho>=.85)&(c.rho<=.92)]
        attained=len(mid)>1 and mid.rho.max()>=.918
        span=longest_span(mid.rho.to_numpy(),mid.ratio.to_numpy(),1.5) if attained else 0.
        row=dict(screen_id=label,reference_rate_C_min=1.,fast_rate_C_min=20.,rho_lo=.85,rho_hi=.92,
                 both_paths_attained=attained,max_ratio=float(mid.ratio.max()) if len(mid) else math.nan,
                 median_ratio=float(mid.ratio.median()) if len(mid) else math.nan,
                 meaningful_span=span,meaningful_trajectory=span>=.03,
                 chen_window_preserved="not_retested_no_meaningful_candidate",
                 target_or_budget_changed=False,conservative_PR_flux=True)
        rows.append(row)
        if not row["meaningful_trajectory"]: rejected.append({**row,"rejection_reason":"ratio_below_1.5_over_Drho_0.03" if attained else "joint_target_nonattainment"})
    return rows,rejected


def plots(fast_all,hr,ts,windows,attainment,slow,fast,high,two):
    ps.apply_style();OUT.mkdir(parents=True,exist_ok=True)
    # 1: honest common-axis trajectories and ratio inset
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),constrained_layout=True)
    for label,q,col in (("slow",slow,ps.COLORS["slow"]),("fast",fast,ps.COLORS["fast"]),("high-T",high,ps.COLORS["highT"]),("two-step",two,ps.COLORS["two_step"])):
        axs[0].plot(q.rho,q.G_nm,label=label,color=col)
    axs[0].set(xlim=(.7,.95),xlabel="Relative density, $\\rho$",ylabel="Mean grain size, $G$ [nm]");axs[0].legend(fontsize=7);ps.clean(axs[0])
    axs[1].plot(hr.rho,hr.ratio,label="slow / fast",color=ps.COLORS["fast"]);axs[1].plot(ts.rho,ts.ratio,label="high-T / two-step",color=ps.COLORS["two_step"])
    for v,ls in zip(THRESHOLDS,(":","--","-")):axs[1].axhline(v,color="#666",ls=ls,lw=.8)
    axs[1].set(xlim=(.7,.95),ylim=(.8,2.1),xlabel="Relative density, $\\rho$",ylabel="Matched-density grain-size ratio");axs[1].legend(fontsize=7);ps.clean(axs[1]);ps.panel_labels(axs);ps.finish(fig,OUT/"rho_G_trajectories_representative")
    # 2
    fig,ax=plt.subplots(figsize=(5.2,3.6));ax.plot(hr.rho,hr.ratio,label="slow / fast");ax.plot(ts.rho,ts.ratio,label="high-T / two-step")
    for v in THRESHOLDS:ax.axhline(v,color="#777",ls="--",lw=.8);ax.text(.705,v+.02,f"{v:g}×",fontsize=7)
    ax.set(xlim=(.7,.95),ylim=(.8,2.1),xlabel="Relative density, $\\rho$",ylabel="Grain-size ratio");ax.legend();ps.clean(ax);ps.finish(fig,OUT/"grain_size_ratio_vs_density")
    # 3
    fig,ax=plt.subplots(figsize=(5.2,3.6));ax.axhspan(0,20,color="#EEEEEE");ax.axhspan(20,33,color="#FDE7B2");ax.axhspan(33,70,color="#CDE8D2");ax.plot(hr.rho,hr.pct,label="HR");ax.plot(ts.rho,ts.pct,label="TS");ax.axhline(33,color="#555",ls="--",lw=.8);ax.set(xlim=(.7,.95),ylim=(-10,70),xlabel="Relative density, $\\rho$",ylabel="Matched-density improvement [%]");ax.legend();ps.clean(ax);ps.finish(fig,OUT/"HR_TS_pct_vs_density")
    # 4
    w=pd.DataFrame(windows);p=w.pivot_table(index="comparison_type",columns="density_window",values="median_ratio",aggfunc="median").reindex(columns=WINDOWS)
    fig,ax=plt.subplots(figsize=(6.2,2.6));masked=np.ma.masked_invalid(p.to_numpy());im=ax.imshow(masked,aspect="auto",vmin=1,vmax=1.5,cmap="viridis");ax.set_xticks(range(len(p.columns)),p.columns,rotation=20);ax.set_yticks(range(len(p.index)),p.index);fig.colorbar(im,ax=ax,label="Median grain-size ratio");ps.finish(fig,OUT/"density_window_effect_heatmap")
    # 5
    fig,ax=plt.subplots(figsize=(5.6,3.2));a=attainment.pivot_table(index="rho_target",columns="comparison_attained",values="n_cases",aggfunc="sum",fill_value=0);a.plot.bar(stacked=True,ax=ax,color=["#CC6677","#44AA99"]);ax.axvspan(2.5,6.5,color="#CC6677",alpha=.12,label="unsupported closed-pore regime");ax.set(xlabel="Target density",ylabel="Production comparisons",title="Attainment, not extrapolation");ps.clean(ax);ps.finish(fig,OUT/"high_density_attainment_map")
    # 6
    fig,axs=plt.subplots(2,2,figsize=(7.2,5.2),constrained_layout=True)
    for ax,field,label in zip(axs.flat,("connected_fine_pore_fraction","connected_mean_radius_nm","cumulative_PR_desintering_work","G_nm"),("Connected fine-pore fraction","Connected radius [nm]","Cumulative PR work","Grain size [nm]")):
        for name,q,col in (("slow",slow,ps.COLORS["slow"]),("fast",fast,ps.COLORS["fast"])):ax.plot(q.rho,q[field],label=name,color=col)
        ax.set(xlabel="Relative density, $\\rho$",ylabel=label);ps.clean(ax)
    axs[0,0].legend();ps.panel_labels(axs);ps.finish(fig,OUT/"PR_memory_vs_visible_trajectory_effect")
    # 7: explicit negative control--internal memory exists, observable threshold fails.
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.2),constrained_layout=True)
    axs[0].plot(slow.rho,slow.cumulative_PR_desintering_work,label="slow",color=ps.COLORS["slow"])
    axs[0].plot(fast.rho,fast.cumulative_PR_desintering_work,label="fast",color=ps.COLORS["fast"])
    axs[0].set(xlabel="Relative density, $\\rho$",ylabel="Cumulative PR work [model units]",title="Internal diagnostic separates")
    axs[1].plot(hr.rho,hr.ratio,color=ps.COLORS["fast"]);axs[1].axhline(1.5,color="#555",ls="--",label="meaningful threshold")
    axs[1].set(xlabel="Relative density, $\\rho$",ylabel="$G_{slow}/G_{fast}$",ylim=(.9,1.55),title="Observable trajectory does not")
    for ax in axs:ps.clean(ax);ax.legend(fontsize=7)
    ps.panel_labels(axs);ps.finish(fig,OUT/"negative_control_internal_vs_observable")


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--skip-rescue",action="store_true");args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    fast_all,fw,summary=fast_curves();hr,ts,rw=representative_curves();windows=fw+rw
    write_csv(OUT/"fast_firing_ratio_curves.csv",fast_all.to_dict("records"))
    write_csv(OUT/"two_step_ratio_curves.csv",ts.to_dict("records"))
    write_csv(OUT/"density_window_effects.csv",windows)
    write_csv(OUT/"trajectory_effect_summary.csv",summary)
    attainment=(pd.read_csv(FAST_SOURCE).groupby(["candidate_id","rho_target","comparison_attained"],dropna=False).size().rename("n_cases").reset_index())
    extra=[]
    for target in (.95,.98,.99):extra.append(dict(candidate_id="current_open_pore_model",rho_target=target,comparison_attained=False,n_cases=0,physics_support=False,reason="not sampled; closed-pore treatment absent"))
    attrows=attainment.assign(physics_support=attainment.rho_target<=.92,reason="production table").to_dict("records")+extra
    write_csv(OUT/"high_density_attainment.csv",attrows)
    weak=[r for r in summary if r["classification"]!="trajectory_meaningful"]
    meaningful=[r for r in summary if r["classification"]=="trajectory_meaningful"]
    summary_fields=list(summary[0]) if summary else ["comparison_id","classification"]
    write_csv(OUT/"weak_or_negligible_cases.csv",weak,summary_fields);write_csv(OUT/"meaningful_trajectory_cases.csv",meaningful,summary_fields)
    # Baseline fails unless a finite meaningful interval exists; only then skip rescue.
    baseline_pass=any(r["trajectory_class"]=="trajectory_meaningful" for r in rw)
    rescue,rejected=rescue_screen(not args.skip_rescue and not baseline_pass)
    write_csv(OUT/"PR_rescue_screen_summary.csv",rescue);write_csv(OUT/"rejected_PR_rescue_cases.csv",rejected)
    slow,fast,high,two=representative_frames();plots(fast_all,hr,ts,windows,pd.DataFrame(attrows),slow,fast,high,two)
    print(f"observable audit: {len(summary)} fast comparisons; representative pass={baseline_pass}; rescue={len(rescue)}")


if __name__=="__main__":main()
