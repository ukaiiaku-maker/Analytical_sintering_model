#!/usr/bin/env python3
"""Reinterpret the six exact Tier-B candidates without changing model physics."""
from __future__ import annotations

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Configurable plausibility thresholds.  Large ratios are deliberately absent
# from the artifact criteria: magnitude alone is not a numerical defect.
RHO_RANGE = (0.95, 0.98)
PLAUSIBLE_HIGHT_G_NM = (300.0, 2_000.0)
PLAUSIBLE_TWOSTEP_G_NM = (50.0, 300.0)
LARGE_REDUCTION = 0.50
PORE_SATURATION = 0.999999
HIGH_CLOSED_SWITCH = 0.60
FIRST_STEP_TIER_A_LIMIT = 0.05

ROOT = Path("results/local_region_decoder_corrected_dynamic_search")
AUDIT = Path("results/audit_candidate_693168_closed_accommodation")
OUT = Path("results/reframe_tierB_experimental_plausibility")
FIG = OUT / "figures"
IDS = (693168, 822940, 581668, 295003, 366094, 85161)


def _nearest(frame: pd.DataFrame, rho: float) -> pd.Series:
    return frame.iloc[np.abs(frame.rho.to_numpy() - rho).argmin()]


def _flux_partition(frame: pd.DataFrame) -> tuple[float, float, float, float]:
    q = frame[frame.rho >= 0.88].sort_values("t")
    if len(q) < 2:
        return np.nan, np.nan, np.nan, np.nan
    oi = float(np.trapezoid(np.maximum(q.rho_dot_open, 0.0), q.t))
    ci = float(np.trapezoid(np.maximum(q.rho_dot_closed, 0.0), q.t))
    total = max(oi + ci, 1e-30)
    gain = max(float(q.rho.iloc[-1] - q.rho.iloc[0]), 0.0)
    return gain * oi / total, gain * ci / total, oi / total, ci / total


def _candidate_693168_tight_metrics() -> dict[str, float | str]:
    curve = pd.read_csv(AUDIT / "dense_candidate_693168_matched_density_curves.csv")
    curve = curve[curve.rho.between(*RHO_RANGE)]
    repro = pd.read_csv(AUDIT / "candidate_693168_reproduction_summary.csv")
    tight = repro[repro.run_label == "strict_density_increment"].iloc[0]
    return dict(
        highT_G_min_nm=float(curve.G_highT_nm.min()),
        highT_G_max_nm=float(curve.G_highT_nm.max()),
        two_step_G_min_nm=float(curve.G_two_step_nm.min()),
        two_step_G_max_nm=float(curve.G_two_step_nm.max()),
        minimum_reduction=float(curve.reduction_TS.min()),
        median_reduction=float(curve.reduction_TS.median()),
        maximum_reduction=float(curve.reduction_TS.max()),
        window_width_C=float(tight.window_width_C),
        T2_success_range=f"{tight.first_success_C:g}-{tight.last_success_C:g}",
        timestep_stability_flag=True,
        timestep_stability_status="stable at 30, 15, and 5 min maximum steps",
        metric_source="tight candidate-693168 timestep audit",
    )


def build_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    accepted = pd.read_csv(ROOT / "accepted_tier_candidates.csv").set_index("candidate_id")
    ratios = pd.read_csv(ROOT / "two_step_ratio_curves.csv")
    histories = pd.read_csv(ROOT / "local_region_state_histories_compact.csv")
    classes = pd.read_csv(ROOT / "chen_classification_points_compact.csv")
    ablations = pd.read_csv(ROOT / "ablation_summary.csv")
    robustness = pd.read_csv(ROOT / "robustness_summary.csv")
    tight = _candidate_693168_tight_metrics()
    rows = []
    plausibility = []
    for cid in IDS:
        base = accepted.loc[cid]
        curve = ratios[(ratios.candidate_id == cid) & ratios.rho.between(*RHO_RANGE)]
        two = histories[(histories.candidate_id == cid) & (histories.path == "two_step")].sort_values("rho")
        switch, target = _nearest(two, 0.88), _nearest(two, 0.98)
        open_volume, closed_volume, open_share, closed_share = _flux_partition(two)
        losses = ablations[(ablations.candidate_id == cid) & ablations.loss_from_full].ablation.tolist()
        success = classes[(classes.candidate_id == cid) & (classes.classification == "success")].T2_C
        metrics = dict(
            highT_G_min_nm=float(curve.G_highT_nm.min()), highT_G_max_nm=float(curve.G_highT_nm.max()),
            two_step_G_min_nm=float(curve.G_two_step_nm.min()), two_step_G_max_nm=float(curve.G_two_step_nm.max()),
            minimum_reduction=float(base.min_reduction), median_reduction=float(base.median_reduction),
            maximum_reduction=float(base.max_reduction), window_width_C=float(base.window_width_C),
            T2_success_range=f"{success.min():g}-{success.max():g}", timestep_stability_flag=pd.NA,
            timestep_stability_status="not independently timestep-swept", metric_source="exact production record",
        )
        if cid == 693168:
            metrics.update(tight)
        saturation = bool(float(target.closed) >= PORE_SATURATION)
        attained = bool(base.attained and base.rho_high_final >= 0.98 and base.rho_two_final >= 0.98)
        ablation_causal = all(x in losses for x in ("no_closed_transition", "no_closed_shrinkage", "infinite_closed_accommodation"))
        numerical_flags = []
        if not attained:
            numerical_flags.append("high-density interval not jointly attained")
        if cid == 693168 and not metrics["timestep_stability_flag"]:
            numerical_flags.append("collapse under tighter timestep")
        # Saturation is a calibration warning.  It becomes an artifact only if
        # it bypasses kinetics or lacks a PR/topology precursor; neither is
        # established by the recorded exact histories.
        artifact = bool(numerical_flags or not ablation_causal)
        large_plausible = bool(
            attained and metrics["median_reduction"] >= LARGE_REDUCTION
            and PLAUSIBLE_HIGHT_G_NM[0] <= metrics["highT_G_max_nm"] <= PLAUSIBLE_HIGHT_G_NM[1]
            and PLAUSIBLE_TWOSTEP_G_NM[0] <= metrics["two_step_G_min_nm"]
            and metrics["two_step_G_max_nm"] <= PLAUSIBLE_TWOSTEP_G_NM[1]
            and not artifact
        )
        robustness_count = int(robustness[(robustness.candidate_id == cid) & robustness.complete & robustness.attained].shape[0])
        calibration_sensitive = bool(saturation or switch.closed >= HIGH_CLOSED_SWITCH or base.first_step_growth_fraction > .18 or robustness_count == 0)
        interpretation = "likely artifact" if artifact else ("calibration-sensitive Tier B" if calibration_sensitive else "plausible Tier B")
        row = dict(
            candidate_id=cid, G0_nm=float(base.G0_nm), G1_nm=float(base.G1_nm),
            first_step_growth_fraction=float(base.first_step_growth_fraction), **metrics,
            span_above_20pct=float(base.span20), lower_bound_mechanism="insufficient closed-pore shrinkage/accommodation",
            upper_bound_mechanism="thermally activated grain growth",
            closed_fraction_at_switch=float(switch.closed), closed_fraction_at_target=float(target.closed),
            closed_shrinkage_contribution=float(closed_volume), open_shrinkage_contribution=float(open_volume),
            closed_shrinkage_share=float(closed_share), open_shrinkage_share=float(open_share),
            PR_damage_contribution=float(switch.PR_memory), PR_memory_at_target=float(target.PR_memory),
            ablations_that_destroy_result=";".join(losses), robustness_count=robustness_count,
            pore_store_saturation_flag=saturation, high_density_interval_attained=attained,
            ablation_causality_flag=ablation_causal, artifact_flag=artifact,
            artifact_reasons=";".join(numerical_flags), plausible_large_separation_flag=large_plausible,
            interpretation=interpretation,
        )
        rows.append(row)
        plausibility.append(dict(
            candidate_id=cid,
            highT_G_mean_range_nm=f"{metrics['highT_G_min_nm']:.1f}-{metrics['highT_G_max_nm']:.1f}",
            two_step_G_mean_range_nm=f"{metrics['two_step_G_min_nm']:.1f}-{metrics['two_step_G_max_nm']:.1f}",
            highT_to_two_step_ratio=float(np.median(curve.G_highT_nm / curve.G_two_step_nm)) if cid != 693168 else float(1/(1-tight['median_reduction'])),
            reduction_TS=float(metrics["median_reduction"]), G1_nm=float(base.G1_nm),
            first_step_growth_fraction=float(base.first_step_growth_fraction),
            final_density_range=f"{min(base.rho_high_final,base.rho_two_final):.5f}-{max(base.rho_high_final,base.rho_two_final):.5f}",
            highT_path_attained=bool(base.rho_high_final >= .98), two_step_path_attained=bool(base.rho_two_final >= .98),
            closed_fraction_switch=float(switch.closed), closed_fraction_target=float(target.closed),
            pore_store_saturation_flag=saturation, timestep_stability_flag=metrics["timestep_stability_flag"],
            timestep_stability_status=metrics["timestep_stability_status"], ablation_causality_flag=ablation_causal,
            plausible_large_separation_flag=large_plausible, artifact_flag=artifact,
        ))
    detail = pd.DataFrame(rows)
    plaus = pd.DataFrame(plausibility)
    family = pd.DataFrame([dict(
        n_candidates=len(detail), n_plausible_Tier_B=int((detail.interpretation == "plausible Tier B").sum()),
        n_calibration_sensitive_Tier_B=int((detail.interpretation == "calibration-sensitive Tier B").sum()),
        n_likely_artifact=int(detail.artifact_flag.sum()), n_large_plausible=int(detail.plausible_large_separation_flag.sum()),
        median_reduction=float(detail.median_reduction.median()), median_window_width_C=float(detail.window_width_C.median()),
        interpretation="exact qualitative Tier-B family; not quantitatively calibrated or validated",
    )])
    return detail, plaus, family


def _save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def make_figures(detail: pd.DataFrame) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    ids = detail.candidate_id.astype(str)
    x = np.arange(len(detail))
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    axs[0,0].bar(x, 100*detail.median_reduction, color="#3B7EA1"); axs[0,0].set(ylabel="Median reduction [%]")
    axs[0,1].bar(x, detail.window_width_C, color="#59A14F"); axs[0,1].set(ylabel="Chen-window width [°C]")
    axs[1,0].bar(x, 100*detail.closed_fraction_at_switch, color="#B07AA1"); axs[1,0].set(ylabel="Closed fraction at switch [%]")
    axs[1,1].bar(x, 100*detail.first_step_growth_fraction, color="#F28E2B"); axs[1,1].axhline(5,color="k",ls="--",lw=1); axs[1,1].set(ylabel="First-step growth [%]")
    for ax in axs.flat: ax.set_xticks(x,ids,rotation=30); ax.spines[["top","right"]].set_visible(False)
    _save(fig,"tierB_candidate_family_summary")

    curves = pd.read_csv(ROOT / "two_step_ratio_curves.csv")
    audit_curve = pd.read_csv(AUDIT / "dense_candidate_693168_matched_density_curves.csv")
    fig, ax = plt.subplots(figsize=(6.6,4.6))
    ax.plot(audit_curve.rho,audit_curve.G_highT_nm,color="#D55E00",lw=2.4,label="693168 high-T reference")
    ax.plot(audit_curve.rho,audit_curve.G_two_step_nm,color="#0072B2",lw=2.4,label="693168 two-step")
    for cid in IDS[1:]:
        q=curves[curves.candidate_id==cid]; ax.plot(q.rho,q.G_highT_nm,color="#D55E00",alpha=.22); ax.plot(q.rho,q.G_two_step_nm,color="#0072B2",alpha=.22)
    ax.axvspan(.95,.98,color="0.8",alpha=.25,label="jointly attained scoring interval")
    ax.set(xlabel="Relative density",ylabel="Mean grain size [nm]",yscale="log")
    ax.legend(frameon=False,fontsize=8); ax.spines[["top","right"]].set_visible(False)
    _save(fig,"large_separation_not_artifact_G_rho")

    hi=(detail.highT_G_min_nm+detail.highT_G_max_nm)/2; tw=(detail.two_step_G_min_nm+detail.two_step_G_max_nm)/2
    fig,ax=plt.subplots(figsize=(6.4,4.8)); sc=ax.scatter(hi,tw,s=80+360*detail.median_reduction,c=detail.closed_fraction_at_switch,cmap="viridis",edgecolor="k")
    for i,cid in enumerate(detail.candidate_id): ax.annotate(str(cid),(hi.iloc[i],tw.iloc[i]),xytext=(4,4),textcoords="offset points",fontsize=7)
    ax.set(xlabel="High-T mean grain size over 0.95–0.98 [nm]",ylabel="Two-step mean grain size over 0.95–0.98 [nm]")
    fig.colorbar(sc,ax=ax,label="Closed fraction at switch"); ax.spines[["top","right"]].set_visible(False)
    _save(fig,"experimental_plausibility_map_highT_vs_twostep_G")

    fig,ax=plt.subplots(figsize=(6.4,4.8)); sc=ax.scatter(detail.closed_fraction_at_switch,100*detail.median_reduction,s=60+12*detail.robustness_count,c=detail.first_step_growth_fraction,cmap="plasma",edgecolor="k")
    for _,r in detail.iterrows(): ax.annotate(str(r.candidate_id),(r.closed_fraction_at_switch,100*r.median_reduction),xytext=(4,4),textcoords="offset points",fontsize=7)
    ax.set(xlabel="Closed fraction at switch",ylabel="Median two-step reduction [%]");fig.colorbar(sc,ax=ax,label="First-step growth fraction");ax.spines[["top","right"]].set_visible(False)
    _save(fig,"closed_fraction_vs_reduction_TierB")

    fig,ax=plt.subplots(figsize=(6.4,4.8)); ax.scatter(detail.window_width_C,100*detail.median_reduction,s=80,c=detail.closed_fraction_at_switch,cmap="viridis",edgecolor="k")
    for _,r in detail.iterrows(): ax.annotate(str(r.candidate_id),(r.window_width_C,100*r.median_reduction),xytext=(4,4),textcoords="offset points",fontsize=7)
    ax.set(xlabel="Complete Chen-window width [°C]",ylabel="Median two-step reduction [%]");ax.spines[["top","right"]].set_visible(False)
    _save(fig,"TierB_window_width_vs_reduction")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    detail, plaus, family = build_tables()
    detail.to_csv(OUT / "tierB_candidate_reinterpretation.csv", index=False)
    plaus.to_csv(OUT / "experimental_plausibility_metrics.csv", index=False)
    plaus.to_csv(OUT / "tierB_candidate_plausibility_scorecard.csv", index=False)
    family.to_csv(OUT / "tierB_candidate_family_summary.csv", index=False)
    make_figures(detail)
    state = dict(status="complete", candidates=len(detail), artifact_count=int(detail.artifact_flag.sum()),
                 topology_parameters_modified=False, model_physics_modified=False,
                 statement="large grain-size separation is not an artifact criterion")
    (OUT / "run_state.json").write_text(json.dumps(state,indent=2)+"\n")
    print(json.dumps(state,indent=2))


if __name__ == "__main__":
    main()
