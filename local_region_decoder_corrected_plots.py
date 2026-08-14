#!/usr/bin/env python3
"""Production diagnostics for the decoder-corrected local-region campaign."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plot_style as ps

OUT = Path("results/local_region_decoder_corrected_dynamic_search")
AUDIT = Path("results/local_region_decoder_audit")
FIG = OUT / "figures"
COLORS = {
    "success": "#009E73", "density_exhaustion": "#0072B2",
    "grain_growth": "#D55E00", "mixed": "#CC79A7",
    "unattainable_first_step": "#777777",
}
INVENTORY = []


def save(fig, stem, title, source, purpose):
    ps.finish(fig, FIG / stem)
    INVENTORY.append(ps.inventory_row(len(INVENTORY) + 1, stem, title, source, purpose, "Results"))


def line_by_path(ax, frame, x, y, ylabel):
    for path, group in frame.groupby("path"):
        ax.plot(group[x], group[y], label=path, color=ps.COLORS.get(path, None))
    ax.set(xlabel=x.replace("_", " "), ylabel=ylabel)
    ax.legend()
    ps.clean(ax)


def main():
    ps.apply_style()
    FIG.mkdir(parents=True, exist_ok=True)
    stage0 = pd.read_csv(OUT / "stage0_massive_screen.csv.gz")
    stage2 = pd.read_csv(OUT / "stage2_exact_dynamic_summary.csv")
    accepted = pd.read_csv(OUT / "accepted_tier_candidates.csv")
    histories = pd.read_csv(OUT / "local_region_state_histories_compact.csv")
    ratios = pd.read_csv(OUT / "two_step_ratio_curves.csv")
    points = pd.read_csv(OUT / "chen_classification_points_compact.csv")
    ablations = pd.read_csv(OUT / "ablation_summary.csv")
    robust = pd.read_csv(OUT / "robustness_summary.csv")
    fast = pd.read_csv(OUT / "fast_firing_preservation.csv")
    best = int(accepted.sort_values("median_reduction", ascending=False).iloc[0].candidate_id)
    h = histories[histories.candidate_id == best]
    ratio = ratios[ratios.candidate_id == best]
    point = points[points.candidate_id == best]

    summary = pd.read_csv(AUDIT / "decoder_fingerprint_summary.csv").iloc[0]
    fig, ax = plt.subplots(figsize=(6, 4))
    labels = ["sampled", "unique vectors", "unique fingerprints"]
    values = [summary.sampled_rows, summary.unique_parameter_vectors, summary.unique_dynamic_fingerprints]
    ax.bar(labels, values, color=["#999999", "#56B4E9", "#009E73"])
    ax.set(ylabel="Candidate count", title="Decoder diversity preflight")
    ax.tick_params(axis="x", rotation=15); ps.clean(ax)
    save(fig, "decoder_fingerprint_diversity", "Decoder fingerprint diversity", "decoder_fingerprint_summary.csv", "Verify true decoder diversity")

    ordered = stage0.sort_values("score", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6, 4)); ax.plot(np.arange(1, len(ordered) + 1), ordered.score)
    ax.set(xlabel="Promoted rank", ylabel="Stage-0 surrogate score", xscale="log", title="Screen convergence")
    ps.clean(ax); save(fig, "optimizer_or_screen_convergence", "Screen convergence", "stage0_massive_screen.csv.gz", "Show retained-score decay")

    fig, ax = plt.subplots(figsize=(6, 4)); sc = ax.scatter(stage0.lower_boundary_likelihood, stage0.upper_boundary_likelihood, c=stage0.approximate_high_density_reduction, s=8, cmap="viridis")
    ax.set(xlabel="Lower-bound likelihood", ylabel="Upper-bound likelihood", title="Stage-0 physical phase map"); fig.colorbar(sc, ax=ax, label="Projected reduction")
    ps.clean(ax); save(fig, "stage0_phase_map", "Stage-0 phase map", "stage0_massive_screen.csv.gz", "Display surrogate boundary trade-off")

    fig, ax = plt.subplots(figsize=(6, 4)); ax.scatter(stage2.window_width_C, stage2.median_reduction, c=stage2.span20, cmap="plasma", s=12, alpha=.65)
    ax.scatter(accepted.window_width_C, accepted.median_reduction, facecolors="none", edgecolors="black", s=45, label="production eligible")
    ax.set(xlabel="Chen window width [°C]", ylabel="Median two-step reduction", title="Exact Pareto front"); ax.legend(); ps.clean(ax)
    save(fig, "pareto_front", "Exact Pareto front", "stage2_exact_dynamic_summary.csv", "Compare window width and matched-density benefit")

    fig, ax = plt.subplots(figsize=(6, 4)); line_by_path(ax, h, "rho", "G_mean", "Mean grain size, G [nm]")
    ax.set_title(f"Candidate {best}: high-T vs two-step")
    save(fig, "best_highT_vs_twostep_G_rho", "Best G-rho trajectories", "local_region_state_histories_compact.csv", "Show matched-density separation")

    fig, ax = plt.subplots(figsize=(6, 4)); ax.plot(ratio.rho, 100 * ratio.reduction_TS, color="#009E73")
    ax.axhline(20, color="#555555", ls="--"); ax.set(xlabel="Relative density, ρ", ylabel="Two-step grain reduction [%]", title=f"Candidate {best}")
    ps.clean(ax); save(fig, "reduction_TS_vs_density_best", "Two-step reduction versus density", "two_step_ratio_curves.csv", "Show attained density span")

    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    for path, group in h.groupby("path"):
        axes[0].plot(group.t / 3600, group.rho, label=path); axes[1].plot(group.t / 3600, group.G_mean, label=path)
    axes[0].set(ylabel="Relative density, ρ"); axes[1].set(xlabel="Physical time [h]", ylabel="Mean G [nm]")
    axes[0].legend(); [ps.clean(ax) for ax in axes]; ps.panel_labels(axes)
    save(fig, "best_physical_time_histories", "Physical-time histories", "local_region_state_histories_compact.csv", "Show continuous path timing")

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
    for path, group in h.groupby("path"):
        for ax, key, label in zip(axes.flat, ("connected", "PR_memory", "sweep_memory", "X_J"), ("Connected removable fraction", "PR memory", "Sweep memory", "Persistent junction state")):
            ax.plot(group.rho, group[key], label=path); ax.set(ylabel=label)
    for ax in axes[-1]: ax.set(xlabel="Relative density, ρ")
    axes[0, 0].legend(); [ps.clean(ax) for ax in axes.flat]; ps.panel_labels(axes)
    save(fig, "local_region_topology_histories", "Local topology histories", "local_region_state_histories_compact.csv", "Show evolving topology memory")

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for path, group in h.groupby("path"):
        axes[0].plot(group.rho, group.connected, label=path); axes[1].plot(group.rho, group.topology_variance, label=path)
    axes[0].set(xlabel="Relative density, ρ", ylabel="Connected fraction"); axes[1].set(xlabel="Relative density, ρ", ylabel="Connectivity variance")
    axes[0].legend(); [ps.clean(ax) for ax in axes]; ps.panel_labels(axes)
    save(fig, "pore_connectivity_distribution_histories", "Pore connectivity distribution", "local_region_state_histories_compact.csv", "Show network heterogeneity")

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for path, group in h.groupby("path"):
        axes[0].plot(group.rho, group.N_total, label=path); axes[1].plot(group.rho, group.pore_radius_proxy, label=path); axes[2].plot(group.rho, group.f_GBseg, label=f"GB {path}"); axes[2].plot(group.rho, group.f_TJ, ls="--", label=f"TJ {path}")
    axes[0].set(ylabel="Pore-number proxy"); axes[1].set(ylabel="Mean-radius proxy"); axes[2].set(ylabel="Pore location fraction")
    for ax in axes: ax.set(xlabel="Relative density, ρ"); ps.clean(ax)
    axes[0].legend(); axes[2].legend(fontsize=7); ps.panel_labels(axes)
    save(fig, "pore_number_D90_location_histories", "Pore number, radius, and location", "local_region_state_histories_compact.csv", "Show observable pore evolution")

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for path, group in h.groupby("path"):
        axes[0].plot(group.rho, group.closed, label=path); axes[1].plot(group.rho, group.closed_accommodation, label=path)
    axes[0].set(xlabel="Relative density, ρ", ylabel="Closed-pore fraction"); axes[1].set(xlabel="Relative density, ρ", ylabel="Closed accommodation")
    axes[0].legend(); [ps.clean(ax) for ax in axes]; ps.panel_labels(axes)
    save(fig, "closed_pore_and_accommodation_histories", "Closed-pore support", "closed_pore_histories_compact.csv", "Show high-density support state")

    success = point[point.classification == "success"]
    fig, ax = plt.subplots(figsize=(6, 4)); ax.axvspan(success.T2_C.min(), success.T2_C.max(), color="#009E73", alpha=.25, label="success band")
    ax.scatter(point.T2_C, point.rho_final, c=[COLORS[x] for x in point.classification], s=25)
    ax.axhline(.98, color="#555555", ls="--"); ax.set(xlabel="Second-step T₂ [°C]", ylabel="Final density, ρ", title=f"Candidate {best}: filled window"); ax.legend(); ps.clean(ax)
    save(fig, "Chen_filled_window_best", "Chen filled window", "chen_classification_points_compact.csv", "Show finite success interval")

    fig, ax = plt.subplots(figsize=(6, 4))
    for classification, group in point.groupby("classification"):
        ax.scatter(group.T2_C, np.ones(len(group)), label=classification.replace("_", " "), c=COLORS[classification], s=45)
    ax.set(xlabel="Second-step T₂ [°C]", ylabel="Classification row", yticks=[], title=f"Candidate {best}: lower/success/upper topology"); ax.legend(fontsize=7); ps.clean(ax)
    save(fig, "Chen_classification_map_best", "Chen classification map", "chen_classification_points_compact.csv", "Distinguish boundary failure modes")

    top_robust = robust[robust.candidate_id == best].pivot(index="rho0", columns="G0_nm", values="median_reduction")
    fig, ax = plt.subplots(figsize=(6, 4)); im = ax.imshow(top_robust.values, origin="lower", aspect="auto", cmap="viridis", extent=(top_robust.columns.min(), top_robust.columns.max(), top_robust.index.min(), top_robust.index.max()))
    ax.set(xlabel="Initial grain size G₀ [nm]", ylabel="Initial density ρ₀", title=f"Candidate {best}: robustness"); fig.colorbar(im, ax=ax, label="Median reduction")
    save(fig, "robustness_heatmap_rho0_G0", "Initial-condition robustness", "robustness_summary.csv", "Check bounded robustness")

    abl = ablations[ablations.candidate_id == best].copy(); abl["median_reduction"] = abl.median_reduction.fillna(0); abl = abl.sort_values("median_reduction")
    fig, ax = plt.subplots(figsize=(8, 6)); ax.barh(abl.ablation.str.replace("_", " "), abl.median_reduction, color=np.where(abl.tier.isin(["Tier_A", "Tier_B"]), "#009E73", "#D55E00"))
    ax.set(xlabel="Median two-step reduction", ylabel="Ablation", title=f"Candidate {best}: causal ablations"); ps.clean(ax)
    save(fig, "ablation_waterfall_best", "Causal ablation waterfall", "ablation_summary.csv", "Identify named controlling channels")

    fig, ax = plt.subplots(figsize=(6, 4)); x = np.arange(len(fast)); width = .35
    ax.bar(x - width / 2, fast.full_max_ratio, width, label="full material"); ax.bar(x + width / 2, np.where(fast.nucleation_facile_meaningful, fast.full_max_ratio, 1.0), width, label="nucleation facile")
    ax.axhline(1.5, color="#555555", ls="--"); ax.set(xticks=x, xticklabels=fast.material_id, ylabel="G reference / G fast", title="Frozen fast-firing envelope"); ax.legend(); ps.clean(ax)
    save(fig, "fast_firing_preservation", "Fast-firing preservation", "fast_firing_preservation.csv", "Verify frozen material causality")

    ps.write_inventory(OUT / "figure_inventory.csv", INVENTORY)


if __name__ == "__main__":
    main()
