#!/usr/bin/env python3
"""Publication-style figures for the exact, no-new-search mechanism synthesis."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd


OUT = Path("results/final_mechanism_synthesis_and_property_windows")
SRC = OUT / "source_tables"
FIG = OUT / "figures"
ATTR = Path("results/relative_material_property_window_attribution/source_tables")
PUB = Path("results/publication_style_sintering_figures_693168/source_tables")
AUDIT = Path("results/audit_candidate_693168_closed_accommodation")

BLUE, ORANGE, GREEN, GRAY, PURPLE = "#0072B2", "#D55E00", "#009E73", "#888888", "#CC79A7"
CLASS_COLORS = {"neither": GRAY, "fast_only": BLUE, "two_step_only": ORANGE, "both_pass": GREEN}


def setup() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
        "axes.labelsize": 9, "legend.fontsize": 8, "axes.spines.top": False,
        "axes.spines.right": False, "figure.constrained_layout.use": True,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def save(fig, stem: str, description: str, inventory: list[dict]) -> None:
    pdf, png = FIG / f"{stem}.pdf", FIG / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    inventory.append({"figure": stem, "pdf": str(pdf), "png": str(png),
                      "description": description, "evidence_level": "exact synthesis"})


def exact_merged() -> pd.DataFrame:
    exact = pd.read_csv(ATTR / "material_property_window_exact_promotions.csv")
    score = pd.read_csv(ATTR / "material_property_window_scorecard.csv", low_memory=False)
    return exact.merge(score, on="property_id", how="left", suffixes=("", "_screen"))


def bools(s: pd.Series) -> np.ndarray:
    return s.astype("boolean").fillna(False).to_numpy(bool)


def mechanism_chain(inventory: list[dict]) -> None:
    fig, axs = plt.subplots(1, 3, figsize=(11.0, 4.4))
    panels = [
        ("A  Fast firing", BLUE, ["Nucleation waiting", "Slow-ramp growth exposure", "Fast ramp crosses quickly"],
         "Nucleation-facile ablation weakens the effect"),
        ("B  Two-step", ORANGE, ["PR prepares closed topology", "Low T₂: density exhaustion", "Mid T₂: success", "High T₂: grain growth"],
         "Closed support and growth create two boundaries"),
        ("C  Shared framework", GREEN, ["Activated source/sink renewal", "Exchange + transport completion", "Matched-density attainment", "Bounded physical stores"],
         "Same perturbation vector; separate model layers"),
    ]
    for ax, (title, color, steps, note) in zip(axs, panels):
        ax.axis("off"); ax.set_title(title, color=color, weight="bold", fontsize=12)
        ys = np.linspace(.82, .25, len(steps))
        for i, (y, text) in enumerate(zip(ys, steps)):
            ax.add_patch(mpl.patches.FancyBboxPatch((.08, y-.07), .84, .13,
                         boxstyle="round,pad=.015", edgecolor=color,
                         facecolor=mpl.colors.to_rgba(color, .10), transform=ax.transAxes))
            ax.text(.50, y, text, ha="center", va="center", transform=ax.transAxes)
            if i < len(steps)-1:
                ax.annotate("", (.50, ys[i+1]+.07), (.50, y-.07), xycoords=ax.transAxes,
                            arrowprops=dict(arrowstyle="->", color="#444444"))
        ax.text(.5, .05, note, ha="center", va="center", fontsize=8, style="italic",
                wrap=True, transform=ax.transAxes)
    fig.suptitle("Compatible outcomes, different dominant mechanisms", fontsize=14)
    save(fig, "mechanism_chain_fast_vs_twostep", "Mechanistic chains for fast firing, two step, and shared constraints.", inventory)


def exact_phase_map(exact: pd.DataFrame, inventory: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    for c in ["neither", "fast_only", "two_step_only", "both_pass"]:
        d = exact[exact.classification_exact.eq(c)]
        label = c.replace("_", " ").title()
        ax.scatter(d.Theta_nuc, d.S_closed_growth, s=18, alpha=.63,
                   color=CLASS_COLORS[c], label=f"{label} (n={len(d)})")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.legend(ncol=2)
    ax.set(xlabel="Nucleation dominance, $\\Theta_{nuc}$",
           ylabel="Closed-shrinkage / growth selectivity, $S_{closed/growth}$",
           title="Exact-promoted relative-property phase map")
    ax.text(.99, .01, "Exact union only; no surrogate rows", ha="right", va="bottom",
            transform=ax.transAxes, fontsize=8, style="italic")
    save(fig, "relative_property_phase_map_exact", "Exact-only behavior phase map.", inventory)


def fast_window(exact: pd.DataFrame, inventory: list[dict]) -> None:
    ratio = pd.read_csv(PUB / "fast_firing_ratio_curves.csv")
    oat = pd.read_csv(SRC / "fast_firing_OAT_window.csv")
    q = oat[oat.perturbed_parameter.eq("Q_nuc_delta_kJ")].sort_values("perturbed_value")
    promoted = exact[exact.fast_firing_pass_exact.notna()]
    fig, axs = plt.subplots(2, 2, figsize=(9.3, 7.0))
    r = ratio[ratio.ablation_mode.eq("full_material_model")]
    axs[0,0].plot(r.rho, r.ratio, color=BLUE, lw=2.3, label="Base exact pass")
    axs[0,0].axhline(1.5, color="k", ls="--", lw=1, label="Pass threshold")
    axs[0,0].fill_between(r.rho, 1.5, r.ratio, where=r.ratio>=1.5, color=BLUE, alpha=.15)
    axs[0,0].legend(); axs[0,0].set(xlabel="Density, $\\rho$", ylabel="$G_{ref}/G_{fast}$",
                                     title="Representative exact pass trajectory")
    passed = bools(q.fast_firing_pass_exact)
    axs[0,1].plot(q.perturbed_value, q.R_fast_exact, color="#BBBBBB", lw=1)
    axs[0,1].scatter(q.perturbed_value[~passed], q.R_fast_exact[~passed], color=GRAY, label="Fail", zorder=3)
    axs[0,1].scatter(q.perturbed_value[passed], q.R_fast_exact[passed], color=BLUE, label="Pass", zorder=3)
    axs[0,1].axhline(1.5, color="k", ls="--", lw=1); axs[0,1].axvspan(0, 50, color=BLUE, alpha=.10)
    axs[0,1].legend(); axs[0,1].set(xlabel="$\\Delta Q_{nuc}$ (kJ mol$^{-1}$)", ylabel="Maximum $R_{fast}$",
                                     title="Finite nucleation-barrier OAT window")
    values = [bools(promoted.fast_firing_pass_exact).mean(),
              bools(promoted.nucleation_facile_pass_exact).mean(),
              bools(promoted.PR_off_pass_exact).mean()]
    bars = axs[1,0].bar(["Full", "Nucleation\nfacile", "PR off"], values,
                        color=[BLUE, "#E69F00", "#56B4E9"])
    axs[1,0].set_ylim(0, 1); axs[1,0].set(ylabel="Passing fraction of exact fast promotions",
                                          title="Ablation attribution")
    for bar, value in zip(bars, values):
        axs[1,0].text(bar.get_x()+bar.get_width()/2, value+.025, f"{100*value:.1f}%", ha="center")
    sc = axs[1,1].scatter(promoted.Theta_nuc, promoted.R_fast_exact,
                          c=bools(promoted.fast_firing_pass_exact),
                          cmap=ListedColormap(["#BBBBBB", BLUE]), s=14, alpha=.7)
    axs[1,1].set_xscale("log"); axs[1,1].axhline(1.5, color="k", ls="--", lw=1)
    axs[1,1].set(xlabel="Nucleation dominance, $\\Theta_{nuc}$", ylabel="$R_{fast}$",
                  title="Exact dominance envelope")
    save(fig, "fast_firing_property_window", "Fast-firing exact trajectory, OAT window, ablations, and dimensionless group.", inventory)


def two_step_window(exact: pd.DataFrame, inventory: list[dict]) -> None:
    matched = pd.read_csv(AUDIT / "dense_candidate_693168_matched_density_curves.csv")
    oat = pd.read_csv(SRC / "two_step_OAT_window.csv")
    q = oat[oat.perturbed_parameter.eq("Q_closed_delta_kJ")].sort_values("perturbed_value")
    kp = oat[oat.perturbed_parameter.eq("k_PR_factor")].sort_values("perturbed_value")
    et = exact[exact.two_step_pass_exact.notna()]
    fig, axs = plt.subplots(2, 3, figsize=(12.0, 7.0))
    m = matched[matched.both_paths_attained]
    axs[0,0].plot(m.rho, m.G_highT_nm, color="#333333", lw=2, label="High-T isothermal")
    axs[0,0].plot(m.rho, m.G_two_step_nm, color=ORANGE, lw=2, label="Two step")
    axs[0,0].set_yscale("log"); axs[0,0].legend(); axs[0,0].set(xlabel="Density, $\\rho$", ylabel="Grain size (nm)", title="Candidate 693168")
    axs[0,1].plot(m.rho, 100*m.reduction_TS, color=ORANGE, lw=2)
    axs[0,1].axhline(20, color="k", ls="--", lw=1); axs[0,1].set(xlabel="Density, $\\rho$", ylabel="Reduction (%)", title="Matched-density benefit")
    qp = bools(q.two_step_pass_exact)
    axs[0,2].plot(q.perturbed_value, 100*q.reduction_TS_exact, color="#BBBBBB")
    axs[0,2].scatter(q.perturbed_value[~qp], 100*q.reduction_TS_exact[~qp], color=GRAY, label="Fail")
    axs[0,2].scatter(q.perturbed_value[qp], 100*q.reduction_TS_exact[qp], color=ORANGE, label="Pass")
    axs[0,2].axvspan(-25, 100, color=ORANGE, alpha=.10); axs[0,2].legend()
    axs[0,2].set(xlabel="$\\Delta Q_{closed}$ (kJ mol$^{-1}$)", ylabel="Reduction (%)", title="Closed-barrier OAT window")
    kpp = bools(kp.two_step_pass_exact)
    axs[1,0].semilogx(kp.perturbed_value, kp.Chen_window_width_C_exact, color="#BBBBBB")
    axs[1,0].scatter(kp.perturbed_value[~kpp], kp.Chen_window_width_C_exact[~kpp], color=GRAY, label="Fail")
    axs[1,0].scatter(kp.perturbed_value[kpp], kp.Chen_window_width_C_exact[kpp], color=ORANGE, label="Pass")
    axs[1,0].axvline(.3, color="k", ls="--", lw=1); axs[1,0].legend()
    axs[1,0].set(xlabel="PR prefactor / base", ylabel="Chen-window width (°C)", title="PR preparation threshold")
    sc = axs[1,1].scatter(et.S_closed_growth, et.Chen_window_width_C_exact,
                          c=et.reduction_TS_exact, cmap="magma", s=18, alpha=.7)
    axs[1,1].set_xscale("log"); axs[1,1].set(xlabel="$S_{closed/growth}$", ylabel="Chen-window width (°C)", title="Selectivity and finite window")
    fig.colorbar(sc, ax=axs[1,1], label="Reduction fraction")
    axs[1,2].axis("off")
    text = ("LOW T₂\nclosed shrinkage / accommodation\ninsufficient to attain density\n\n"
            "INTERMEDIATE T₂\ndensity attained with bounded growth\n\n"
            "HIGH T₂\nthermally activated grain growth")
    axs[1,2].text(.05, .95, text, va="top", linespacing=1.45, transform=axs[1,2].transAxes)
    axs[1,2].set_title("Boundary interpretation")
    save(fig, "two_step_property_window", "Two-step exact trajectories, OAT windows, PR threshold, and selectivity.", inventory)


def chen_boundaries(inventory: list[dict]) -> None:
    d = pd.read_csv(AUDIT / "final_tables/candidate_693168_T2_classification_fine.csv")
    colors = {"DENSIFICATION_EXHAUSTION_FAILURE": "#56B4E9", "SUCCESS": GREEN, "GRAIN_GROWTH_FAILURE": PURPLE}
    code = {name: i for i, name in enumerate(colors)}
    fig, axs = plt.subplots(2, 2, figsize=(9.5, 7.0))
    for name, g in d.groupby("classification"):
        axs[0,0].scatter(g.T2_C, np.zeros(len(g)), marker="s", s=65, color=colors[name],
                         label=name.replace("_", " ").title())
    axs[0,0].axvline(925, color="k", ls="--", lw=1); axs[0,0].axvline(1205, color="k", ls="--", lw=1)
    axs[0,0].set_yticks([]); axs[0,0].legend(fontsize=7); axs[0,0].set(xlabel="Second-step temperature, T₂ (°C)", title="Fine exact Chen classification")
    for x, label in [(925, "Lower boundary"), (1205, "Upper boundary")]:
        axs[0,0].annotate(label, (x, 0), xytext=(0, 22), textcoords="offset points", ha="center", fontsize=8)
    axs[0,1].plot(d.T2_C, d.rho2, color=BLUE, lw=2); axs[0,1].axhline(.98, color="k", ls="--", lw=1, label="Target density")
    axs[0,1].set(xlabel="T₂ (°C)", ylabel="Final density", title="Lower boundary: density exhaustion"); axs[0,1].legend()
    axs[1,0].plot(d.T2_C, 100*d.growth_fraction, color=PURPLE, lw=2); axs[1,0].axhline(20, color="k", ls="--", lw=1, label="Tier-B growth bound")
    axs[1,0].set(xlabel="T₂ (°C)", ylabel="Second-step growth (%)", title="Upper boundary: growth activation"); axs[1,0].legend()
    axs[1,1].plot(d.T2_C, d.closed_shrinkage_contribution, color=ORANGE, lw=2, label="Closed shrinkage contribution")
    axs[1,1].plot(d.T2_C, d.closed_accommodation_factor, color=GREEN, lw=2, label="Accommodation factor")
    axs[1,1].set(xlabel="T₂ (°C)", ylabel="Fraction", title="Closed-store diagnostic"); axs[1,1].legend()
    save(fig, "chen_window_mechanism_boundaries", "Candidate 693168 fine Chen map with lower and upper diagnostics.", inventory)


def family_summary(inventory: list[dict]) -> None:
    f = pd.read_csv(SRC / "six_TierB_family_mechanism_summary.csv").sort_values("candidate_id")
    labels = f.candidate_id.astype(str); x = np.arange(len(f))
    fig, axs = plt.subplots(2, 2, figsize=(10.0, 7.0))
    axs[0,0].bar(x, 100*f.median_reduction, color=ORANGE); axs[0,0].set(ylabel="Median reduction (%)", title="High-density trajectory separation")
    axs[0,1].bar(x, f.window_width_C, color=GREEN); axs[0,1].set(ylabel="Window width (°C)", title="Finite Chen window")
    axs[1,0].bar(x-.18, 100*f.first_step_growth_fraction, .36, color=PURPLE, label="First-step growth")
    axs[1,0].bar(x+.18, 100*f.closed_fraction_at_switch, .36, color=BLUE, label="Closed fraction at switch")
    axs[1,0].set(ylabel="Percent", title="Preparation state"); axs[1,0].legend()
    axs[1,1].bar(x, f.destructive_ablation_count, color=GRAY)
    axs[1,1].set(ylabel="Number of destructive ablations", title="Mechanism dependence")
    for ax in axs.flat:
        ax.set_xticks(x, labels, rotation=30); ax.set_xlabel("Tier-B candidate")
    fig.suptitle("Six exact Tier-B base candidates; material robustness is a reduced transfer", fontsize=13)
    save(fig, "six_TierB_family_property_summary", "Six-candidate Tier-B mechanism and property summary.", inventory)


def surrogate_warning(inventory: list[dict]) -> None:
    d = pd.read_csv(SRC / "surrogate_vs_exact_comparison.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.3))
    x = np.arange(len(d)); width=.36
    axs[0].bar(x-width/2, d.surrogate_count, width, color="#BBBBBB", label="Surrogate screen")
    axs[0].bar(x+width/2, d.exact_count, width, color=[CLASS_COLORS[c] for c in d.classification], label="Exact promoted union")
    axs[0].set_yscale("log"); axs[0].set_xticks(x, d.classification.str.replace("_", " "), rotation=25)
    axs[0].set(ylabel="Case count (log scale)", title="Screening classifications are not final evidence"); axs[0].legend()
    both = d[d.classification.eq("both_pass")].iloc[0]
    bars = axs[1].bar(["Surrogate both", "Exact both"], [both.surrogate_count, both.exact_count], color=["#BBBBBB", GREEN])
    axs[1].set_yscale("log"); axs[1].set(ylabel="Both-pass count (log scale)", title=f"Overlap overprediction: {both.surrogate_to_exact_ratio:.1f}×")
    for bar, value in zip(bars, [both.surrogate_count, both.exact_count]):
        axs[1].text(bar.get_x()+bar.get_width()/2, value*1.18, f"{int(value):,}", ha="center")
    save(fig, "surrogate_vs_exact_warning", "Surrogate-to-exact discrepancy and overlap warning.", inventory)


def falsification_figure(inventory: list[dict]) -> None:
    d = pd.read_csv(SRC / "experimental_falsification_targets.csv")
    fig, ax = plt.subplots(figsize=(12.0, 6.5)); ax.axis("off")
    cell = [[r.measurement, r.mechanism_constrained, r.falsifiable_expected_signature] for _, r in d.iterrows()]
    table = ax.table(cellText=cell, colLabels=["Measurement", "Mechanism constrained", "Expected falsifiable signature"],
                     colWidths=[.26, .25, .49], cellLoc="left", loc="center")
    table.auto_set_font_size(False); table.set_fontsize(8); table.scale(1, 1.65)
    for (row, col), c in table.get_celld().items():
        c.set_edgecolor("#D0D0D0")
        if row == 0:
            c.set_facecolor("#E8F1F8"); c.set_text_props(weight="bold", color="#222222")
        elif row % 2 == 0:
            c.set_facecolor("#F7F7F7")
    ax.set_title("Measurements that can falsify the proposed mechanism attribution", fontsize=13, pad=18)
    ax.text(.5, .01, "Primary calibration targets: closed-pore fraction and accommodation trajectory",
            ha="center", weight="bold", color=ORANGE, transform=ax.transAxes)
    save(fig, "experimental_falsification_targets", "Experimental measurements, constrained mechanisms, and expected signatures.", inventory)


def main() -> None:
    setup(); inventory=[]; exact=exact_merged()
    mechanism_chain(inventory)
    exact_phase_map(exact, inventory)
    fast_window(exact, inventory)
    two_step_window(exact, inventory)
    chen_boundaries(inventory)
    family_summary(inventory)
    surrogate_warning(inventory)
    falsification_figure(inventory)
    pd.DataFrame(inventory).to_csv(SRC / "final_figure_inventory.csv", index=False)
    print(f"wrote {len(inventory)} synthesis figures as PDF and 600 dpi PNG")


if __name__ == "__main__":
    main()
