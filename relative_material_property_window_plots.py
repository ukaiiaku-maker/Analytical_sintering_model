#!/usr/bin/env python3
"""Create attribution figures for the frozen-topology material-window audit."""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd


OUT = Path("results/relative_material_property_window_attribution")
SRC = OUT / "source_tables"
FIG = OUT / "figures"
PUB = Path("results/publication_style_sintering_figures_693168/source_tables")

COLORS = {
    "fast": "#0072B2", "two": "#D55E00", "both": "#009E73",
    "neither": "#8A8A8A", "lower": "#56B4E9", "upper": "#CC79A7",
}


def setup():
    FIG.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
        "axes.labelsize": 9, "legend.fontsize": 8, "axes.spines.top": False,
        "axes.spines.right": False, "figure.constrained_layout.use": True,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def save(fig, stem, description, inventory):
    pdf = FIG / f"{stem}.pdf"
    png = FIG / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    inventory.append({"figure": stem, "pdf": str(pdf), "png": str(png),
                      "description": description, "evidence_level": "exact-first attribution"})


def load():
    score = pd.read_csv(SRC / "material_property_window_scorecard.csv", low_memory=False)
    ep = pd.read_csv(SRC / "material_property_window_exact_promotions.csv")
    exact = ep.merge(score, on="property_id", how="left", suffixes=("", "_screen"))
    rank = pd.read_csv(SRC / "material_property_sensitivity_rankings.csv")
    thresh = pd.read_csv(SRC / "dimensionless_thresholds.csv")
    fam = pd.read_csv(SRC / "tierB_family_material_window_comparison.csv")
    attr = pd.read_csv(SRC / "mechanism_attribution_summary.csv")
    ratio = pd.read_csv(PUB / "fast_firing_ratio_curves.csv")
    histories = pd.read_csv(PUB / "dense_time_histories.csv")
    matched = pd.read_csv(PUB / "dense_matched_density_curves.csv")
    return score, exact, rank, thresh, fam, attr, ratio, histories, matched


def numeric_bool(series):
    return series.astype("boolean").fillna(False).to_numpy(bool)


def binned_mean(ax, x, y, z, xlabel, ylabel, title, cmap="viridis"):
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = np.asarray(x)[mask], np.asarray(y)[mask], np.asarray(z)[mask]
    xb = np.linspace(x.min(), x.max(), 24)
    yb = np.linspace(y.min(), y.max(), 24)
    count, _, _ = np.histogram2d(x, y, bins=(xb, yb))
    total, _, _ = np.histogram2d(x, y, bins=(xb, yb), weights=z)
    mean = np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)
    mesh = ax.pcolormesh(xb, yb, mean.T, shading="auto", cmap=cmap)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    return mesh


def figure_ingredient(attr, inventory):
    fig, ax = plt.subplots(figsize=(9.0, 4.7)); ax.axis("off")
    ax.text(.02, .94, "Distinct material-property chains, shared attainment constraint",
            fontsize=13, weight="bold", transform=ax.transAxes)
    boxes = [
        (.04, .62, .24, .20, "FAST FIRING", "Nucleation-limited waiting\n+ ramp duration through low activity", COLORS["fast"]),
        (.38, .62, .24, .20, "SHARED", "Matched-density attainment\n bounded physical states", COLORS["both"]),
        (.72, .62, .24, .20, "TWO STEP", "PR-prepared closed store\n closed shrinkage vs migration", COLORS["two"]),
        (.04, .16, .24, .22, "FAST ORDERING", "nucleation onset relative to\nexchange / transport / growth", COLORS["fast"]),
        (.38, .16, .24, .22, "JOINT RESULT", "same perturbation vector passes\ntwo separate model layers", COLORS["both"]),
        (.72, .16, .24, .22, "TWO-STEP ORDERING", "closed accommodation support\nrelative to growth recovery", COLORS["two"]),
    ]
    for x, y, w, h, head, body, c in boxes:
        ax.add_patch(mpl.patches.FancyBboxPatch((x,y),w,h,boxstyle="round,pad=.015",
                     edgecolor=c, facecolor=mpl.colors.to_rgba(c,.10), linewidth=1.8,
                     transform=ax.transAxes))
        ax.text(x+.015,y+h-.045,head,weight="bold",color=c,transform=ax.transAxes)
        ax.text(x+.015,y+.045,body,va="bottom",transform=ax.transAxes)
    for x0, x1 in [(.28,.38),(.62,.72)]:
        ax.annotate("",(x1,.72),(x0,.72),xycoords=ax.transAxes,
                    arrowprops=dict(arrowstyle="->",lw=1.4,color="#444444"))
    ax.text(.5,.04,"Topology frozen in the primary scan; no shared dynamic coupling was introduced.",
            ha="center",style="italic",transform=ax.transAxes)
    save(fig,"mechanism_ingredient_summary","Conceptual separation of fast-firing and two-step ingredients.",inventory)


def figure_fast(exact, ratio, inventory):
    fig, axs = plt.subplots(2, 2, figsize=(9.2, 7.0))
    d = ratio[ratio.ablation_mode.eq("full_material_model")]
    axs[0,0].plot(d.rho,d.ratio,color=COLORS["fast"],lw=2)
    axs[0,0].axhline(1.5,color="k",ls="--",lw=1); axs[0,0].fill_between(d.rho,1.5,d.ratio,where=d.ratio>=1.5,color=COLORS["fast"],alpha=.15)
    axs[0,0].set(xlabel="Density, $\\rho$",ylabel="$G_{ref}/G_{fast}$",title="Matched-density fast-firing separation")
    fmask = exact.fast_firing_pass_exact.notna()
    p = exact[fmask]
    sc=axs[0,1].scatter(p.Theta_nuc,p.R_fast_exact,c=numeric_bool(p.fast_firing_pass_exact),cmap=ListedColormap(["#BBBBBB",COLORS["fast"]]),s=12,alpha=.7)
    axs[0,1].set_xscale("log"); axs[0,1].axhline(1.5,color="k",ls="--",lw=1)
    axs[0,1].set(xlabel="$\\Theta_{nuc}=\\tau_{nuc}/(\\tau_{ex}+\\tau_{tr})$",ylabel="$R_{fast}$",title="Exact promoted material perturbations")
    axs[1,0].scatter(p.I_low_slow,p.R_fast_exact,c=p.f_nuc,cmap="viridis",s=13,alpha=.7)
    axs[1,0].set_xscale("log"); axs[1,0].set(xlabel="Slow-path low-activity exposure, $I_{low}$ (s)",ylabel="$R_{fast}$",title="Waiting exposure and nucleation fraction")
    vals=[numeric_bool(p.fast_firing_pass_exact).mean(),numeric_bool(p.nucleation_facile_pass_exact).mean(),numeric_bool(p.PR_off_pass_exact).mean()]
    bars=axs[1,1].bar(["Full","Nucleation\nfacile","PR off"],vals,color=[COLORS["fast"],"#E69F00","#56B4E9"])
    axs[1,1].set_ylim(0,1);axs[1,1].set(ylabel="Fraction passing among exact fast promotions",title="Required ablation direction")
    for b,v in zip(bars,vals):axs[1,1].text(b.get_x()+b.get_width()/2,v+.025,f"{100*v:.1f}%",ha="center")
    save(fig,"fast_firing_mechanism_attribution","Exact fast-firing trajectory and nucleation/PR attribution.",inventory)


def figure_two(hist, matched, inventory):
    fig, axs = plt.subplots(2,2,figsize=(9.2,7.0))
    s=hist[(hist.path_label=="success") & (hist.record_kind=="solver_state")]
    h=hist[(hist.path_label=="highT_reference") & (hist.record_kind=="solver_state")]
    axs[0,0].plot(h.rho,h.G_mean_nm,label="1400 °C isothermal",color="#333333",lw=2)
    axs[0,0].plot(s.rho,s.G_mean_nm,label="1400→1100 °C",color=COLORS["two"],lw=2)
    axs[0,0].set_yscale("log");axs[0,0].legend();axs[0,0].set(xlabel="Density, $\\rho$",ylabel="Mean grain size (nm)",title="Two-step trajectory separation")
    q=matched[matched.both_paths_attained]
    axs[0,1].plot(q.rho,100*q.reduction_TS,color=COLORS["two"],lw=2);axs[0,1].axhline(20,color="k",ls="--",lw=1)
    axs[0,1].set(xlabel="Density, $\\rho$",ylabel="Grain-size reduction (%)",title="Attained high-density interval")
    axs[1,0].plot(s.rho,s.closed_fraction,label="Closed pore fraction",color=COLORS["two"])
    axs[1,0].plot(s.rho,s.closed_accommodation_factor,label="Accommodation available",color=COLORS["both"])
    axs[1,0].plot(s.rho,s.closed_pore_contribution_to_rho_dot,label="Closed $\\dot\\rho$ share",color=COLORS["fast"])
    axs[1,0].legend();axs[1,0].set(xlabel="Density, $\\rho$",ylabel="Fraction",title="Closed-store support during successful path")
    labels=["Lower\n900 °C","Success\n1100 °C","Upper\n1220 °C"]
    growth=[];finalrho=[]
    for path in ["lower_failure","success","upper_failure"]:
        d=hist[(hist.path_label==path)&(hist.stage=="second_step")]
        growth.append(float(d.growth_fraction_from_switch.max()))
        finalrho.append(float(d.rho.max()))
    x=np.arange(3);axs[1,1].bar(x-.18,finalrho,.36,label="$\\rho_{final}$",color="#56B4E9")
    axs[1,1].bar(x+.18,growth,.36,label="Second-step growth",color="#CC79A7")
    axs[1,1].axhline(.98,color="#56B4E9",ls="--",lw=1);axs[1,1].set_xticks(x,labels);axs[1,1].legend();axs[1,1].set(title="Bracketed lower and upper boundaries",ylabel="Fraction")
    save(fig,"two_step_mechanism_attribution","Candidate 693168 two-step trajectory, closed support, and finite boundaries.",inventory)


def figure_barriers(exact, inventory):
    ef=exact[exact.fast_firing_pass_exact.notna()]
    fig,axs=plt.subplots(1,2,figsize=(9.3,3.9))
    m=binned_mean(axs[0],ef.Q_nuc_minus_Q_growth_kJ,ef.Q_nuc_minus_Q_transport_kJ,ef.span_fast_1p5_exact,
                  "$Q_{nuc}-Q_{growth}$ (kJ mol$^{-1}$)","$Q_{nuc}-Q_{transport}$ (kJ mol$^{-1}$)","Exact fast-firing density span")
    fig.colorbar(m,ax=axs[0],label="$\\Delta\\rho$ with ratio ≥1.5")
    m=binned_mean(axs[1],ef.Q_nuc_minus_Q_growth_kJ,ef.Q_nuc_minus_Q_PR_kJ,ef.R_fast_exact,
                  "$Q_{nuc}-Q_{growth}$ (kJ mol$^{-1}$)","$Q_{nuc}-Q_{PR}$ (kJ mol$^{-1}$)","Maximum matched-density ratio")
    fig.colorbar(m,ax=axs[1],label="$R_{fast}$")
    save(fig,"relative_barrier_window_fast_firing","Exact promoted fast-firing relative-barrier maps.",inventory)
    et=exact[exact.two_step_pass_exact.notna()]
    fig,axs=plt.subplots(1,2,figsize=(9.3,3.9))
    m=binned_mean(axs[0],et.Q_closed_minus_Q_growth_kJ,et.Q_PR_minus_Q_closed_kJ,et.reduction_TS_exact,
                  "$Q_{closed}-Q_{growth}$ (kJ mol$^{-1}$)","$Q_{PR}-Q_{closed}$ (kJ mol$^{-1}$)","Exact high-density reduction")
    fig.colorbar(m,ax=axs[0],label="Reduction fraction")
    m=binned_mean(axs[1],et.log10_kclosed_over_kgrowth,et.log10_kPR_over_kgrowth,et.Chen_window_width_C_exact,
                  "$\\log_{10}(k_{closed}/k_{growth})$","$\\log_{10}(k_{PR}/k_{growth})$","Exact Chen-window width")
    fig.colorbar(m,ax=axs[1],label="Window width (°C)")
    save(fig,"relative_barrier_window_two_step","Exact promoted two-step relative-barrier and prefactor maps.",inventory)


def figure_joint(exact, inventory):
    fig,ax=plt.subplots(figsize=(7.3,5.2))
    classes=["neither","fast_only","two_step_only","both_pass"]
    labels={"both_pass":"Both","fast_only":"Fast only","two_step_only":"Two-step only","neither":"Neither"}
    class_colors={"neither":COLORS["neither"],"fast_only":COLORS["fast"],
                  "two_step_only":COLORS["two"],"both_pass":COLORS["both"]}
    for c in classes:
        d=exact[exact.classification_exact.eq(c)]
        ax.scatter(d.Theta_nuc,d.S_closed_growth,s=18,alpha=.65,label=f"{labels[c]} (n={len(d)})",color=class_colors[c])
    ax.set_xscale("log");ax.set_yscale("log");ax.legend(ncol=2)
    ax.set(xlabel="$\\Theta_{nuc}$ at fast-firing reference state",ylabel="$S_{closed/growth}$ at second step",title="Exact joint behavior phase map (separate model layers)")
    save(fig,"joint_behavior_phase_map","Exact classifications in fast and two-step dimensionless coordinates.",inventory)


def figure_tornado(rank, inventory):
    for target,stem,title,color in [
        ("R_fast_exact","sensitivity_tornado_fast_firing","Fast-firing exact rank sensitivity",COLORS["fast"]),
        ("Chen_window_width_C_exact","sensitivity_tornado_two_step","Two-step exact rank sensitivity",COLORS["two"]),
    ]:
        d=rank[(rank.target==target)&(rank.evidence_level=="exact")].dropna().head(12).sort_values("absolute_rank_score")
        fig,ax=plt.subplots(figsize=(7.4,4.8));colors=[color if v>=0 else "#999999" for v in d.spearman_r]
        ax.barh(d.feature,d.spearman_r,color=colors);ax.axvline(0,color="k",lw=.8)
        ax.set(xlabel="Spearman rank correlation",title=title)
        ax.text(.01,.01,"Exact promoted subset; association is not sufficiency",transform=ax.transAxes,fontsize=8,style="italic")
        save(fig,stem,title+" across exact promotions.",inventory)


def figure_closed(exact, inventory):
    d=exact[exact.two_step_pass_exact.notna()].copy()
    fig,axs=plt.subplots(1,2,figsize=(9.2,4.0))
    sc=axs[0].scatter(d.closed_fraction_switch_exact,d.closed_accommodation_fraction_exact,
                     c=d.reduction_TS_exact,cmap="magma",s=20,alpha=.7)
    axs[0].set(xlabel="Closed fraction at switch",ylabel="Accommodation available / capacity",title="High-density reduction support")
    fig.colorbar(sc,ax=axs[0],label="Reduction fraction")
    sc=axs[1].scatter(d.A_closed_fraction,d.S_closed_growth,c=d.Chen_window_width_C_exact,cmap="viridis",s=20,alpha=.7)
    axs[1].set_yscale("log");axs[1].set(xlabel="Reference accommodation fraction",ylabel="$S_{closed/growth}$",title="Finite Chen-window response")
    fig.colorbar(sc,ax=axs[1],label="Window width (°C)")
    save(fig,"closed_accommodation_plausibility_window","Closed-store/accommodation material-property plausibility window.",inventory)


def figure_family(fam, inventory):
    f=fam.sort_values("transferred_material_margin",ascending=False)
    fig,axs=plt.subplots(1,2,figsize=(9.2,4.3));x=np.arange(len(f));labels=f.candidate_id.astype(str)
    axs[0].bar(x-.2,f.predicted_Qclosed_window_kJ_mol,.4,label="$Q_{closed}$",color=COLORS["two"])
    axs[0].bar(x+.2,f.predicted_Qgrowth_window_kJ_mol,.4,label="$Q_{growth}$",color=COLORS["fast"])
    axs[0].set_xticks(x,labels,rotation=35);axs[0].set(ylabel="Transferred half-range proxy (kJ mol$^{-1}$)",title="Reduced family transfer",xlabel="Tier-B candidate");axs[0].legend()
    sc=axs[1].scatter(f.closed_fraction_at_switch,f.base_reduction,s=50+20*f.robustness_count,
                     c=f.transferred_material_margin,cmap="viridis",vmin=0,vmax=1)
    for _,r in f.iterrows():axs[1].annotate(str(int(r.candidate_id)),(r.closed_fraction_at_switch,r.base_reduction),xytext=(3,3),textcoords="offset points",fontsize=7)
    axs[1].set(xlabel="Closed fraction at switch",ylabel="Base median reduction",title="Family topology and transferred margin")
    fig.colorbar(sc,ax=axs[1],label="Transferred material margin")
    fig.text(.5,.005,"Reduced transfer from exact 693168 OAT; not an exact six-candidate material sweep.",ha="center",fontsize=8,style="italic")
    save(fig,"six_TierB_material_window_comparison","Reduced material-window transfer across six exact Tier-B candidates.",inventory)


def figure_trends(inventory):
    fig,axs=plt.subplots(1,3,figsize=(10.0,4.2))
    panels=[
      ("Fast firing",["Delayed densification onset","Short low-activity exposure on fast ramp","Exchange/transport still attain density"],COLORS["fast"]),
      ("Two step",["First step prepares closed store","Closed shrinkage survives at low T₂","Migration recovers at upper T₂"],COLORS["two"]),
      ("Shared tests",["Matched-density support","Finite bounded states","Independent kinetic calibration"],COLORS["both"]),
    ]
    for ax,(head,items,c) in zip(axs,panels):
        ax.axis("off");ax.set_title(head,color=c,weight="bold",fontsize=12)
        for i,item in enumerate(items):
            y=.78-i*.27
            ax.add_patch(mpl.patches.FancyBboxPatch((.05,y-.1),.90,.16,boxstyle="round,pad=.02",facecolor=mpl.colors.to_rgba(c,.10),edgecolor=c,transform=ax.transAxes))
            ax.text(.50,y-.02,item,ha="center",va="center",wrap=True,transform=ax.transAxes)
    fig.suptitle("Qualitative relative-property trends for crystalline particle systems",fontsize=13)
    fig.text(.5,.01,"These are falsifiable orderings and timescale separations, not universal numerical constants.",ha="center",style="italic")
    save(fig,"common_particle_system_trends","Qualitative cross-system trends supported by the attribution audit.",inventory)


def main():
    setup(); inventory=[]
    score, exact, rank, thresh, fam, attr, ratio, histories, matched = load()
    figure_ingredient(attr,inventory)
    figure_fast(exact,ratio,inventory)
    figure_two(histories,matched,inventory)
    figure_barriers(exact,inventory)
    figure_joint(exact,inventory)
    figure_tornado(rank,inventory)
    figure_closed(exact,inventory)
    figure_family(fam,inventory)
    figure_trends(inventory)
    pd.DataFrame(inventory).to_csv(SRC/"figure_inventory.csv",index=False)
    print(f"wrote {len(inventory)} figures as PDF and 600 dpi PNG")


if __name__ == "__main__":
    main()
