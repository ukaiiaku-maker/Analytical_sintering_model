#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import plot_style as ps

ROOT=Path("results/separated_mechanism_production_search");FIG=ROOT/"production_figures"

def main():
    ps.apply_style();FIG.mkdir(parents=True,exist_ok=True);inventory=[]
    screen=pd.read_csv(ROOT/"fast_firing_material_screen.csv");registry=pd.read_csv(ROOT/"material_parameter_registry.csv");ab=pd.read_csv(ROOT/"fast_ablation_summary.csv");curves=pd.read_csv(ROOT/"fast_ablation_ratio_curves.csv");hist=pd.read_csv(ROOT/"fast_ablation_state_histories.csv")
    best=screen.sort_values("span_ge_1p5").groupby("material_id").tail(1).merge(registry,on="material_id",suffixes=("","_p"))
    fig,ax=plt.subplots(figsize=(6.2,4.6));sc=ax.scatter(best.Q_disconnection_nucleation/1e3,best.Q_surface_diffusion/1e3,c=best.span_ge_1p5,cmap="viridis",s=30,edgecolor="none");ax.set(xlabel=r"$Q_{\rm nuc}$ [kJ mol$^{-1}$]",ylabel=r"$Q_{\rm surface}$ [kJ mol$^{-1}$]");plt.colorbar(sc,ax=ax,label=r"Longest $G_{ref}/G_{fast}\geq1.5$ span");ps.clean(ax);ps.finish(fig,FIG/"Figure2_fast_material_phase_map");inventory.append(ps.inventory_row("2","Figure2_fast_material_phase_map","Fast-firing material phase map","fast_firing_material_screen.csv","Shows strict span before causal rejection","Results"))
    piv=ab.pivot(index="material_id",columns="ablation_mode",values="span_ge_1p5").fillna(0);cols=[c for c in ("full_material_model","no_PR_redistribution","no_nucleation_limitation","transport_only","exchange_limited_variant") if c in piv];fig,ax=plt.subplots(figsize=(7,3.8));im=ax.imshow(piv[cols],aspect="auto",cmap="magma",vmin=0,vmax=.17);ax.set_xticks(range(len(cols)),[c.replace("_","\n") for c in cols]);ax.set_yticks(range(len(piv)),piv.index);plt.colorbar(im,ax=ax,label=r"Span with ratio $\geq1.5$");ps.finish(fig,FIG/"Figure3_causal_ablation_matrix");inventory.append(ps.inventory_row("3","Figure3_causal_ablation_matrix","Causal ablation matrix","fast_ablation_summary.csv","Demonstrates failure of PR causal gate","Results"))
    full=ab[ab.ablation_mode=="full_material_model"].sort_values("max_ratio");mid=full.iloc[len(full)//2].material_id;q=curves[(curves.material_id==mid)&curves.ablation_mode.isin(("full_material_model","no_PR_redistribution","no_nucleation_limitation"))];fig,ax=plt.subplots(figsize=(6.2,4.2));
    for mode,g in q.groupby("ablation_mode"):ax.plot(g.rho,g.ratio,label=mode.replace("_"," "))
    for y in (1.2,1.5,2):ax.axhline(y,color="#777777",ls="--",lw=.8)
    ax.set(xlabel=r"Relative density, $\rho$",ylabel=r"$G_{reference}/G_{fast}$");ax.legend(fontsize=7);ps.clean(ax);ps.finish(fig,FIG/"Figure4_fast_trajectory_ablation");inventory.append(ps.inventory_row("4","Figure4_fast_trajectory_ablation","Representative fast-firing ablation","fast_ablation_ratio_curves.csv","Shows PR-off invariance and nucleation sensitivity","Results"))
    q=hist[(hist.material_id==mid)&(hist.ablation_mode=="full_material_model")];fig,axs=plt.subplots(1,2,figsize=(9,3.8));
    for path,g in q.groupby("path"):
        axs[0].semilogy(g.T_C,g.tau_nuc,label=f"{path} nucleation");axs[0].semilogy(g.T_C,g.tau_exchange,ls="--",label=f"{path} exchange");axs[0].semilogy(g.T_C,g.tau_transport,ls=":",label=f"{path} transport");axs[1].plot(g.rho,g.PR_exposure,label=path)
    axs[0].set(xlabel=r"Temperature [$^\circ$C]",ylabel="Serial time [s]");axs[1].set(xlabel=r"Relative density, $\rho$",ylabel="Cumulative PR exposure");axs[0].legend(fontsize=6);axs[1].legend();[ps.clean(a) for a in axs];ps.panel_labels(axs);ps.finish(fig,FIG/"Figure5_nucleation_PR_timing");inventory.append(ps.inventory_row("5","Figure5_nucleation_PR_timing","Nucleation and PR timing","fast_ablation_state_histories.csv","Diagnoses timing without asserting PR causality","Results"))
    ps.write_inventory(ROOT/"production_figure_inventory.csv",inventory)
if __name__=="__main__":main()
