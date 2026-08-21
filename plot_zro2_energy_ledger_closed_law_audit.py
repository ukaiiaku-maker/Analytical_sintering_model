"""Render diagnostic figures for the energy-ledger/closed-law audit."""
from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd

from derive_zro2_energy_ledger_and_closed_pore_laws import OUT, FIG

plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "savefig.bbox":"tight","figure.dpi":150})


def save(fig, name):
    FIG.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf","png"):
        fig.savefig(FIG/f"{name}.{ext}",dpi=220)
    plt.close(fig)


def schematic():
    fig,ax=plt.subplots(figsize=(9,5)); ax.axis("off")
    boxes={
        "available":(.08,.72,"Stored interfacial energy\n$E_{surf}+E_{GB}+E_{pore}$"),
        "gas":(.72,.72,"Gas free energy / work\n$E_{gas}$"),
        "dens":(.05,.16,"Open + closed\ndensification work"),
        "prep":(.29,.16,"Surface smoothing +\nPR/coarsening"),
        "growth":(.55,.16,"GB migration +\npore/junction drag"),
        "other":(.79,.16,"Named unresolved\nresidual"),
    }
    for x,y,label in boxes.values():
        ax.add_patch(FancyBboxPatch((x,y),.17,.12,boxstyle="round,pad=.02",fc="#eef4f8",ec="#294c60"))
        ax.text(x+.085,y+.06,label,ha="center",va="center")
    for key in ("dens","prep","growth","other"):
        x,y,_=boxes[key]
        ax.add_patch(FancyArrowPatch((.165,.72),(x+.085,y+.14),arrowstyle="->",mutation_scale=12,color="#555"))
    ax.add_patch(FancyArrowPatch((.72,.76),(.25,.76),arrowstyle="-|>",mutation_scale=12,color="#b24c3d"))
    ax.text(.5,.80,"gas counter-work",ha="center",color="#b24c3d")
    ax.text(.5,.94,"Diagnostic signed energy and dissipation ledger",ha="center",fontsize=14,weight="bold")
    ax.text(.5,.04,"No channel is silently rescaled; budget violations remain visible.",ha="center",style="italic")
    save(fig,"energy_ledger_schematic")


def law_comparison(scan):
    q=scan[(scan.phi_closed==.05)&(scan.accommodation==.5)&(scan.gas_fraction==.25)&(scan.hold_h==96)&(scan.exponent==4)&(scan.radius_nm.isin([10,25,100]))]
    modes=[m for m in q["mode"].unique() if m!="empirical_reduced_closure"]
    fig,axs=plt.subplots(1,2,figsize=(10,6))
    for mode in modes:
      for r,ls in zip((10,25,100),("-","--",":")):
        z=q[(q["mode"]==mode)&(q.radius_nm==r)].sort_values("T_C")
        label=f"{mode.replace('_',' ')}; {r} nm"
        axs[0].semilogy(z.T_C,np.maximum(z.tau_closed_s,1e-12),ls,label=label)
        axs[1].semilogy(z.T_C,np.maximum(z.density_gain,1e-16),ls,label=label)
    axs[0].set(ylabel=r"$\tau_{closed}$ (s)",xlabel="Temperature (°C)")
    axs[1].set(ylabel="96 h density gain",xlabel="Temperature (°C)")
    handles,labels=axs[1].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=6,ncol=3,loc="lower center",bbox_to_anchor=(.5,.01))
    fig.suptitle("Fixed-input closed-pore candidate comparison")
    fig.tight_layout(rect=(0,.22,1,.95)); save(fig,"closed_pore_law_comparison")


def reduced_overlay(scan):
    q=scan[(scan.radius_nm==25)&(scan.phi_closed==.05)&(scan.accommodation==.5)&(scan.gas_fraction==.25)&(scan.hold_h==96)&(scan.exponent==4)]
    fig,ax=plt.subplots(figsize=(8,4.5))
    ax.axhspan(1e-4,.15,color="#70ad47",alpha=.16,label="diagnostic finite-gain band")
    for mode,z in q.groupby("mode"):
        ax.semilogy(z.T_C,np.maximum(z.density_gain,1e-16),marker="o",ms=3,label=mode.replace("_"," "))
    ax.set(xlabel="Temperature (°C)",ylabel="96 h density gain",title="Physical candidates versus reduced works target")
    ax.legend(fontsize=6,loc="upper left"); fig.tight_layout()
    save(fig,"reduced_property_window_overlay")


def boundary_plot(boundary):
    order=["DENSIFICATION_EXHAUSTION_FAILURE","SUCCESS","GRAIN_GROWTH_FAILURE","MIXED_FAILURE"]
    colors={order[0]:"#4c78a8",order[1]:"#59a14f",order[2]:"#e15759",order[3]:"#b07aa1"}
    groups=list(boundary.groupby(["law_id","exponent"]))
    fig,axs=plt.subplots(len(groups),1,figsize=(8,max(5,1.35*len(groups))),sharex=True,squeeze=False)
    for ax,((mode,m),z) in zip(axs[:,0],groups):
        for cls in order:
            x=z[z.classification==cls]
            ax.scatter(x.T2_C,np.zeros(len(x)),s=45,marker="s",color=colors[cls],label=cls)
        ax.set_yticks([]); ax.set_ylabel(f"{mode}\nm={m}",rotation=0,ha="right",va="center",fontsize=7)
    axs[-1,0].set_xlabel("Second-step temperature (°C)")
    handles=[plt.Line2D([],[],marker="s",ls="",color=colors[c],label=c.replace("_"," ")) for c in order]
    fig.suptitle("Boundary preservation at the fixed selected first-step state",y=1.025)
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,1.005),ncol=2,fontsize=7)
    fig.tight_layout(rect=(0,0,1,.96)); save(fig,"boundary_preservation_fixed_state")


def ledger_histories(hist):
    fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=False)
    channels=["P_available_W_m3","P_open_dens_W_m3","P_closed_dens_W_m3","P_pore_coarsen_W_m3","P_GB_growth_W_m3","P_drag_W_m3","P_gas_W_m3"]
    for ax,((src,run),z) in zip(axs,hist.groupby(["source_history","run_id"])):
        x=z.physical_time_s/3600
        for c in channels:
            ax.semilogy(x,np.maximum(np.abs(z[c]),1e-12),label=c.replace("_W_m3",""))
        bad=z[z.budget_violation]
        if len(bad): ax.scatter(bad.physical_time_s/3600,np.maximum(np.abs(bad.P_available_W_m3),1e-12),marker="x",c="k",s=12,label="budget violation")
        ax.set(ylabel="|power| (W m$^{-3}$)",title=f"{run} — reconstructed ledger")
    axs[-1].set_xlabel("Physical time (h)")
    axs[0].legend(fontsize=6,ncol=3); fig.tight_layout(); save(fig,"energy_ledger_histories")


def inventory():
    rows=[]
    sources={
        "energy_ledger_schematic":"energy_ledger_channel_registry.csv",
        "closed_pore_law_comparison":"analytical_closed_law_rate_scan.csv",
        "reduced_property_window_overlay":"analytical_closed_law_rate_scan.csv",
        "boundary_preservation_fixed_state":"boundary_preservation_test.csv",
        "energy_ledger_histories":"energy_ledger_diagnostic_histories.csv",
    }
    for name,source in sources.items():
        rows.append({"figure_id":name,"pdf":f"figures/{name}.pdf","png":f"figures/{name}.png","source_table":source,
                     "pdf_nonempty":(FIG/f"{name}.pdf").stat().st_size>1000,"png_nonempty":(FIG/f"{name}.png").stat().st_size>1000,
                     "placeholder":False,"validation_claim":False})
    pd.DataFrame(rows).to_csv(OUT/"figure_inventory.csv",index=False)


def main():
    scan=pd.read_csv(OUT/"analytical_closed_law_rate_scan.csv")
    boundary=pd.read_csv(OUT/"boundary_preservation_test.csv")
    hist=pd.read_csv(OUT/"energy_ledger_diagnostic_histories.csv")
    schematic(); law_comparison(scan); reduced_overlay(scan); boundary_plot(boundary); ledger_histories(hist); inventory()
    print(f"wrote diagnostic figures to {FIG}")


if __name__=="__main__": main()
