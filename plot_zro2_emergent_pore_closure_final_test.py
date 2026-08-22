"""Diagnostic figures for the emergent pore-closure final promotion audit."""
from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results/zro2_forward_emergent_pore_closure_final_test"
FIG=OUT/"figures"
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150})


def save(fig,name):
    FIG.mkdir(parents=True,exist_ok=True)
    for ext in ("pdf","png"): fig.savefig(FIG/f"{name}.{ext}",dpi=220,bbox_inches="tight")
    plt.close(fig)


def schematic():
    fig,ax=plt.subplots(figsize=(10,5.2)); ax.axis("off")
    nodes=[(.05,.62,"Open pore bins"),(.30,.62,"Precursor / isolated"),(.56,.62,"Closed pore bins"),(.78,.62,"Named closed\nshrinkage"),
           (.30,.18,"PR / surface\npreparation"),(.56,.18,"Energy ledger\nand gas work"),(.78,.18,"GB growth, pinning\nand drag")]
    for x,y,s in nodes:
        ax.add_patch(FancyBboxPatch((x,y),.17,.14,boxstyle="round,pad=.02",fc="#edf3f7",ec="#31576b")); ax.text(x+.085,y+.07,s,ha="center",va="center")
    arrows=[((.22,.69),(.30,.69)),((.47,.69),(.56,.69)),((.73,.69),(.78,.69)),((.135,.62),(.36,.32)),((.385,.32),(.64,.62)),
            ((.645,.62),(.645,.32)),((.73,.25),(.78,.25))]
    for a,b in arrows: ax.add_patch(FancyArrowPatch(a,b,arrowstyle="->",mutation_scale=12,color="#555"))
    ax.text(.5,.91,"emergent_pore_closure_v1 — diagnostic promotion architecture",ha="center",fontsize=14,weight="bold")
    ax.text(.5,.05,"Only open and closed shrinkage remove pore volume; preparation and smoothing are conservative.",ha="center",style="italic")
    save(fig,"emergent_closure_model_schematic")


def rate_scan(scan):
    q=scan[(scan.phi_closed==.05)&(scan.A_closed==.5)&(scan.gas_fraction==.25)&(scan.hold_h==96)&(scan.radius_exponent==3)]
    fig,axs=plt.subplots(1,2,figsize=(10,4.5),sharey=True)
    for ax,kernel in zip(axs,("renewal","GB_diffusion")):
      for r in (5,10,25,50,100,250):
        z=q[(q.kernel==kernel)&(q.radius_nm==r)].sort_values("T_C")
        ax.semilogy(z.T_C,np.maximum(z.density_gain,1e-16),label=f"{r} nm")
      ax.set(xlabel="Temperature (°C)",title=kernel.replace("_"," "))
    axs[0].set_ylabel("96 h closed-density gain"); axs[1].legend(fontsize=7,ncol=2)
    fig.suptitle("Closed-law rate scan (fixed inventory and accommodation)"); fig.tight_layout(); save(fig,"closed_law_rate_scan")


def gas_radius(scan):
    q=scan[(scan.T_C==1100)&(scan.phi_closed==.05)&(scan.A_closed==.5)&(scan.hold_h==96)&(scan.radius_exponent==3)&(scan.kernel=="renewal")]
    fig,axs=plt.subplots(1,2,figsize=(9,4))
    for r in (5,10,25,50,100,250):
      z=q[q.radius_nm==r].sort_values("gas_fraction")
      axs[0].plot(z.gas_fraction,z.sigma_c_Pa/1e6,marker="o",label=f"{r} nm")
      axs[1].semilogy(z.gas_fraction,np.maximum(z.rho_dot_closed_sinv,1e-30),marker="o",label=f"{r} nm")
    axs[0].set(xlabel=r"$P_{gas}/P_{cap}$",ylabel=r"$\sigma_c$ (MPa)")
    axs[1].set(xlabel=r"$P_{gas}/P_{cap}$",ylabel=r"$\dot\rho_c$ (s$^{-1}$)")
    axs[1].legend(fontsize=7,ncol=2); fig.suptitle("Gas counterpressure and pore-radius limits"); fig.tight_layout(); save(fig,"gas_pressure_and_radius_limits")


def ledger_plot(ledger):
    groups=list(ledger.groupby(["source_history","run_id"]))[:3]
    fig,axs=plt.subplots(len(groups),1,figsize=(9,7),squeeze=False)
    cols=["P_available_W_m3","P_open_dens_W_m3","P_closed_dens_W_m3","P_PR_W_m3","P_GB_growth_W_m3","P_drag_W_m3","P_gas_W_m3"]
    for ax,((src,run),z) in zip(axs[:,0],groups):
      for c in cols: ax.semilogy(z.physical_time_s/3600,np.maximum(abs(z[c]),1e-12),label=c.replace("_W_m3",""))
      bad=z[z.budget_violation]; ax.scatter(bad.physical_time_s/3600,np.maximum(abs(bad.P_available_W_m3),1e-12),s=10,c="k",marker="x")
      ax.set(ylabel="|power| (W m$^{-3}$)",title=str(run))
    axs[-1,0].set_xlabel("Physical time (h)"); axs[0,0].legend(fontsize=6,ncol=4)
    fig.suptitle("Selected-path energy ledger (violations retained)"); fig.tight_layout(); save(fig,"energy_ledger_for_selected_paths")


COLORS={"DENSIFICATION_EXHAUSTION_FAILURE":"#4c78a8","SUCCESS":"#59a14f","GRAIN_GROWTH_FAILURE":"#e15759","MIXED_FAILURE":"#b07aa1"}


def boundary_plot(boundary):
    combos=list(boundary.groupby(["state_id","kernel","m"]))
    fig,axs=plt.subplots(3,1,figsize=(10,8),sharex=True)
    for ax,(state_id,s) in zip(axs,boundary.groupby("state_id",sort=False)):
      for j,((kernel,m),z) in enumerate(s.groupby(["kernel","m"])):
       for cls,c in COLORS.items():
        x=z[z.classification==cls]; ax.scatter(x.T2_C,np.full(len(x),j),s=24,marker="s",color=c)
      ax.set_yticks(range(4)); ax.set_yticklabels([f"{k}, m={m}" for k,m in s.groupby(["kernel","m"]).groups],fontsize=7); ax.set_title(state_id)
    axs[-1].set_xlabel("Second-step temperature (°C)")
    handles=[plt.Line2D([],[],marker="s",ls="",color=c,label=k.replace("_"," ")) for k,c in COLORS.items()]
    fig.legend(handles=handles,ncol=2,loc="upper center",bbox_to_anchor=(.5,.985),fontsize=7)
    fig.suptitle("Fixed-state boundary preservation",y=1.02); fig.tight_layout(rect=(0,0,1,.94)); save(fig,"boundary_preservation_test")


def failure_maps(boundary):
    q=boundary[boundary.state_id=="actual_selected"].copy()
    codes={"DENSIFICATION_EXHAUSTION_FAILURE":0,"MIXED_FAILURE":1,"GRAIN_GROWTH_FAILURE":2,"SUCCESS":3}
    q["code"]=q.classification.map(codes)
    for name,ylabel in (("final_process_map_G1_T2_failure_modes","First-step grain size (nm)"),("final_process_map_rho_switch_T2_failure_modes","Switch density")):
      fig,axs=plt.subplots(2,2,figsize=(10,6),sharex=True,sharey=True)
      for ax,((kernel,m),z) in zip(axs.flat,q.groupby(["kernel","m"])):
        y=np.full(len(z),z.G_nm.iloc[0] if "G1" in name else z.rho.iloc[0])
        ax.scatter(z.T2_C,y,c=[COLORS[x] for x in z.classification],marker="s",s=35)
        ax.set_title(f"{kernel}, m={m}"); ax.set_xlabel("T2 (°C)"); ax.set_ylabel(ylabel)
      fig.suptitle("Naturally prepared state: failure modes only (no process map run)")
      fig.tight_layout(); save(fig,name)


def ablation_plot(ablation):
    z=ablation.sort_values("density_gain")
    fig,ax=plt.subplots(figsize=(9,6)); c=["#59a14f" if x else "#e15759" for x in z.mechanism_signature_supported]
    ax.barh(z.ablation,z.density_gain,color=c); ax.axvline(.08,color="k",ls="--",lw=1,label="gain needed from rho=0.82 to 0.90")
    ax.set(xlabel="96 h density gain at 1100 °C",title="Emergent-closure ablation matrix")
    ax.legend(fontsize=7); fig.tight_layout(); save(fig,"ablation_matrix")


def inventory():
    items={"emergent_closure_model_schematic":"PR_preparation_flux_registry.csv","closed_law_rate_scan":"final_closed_law_rate_scan.csv",
           "gas_pressure_and_radius_limits":"final_closed_law_rate_scan.csv","energy_ledger_for_selected_paths":"final_energy_ledger_selected_paths.csv",
           "boundary_preservation_test":"final_boundary_preservation_test.csv","final_process_map_G1_T2_failure_modes":"final_process_map_failure_modes.csv",
           "final_process_map_rho_switch_T2_failure_modes":"final_process_map_failure_modes.csv","ablation_matrix":"final_emergent_closure_ablation_matrix.csv"}
    rows=[]
    for name,source in items.items():
      rows.append({"figure_id":name,"pdf":f"figures/{name}.pdf","png":f"figures/{name}.png","source_table":source,
                   "pdf_nonempty":(FIG/f"{name}.pdf").stat().st_size>1000,"png_nonempty":(FIG/f"{name}.png").stat().st_size>1000,
                   "placeholder":False,"success_colored_process_map":False,"validation_claim":False})
    pd.DataFrame(rows).to_csv(OUT/"figure_inventory.csv",index=False)


def main():
    scan=pd.read_csv(OUT/"final_closed_law_rate_scan.csv"); ledger=pd.read_csv(OUT/"final_energy_ledger_selected_paths.csv")
    boundary=pd.read_csv(OUT/"final_boundary_preservation_test.csv"); ablation=pd.read_csv(OUT/"final_emergent_closure_ablation_matrix.csv")
    schematic(); rate_scan(scan); gas_radius(scan); ledger_plot(ledger); boundary_plot(boundary); failure_maps(boundary); ablation_plot(ablation); inventory()
    print(f"wrote diagnostic figures to {FIG}")


if __name__=="__main__": main()
