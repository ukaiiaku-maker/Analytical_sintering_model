#!/usr/bin/env python3
"""Plot the frozen relative-property evidence and its explicit coverage gaps."""
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

OUT=Path("results/reframe_tierB_experimental_plausibility");FIG=OUT/"figures"


def main():
    d=pd.read_csv(OUT/"relative_material_property_window_reframed.csv")
    q=d[d.material_id.isin(("E0021","E0142"))].copy();FIG.mkdir(parents=True,exist_ok=True)
    fig,axs=plt.subplots(1,2,figsize=(9.2,4.2))
    colors=q.two_step_pass.map({True:"#009E73",False:"#D55E00"})
    axs[0].scatter(q.Q_nuc_minus_Q_GB_kJ_mol,q.Q_surface_minus_Q_GB_kJ_mol,c=colors,s=110,edgecolor="k")
    for _,r in q.iterrows():axs[0].annotate(r.material_id,(r.Q_nuc_minus_Q_GB_kJ_mol,r.Q_surface_minus_Q_GB_kJ_mol),xytext=(5,5),textcoords="offset points")
    axs[0].set(xlabel=r"$Q_{nuc}-Q_{GB}$ [kJ mol$^{-1}$]",ylabel=r"$Q_{surface}-Q_{GB}$ [kJ mol$^{-1}$]",title="Causally audited fast-firing materials")
    axs[1].bar(q.material_id,q.causal_fast_ratio,color="#0072B2",label="causally audited fast-firing ratio")
    axs[1].axhline(1.5,color="k",ls="--",lw=1,label="fast threshold")
    for i,r in q.reset_index(drop=True).iterrows():axs[1].text(i,r.causal_fast_ratio+.02,r.two_step_tier.replace("_"," "),ha="center",fontsize=8)
    axs[1].set(ylabel=r"Maximum $G_{ref}/G_{fast}$",title="Fast-firing and Chen evidence remain separate");axs[1].legend(frameon=False,fontsize=8)
    for ax in axs:ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Relative material-property window: observed ordering, not calibrated bounds")
    fig.tight_layout()
    for ext in ("pdf","png"):fig.savefig(FIG/f"relative_material_property_window_summary.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight")
    plt.close(fig)
    print("wrote relative material-property window figure")


if __name__=="__main__":main()
