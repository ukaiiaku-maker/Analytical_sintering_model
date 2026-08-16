#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT=Path("results/zro2_forward_closed_channel_physical_law_comparison");FIG=OUT/"figures"
def save(fig,name):
    fig.tight_layout();fig.savefig(FIG/f"{name}.png",dpi=180);fig.savefig(FIG/f"{name}.pdf");plt.close(fig)
def main():
    FIG.mkdir(parents=True,exist_ok=True)
    reg=pd.read_csv(OUT/"closed_law_registry.csv");fig,ax=plt.subplots(figsize=(11,4));ax.axis("off");ax.table(cellText=reg[["closed_channel_law","direct_density_change","empirical"]].values,colLabels=["law","density removal","empirical"],loc="center");save(fig,"closed_law_schematic")
    f=pd.read_csv(OUT/"fixed_path_closed_law_summary.csv");g=f.groupby(["law","path"]).final_rho.median().unstack();fig,ax=plt.subplots(figsize=(10,5));g.plot.bar(ax=ax);ax.set_ylabel("final density");ax.axhline(.976,color="k",ls="--");save(fig,"fixed_path_density_fluxes_by_law")
    for file,name in (("candidate_state_closed_law_T2_scan.csv","candidate_state_T2_classification_by_law"),("natural_state_closed_law_T2_scan.csv","natural_state_T2_classification_by_law")):
        x=pd.read_csv(OUT/file);x=x[x.prefactor_factor.eq(1)];fig,ax=plt.subplots(figsize=(10,5));
        for law,q in x.groupby("law"):ax.plot(q.T2_C,q.final_rho,label=law)
        ax.axhline(.976,color="k",ls="--");ax.set(xlabel="T2 (°C)",ylabel="final density");ax.legend(fontsize=6);save(fig,name)
    x=pd.read_csv(OUT/"closed_law_lower_boundary_preservation.csv");fig,ax=plt.subplots(figsize=(11,5));z=pd.crosstab(x.law,x.preservation_class);z.plot.bar(stacked=True,ax=ax);ax.set_ylabel("configuration count");save(fig,"closed_law_lower_boundary_preservation")
    x=pd.read_csv(OUT/"closed_law_boundary_gap_summary.csv");fig,ax=plt.subplots(figsize=(10,4));ax.barh(x.config_id,x.finite_window_count);ax.set_xlabel("finite strict window count");save(fig,"closed_law_boundary_gap_summary")
    x=pd.read_csv(OUT/"candidate693168_vs_physical_closed_laws.csv");x=x[(x.model.eq("candidate_693168"))|x.model.str.endswith("__x1")];fig,ax=plt.subplots(figsize=(11,5));ax.bar(x.model,x.closed_density_contribution);ax.set_ylabel("closed density contribution");plt.setp(ax.get_xticklabels(),rotation=35,ha="right");save(fig,"candidate693168_vs_closed_laws")
    x=pd.read_csv(OUT/"candidate_state_closed_law_T2_scan.csv");x=x[x.law.isin(["empirical_closed_rate_scale","GB_diffusion_closed_shrinkage","renewal_limited_closed_shrinkage"])];fig,ax=plt.subplots(figsize=(10,5));q=x.groupby(["law","prefactor_factor"]).Delta_rho_closed.max().reset_index();
    for law,g in q.groupby("law"):ax.semilogx(g.prefactor_factor,g.Delta_rho_closed,"o-",label=law)
    ax.set(xlabel="prefactor uncertainty / empirical scale",ylabel="max inventory-bounded closed density contribution");ax.legend(fontsize=7);save(fig,"empirical_vs_physical_law_comparison")
if __name__=="__main__":main()
