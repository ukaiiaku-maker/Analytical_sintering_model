#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
OUT=Path("results/zro2_forward_defined_law_parameter_mapping");FIG=OUT/"figures"
def save(fig,name):fig.tight_layout();fig.savefig(FIG/f"{name}.png",dpi=180);fig.savefig(FIG/f"{name}.pdf");plt.close(fig)
def main():
 FIG.mkdir(parents=True,exist_ok=True);reg=pd.read_csv(OUT/"defined_law_registry.csv")
 fig,ax=plt.subplots(figsize=(12,6));ax.axis("off");ax.table(cellText=reg[["law_id","affected_process","implementation_status"]].values,colLabels=["law","process","status"],loc="center",fontsize=6);save(fig,"defined_law_architecture_schematic")
 x=pd.read_csv(OUT/"fixed_path_flux_integrals.csv");g=x.groupby("mode")[["Delta_rho_open","Delta_rho_closed","PR_topology_transfer"]].max();fig,ax=plt.subplots(figsize=(11,5));g.plot.bar(ax=ax);ax.set_ylabel("integrated channel");save(fig,"renewal_densification_channels")
 fig,ax=plt.subplots(figsize=(11,5));ax.bar(g.index,g.PR_topology_transfer);plt.setp(ax.get_xticklabels(),rotation=35,ha="right");ax.set_ylabel("conservative PR transfer integral");save(fig,"PR_topology_transfer_vs_density")
 x=pd.read_csv(OUT/"natural_state_T2_scan_by_mode.csv");fig,ax=plt.subplots(figsize=(11,5));
 for mode,q in x.groupby("mode"):ax.plot(q.T2_C,q.closed_inventory_formed,label=mode)
 ax.set(xlabel="T2 (°C)",ylabel="maximum closed inventory");ax.legend(fontsize=6);save(fig,"closed_inventory_accommodation_history")
 r=np.logspace(-9,-6);fv=.05;fig,ax=plt.subplots();ax.loglog(r*1e9,6*r/fv,label="R_Z=6r/f_v");ax.loglog(r*1e9,fv/r/1e6,label="pinning pressure proxy (MPa)");ax.set(xlabel="pore radius (nm)",ylabel="geometric pinning quantity");ax.legend();save(fig,"Zener_pinning_vs_pore_size")
 for file,name in (("candidate_state_T2_scan_by_mode.csv","candidate_state_T2_classification_by_mode"),("natural_state_T2_scan_by_mode.csv","natural_state_T2_classification_by_mode")):
  x=pd.read_csv(OUT/file);fig,ax=plt.subplots(figsize=(11,5));
  for mode,q in x.groupby("mode"):ax.plot(q.T2_C,q.final_rho,label=mode)
  ax.axhline(.976,color="k",ls="--");ax.set(xlabel="T2 (°C)",ylabel="final density");ax.legend(fontsize=6);save(fig,name)
 x=pd.read_csv(OUT/"mini_map_boundary_gap_summary.csv");fig,ax=plt.subplots(figsize=(9,4));ax.text(.5,.5,"No mode passed the fixed-path topology gate\nMini-map not run",ha="center",va="center");ax.axis("off");save(fig,"mini_map_window_by_mode")
 x=pd.read_csv(OUT/"candidate693168_defined_law_comparison.csv");fig,ax=plt.subplots(figsize=(11,5));ax.bar(x.model,x.closed_density_contribution);plt.setp(ax.get_xticklabels(),rotation=35,ha="right");ax.set_ylabel("closed density contribution");save(fig,"candidate693168_vs_forward_defined_laws")
 x=pd.read_csv(OUT/"defined_law_parameter_mapping.csv");z=pd.crosstab(x.parameter_class,x.requires_calibration);fig,ax=plt.subplots(figsize=(8,6));ax.imshow(z.values,aspect="auto",cmap="Blues");ax.set(yticks=np.arange(len(z)),yticklabels=z.index,xticks=np.arange(len(z.columns)),xticklabels=z.columns,xlabel="requires calibration");save(fig,"parameter_mapping_status_matrix")
if __name__=="__main__":main()
