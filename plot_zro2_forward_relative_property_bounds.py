#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
OUT=Path("results/zro2_forward_relative_property_bound_audit")
def save(fig,n):fig.suptitle("ZrO2 relative-property bound audit — diagnostic, not validation",fontsize=10);fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(OUT/f"{n}.png",dpi=180);fig.savefig(OUT/f"{n}.pdf");plt.close(fig)
def main():
 b=pd.read_csv(OUT/"effective_barrier_values_vs_T.csv");r=b[b.source.eq("representative_state")];fig,ax=plt.subplots(1,2,figsize=(10,4))
 for k,g in r.groupby("state_id"):ax[0].plot(g.T_C,g.Gstar_eV,label=k);ax[1].plot(g.T_C,g.sigma_eff_Pa/1e6,label=k)
 ax[0].set(xlabel="temperature (°C)",ylabel="G* (eV)");ax[1].set(xlabel="temperature (°C)",ylabel="effective stress (MPa)");ax[0].legend(fontsize=6);save(fig,"effective_Gstar_vs_T_and_sigma")
 q=pd.read_csv(OUT/"relative_barrier_group_audit.csv");q=q[q.source.eq("representative_state")];fig,ax=plt.subplots(figsize=(9,5));
 for k,g in q.groupby("state_id"):ax.plot(g.T_C,g.Q_nuc_eff_minus_Q_growth_kJ,label=k)
 ax.axhspan(-50,96,alpha=.2,color="green",label="promoted envelope");ax.set(xlabel="temperature (°C)",ylabel="Qnuc,eff − Qgrowth (kJ/mol)");ax.legend(fontsize=6);save(fig,"relative_barrier_groups_vs_T")
 c=q.dropna(subset=["Q_closed_eff_minus_Q_growth_kJ"]);fig,ax=plt.subplots();
 for k,g in c.groupby("state_id"):ax.plot(g.T_C,g.Q_closed_eff_minus_Q_growth_kJ,label=k)
 ax.axhspan(-252,-127,alpha=.2,color="green",label="promoted envelope");ax.set(xlabel="temperature (°C)",ylabel="Qclosed,eff − Qgrowth (kJ/mol)");ax.legend(fontsize=7);save(fig,"closed_growth_gap_vs_success_envelope")
 p=pd.read_csv(OUT/"relative_prefactor_group_audit.csv");p=p[p.source.eq("representative_state")];fig,ax=plt.subplots(figsize=(9,5));
 for k,g in p.groupby("state_id"):ax.plot(g.T_C,g.log10_kclosed_over_kgrowth,label=k)
 ax.axhspan(-1.52,1.50,alpha=.2,color="green");ax.set(xlabel="temperature (°C)",ylabel="log10(kclosed/kgrowth)");ax.legend(fontsize=7);save(fig,"prefactor_ratio_envelope")
 f=pd.read_csv(OUT/"statewise_property_window_classification.csv");cols=[c for c in f if c.startswith("inside_")];z=f.groupby("state_id")[cols].mean();fig,ax=plt.subplots(figsize=(11,5));ax.imshow(z,aspect="auto",cmap="RdYlGn",vmin=0,vmax=1);ax.set(yticks=range(len(z)),yticklabels=z.index,xticks=range(len(cols)),xticklabels=cols);plt.setp(ax.get_xticklabels(),rotation=40,ha="right");save(fig,"property_bound_pass_fail_matrix")
 s=pd.read_csv(OUT/"property_outside_bounds_summary.csv");fig,ax=plt.subplots(figsize=(10,5));ax.barh(s.property_group,s.outside_fraction,color=np.where(s.primary_unknown,"#d95f02","#1b9e77"));ax.set_xlabel("outside fraction");save(fig,"statewise_outside_bound_summary")
if __name__=="__main__":main()
