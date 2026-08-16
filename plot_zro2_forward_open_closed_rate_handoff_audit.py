#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
OUT=Path("results/zro2_forward_open_closed_rate_handoff_audit")
def save(fig,n):fig.suptitle("ZrO2 open/closed handoff audit — conditional, not validation",fontsize=10);fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(OUT/f"{n}.png",dpi=180);fig.savefig(OUT/f"{n}.pdf");plt.close(fig)
def main():
 d=pd.read_csv(OUT/"handoff_mode_design.csv");fig,ax=plt.subplots(figsize=(12,3));ax.axis("off");xs=np.linspace(.05,.95,5);labs=["open shrinkage","conservative PR/transfer","closed availability","bounded closed shrinkage","density = open + closed"]
 for x,l in zip(xs,labs):ax.text(x,.5,l,ha="center",bbox=dict(boxstyle="round",fc="#e8f1fa"))
 for a,b in zip(xs[:-1],xs[1:]):ax.annotate("",(b-.07,.5),(a+.07,.5),arrowprops=dict(arrowstyle="->"))
 save(fig,"open_closed_handoff_schematic")
 f=pd.read_csv(OUT/"fast_rate_handoff_summary.csv").query("rate_C_min==50");fig,ax=plt.subplots(figsize=(11,5));ax.bar(f.mode_id,f.final_rho);ax.axhline(.959216,ls="--",c="k",label="previous PDF-conditioned");ax.set_ylabel("50 °C/min final density");plt.setp(ax.get_xticklabels(),rotation=45,ha="right");ax.legend();save(fig,"fast_rate_density_recovery_by_handoff")
 x=pd.read_csv(OUT/"fast_rate_flux_integrals.csv").query("rate_C_min==50");fig,ax=plt.subplots(figsize=(11,5));q=np.arange(len(x));ax.bar(q-.2,x.Delta_rho_open,.4,label="open");ax.bar(q+.2,x.Delta_rho_closed,.4,label="closed");ax.set(xticks=q,xticklabels=x.mode_id,ylabel="integrated density flux");plt.setp(ax.get_xticklabels(),rotation=45,ha="right");ax.legend();save(fig,"open_vs_closed_shrinkage_integrals")
 h=pd.read_csv(OUT/"handoff_diagnostic_histories.csv");sel=h[(h.path.eq("50C"))&h.mode_id.isin(["resolved_default","diagnostic_open_recovery","balanced_handoff","closed_rate_100x"])]
 fig,ax=plt.subplots();
 for k,g in sel.groupby("mode_id"):ax.plot(g.T_C,g.handoff_readiness,label=k)
 ax.set(xlabel="temperature (°C)",ylabel="handoff readiness");ax.legend(fontsize=7);save(fig,"handoff_readiness_vs_temperature")
 sel=h[(h.path.eq("near_best"))&h.mode_id.isin(["resolved_default","balanced_handoff","closed_rate_100x","candidate_state_injection_diagnostic"])]
 fig,ax=plt.subplots(1,3,figsize=(12,4));
 for k,g in sel.groupby("mode_id"):
  ax[0].plot(g.rho,g.phi_open,label=k);ax[1].plot(g.rho,g.phi_precursor);ax[2].plot(g.rho,g.phi_closed)
 for a,l in zip(ax,["open","precursor","closed"]):a.set(xlabel="density",ylabel=f"{l} fraction")
 ax[0].legend(fontsize=6);save(fig,"pore_store_evolution_by_handoff")
 c=pd.read_csv(OUT/"chen_handoff_classification_points.csv");fig,ax=plt.subplots(figsize=(10,6));codes=pd.Categorical(c.classification).codes;ax.scatter(c.T2_C,c.mode_id,c=codes,cmap="tab10",s=8);ax.set_xlabel("T2 (°C)");save(fig,"Chen_map_by_handoff_mode")
 b=pd.read_csv(OUT/"chen_handoff_boundary_gap_summary.csv");fig,ax=plt.subplots(figsize=(11,5));ax.bar(b.mode_id,b.max_gap_C);ax.axhline(0,c="k",ls="--");ax.set_ylabel("best upper − lower boundary (°C)");plt.setp(ax.get_xticklabels(),rotation=45,ha="right");save(fig,"boundary_gap_by_handoff_mode")
 i=pd.read_csv(OUT/"candidate_state_injection_diagnostic.csv");fig,ax=plt.subplots(1,2,figsize=(9,4));q=i.dropna(subset=["rate_C_min"]);ax[0].bar(q.rate_C_min.astype(str),q.final_rho);ax[0].set(xlabel="rate (°C/min)",ylabel="final density");t=i.dropna(subset=["T2_C"]);ax[1].bar(t.path,t.final_rho);ax[1].axhline(.976,c="k",ls="--");ax[1].set(ylabel="two-step final density");save(fig,"candidate_state_injection_result")
 g=pd.read_csv(OUT/"pathway_gate_summary.csv");cols=["fast_density_substantially_recovered","smaller_grain_sign_preserved","closed_shrinkage_non_negligible","finite_bracketed_window","all_gates"];fig,ax=plt.subplots(figsize=(10,6));im=ax.imshow(g[cols].astype(int),aspect="auto",cmap="RdYlGn",vmin=0,vmax=1);ax.set(yticks=range(len(g)),yticklabels=g.mode_id,xticks=range(len(cols)),xticklabels=cols);plt.setp(ax.get_xticklabels(),rotation=35,ha="right");save(fig,"pathway_gate_matrix")
if __name__=="__main__":main()
