#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT=Path("results/zro2_forward_resolved_rule_rate_balance_audit")
def save(fig,name):
    fig.suptitle("Resolved-rule rate-balance audit — diagnostic only, not validation",fontsize=10)
    fig.tight_layout(rect=[0,0,1,.96]); fig.savefig(OUT/f"{name}.png",dpi=180); fig.savefig(OUT/f"{name}.pdf"); plt.close(fig)
def main():
    d=pd.read_csv(OUT/"density_flux_comparison.csv")
    f=d[d.audit_case.isin(["previous_PDF_conditioned_50C","resolved_50C"])]
    fig,ax=plt.subplots();
    for k,g in f.groupby("audit_case"): ax.plot(g.T_C,g.rho,label=k)
    ax.set(xlabel="temperature (°C)",ylabel="density");ax.legend(fontsize=7);save(fig,"previous_vs_resolved_density_trajectory")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in f.groupby("audit_case"):
        ax[0].plot(g.T_C,g.rho_dot_open_sinv,label=k);ax[1].plot(g.T_C,g.rho_dot_closed_sinv,label=k)
    ax[0].set(ylabel="open density flux (s⁻¹)");ax[1].set(ylabel="closed density flux (s⁻¹)");[a.set_xlabel("temperature (°C)") for a in ax];ax[0].legend(fontsize=7);save(fig,"open_vs_closed_density_fluxes")
    loss=pd.read_csv(OUT/"fast_rate_density_loss_decomposition.csv").iloc[0]
    fig,ax=plt.subplots();names=["open-flux change","closed attained","net density loss"];vals=[loss.open_shrinkage_change,loss.closed_shrinkage_attained,-loss.density_loss];ax.bar(names,vals,color=["#d95f02","#1b9e77","#7570b3"]);ax.axhline(0,color="k",lw=.8);ax.set_ylabel("density contribution");save(fig,"fast_density_loss_decomposition")
    p=pd.read_csv(OUT/"pore_store_comparison.csv");p=p[p.audit_case.isin(["resolved_5C","resolved_50C"])]
    fig,ax=plt.subplots(1,3,figsize=(12,4))
    for k,g in p.groupby("audit_case"):
        ax[0].plot(g.T_C,g.phi_open,label=k);ax[1].plot(g.T_C,g.phi_precursor,label=k);ax[2].plot(g.T_C,g.phi_closed,label=k)
    for a,y in zip(ax,["open","precursor","closed"]):a.set(xlabel="temperature (°C)",ylabel=f"{y} pore fraction")
    ax[0].legend(fontsize=7);save(fig,"pore_store_evolution_5C_50C")
    a=pd.read_csv(OUT/"accommodation_comparison.csv");a=a[a.audit_case.isin(["resolved_5C","resolved_50C"])]
    fig,ax=plt.subplots();
    for k,g in a.groupby("audit_case"):ax.plot(g.T_C,g.closed_accommodation_used,label=k)
    ax.axhline(a.closed_accommodation_capacity.max(),ls="--",c="k",label="capacity");ax.set(xlabel="temperature (°C)",ylabel="accommodation used");ax.legend(fontsize=7);save(fig,"accommodation_usage_vs_temperature")
    lo=pd.read_csv(OUT/"lower_boundary_rate_balance.csv");fig,ax=plt.subplots();x=np.arange(len(lo));ax.bar(x-.2,lo.Delta_rho_open,.4,label="open");ax.bar(x+.2,lo.Delta_rho_closed,.4,label="closed");ax.set(xticks=x,xticklabels=lo.path,ylabel="integrated density contribution");ax.legend();save(fig,"lower_boundary_rate_balance")
    up=pd.read_csv(OUT/"upper_boundary_growth_balance.csv");fig,ax=plt.subplots();ax.bar(up.path,up.G_dot_intrinsic_integrated_m*1e6,label="intrinsic");ax.bar(up.path,up.G_dot_actual_integrated_m*1e6,label="actual");ax.set_ylabel("integrated growth (µm)");ax.legend();save(fig,"upper_boundary_growth_balance")
    c=pd.read_csv(OUT/"resolved_vs_candidate693168_state_comparison.csv");fig,ax=plt.subplots();metrics=["closed_fraction_at_switch","A_closed_at_switch","PR_memory_at_switch","first_step_growth"];x=np.arange(len(metrics));w=.35
    for i,(_,r) in enumerate(c.iterrows()):ax.bar(x+(i-.5)*w,[r[m] for m in metrics],w,label=r.model)
    ax.set(xticks=x,xticklabels=metrics);plt.setp(ax.get_xticklabels(),rotation=25,ha="right");ax.legend(fontsize=7);save(fig,"resolved_vs_candidate693168_state_comparison")
if __name__=="__main__":main()
