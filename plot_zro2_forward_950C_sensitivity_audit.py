#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT=Path("results/zro2_forward_950C_sensitivity_chen_failure_audit");FIG=OUT/"figures"
def save(fig,name,evidence):
    fig.suptitle(f"950 °C conditioned start — {evidence} — no validation claim",fontsize=10);fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(FIG/f"{name}.png",dpi=180);fig.savefig(FIG/f"{name}.pdf");plt.close(fig)
def main():
    FIG.mkdir(parents=True,exist_ok=True);base=pd.read_csv(OUT/"baseline_pdf_conditioned_summary.csv");s=pd.read_csv(OUT/"common_state_fast_rate_summary.csv");q=s.query("state_class=='primary_common_state'");m=pd.read_csv(OUT/"common_state_matched_density_curves.csv");d=pd.read_csv(OUT/"common_state_factorial_design.csv")
    fig,ax=plt.subplots(1,2,figsize=(9,4));ax[0].bar(["full 50","conditioned 50"],[.887561,base.rate50_final_rho.iloc[0]]);ax[1].bar(["full 50","conditioned 50"],[.559961,base.rate50_final_G_um.iloc[0]]);ax[0].set_ylabel("final density");ax[1].set_ylabel("final G (µm)");save(fig,"baseline_conditioned_summary","baseline separation")
    fig,ax=plt.subplots(2,3,figsize=(12,7));
    for i,(a,col) in enumerate(zip(ax.ravel(),["pore_D50_nm","rho_start","G_start_nm"]*2)):
        metric="final_rho_50" if i<3 else "final_G_um_50";g=q.groupby(col)[metric].agg(['mean','min','max']);a.plot(g.index,g['mean'],'o-');a.fill_between(g.index,g['min'],g['max'],alpha=.2);a.set(xlabel=col,ylabel=metric)
    save(fig,"common_state_sensitivity_density_grain","controlled sensitivity")
    ids=d.query("rho_start==.66 and G_start_nm==50 and pore_log_width==.65 and phi_iso_fraction==0 and phi_closed_fraction==0")
    fig,ax=plt.subplots();
    for _,r in ids.iterrows():g=m[m.case_id==r.case_id];ax.plot(g.rho,g.G_5_over_G_50,label=f"D50={r.pore_D50_nm:g} nm")
    ax.axhline(1,c='k',ls='--');ax.set(xlabel="density",ylabel="G(5)/G(50)");ax.legend();save(fig,"common_state_fast_rate_G_vs_rho","controlled sensitivity")
    fig,ax=plt.subplots(1,2,figsize=(9,4));
    for _,r in ids.iterrows():g=m[m.case_id==r.case_id];ax[0].plot(g.rho,g.pore_D90_m_50*1e9,label=f"{r.pore_D50_nm:g} nm");ax[1].plot(g.rho,g.tau_remove_D90_s_50)
    ax[0].set(xlabel="density",ylabel="50 C/min D90 (nm)");ax[1].set(xlabel="density",ylabel="tau D90 (s)",yscale='log');ax[0].legend();save(fig,"common_state_pore_D90_tau_remove","controlled sensitivity")
    fig,ax=plt.subplots(1,3,figsize=(12,4));
    for _,r in ids.iterrows():g=m[m.case_id==r.case_id];ax[0].plot(g.rho,g.R_Z_eff_m_50*1e6);ax[1].plot(g.rho,g.S_Z_50);ax[2].plot(g.rho,g.Gamma_growth_50)
    for a,y in zip(ax,["R_Z (µm)","S_Z","Gamma growth"]):a.set(xlabel="density",ylabel=y)
    save(fig,"common_state_zener_pinning","controlled sensitivity")
    x=pd.read_csv(OUT/"chen_failure_decomposition_full.csv");cats=pd.Categorical(x.failure_mechanism);fig,ax=plt.subplots();sc=ax.scatter(x.T1_C,x.T2_C,c=cats.codes,cmap='tab10');ax.set(xlabel="T1 (°C)",ylabel="T2 (°C)");save(fig,"chen_failure_map_strict","strict 0.976 / 0.29 µm")
    b=pd.read_csv(OUT/"chen_boundary_ordering_table.csv");fig,ax=plt.subplots();sc=ax.scatter(b.T1_C,b.switch_density,c=b.gap_C,cmap='coolwarm',vmin=-200,vmax=200,s=35+2*b.hold_h);fig.colorbar(sc,ax=ax,label="upper − lower boundary (°C)");ax.set(xlabel="T1 (°C)",ylabel="switch density");save(fig,"chen_boundary_gap_map","strict boundary audit")
    r=pd.read_csv(OUT/"chen_target_relaxation_diagnostics.csv");fig,ax=plt.subplots(1,2,figsize=(9,4));
    for density,g in r.groupby('density_target'):ax[0].plot(g.grain_threshold_um,g.success_count,'o-',label=f"rho {density}");ax[1].plot(g.grain_threshold_um,g.finite_window_count,'o-')
    ax[0].set(xlabel="grain threshold (µm)",ylabel="success count");ax[1].set(xlabel="grain threshold (µm)",ylabel="finite windows");ax[0].legend();save(fig,"chen_relaxed_threshold_maps","relaxed diagnostic")
    # Required single-metric aliases.
    fig,ax=plt.subplots();
    for density,g in r.groupby('density_target'):ax.plot(g.grain_threshold_um,g.success_count,'o-',label=f"rho {density}")
    ax.set(xlabel="grain threshold (µm)",ylabel="success count");ax.legend();save(fig,"chen_success_count_vs_grain_threshold","relaxed diagnostic")
    maps=pd.read_csv(OUT/"chen_relaxed_threshold_boundary_maps.csv");fig,ax=plt.subplots();g=maps.groupby('grain_threshold_um').gap_C.mean();ax.plot(g.index,g.values,'o-');ax.axhline(0,c='k',ls='--');ax.set(xlabel="grain threshold (µm)",ylabel="mean boundary gap (°C)");save(fig,"chen_boundary_gap_vs_grain_threshold","relaxed diagnostic")
    h=pd.read_csv(OUT/"representative_chen_path_histories.csv");fig,ax=plt.subplots(1,2,figsize=(9,4));
    for k,g in h.groupby('representative_index'):ax[0].plot(g.t_s/3600,g.rho,label=str(k));ax[1].plot(g.t_s/3600,g.G_um,label=str(k))
    ax[0].set(xlabel="time (h)",ylabel="density");ax[1].set(xlabel="time (h)",ylabel="G (µm)");ax[0].legend();save(fig,"representative_chen_histories","strict path histories")
    o=pd.read_csv(OUT/"chen_failure_OAT_summary.csv");fig,ax=plt.subplots();effect=o.groupby('modified_parameter').boundary_gap_C.mean().fillna(-250).sort_values();ax.barh(effect.index,effect.values);ax.set_xlabel("mean boundary gap (°C; censored=-250)");save(fig,"oat_sensitivity_tornado","bounded OAT sensitivity")
    tri=pd.read_csv(OUT/"candidate_triage_950C_sensitivity.csv");cols=['fast_density_attained','fast_smaller_grain_sign','fast_smaller_D90_sign','strict_Chen_window','pathway_consistent'];fig,ax=plt.subplots(figsize=(8,8));ax.imshow(tri[cols].astype(float),aspect='auto',cmap='RdYlGn',vmin=0,vmax=1);ax.set(xticks=range(len(cols)),xticklabels=cols,yticks=[]);plt.setp(ax.get_xticklabels(),rotation=35,ha='right');save(fig,"pathway_consistency_matrix","bounded OAT sensitivity")
if __name__=="__main__":main()
