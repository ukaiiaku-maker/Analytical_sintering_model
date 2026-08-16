#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from zro2_forward.conditioned_950c import OUT

TARGET=Path("data/targets/mazaheri_8ysz_2008")
def save(fig,name):fig.tight_layout();fig.savefig(OUT/f"{name}.png",dpi=180);fig.savefig(OUT/f"{name}.pdf");plt.close(fig)
def main():
    fig,ax=plt.subplots();ax.plot([25,950],[.484,.66],'-o',label="full-process prediction domain");ax.plot([950,1500],[.66,.975],'-o',label="PDF-conditioned comparison domain");ax.set(xlabel="T (°C)",ylabel="reference density");ax.legend();save(fig,"full_vs_conditioned_schematic")
    full=pd.read_csv(OUT/"full_process_dense_histories.csv");fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in full.groupby("case"):q=g[g.T_C<=950];ax[0].plot(q.T_C,q.rho,label=k);ax[1].plot(q.T_C,q.G_um*1000,label=k)
    ax[0].plot(950,.66,'k*',ms=12,label="common target state");ax[1].plot(950,50,'k*',ms=12);ax[0].set(xlabel="T (°C)",ylabel="density");ax[1].set(xlabel="T (°C)",ylabel="G (nm)");ax[0].legend();save(fig,"full_process_prediction_to_950C")
    hist=pd.read_csv(OUT/"pdf_conditioned_CS_calibration_histories.csv");td=pd.read_csv(TARGET/"density_vs_temperature_digitized.csv").query("method=='CS'");tg=pd.read_csv(TARGET/"grain_size_vs_temperature_digitized.csv").query("method=='CS'")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in hist.groupby("case"):ax[0].plot(g.T_C,g.rho,label=k);ax[1].plot(g.T_C,g.G_um,label=k)
    ax[0].scatter(td.T_C,td.fractional_density,c='k');ax[1].scatter(tg.T_C,tg.G_um,c='k');ax[0].set(xlabel="T (°C)",ylabel="density");ax[1].set(xlabel="T (°C)",ylabel="G (µm)");ax[0].legend(fontsize=6);save(fig,"pdf_conditioned_CS_overlay")
    # Required aliases with single-observable emphasis.
    fig,ax=plt.subplots();
    for k,g in hist.groupby("case"):ax.plot(g.T_C,g.rho,label=k)
    ax.scatter(td.T_C,td.fractional_density,c='k');ax.set(xlabel="T (°C)",ylabel="density");ax.legend(fontsize=7);save(fig,"pdf_conditioned_CS_density_curve")
    fig,ax=plt.subplots();
    for k,g in hist.groupby("case"):ax.plot(g.T_C,g.G_um,label=k)
    ax.scatter(tg.T_C,tg.G_um,c='k');ax.set(xlabel="T (°C)",ylabel="G (µm)");ax.legend(fontsize=7);save(fig,"pdf_conditioned_CS_grain_curve")
    fig,ax=plt.subplots();
    for k,g in hist.groupby("case"):ax.plot(g.rho,g.G_um,label=k)
    ax.set(xlabel="density",ylabel="G (µm)");ax.legend(fontsize=7);save(fig,"pdf_conditioned_CS_G_vs_rho")
    h=pd.read_csv(OUT/"pdf_conditioned_dense_histories.csv");match=pd.read_csv(OUT/"pdf_conditioned_matched_density_curves.csv");fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in h.groupby("case"):ax[0].plot(g.T_C,g.rho,label=k);ax[1].plot(g.rho,g.G_um,label=k)
    ax[0].set(xlabel="T (°C)",ylabel="density");ax[1].set(xlabel="density",ylabel="G (µm)");ax[0].legend();save(fig,"pdf_conditioned_fast_rate_overlay")
    fig,ax=plt.subplots();
    for k,g in h.groupby("case"):ax.plot(g.T_C,g.rho,label=k)
    ax.set(xlabel="T (°C)",ylabel="density");ax.legend();save(fig,"pdf_conditioned_fast_rate_density_vs_temperature")
    fig,ax=plt.subplots();ax.plot(match.rho,match.G_um_5,label="5 C/min");ax.plot(match.rho,match.G_um_50,label="50 C/min");ax.set(xlabel="density",ylabel="G (µm)");ax.legend();save(fig,"pdf_conditioned_fast_rate_G_vs_rho")
    fig,ax=plt.subplots();ax.plot(match.rho,match.pore_D90_m_5*1e9,label="5 C/min");ax.plot(match.rho,match.pore_D90_m_50*1e9,label="50 C/min");ax.set(xlabel="density",ylabel="D90 (nm)");ax.legend();save(fig,"pdf_conditioned_fast_rate_pore_D90_vs_rho")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in h.groupby("case"):ax[0].plot(g.T_C,g.activity,label=k);ax[1].plot(g.T_C,g.sigma_eff_Pa/1e6,label=k)
    ax[0].set(xlabel="T (°C)",ylabel="activity");ax[1].set(xlabel="T (°C)",ylabel="stress (MPa)");ax[0].legend();save(fig,"pdf_conditioned_fast_rate_activity_sigma")
    fig,ax=plt.subplots(2,2,figsize=(10,8));ax=ax.ravel();
    for tag in ['5','50']:
        ax[0].plot(match.rho,match[f'pore_D90_m_{tag}']*1e9,label=f"{tag} C/min");ax[1].plot(match.rho,match[f'fine_pore_fraction_{tag}']);ax[2].plot(match.rho,match[f'tau_remove_D90_s_{tag}']);ax[3].plot(match.rho,match[f'R_Z_eff_m_{tag}']*1e6)
    ax[0].set_ylabel("D90 (nm)");ax[1].set_ylabel("fine fraction");ax[2].set(ylabel="D90 removal time (s)",yscale="log");ax[3].set_ylabel("R_Z (µm)");ax[0].legend();save(fig,"pdf_conditioned_pore_evolution")
    c=pd.read_csv(OUT/"pdf_conditioned_chen_classification_points.csv");
    def scatter(name,x,y):
        fig,ax=plt.subplots();sc=ax.scatter(c[x],c[y],c=c.final_rho,s=25+70*c.chen_success,cmap='viridis');fig.colorbar(sc,ax=ax,label="final density");ax.set(xlabel=x,ylabel=y);save(fig,name)
    scatter("pdf_conditioned_chen_map_T1_T2","T1_C","T2_C");scatter("pdf_conditioned_chen_map_G1_T2","G1_um","T2_C");scatter("pdf_conditioned_chen_map_switch_density_T2","switch_density","T2_C")
    fig,ax=plt.subplots(1,3,figsize=(13,4));
    for a,x,y in zip(ax,["T1_C","G1_um","switch_density"],["T2_C"]*3):sc=a.scatter(c[x],c[y],c=c.final_G_um,cmap='plasma');a.set(xlabel=x,ylabel=y)
    fig.colorbar(sc,ax=ax,label="final G (µm)");save(fig,"pdf_conditioned_chen_maps")
    best=pd.read_csv(OUT/"pdf_conditioned_best_TSS_like_paths.csv");fig,ax=plt.subplots();ax.scatter(c.final_rho,c.final_G_um,c=c.T2_C,cmap='viridis',label="map");ax.scatter([.976],[.29],marker='*',s=180,c='red',label="TSS target");ax.set(xlabel="final density",ylabel="final G (µm)");ax.legend();save(fig,"pdf_conditioned_best_two_step_vs_highT_G_rho")
    fig,ax=plt.subplots(1,2,figsize=(10,4));sc=ax[0].scatter(c.final_rho,c.final_D90_nm,c=c.T2_C);ax[1].scatter(c.final_rho,c.final_fine_fraction,c=c.T2_C);ax[0].set(xlabel="density",ylabel="D90 (nm)");ax[1].set(xlabel="density",ylabel="fine fraction");save(fig,"pdf_conditioned_two_step_pore_D90_fine_fraction")
    old=pd.read_csv("results/zro2_forward_natural_pore_evolution_target_8ysz/chen_classification_points.csv");fig,ax=plt.subplots(1,2,figsize=(10,4));ax[0].scatter(old.final_rho,old.final_G_um,c=old.T2_C);ax[1].scatter(c.final_rho,c.final_G_um,c=c.T2_C);ax[0].set_title("full process");ax[1].set_title("conditioned");
    for a in ax:a.set(xlabel="density",ylabel="G (µm)")
    save(fig,"failure_decomposition_full_vs_conditioned")
    print(f"wrote conditioned figures to {OUT}")
if __name__=="__main__":main()
