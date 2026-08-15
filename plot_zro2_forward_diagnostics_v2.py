#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from zro2_forward.diagnostics_v2 import OUT,MODES,build
from zro2_forward.densification import kinetic_state

def save(fig,name):
    fig.tight_layout(); fig.savefig(OUT/f"{name}.png",dpi=180); fig.savefig(OUT/f"{name}.pdf"); plt.close(fig)
def main():
    temps=np.linspace(900,1600,180); fig,ax=plt.subplots()
    for mode in MODES:
        b=build(mode).barrier
        for sigma in [1e7,2.5e8]: ax.plot(temps,[b.Gstar(sigma,(T+273.15))/1.602176634e-19 for T in temps],label=f"{mode}, {sigma/1e6:g} MPa")
    ax.axvline(1557,color='k',ls='--'); ax.set(xlabel="T (°C)",ylabel="G* (eV)"); ax.legend(fontsize=6,ncol=2); save(fig,"Gstar_vs_T_for_relevant_sigma")
    fig,ax=plt.subplots()
    for mode in MODES:
        m=build(mode); ax.plot(temps,[kinetic_state(2.5e8,T+273.15,100e-9,.8,m.barrier,m.material,1,39.5)["activity"] for T in temps],label=mode)
    ax.axvline(1557,color='k',ls='--'); ax.set(xlabel="T (°C)",ylabel="activity"); ax.legend(); save(fig,"activity_vs_T_barrier_modes")
    fig,ax=plt.subplots(); target=pd.read_csv("data/targets/mazaheri_8ysz_2008/density_vs_temperature_digitized.csv").query("method=='CS'")
    ax.errorbar(target.T_C,target.fractional_density,yerr=target.uncertainty_density,fmt='ko',label="PDF")
    for mode in MODES:
        h=pd.read_csv(OUT/f"history_CS_{mode}.csv"); ax.plot(h.T_C,h.rho,label=mode)
    ax.set(xlabel="T (°C)",ylabel="density"); ax.legend(fontsize=7); save(fig,"CS_density_curve_by_barrier_mode")
    counts=pd.read_csv(OUT/"barrier_mode_chen_counts.csv"); fig,ax=plt.subplots(); counts.set_index("barrier_mode")[["success_count","density_ok_count","growth_ok_count"]].plot.bar(ax=ax); ax.set_ylabel("coarse-map count"); save(fig,"Chen_count_by_barrier_mode")
    fig,ax=plt.subplots()
    for label in ["endpoint_only","density_curve_only","density_plus_grain_endpoint"]:
        h=pd.read_csv(OUT/f"history_calibration_{label}.csv"); ax.plot(h.T_C,h.rho,label=label)
    ax.errorbar(target.T_C,target.fractional_density,yerr=target.uncertainty_density,fmt='ko',label="PDF"); ax.set(xlabel="T (°C)",ylabel="density"); ax.legend(); save(fig,"CS_density_curve_calibration_modes")
    grain=pd.read_csv(OUT/"CS_grain_endpoint_residuals.csv"); fig,ax=plt.subplots(); ax.bar(grain.objective_mode,grain.model_G_um); ax.axhline(2.14,color='k',ls='--'); ax.tick_params(axis='x',rotation=20); ax.set_ylabel("final G (µm)"); save(fig,"CS_grain_endpoint_calibration_modes")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for label in ["endpoint_only","density_curve_only","density_plus_grain_endpoint"]:
        h=pd.read_csv(OUT/f"history_calibration_{label}.csv"); ax[0].plot(h.T_C,h.activity,label=label); ax[1].plot(h.T_C,h.sigma_eff_Pa/1e6,label=label)
    ax[0].set(xlabel="T (°C)",ylabel="activity"); ax[1].set(xlabel="T (°C)",ylabel="effective stress (MPa)"); ax[0].legend(); save(fig,"calibrated_activity_and_sigma_CS")
    chen=pd.read_csv(OUT/"chen_failure_decomposition.csv"); classes={x:i for i,x in enumerate(sorted(chen.failure_class.unique()))}; colors=chen.failure_class.map(classes)
    for name,x,y in [("chen_map_failure_decomposition_T1_T2","T1_C","T2_C"),("chen_map_failure_decomposition_G1_T2","first_step_G1_um","T2_C")]:
        fig,ax=plt.subplots(); sc=ax.scatter(chen[x],chen[y],c=colors,cmap="tab10",alpha=.75); ax.set(xlabel=x,ylabel=y)
        handles=[plt.Line2D([],[],marker='o',ls='',color=plt.cm.tab10(i/max(len(classes)-1,1)),label=k) for k,i in classes.items()]; ax.legend(handles=handles,fontsize=6); save(fig,name)
    fig,ax=plt.subplots();
    for _,g in chen.groupby("switch_density"): ax.scatter(g.T2_C,g.final_rho,c=g.final_G_um,cmap="viridis",label=f"switch {g.switch_density.iloc[0]}")
    ax.axhspan(.966,.986,color='green',alpha=.1); ax.set(xlabel="T2 (°C)",ylabel="final density"); ax.legend(); save(fig,"lower_upper_boundary_overlay")
    reps=pd.read_csv(OUT/"representative_chen_path_histories.csv"); fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in reps.groupby("path_label"): ax[0].plot(g.t_s/3600,g.rho,label=k); ax[1].plot(g.t_s/3600,g.G_um,label=k)
    ax[0].set(xlabel="time (h)",ylabel="density"); ax[1].set(xlabel="time (h)",ylabel="G (µm)"); ax[0].legend(); save(fig,"representative_lowT_midT_highT_two_step_histories")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in reps.groupby("path_label"): ax[0].plot(g.t_s/3600,g.pore_D90_m*1e9,label=k); ax[1].plot(g.t_s/3600,[max([v for v in __import__('json').loads(x) if np.isfinite(v)] or [np.nan]) for x in g.tau_remove_s_json],label=k)
    ax[0].set(xlabel="time (h)",ylabel="D90 (nm)"); ax[1].set(xlabel="time (h)",ylabel="max removal time (s)",yscale="log"); ax[0].legend(); save(fig,"pore_D90_and_tau_remove_along_chen_paths")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in reps.groupby("path_label"): ax[0].plot(g.t_s/3600,g.S_Z,label=k); ax[1].plot(g.t_s/3600,g.Gamma_growth,label=k)
    ax[0].set(xlabel="time (h)",ylabel="Zener factor"); ax[1].set(xlabel="time (h)",ylabel="growth factor"); ax[0].legend(); save(fig,"zener_and_growth_mobility_along_chen_paths")
    print(f"wrote diagnostics to {OUT}")
if __name__=="__main__": main()
