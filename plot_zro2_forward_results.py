#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

OUT=Path("results/zro2_forward_natural_pore_evolution_target_8ysz")
DATA=Path("data/targets/mazaheri_8ysz_2008")
def save(fig,name):
    fig.tight_layout(); fig.savefig(OUT/f"{name}.png",dpi=180); fig.savefig(OUT/f"{name}.pdf"); plt.close(fig)
def main():
    h=pd.read_csv(OUT/"dense_histories.csv"); matched=pd.read_csv(OUT/"matched_density_curves.csv")
    target=pd.read_csv(DATA/"density_vs_temperature_digitized.csv"); gt=pd.read_csv(DATA/"grain_size_vs_temperature_digitized.csv")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for k,g in target.groupby("method"): ax[0].plot(g.T_C,g.fractional_density,'o-',label=k)
    for k,g in gt.groupby("method"): ax[1].plot(g.T_C,g.G_um,'o-',label=k)
    ax[0].set(xlabel="T (°C)",ylabel="fractional density"); ax[1].set(xlabel="T (°C)",ylabel="grain size (µm)"); ax[0].legend(); save(fig,"target_data_summary")
    cs=h[h.case.eq("CS_thermal")]; fig,ax=plt.subplots(1,2,figsize=(10,4)); ax[0].plot(cs.T_C,cs.rho,label="model")
    q=target[target.method.eq("CS")]; ax[0].errorbar(q.T_C,q.fractional_density,yerr=q.uncertainty_density,fmt='o',label="PDF")
    ax[1].plot(cs.T_C,cs.G_um,label="model"); q=gt[gt.method.eq("CS")]; ax[1].errorbar(q.T_C,q.G_um,yerr=q.uncertainty_G_um,fmt='o',label="PDF")
    for a in ax:a.legend(); a.set_xlabel("T (°C)")
    ax[0].set_ylabel("fractional density"); ax[1].set_ylabel("grain size (µm)"); save(fig,"CS_density_grain_fit")
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for case in ["CS_thermal","rate50_thermal"]:
        g=h[h.case.eq(case)]; ax[0].plot(g.t_s/60,g.rho,label=case); ax[1].plot(g.rho,g.G_um,label=case)
    for a in ax:a.legend()
    ax[0].set(xlabel="time (min)",ylabel="density"); ax[1].set(xlabel="density",ylabel="G (µm)"); save(fig,"fast_rate_density_grain_comparison")
    fig,ax=plt.subplots()
    for case,g in matched.groupby("case"): ax.plot(g.rho,g.G_um,label=case)
    q=pd.read_csv(DATA/"grain_size_vs_density_digitized.csv")
    for case,g in q.groupby("method"): ax.scatter(g.fractional_density,g.G_um,s=15,label=f"{case} PDF")
    ax.set(xlabel="density",ylabel="G (µm)"); ax.legend(fontsize=7,ncol=2); save(fig,"G_vs_rho_CS_HMS_LMS_model_vs_target")
    for name,cols,ylabel in [("pore_D50_D90_vs_density",["pore_D50_m","pore_D90_m"],"pore diameter (nm)"),("fine_pore_fraction_vs_density",["fine_pore_fraction"],"fine connected fraction"),("closed_pore_fraction_vs_density",["closed_fraction"],"closed pore fraction")]:
        fig,ax=plt.subplots()
        for case in ["CS_thermal","rate50_thermal","two_level_best_map_case"]:
            g=h[h.case.eq(case)]
            for col in cols: ax.plot(g.rho,g[col]*(1e9 if col.startswith('pore_D') else 1),label=f"{case} {col}")
        ax.set(xlabel="density",ylabel=ylabel); ax.legend(fontsize=7); save(fig,name)
    tau=pd.read_csv(OUT/"tau_remove_by_bin.csv"); fig,ax=plt.subplots()
    for case,g in tau.groupby("case"): ax.plot(g.rho,g.tau_remove_s.clip(upper=1e12),'.',ms=1,label=case)
    ax.set(xlabel="density",ylabel="removal time (s)",yscale="log"); ax.legend(); save(fig,"tau_remove_vs_density")
    z=pd.read_csv(OUT/"zener_pinning_histories.csv"); fig,ax=plt.subplots()
    for case,g in z.groupby("case"): ax.plot(g.rho,g.S_Z,label=case)
    ax.set(xlabel="density",ylabel="smooth Zener factor"); ax.legend(); save(fig,"zener_pinning_vs_density")
    e=pd.read_csv(OUT/"energy_balance_histories.csv"); fig,ax=plt.subplots(); g=e[e.case.eq("CS_thermal")]
    for col in ["P_surf_W_m3","P_dens_W_m3","P_excess_W_m3"]: ax.plot(g.t_s/60,g[col].clip(lower=1e-20),label=col)
    ax.set(xlabel="time (min)",ylabel="power density (W/m³)",yscale="log"); ax.legend(); save(fig,"energy_balance_power_channels")
    chen=pd.read_csv(OUT/"chen_classification_points.csv")
    for name,x,y in [("chen_map_T1_T2","T1_C","T2_C"),("chen_map_G1_T2","final_G_um","T2_C"),("chen_map_switch_density_T2","switch_density","T2_C")]:
        fig,ax=plt.subplots(); sc=ax.scatter(chen[x],chen[y],c=chen.final_rho,s=30+80*chen.chen_success,cmap="viridis"); fig.colorbar(sc,ax=ax,label="final density"); ax.set(xlabel=x,ylabel=y); save(fig,name)
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for case in ["CS_thermal","rate50_thermal","two_level_best_map_case"]:
        g=h[h.case.eq(case)]; ax[0].plot(g.t_s/3600,g.rho,label=case); ax[1].plot(g.rho,g.G_um,label=case)
    for a in ax:a.legend()
    ax[0].set(xlabel="time (h)",ylabel="density"); ax[1].set(xlabel="density",ylabel="G (µm)"); save(fig,"two_step_best_path_vs_CS_HMS")
    print(f"wrote figures to {OUT}"); return 0
if __name__=="__main__": raise SystemExit(main())
