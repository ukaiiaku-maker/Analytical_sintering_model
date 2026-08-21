#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=Path("results/zro2_forward_final_summary_figures");BASE=Path("results/zro2_forward_processing_window_prediction_figures");MAIN=OUT/"figures_main";SUP=OUT/"figures_supplement";SRC=OUT/"source_tables"
plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42,"font.size":9,"axes.labelsize":10.5,"xtick.labelsize":8.5,"ytick.labelsize":8.5,"legend.fontsize":8,"lines.linewidth":2.2,"lines.markersize":5})
C={"slow":"#333333","intermediate":"#0072B2","fast":"#D55E00","comparator":"#CC79A7","low":"#0072B2","success":"#009E73","high":"#D55E00","mixed":"#CC79A7","unattained":"#BDBDBD","ineligible":"#666666"};inv=[]
def labels(ax):
 for a,l in zip(ax.ravel(),"ABCDEF"):a.text(.015,.985,l,transform=a.transAxes,ha="left",va="top",fontweight="bold",fontsize=12)
def style(a):a.grid(color=".9",lw=.45);a.tick_params(direction="out")
def save(fig,name,source,root=MAIN,panels=1):
 root.mkdir(exist_ok=True);fig.tight_layout(h_pad=1.5,w_pad=1.2);pdf=root/f"{name}.pdf";png=root/f"{name}.png";fig.savefig(pdf,bbox_inches="tight");fig.savefig(png,dpi=600,bbox_inches="tight");inv.append(dict(figure_id=name,figure_group=root.name,pdf_file=str(pdf.relative_to(OUT)),png_file=str(png.relative_to(OUT)),source_table=str(source.relative_to(OUT)),panel_count=panels,validation=False));plt.close(fig)
def main():
 MAIN.mkdir(exist_ok=True);SUP.mkdir(exist_ok=True);h=pd.read_csv(OUT/"final_heating_rate_histories.csv");t=pd.read_csv(OUT/"final_twostep_histories.csv");chen=pd.read_csv(OUT/"final_chen_map_source.csv")
 # Final Figure 1.
 fig,ax=plt.subplots(3,2,figsize=(178/25.4,160/25.4));colors={"0.2 C/min":C["slow"],"1 C/min":C["intermediate"],"100 C/min":C["fast"]}
 for path,g in h.groupby("path_label",sort=False):
  co=colors[path];ax[0,0].plot(g.physical_time_h,g.T_C,color=co,label=path);ax[0,1].plot(g.rho,g.G_um,color=co,label=path);ax[1,0].plot(g.physical_time_h,g.G_um,color=co,label=path);ax[1,1].plot(g.physical_time_h,g.rho,color=co,label=path);ax[2,0].plot(g.T_C,g.densification_rate,color=co,label=path);ax[2,1].plot(g.T_C,g.sigma_eff_MPa,color=co,label=path)
 ax[0,0].set(xlabel="Time (h)",ylabel="Temperature (°C)");ax[0,1].set(xlabel=r"Relative density, $\rho$",ylabel=r"Grain size, $G$ (µm)",yscale="log");ax[1,0].set(xlabel="Time (h)",ylabel=r"Grain size, $G$ (µm)",yscale="log");ax[1,1].set(xlabel="Time (h)",ylabel=r"Relative density, $\rho$");ax[2,0].set(xlabel="Temperature (°C)",ylabel=r"Densification rate, $d\rho/dt$ (s$^{-1}$)");ax[2,0].set_yscale("symlog",linthresh=1e-10);ax[2,1].set(xlabel="Temperature (°C)",ylabel=r"Effective stress, $\sigma_{eff}$ (MPa)")
 for a in ax.ravel():style(a);a.legend(loc="best",frameon=False)
 labels(ax);save(fig,"final_fig1_heating_rate_response",SRC/"final_fig1_heating_rate_response_source.csv",panels=6)
 # Final Figure 2.
 fig,ax=plt.subplots(3,2,figsize=(178/25.4,164/25.4));tc={"high-T isothermal comparator":C["comparator"],"low T2 density failure":C["low"],"central strict success":C["success"],"high T2 growth failure":C["high"]};short={"high-T isothermal comparator":"1500 °C comparator","low T2 density failure":"Low $T_2$","central strict success":"Central success","high T2 growth failure":"High $T_2$"}
 for path,g in t.groupby("path_label",sort=False):
  co=tc[path];lab=short[path];ax[0,0].plot(g.physical_time_h,g.T_C,color=co,label=lab);ax[0,1].plot(g.rho,g.G_um,color=co,label=lab);ax[1,0].plot(g.physical_time_h,g.G_um,color=co,label=lab);ax[1,1].plot(g.physical_time_h,g.rho,color=co,label=lab);ax[2,0].plot(g.physical_time_h,g.rho_dot_total,color=co,label=f"{lab}: total")
 central=t[t.path_label=="central strict success"];ax[2,0].plot(central.physical_time_h,central.rho_dot_open,color=C["success"],ls="--",label="Success: open");ax[2,0].plot(central.physical_time_h,central.rho_dot_closed,color=C["success"],ls=":",label="Success: closed")
 switch=t[t.path_label=="central strict success"].iloc[np.argmin(np.abs(t[t.path_label=="central strict success"].rho-.831987))];
 for a in ax[:2].ravel():a.scatter([switch.physical_time_h if a in (ax[0,0],ax[1,0],ax[1,1]) else switch.rho],[switch.T_C if a is ax[0,0] else switch.G_um if a in (ax[0,1],ax[1,0]) else switch.rho],marker="D",s=28,c="white",edgecolor="black",zorder=5)
 ax[1,1].axhline(.90,color=".35",ls="--",lw=1.2,label=r"target $\rho=0.90$")
 order=["DENSIFICATION_EXHAUSTION_FAILURE","SUCCESS","GRAIN_GROWTH_FAILURE","MIXED_FAILURE","UNATTAINABLE_FIRST_STEP","INELIGIBLE_TARGET_ALREADY_REACHED"];cc=[C["low"],C["success"],C["high"],C["mixed"],C["unattained"],C["ineligible"]]
 map_labels=["Density exhaustion","Success","Growth failure","Mixed failure","Unattainable","Target already reached"]
 for cl,co,lab in zip(order,cc,map_labels):
  g=chen[chen.classification==cl];ax[2,1].scatter(g.G1_nm,g.T2_C,c=co,s=13,label=lab,alpha=.8)
 b=chen[chen.finite_window.fillna(False)].drop_duplicates("case_id");ax[2,1].scatter(b.G1_nm,b.lower_boundary_C,marker="v",c=C["low"],s=30);ax[2,1].scatter(b.G1_nm,b.upper_boundary_C,marker="^",c=C["high"],s=30);ax[2,1].scatter([37.816]*3,[850,925,1025],marker="*",s=65,c=[C["low"],C["success"],C["high"]],edgecolors="black",linewidths=.4)
 ax[0,0].set(xlabel="Physical time (h)",ylabel="Temperature (°C)");ax[0,1].set(xlabel=r"Relative density, $\rho$",ylabel=r"Grain size, $G$ (µm)",yscale="log");ax[1,0].set(xlabel="Physical time (h)",ylabel=r"Grain size, $G$ (µm)",yscale="log");ax[1,1].set(xlabel="Physical time (h)",ylabel=r"Relative density, $\rho$");ax[2,0].set(xlabel="Physical time (h)",ylabel=r"Density rate, $d\rho/dt$ (s$^{-1}$)");ax[2,0].set_yscale("symlog",linthresh=1e-10);ax[2,1].set(xlabel=r"First-step grain size, $G_1$ (nm)",ylabel=r"Second-step temperature, $T_2$ (°C)")
 for a in ax.ravel():
  style(a)
  if a is ax[2,0]:a.legend(loc="upper right",frameon=False,fontsize=5.1,ncol=2,columnspacing=.7,handlelength=1.8,handletextpad=.35)
  elif a is ax[2,1]:a.legend(loc="upper left",frameon=True,facecolor="white",edgecolor="none",framealpha=.88,fontsize=4.8,ncol=2,columnspacing=.6,handletextpad=.25,borderpad=.25)
  else:a.legend(loc="best",frameon=False,fontsize=6.2,handlelength=1.8,handletextpad=.4)
 labels(ax);save(fig,"final_fig2_twostep_vs_isothermal_response",SRC/"final_fig2_twostep_vs_isothermal_response_source.csv",panels=6)
 # Standalone Chen map.
 fig,a=plt.subplots(figsize=(178/25.4,115/25.4))
 for cl,co,lab in zip(order,cc,map_labels):g=chen[chen.classification==cl];a.scatter(g.G1_nm,g.T2_C,c=co,s=22,label=lab,alpha=.82)
 a.scatter(b.G1_nm,b.lower_boundary_C,marker="v",c=C["low"],s=45,label="Lower boundary");a.scatter(b.G1_nm,b.upper_boundary_C,marker="^",c=C["high"],s=45,label="Upper boundary");a.scatter([37.816]*3,[850,925,1025],marker="*",s=90,c=[C["low"],C["success"],C["high"]],edgecolors="black",label="Selected paths");a.set(xlabel=r"First-step grain size, $G_1$ (nm)",ylabel=r"Second-step temperature, $T_2$ (°C)");style(a);a.legend(loc="lower center",bbox_to_anchor=(.5,1.01),ncol=3,frameon=False,fontsize=7,columnspacing=.9,handletextpad=.35);save(fig,"final_fig3_standalone_chen_map",SRC/"final_fig3_standalone_chen_map_source.csv")
 # Supplements from available source data.
 full=pd.read_csv(BASE/"heating_rate_histories.csv");q=full[(full.case_id=="P001")&full.schedule_id.str.contains("p1500_h0")];q.to_csv(SRC/"supp_all_heating_rate_family_source.csv",index=False);fig,a=plt.subplots(figsize=(7,4.5));
 for path,g in q.groupby("schedule_id"):a.plot(g.rho,g.G_um,label=path.split("_")[0][1:]+" °C/min")
 a.set(xlabel=r"Relative density, $\rho$",ylabel=r"Grain size, $G$ (µm)",yscale="log");style(a);a.legend(ncol=2,frameon=False);save(fig,"all_heating_rate_family",SRC/"supp_all_heating_rate_family_source.csv",SUP)
 ratio=pd.read_csv(OUT/"final_heating_rate_summary.csv");ratio.to_csv(SRC/"supp_heating_rate_ratio_source.csv",index=False);fig,a=plt.subplots(figsize=(7,4.5));
 for fast,g in ratio.groupby("fast_rate_C_min"):a.plot(g.rho,g.G_reference_over_G_fast,label=f"0.2/{fast:g} °C min⁻¹")
 for y in (1.2,1.5,2):a.axhline(y,color=".6",ls="--",lw=1)
 a.set(xlabel=r"Relative density, $\rho$",ylabel=r"$G_{ref}/G_{fast}$");style(a);a.legend(frameon=False);save(fig,"heating_rate_matched_density_ratio",SRC/"supp_heating_rate_ratio_source.csv",SUP)
 t.to_csv(SRC/"supp_twostep_state_variables_source.csv",index=False);fig,aa=plt.subplots(2,2,figsize=(7,5.5))
 for path,g in t.groupby("path_label"):
  co=tc[path];aa[0,0].plot(g.physical_time_h,g.phi_closed_total,color=co,label=path);aa[0,1].plot(g.physical_time_h,g.A_closed_available,color=co,label=path);aa[1,0].plot(g.physical_time_h,g.PR_memory,color=co,label=path);aa[1,1].plot(g.physical_time_h,g.pore_D50_nm,color=co,label=path)
 for a,y in zip(aa.ravel(),["Closed pore fraction","Accommodation available","PR memory","Pore D50 (nm)"]):a.set(xlabel="Physical time (h)",ylabel=y);style(a);a.legend(fontsize=6,frameon=False)
 fig.tight_layout();save(fig,"twostep_state_variables",SRC/"supp_twostep_state_variables_source.csv",SUP,panels=4)
 chen.to_csv(SRC/"supp_twostep_full_map_source.csv",index=False);fig,a=plt.subplots(figsize=(7,4.5));
 for cl,co in zip(order,cc):g=chen[chen.classification==cl];a.scatter(g.G1_nm,g.T2_C,c=co,s=15,label=cl.replace("_"," ").title())
 a.set(xlabel=r"First-step grain size, $G_1$ (nm)",ylabel=r"Second-step temperature, $T_2$ (°C)");style(a);a.legend(ncol=2,fontsize=6,frameon=False);save(fig,"twostep_window_full_map",SRC/"supp_twostep_full_map_source.csv",SUP)
 pd.DataFrame(inv).to_csv(OUT/"final_figure_inventory.csv",index=False)
 (OUT/"final_caption_drafts.md").write_text("# Final caption drafts\n\n**Figure 1.** Fixed-parameter heating-rate prediction for P001 at 1500 °C with zero hold. Heating rate alone varies (0.2, 1, 100 °C min⁻¹); density is ρ and grain size is G. Densification uses model `rho_dot_total`; effective stress is in MPa. Matched-density separation is read from panel B and the ratio table. Not a fit or validation.\n\n**Figure 2.** Fixed-parameter P014 trajectories from the identical ρ≈0.832, G₁≈37.82 nm switch state: 850 °C density failure, 925 °C strict success, 1025 °C growth failure, and 1500 °C comparator. Physical time is continuous. Chen colors are blue/green/red/purple/gray for density failure/success/growth failure/mixed/ineligible. Finite windows require both lower and upper boundaries. Not validation.\n")
if __name__=="__main__":main()
