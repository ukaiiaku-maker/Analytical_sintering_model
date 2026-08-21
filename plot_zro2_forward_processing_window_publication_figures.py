#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
OUT=Path("results/zro2_forward_processing_window_prediction_figures");MAIN=OUT/"figures_main";SUP=OUT/"figures_supplement";SRC=OUT/"source_tables"
plt.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"ps.fonttype":42,"font.size":9,"axes.labelsize":10,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"lines.linewidth":2})
C={"slow":"#333333","mid":"#0072B2","fast":"#D55E00","two":"#009E73","high":"#CC79A7","density":"#0072B2","success":"#009E73","growth":"#D55E00","mixed":"#CC79A7","un":"#BDBDBD"};inventory=[]
def panels(ax):
 for a,l in zip(np.ravel(ax),"ABCD"):a.text(.01,.99,l,transform=a.transAxes,va="top",fontweight="bold",fontsize=12)
def save(fig,name,source,root=MAIN):
 root.mkdir(exist_ok=True);fig.tight_layout()
 for ext in ("pdf","png"):fig.savefig(root/f"{name}.{ext}",dpi=600 if ext=="png" else None,bbox_inches="tight")
 inventory.append(dict(figure_id=name,figure_group=root.name,pdf_file=str((root/f'{name}.pdf').relative_to(OUT)),png_file=str((root/f'{name}.png').relative_to(OUT)),source_table=str(source.relative_to(OUT)),status="publication_candidate",validation=False));plt.close(fig)
def source(df,name):p=SRC/f"{name}.csv";df.to_csv(p,index=False);return p
def main():
 MAIN.mkdir(exist_ok=True);SUP.mkdir(exist_ok=True);SRC.mkdir(exist_ok=True)
 scr=pd.read_csv(OUT/"first_step_processing_state_screen.csv");prom=pd.read_csv(OUT/"promoted_first_step_states.csv");pts=pd.read_csv(OUT/"twostep_second_step_classification_points.csv");win=pd.read_csv(OUT/"twostep_window_boundaries.csv");hh=pd.read_csv(OUT/"heating_rate_histories.csv");hr=pd.read_csv(OUT/"heating_rate_matched_density_ratios.csv");rep=pd.read_csv(OUT/"twostep_representative_histories.csv");iso=pd.read_csv(OUT/"isothermal_comparator_histories.csv")
 # Figure 1: explicit design rather than fitted values.
 design=pd.DataFrame({"panel":["A","A","A","B","B","B","C","C","D","D","D"],"item":["barrier G*(sigma,T)","fixed diffusivities","provisional default mobility","rho0/G0/pore state","T1/switch/G1","T2/hold","693168 response form","ZrO2 independent state","coarse screen","exact cloned T2 scan","boundary/QC figures"]});sp=source(design,"main_figure_1_design")
 fig,ax=plt.subplots(2,2,figsize=(7.2,5.6));texts=[("Fixed ZrO2 physics",["stress-resolved barrier","GB/surface diffusion fixed","defined open/closed laws","default mobility provisional"]),("Processing variables",["ρ₀, G₀, pore D50/width","T₁, heating rate","ρswitch and resulting G₁","T₂ and hold time"]),("Prior reference",["693168: response form only","not ZrO2 parameterization","no state values copied","not calibration"]),("Prediction workflow",["5,020-state screen","24 promoted switch states","exact cloned-state T₂ grid","classification → figures → QC"])]
 for a,(title,lines) in zip(ax.ravel(),texts):a.axis("off");a.set_title(title);a.text(.08,.82,"\n".join(f"• {x}" for x in lines),va="top",linespacing=1.7)
 panels(ax);save(fig,"main_figure_1_fixed_model_search_design",sp)
 # Figure 2 heating response.
 bestrow=hr.loc[hr.G_reference_over_G_fast.idxmax()];case=bestrow.case_id;z=hh[hh.case_id.eq(case) & hh.schedule_id.str.contains(f"p{int(bestrow.T_peak_C)}_h{int(bestrow.hold_h)}")];sp=source(z,"main_figure_2_heating_response")
 fig,ax=plt.subplots(2,2,figsize=(7.2,5.8));
 for run,g in z.groupby("run_id"):
  rate=float(run.split("_")[-3]);color=C["slow"] if rate<=1 else C["fast"] if rate>=20 else C["mid"];lab=f"{rate:g} °C min⁻¹";ax[0,0].plot(g.physical_time_h,g.T_C,color=color,label=lab);ax[0,1].plot(g.T_C,g.rho,color=color,label=lab);ax[1,0].plot(g.T_C,g.G_nm,color=color,label=lab);ax[1,1].plot(g.rho,g.G_nm,color=color,label=lab)
 ax[0,0].set(xlabel="Physical time (h)",ylabel="Temperature (°C)");ax[0,1].set(xlabel="Temperature (°C)",ylabel="Relative density");ax[1,0].set(xlabel="Temperature (°C)",ylabel="Grain size (nm)");ax[1,1].set(xlabel="Relative density",ylabel="Grain size (nm)");ax[0,0].legend(ncol=2);panels(ax);save(fig,"main_figure_2_heating_rate_response",sp)
 # Figure 3 fast map.
 sp=source(hr,"main_figure_3_fast_firing_map");fig,ax=plt.subplots(2,2,figsize=(7.2,5.8));best=hr[(hr.case_id==case)&(hr.T_peak_C==bestrow.T_peak_C)&(hr.hold_h==bestrow.hold_h)]
 for rate,g in best.groupby("fast_rate_C_min"):ax[0,0].plot(g.rho,g.G_reference_over_G_fast,label=f"{rate:g}")
 pivot=hr.groupby(["case_id","T_peak_C"]).G_reference_over_G_fast.max().unstack();im=ax[0,1].imshow(pivot,aspect="auto",cmap="cividis",vmin=1);fig.colorbar(im,ax=ax[0,1],label="max Gref/Gfast");ax[0,1].set(yticks=np.arange(len(pivot)),yticklabels=pivot.index,xlabel="Peak T column",ylabel="Initial case")
 spans=[]
 for t in (1.2,1.5,2):spans.append((t,float(best.loc[best.G_reference_over_G_fast>=t,"rho"].max()-best.loc[best.G_reference_over_G_fast>=t,"rho"].min()) if (best.G_reference_over_G_fast>=t).any() else 0))
 ax[1,0].bar([str(x[0]) for x in spans],[x[1] for x in spans],color=C["mid"]);ax[1,0].set(xlabel="Ratio threshold",ylabel="Density span")
 diag=z.groupby("run_id").agg(activity=("activity","median"),Gstar=("Gstar_eV","median"),PR=("PR_memory","max")).reset_index();u=np.arange(len(diag));ax[1,1].plot(u,diag.activity,"o-",label="activity");ax[1,1].plot(u,diag.PR,"s-",label="PR memory");ax[1,1].set(xlabel="Schedule index",ylabel="State diagnostic");ax[1,1].legend();ax[0,0].set(xlabel="Relative density",ylabel="Gref/Gfast");ax[0,0].legend(title="Fast rate");panels(ax);save(fig,"main_figure_3_matched_density_fast_firing_map",sp)
 # Figure 4 preparation map.
 sample=scr.sample(min(2000,len(scr)),random_state=1);sp=source(sample,"main_figure_4_processing_state_map");fig,ax=plt.subplots(2,2,figsize=(7.2,5.8));col=np.where(sample.promotable,C["success"],C["un"]);ax[0,0].scatter(sample.actual_rho_switch,sample.G1_nm,c=col,s=5);ax[0,1].scatter(sample.G1_nm,sample.phi_closed_total,c=sample.A_closed,cmap="cividis",s=6);p=ax[1,0].scatter(sample.rho_switch,sample.G0_nm,c=sample.first_step_growth_fraction,cmap="viridis",s=6);fig.colorbar(p,ax=ax[1,0],label="First-step growth");ax[1,1].scatter(sample.G1_nm,sample.PR_memory,c=C["un"],s=3);ax[1,1].scatter(prom.G1_nm,prom.PR_memory,c=C["success"],s=14,label="promoted");ax[0,0].set(xlabel="Actual switch density",ylabel="G₁ (nm)");ax[0,1].set(xlabel="G₁ (nm)",ylabel="Closed fraction");ax[1,0].set(xlabel="Requested switch density",ylabel="G₀ (nm)");ax[1,1].set(xlabel="G₁ (nm)",ylabel="PR memory");ax[1,1].legend();panels(ax);save(fig,"main_figure_4_twostep_preparation_state_map",sp)
 # Figure 5 Chen axes.
 q=pts[(pts.density_target==.90)&(pts.hold_h==96)&(pts.tolerance_tier=="Tier_like_B")];sp=source(q,"main_figure_5_chen_map");fig,ax=plt.subplots(2,2,figsize=(7.2,5.8));order=["DENSIFICATION_EXHAUSTION_FAILURE","SUCCESS","GRAIN_GROWTH_FAILURE","MIXED_FAILURE"];colors=[C["density"],C["success"],C["growth"],C["mixed"]]
 for cl,co in zip(order,colors):g=q[q.classification==cl];ax[0,0].scatter(g.G1_nm,g.T2_C,c=co,s=12,label=cl.replace("_"," ").title())
 f=win[(win.density_target==.90)&(win.hold_h==96)&(win.tolerance_tier=="Tier_like_B")&win.finite_window];ax[0,0].scatter(f.G1_nm,f.lower_boundary_C,marker="v",c=C["density"]);ax[0,0].scatter(f.G1_nm,f.upper_boundary_C,marker="^",c=C["growth"]);ax[0,0].set(xlabel="First-step grain size G₁ (nm)",ylabel="Second-step T₂ (°C)");ax[0,0].legend(fontsize=5)
 for a,y,label in [(ax[0,1],q.final_rho,"Final density"),(ax[1,0],q.second_step_growth_fraction,"Second-step growth fraction")]:s=a.scatter(q.G1_nm,q.T2_C,c=y,cmap="cividis",s=10);fig.colorbar(s,ax=a,label=label);a.set(xlabel="G₁ (nm)",ylabel="T₂ (°C)")
 ax[1,1].bar(f.case_id,f.window_width_C,color=C["success"]);ax[1,1].set(xlabel="Case",ylabel="Finite window width (°C)");ax[1,1].tick_params(axis="x",rotation=45);panels(ax);save(fig,"main_figure_5_chen_style_G1_T2_map",sp)
 # Figure 6 comparators.
 sp=source(iso,"main_figure_6_twostep_isothermal");fig,ax=plt.subplots(2,2,figsize=(7.2,5.8));cm={"two_step":C["two"],"high_T_isothermal":C["high"],"low_T_isothermal":C["mid"]}
 for typ,g in iso.groupby("path_type"):
  co=cm.get(typ,C["un"]);ax[0,0].plot(g.physical_time_h,g.T_C,color=co,label=typ);ax[0,1].plot(g.physical_time_h,g.rho,color=co,label=typ);ax[1,0].plot(g.physical_time_h,g.G_nm,color=co,label=typ);ax[1,1].plot(g.rho,g.G_nm,color=co,label=typ)
 ax[0,0].set(xlabel="Physical time (h)",ylabel="Temperature (°C)");ax[0,1].set(xlabel="Physical time (h)",ylabel="Density");ax[1,0].set(xlabel="Physical time (h)",ylabel="G (nm)");ax[1,1].set(xlabel="Density",ylabel="G (nm)");ax[0,0].legend();panels(ax);save(fig,"main_figure_6_twostep_vs_isothermal",sp)
 # Figure 7 state evolution.
 z=rep[rep.path_type.eq("two_step")];sp=source(z,"main_figure_7_state_evolution");fig,ax=plt.subplots(2,2,figsize=(7.2,5.8))
 for run,g in z.groupby("run_id"):
  ax[0,0].plot(g.physical_time_h,g.rho_dot_open,label=run);ax[0,0].plot(g.physical_time_h,g.rho_dot_closed,ls="--");ax[0,1].plot(g.physical_time_h,g.phi_closed_total);ax[0,1].plot(g.physical_time_h,g.A_closed_available,ls="--");ax[1,0].plot(g.physical_time_h,g.pore_D50_nm);ax[1,0].plot(g.physical_time_h,g.pore_D90_nm,ls="--");ax[1,1].plot(g.physical_time_h,g.Gamma_migration);ax[1,1].plot(g.physical_time_h,g.G_dot_intrinsic/g.G_dot_intrinsic.max(),ls="--")
 ax[0,0].set(xlabel="Physical time (h)",ylabel="Density flux (s⁻¹)");ax[0,1].set(xlabel="Physical time (h)",ylabel="Closed fraction / A");ax[1,0].set(xlabel="Physical time (h)",ylabel="Pore diameter (nm)");ax[1,1].set(xlabel="Physical time (h)",ylabel="Activity / normalized intrinsic rate");ax[0,0].legend(fontsize=5);panels(ax);save(fig,"main_figure_7_mechanistic_state_evolution",sp)
 # Figure 8 sensitivity.
 q=pts[(pts.density_target==.90)&(pts.hold_h==96)&(pts.tolerance_tier=="Tier_like_B")];sp=source(q,"main_figure_8_initial_switch_sensitivity");fig,ax=plt.subplots(2,2,figsize=(7.2,5.8));code={x:i for i,x in enumerate(order)};s=ax[0,0].scatter(q.rho0,q.G0_nm,c=q.classification.map(code),cmap="viridis",s=8);ax[0,1].scatter(q.rho_switch,q.G1_nm,c=q.classification.map(code),cmap="viridis",s=8);ax[1,0].scatter(q.pore_D50_nm,q.final_rho,c=q.T2_C,cmap="cividis",s=8);ax[1,1].scatter(q.pore_D50_nm,q.final_G_nm,c=q.T2_C,cmap="cividis",s=8);ax[0,0].set(xlabel="ρ₀",ylabel="G₀ (nm)");ax[0,1].set(xlabel="ρswitch",ylabel="G₁ (nm)");ax[1,0].set(xlabel="Initial pore D50 (nm)",ylabel="Final density");ax[1,1].set(xlabel="Initial pore D50 (nm)",ylabel="Final G (nm)");panels(ax);save(fig,"main_figure_8_initial_switch_state_sensitivity",sp)
 # Figure 9 prior vs ZrO2.
 ref=pd.DataFrame([dict(source="prior reduced mechanism reference",G1_nm=117,lower_C=925,upper_C=1200,parameterization="not ZrO2"),dict(source="ZrO2 fixed forward",G1_nm=float(f.G1_nm.min()) if len(f) else np.nan,lower_C=float(f.lower_boundary_C.min()) if len(f) else np.nan,upper_C=float(f.upper_boundary_C.max()) if len(f) else np.nan,parameterization="independently searched")]);sp=source(ref,"main_figure_9_prior_reference")
 fig,ax=plt.subplots(1,2,figsize=(7.2,3.2));ax[0].errorbar(ref.G1_nm,(ref.lower_C+ref.upper_C)/2,yerr=(ref.upper_C-ref.lower_C)/2,fmt="o",capsize=5);ax[0].set(xticks=[0,1],xticklabels=["prior reference","ZrO₂ forward"],ylabel="T₂ interval (°C)");ax[1].axis("off");ax[1].text(.05,.9,"693168\n• response form only\n• not ZrO₂ parameterization\n• not calibration\n\nZrO₂ result\n• fixed physics\n• state searched independently\n• not validated",va="top");save(fig,"main_figure_9_prior_reference_vs_zro2",sp)
 # Supplements.
 sp=source(win,"supplement_window_boundaries");fig,ax=plt.subplots();
 for tier,g in win.groupby("tolerance_tier"):ax.scatter(g.G1_nm,g.window_width_C,label=tier,s=10)
 ax.set(xlabel="G₁ (nm)",ylabel="Window width (°C)");ax.legend();save(fig,"supplement_figure_1_all_window_boundaries",sp,SUP)
 sp=source(scr.sample(min(2000,len(scr)),random_state=2),"supplement_screen");fig,ax=plt.subplots();ax.hist(scr.first_step_growth_fraction.clip(-.1,1),bins=50,color=C["mid"]);ax.axvline(.2,color=C["growth"]);ax.set(xlabel="First-step growth fraction",ylabel="State count");save(fig,"supplement_figure_2_screen_distribution",sp,SUP)
 pd.DataFrame(inventory).to_csv(OUT/"figure_inventory.csv",index=False)
 (OUT/"figure_caption_drafts.md").write_text("# Figure caption drafts\n\n"+"\n".join(f"- **{r['figure_id']}**: Fixed-parameter ZrO2 forward prediction; no fitting, no validation. Source: `{r['source_table']}`." for r in inventory)+"\n")
if __name__=="__main__":main()
