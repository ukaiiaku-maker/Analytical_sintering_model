#!/usr/bin/env python3
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT=Path("results/zro2_forward_required_chen_physics_gap_analysis");FIG=OUT/"figures"
def save(fig,name):fig.suptitle("Read-only diagnostic — strict result preserved — no validation claim",fontsize=10);fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(FIG/f"{name}.png",dpi=180);fig.savefig(FIG/f"{name}.pdf");plt.close(fig)
def main():
    FIG.mkdir(parents=True,exist_ok=True);b=pd.read_csv(OUT/"boundary_gap_by_state_ranked.csv");finite=b.dropna(subset=['gap_C'])
    fig,ax=plt.subplots();ax.hist(finite.gap_C,bins=np.arange(-225,-50,25),edgecolor='k');ax.axvline(0,c='k',ls='--');ax.set(xlabel="strict gap: upper growth − lower density (°C)",ylabel="first-step groups");save(fig,"boundary_gap_histogram")
    fig,ax=plt.subplots();sc=ax.scatter(finite.T1_C,finite.switch_density,c=finite.gap_C,s=30+2*finite.hold_h,cmap='coolwarm',vmin=-200,vmax=200);fig.colorbar(sc,ax=ax,label="gap (°C)");ax.set(xlabel="T1 (°C)",ylabel="switch density");save(fig,"boundary_gap_map_T1_switch")
    fig,ax=plt.subplots();g=finite.sort_values(['T1_C','switch_density','hold_h']).reset_index();x=np.arange(len(g));ax.plot(x,g.T_lower_density_C,'o',label='lower density boundary');ax.plot(x,g.T_upper_growth_C,'o',label='upper growth boundary');ax.fill_between(x,g.T_upper_growth_C,g.T_lower_density_C,alpha=.2);ax.set(xlabel="ranked bracketed first-step state",ylabel="T2 boundary (°C)");ax.legend();save(fig,"lower_upper_boundary_overlay")
    fig,ax=plt.subplots();g=finite.sort_values('required_shift_C');ax.bar(np.arange(len(g)),g.required_shift_C);ax.axhline(25,c='k',ls='--',label='desired width');ax.set(xlabel="ranked first-step state",ylabel="required boundary shift (°C)");ax.legend();save(fig,"required_shift_ranked_states")
    c=pd.read_csv(OUT/"boundary_sensitivity_coefficients.csv");fig,ax=plt.subplots(figsize=(9,5));x=np.arange(len(c));ax.barh(x-.18,c.dT_lower_density.fillna(0),height=.35,label='lower boundary response');ax.barh(x+.18,c.dT_upper_growth.fillna(0),height=.35,label='upper boundary response');ax.set(yticks=x,yticklabels=c.parameter,xlabel="coarse finite-difference response (°C / stated coordinate)");ax.legend();ax.text(.99,.02,"gap slopes not ranked: ≤2 finite gap points per lever",transform=ax.transAxes,ha='right');save(fig,"OAT_boundary_sensitivity_tornado")
    t=pd.read_csv(OUT/"threshold_relaxation_transition_table.csv");fig,ax=plt.subplots(1,2,figsize=(10,4));p=t.pivot(index='density_target',columns='grain_threshold_um',values='success_count');im=ax[0].imshow(p,aspect='auto',cmap='Blues');fig.colorbar(im,ax=ax[0],label='successes');w=t.pivot(index='density_target',columns='grain_threshold_um',values='finite_window_count');im=ax[1].imshow(w,aspect='auto',cmap='Oranges');fig.colorbar(im,ax=ax[1],label='finite windows')
    for a in ax:a.set(xticks=range(len(p.columns)),xticklabels=p.columns,yticks=range(len(p.index)),yticklabels=p.index,xlabel="grain threshold (µm)",ylabel="density target")
    save(fig,"threshold_relaxation_phase_diagram")
    m=pd.read_csv(OUT/"common_state_sensitivity_main_effects.csv");fig,ax=plt.subplots(1,3,figsize=(12,4));
    for a,f in zip(ax,["pore_D50_nm","rho_start","G_start_nm"]):g=m[m.factor==f];a.plot(g.level,g.mean_final_rho_50,'o-');a.set(xlabel=f,ylabel="mean final rho, 50 C/min")
    save(fig,"common_state_main_effects_fast_density")
    fig,ax=plt.subplots(figsize=(9,4));factors=m.factor.unique();ax.bar(factors,[np.nan]*len(factors));ax.set(ylim=(0,1),ylabel="boundary-gap effect");ax.text(.5,.55,"Not identifiable: factorial stored one representative\ntwo-step schedule per state, not boundary maps.",transform=ax.transAxes,ha='center',va='center',fontsize=13);plt.setp(ax.get_xticklabels(),rotation=25,ha='right');save(fig,"common_state_main_effects_boundary_gap")
    p=pd.read_csv(OUT/"common_state_pathway_consistency_summary.csv");cols=['common_interval_attained','fast_smaller_grain_gate','fast_smaller_D90_gate','finite_strict_Chen_window','both_boundaries_present','highT_comparator_gates_available','density_ok','grain_ok','all_pathway_gates_pass'];fig,ax=plt.subplots(figsize=(10,8));ax.imshow(p[cols].astype(float),aspect='auto',cmap='RdYlGn',vmin=0,vmax=1);ax.set(xticks=range(len(cols)),xticklabels=cols,yticks=[]);plt.setp(ax.get_xticklabels(),rotation=35,ha='right');save(fig,"pathway_consistency_failure_matrix")
    ph=pd.read_csv(OUT/"required_physics_gap_interpretation.csv");fig,ax=plt.subplots(figsize=(13,7));ax.axis('off');table=ax.table(cellText=[[textwrap.fill(a,24),textwrap.fill(b,30)] for a,b in zip(ph.failure_mode,ph.experimental_measurement_needed)],colLabels=['candidate failure mode','measurement needed before model change'],loc='center',cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(7);table.scale(1,1.7);save(fig,"recommended_next_data_targets")
if __name__=="__main__":main()
