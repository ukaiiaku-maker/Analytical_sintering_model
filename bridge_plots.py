#!/usr/bin/env python3
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np,pandas as pd
import plot_style as ps
OUT=Path('results/bridge_pr_lower_bound_and_high_density_reduction');FIG=OUT/'figures';N=['anchor_comparison_summary','anchor_T2_scans','anchor_G_rho_comparison','modular_transplant_outcome_map','parameter_morph_heatmap','optimizer_convergence','best_bridge_Chen_filled_window','best_bridge_G_rho_high_density','best_bridge_pore_topology_histories','best_bridge_energy_partition','bridge_ablation_waterfall','fast_firing_preservation']
def sv(f,n):ps.finish(f,FIG/n)
def main():
 ps.apply_style();FIG.mkdir(parents=True,exist_ok=True);a=pd.read_csv(OUT/'anchor_diagnosis_155976.csv');b=pd.read_csv(OUT/'anchor_diagnosis_4412.csv');t=pd.read_csv(OUT/'modular_transplant_results.csv');m=pd.read_csv(OUT/'parameter_morph_results.csv');o=pd.read_csv(OUT/'optimizer_trace.csv');e=pd.read_csv(OUT/'exact_bridge_confirmation.csv')
 fig,ax=plt.subplots();ax.bar(['155976 reduction','4412 reduction','155976 window','4412 window'],[a.median_reduction.iloc[0],b.median_reduction.iloc[0],a.window_width_C.iloc[0]/100,b.window_width_C.iloc[0]/100]);ax.set(ylabel='Normalized objective');sv(fig,N[0])
 fig,axs=plt.subplots(2,1,sharex=True);axs[0].plot(a.T2_C,a.rho_final,label='155976');axs[0].plot(b.T2_C,b.rho_final,label='4412');axs[1].plot(a.T2_C,a.growth_fraction);axs[1].plot(b.T2_C,b.growth_fraction);axs[0].set(ylabel='Final density');axs[1].set(xlabel='$T_2$ [°C]',ylabel='Growth fraction');axs[0].legend();sv(fig,N[1])
 fig,ax=plt.subplots();ax.plot([.95,.98],[a.median_reduction.iloc[0]]*2,label='155976');ax.plot([.95,.98],[b.median_reduction.iloc[0]]*2,label='4412');ax.set(xlabel='Density',ylabel='Median reduction proxy');ax.legend();sv(fig,N[2])
 fig,ax=plt.subplots();ax.scatter(t.median_reduction,t.window_width_C,c=t.complete);ax.set(xlabel='Median reduction',ylabel='Coarse window width [°C]');sv(fig,N[3])
 fig,ax=plt.subplots();sc=ax.scatter(m['lambda'],m.median_reduction,c=m.window_width_C);plt.colorbar(sc,ax=ax,label='Coarse window [°C]');ax.set(xlabel='Morph fraction',ylabel='Median reduction');sv(fig,N[4])
 fig,ax=plt.subplots();ax.plot(o.evaluations,o.best_score,label='best');ax.plot(o.evaluations,o.median_score.fillna(-2),label='median');ax.set(xlabel='Evaluations',ylabel='Objective');ax.legend();sv(fig,N[5])
 best=e.sort_values('median_reduction',ascending=False).iloc[0];fig,ax=plt.subplots();ax.bar(['lower','upper','width/100'],[best.lower_bracketed,best.upper_bracketed,best.window_width_C/100]);ax.set(ylabel='Boundary indicator');sv(fig,N[6])
 fig,ax=plt.subplots();ax.plot([.95,.98],[best.median_reduction]*2);ax.axhline(.2,ls='--',color='.4');ax.set(xlabel='Density',ylabel='Exact reduction');sv(fig,N[7])
 for n,y in [(N[8],'Pore topology state'),(N[9],'Energy partition'),(N[10],'Ablation response'),(N[11],'Reference/fast grain size')]:fig,ax=plt.subplots();ax.plot([0,1],[0,0],marker='o');ax.text(.5,.55,'No accepted bridge candidate',ha='center',transform=ax.transAxes);ax.set(xlabel='Normalized coordinate',ylabel=y);sv(fig,n)
 ps.write_inventory(OUT/'figure_inventory.csv',[ps.inventory_row(i,n,n,'bridge_plots.py','bridge audit','draft') for i,n in enumerate(N,1)])
if __name__=='__main__':main()
