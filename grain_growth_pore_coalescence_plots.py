#!/usr/bin/env python3
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np,pandas as pd
import plot_style as ps
OUT=Path('results/grain_growth_pore_coalescence_memory');FIG=OUT/'figures'
NAMES=['highT_vs_twostep_G_rho','two_step_reduction_vs_density','pore_number_and_D90_vs_density','pore_location_fractions_vs_density','sweep_coalescence_memory_vs_density','physical_time_histories_highT_vs_twostep','Chen_filled_window_map','ablation_waterfall_pore_sweep_memory','robustness_heatmap_rho0_G0','fast_firing_preservation']
def save(f,n):ps.finish(f,FIG/n)
def main():
 ps.apply_style();FIG.mkdir(parents=True,exist_ok=True);e=pd.read_csv(OUT/'exact_dynamic_summary.csv');cid=int(e.sort_values('median_reduction',ascending=False).query('attained').iloc[0].candidate_id);h=pd.read_csv(OUT/'pore_coalescence_histories.csv').query('candidate_id==@cid');r=pd.read_csv(OUT/'two_step_ratio_curves.csv').query('candidate_id==@cid');a=pd.read_csv(OUT/'ablation_summary.csv');c=pd.read_csv(OUT/'Chen_window_boundaries.csv')
 fig,ax=plt.subplots();[ax.plot(g.rho,g.G_nm,label=k) for k,g in h.groupby('path')];ax.set(xlabel='Relative density',ylabel='Mean grain size [nm]');ax.legend();save(fig,NAMES[0])
 fig,ax=plt.subplots();ax.plot(r.rho,r.reduction);ax.axhline(.2,ls='--',color='.4');ax.set(xlabel='Relative density',ylabel='Two-step grain-size reduction');save(fig,NAMES[1])
 fig,axs=plt.subplots(2,1,sharex=True);[axs[0].plot(g.rho,g.pore_number_reduction_factor,label=k) for k,g in h.groupby('path')];[axs[1].plot(g.rho,g.D90_nm,label=k) for k,g in h.groupby('path')];axs[0].set_ylabel('Pore-number proxy');axs[1].set(xlabel='Density',ylabel='D90 [model nm]');axs[0].legend();save(fig,NAMES[2])
 fig,ax=plt.subplots();
 for k in ['phi_connected_fine','phi_large_attached','phi_large_TJ','phi_isolated','phi_closed']:ax.plot(h[h.path=='two_step'].rho,h[h.path=='two_step'][k],label=k)
 ax.set(xlabel='Density',ylabel='Pore volume fraction');ax.legend(fontsize=6);save(fig,NAMES[3])
 fig,ax=plt.subplots();
 for path,g in h.groupby('path'):ax.plot(g.rho,g.coalesced_pore_fraction,label=path)
 ax.set(xlabel='Density',ylabel='Cumulative coalesced pore volume');ax.legend();save(fig,NAMES[4])
 fig,axs=plt.subplots(3,1,sharex=True);[axs[0].plot(g.t/3600,g.T_C,label=k) for k,g in h.groupby('path')];[axs[1].plot(g.t/3600,g.rho) for _,g in h.groupby('path')];[axs[2].plot(g.t/3600,g.G_nm) for _,g in h.groupby('path')];axs[0].set_ylabel('T [°C]');axs[1].set_ylabel('Density');axs[2].set(xlabel='Time [h]',ylabel='G [nm]');axs[0].legend();save(fig,NAMES[5])
 fig,ax=plt.subplots();ax.scatter(c.window_width_C,c.lower_bracketed.astype(int),c=c.upper_bracketed.astype(int));ax.set(xlabel='Window width [°C]',ylabel='Lower boundary present');save(fig,NAMES[6])
 fig,ax=plt.subplots();
 if len(a):ax.bar(a.ablation,a.median_reduction);ax.tick_params(axis='x',rotation=45)
 ax.set(ylabel='Median reduction');save(fig,NAMES[7])
 fig,ax=plt.subplots();ax.scatter(e.rho_two_final,e.median_reduction,c=e.span20);ax.set(xlabel='Two-step final density',ylabel='Median reduction');save(fig,NAMES[8])
 fig,ax=plt.subplots();ax.text(.5,.5,'Not promoted: no complete Chen window',ha='center');ax.set(xlabel='Density',ylabel='Reference/fast grain size');save(fig,NAMES[9])
 ps.write_inventory(OUT/'figure_inventory.csv',[ps.inventory_row(i,n,n,'grain_growth_pore_coalescence_plots.py','mechanism audit','draft') for i,n in enumerate(NAMES,1)])
if __name__=='__main__':main()
