#!/usr/bin/env python3
"""Production figure generation for the massive topology search."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import numpy as np,pandas as pd,matplotlib.pyplot as plt
import plot_style as ps

OUT=Path('results/massive_latent_topology_search_high_density_twostep');FIG=OUT/'figures'
FIGURES=['optimizer_convergence','massive_screen_phase_map','pareto_front','best_two_step_G_rho','best_two_step_time_histories','best_two_step_topology_histories','best_Chen_filled_window_map','high_density_attainment_map','ablation_waterfall_best','robustness_heatmap_rho0_G0','qTJ_comparison','fast_firing_preservation']
def save(fig,name):ps.finish(fig,FIG/name)
def main():
 ps.apply_style();FIG.mkdir(parents=True,exist_ok=True);s0=pd.read_csv(OUT/'massive_stage0_screen.csv');s1=pd.read_csv(OUT/'stage1_reduced_dynamic_summary.csv');s2=pd.read_csv(OUT/'stage2_exact_dynamic_summary.csv');h=pd.read_csv(OUT/'topology_state_histories.csv');c=pd.read_csv(OUT/'chen_classification_points_compact.csv');best=int(s2.sort_values('median_reduction',ascending=False).iloc[0].candidate_id) if len(s2) else -1
 fig,ax=plt.subplots();z=s0.stage0_score.cummax();ax.plot(np.arange(len(z)),z);ax.set(xlabel='Retained Stage-0 rank',ylabel='Best surrogate score');save(fig,'optimizer_convergence')
 fig,ax=plt.subplots();q=s0.iloc[::max(1,len(s0)//5000)];sc=ax.scatter(q.topology_divergence,q.topology_persistence,c=q.projected_reduction,s=8);plt.colorbar(sc,ax=ax,label='Projected reduction');ax.set(xlabel='Topology divergence',ylabel='Persistence');save(fig,'massive_screen_phase_map')
 fig,ax=plt.subplots();sc=ax.scatter(s2.median_reduction,s2.rho_two_final,c=s2.span20);plt.colorbar(sc,ax=ax,label=r'$\Delta\rho$ at reduction ≥20%');ax.set(xlabel='Median two-step reduction',ylabel='Two-step final density');save(fig,'pareto_front')
 hb=h[h.candidate_id==best]
 fig,ax=plt.subplots();
 for name,g in hb.groupby('path'):ax.plot(g.rho,g.G_nm,label=name)
 ax.axvspan(.95,.98,color='.8',alpha=.4);ax.set(xlabel='Relative density',ylabel='Mean grain size [nm]');ax.legend();save(fig,'best_two_step_G_rho')
 fig,axs=plt.subplots(3,1,sharex=True,figsize=(6,7));
 for name,g in hb.groupby('path'):
  axs[0].plot(g.t/3600,g.T_C,label=name);axs[1].plot(g.t/3600,g.rho);axs[2].plot(g.t/3600,g.G_nm)
 axs[0].set(ylabel='Temperature [°C]');axs[1].set(ylabel='Density');axs[2].set(xlabel='Time [h]',ylabel='Grain size [nm]');axs[0].legend();save(fig,'best_two_step_time_histories')
 fig,axs=plt.subplots(2,2,sharex=True,figsize=(7,6));
 for name,g in hb.groupby('path'):
  for ax,k,l in zip(axs.flat,['phi_connected','phi_closed','XJ','stress'],['Connected pore volume','Closed pore volume','$X_J$','Residual stress state']):ax.plot(g.rho,g[k],label=name);ax.set_ylabel(l)
 axs[-1,0].set_xlabel('Density');axs[-1,1].set_xlabel('Density');axs[0,0].legend();save(fig,'best_two_step_topology_histories')
 cb=c[c.candidate_id==best];colors={'success':'#009E73','density_exhaustion':'#6A3D9A','grain_growth':'#D55E00','mixed':'#999999'}
 fig,ax=plt.subplots();
 for k,g in cb.groupby('classification'):ax.scatter(np.full(len(g),100),g.T2_C,label=k,color=colors[k])
 ax.set(xlabel='$G_1$ [nm]',ylabel='$T_2$ [°C]');ax.legend();save(fig,'best_Chen_filled_window_map')
 a=pd.read_csv(OUT/'high_density_attainment.csv');fig,ax=plt.subplots();ax.bar(np.arange(len(a)),a.attain_095.astype(int)+a.attain_098.astype(int));ax.set(xlabel='Exact candidate rank',ylabel='Targets attained (0–2)');save(fig,'high_density_attainment_map')
 for name,xlabel,ylabel in [('ablation_waterfall_best','Ablation','Median reduction'),('robustness_heatmap_rho0_G0','Initial grain size [nm]','Initial density'),('qTJ_comparison','TJ size exponent','Median reduction'),('fast_firing_preservation','Density','Reference/fast grain size')]:
  fig,ax=plt.subplots();ax.plot([0,1],[0,0],marker='o');ax.text(.5,.55,'Production confirmation pending',ha='center',transform=ax.transAxes);ax.set(xlabel=xlabel,ylabel=ylabel);save(fig,name)
 rows=[]
 for i,n in enumerate(FIGURES,1):rows.append(dict(figure_id=i,filename_pdf=f'figures/{n}.pdf',filename_png=f'figures/{n}.png',short_title=n,source_table_or_script='massive_latent_topology_plots.py',purpose='search diagnostic',manuscript_location='audit'))
 ps.write_inventory(OUT/'figure_inventory.csv',rows)
if __name__=='__main__':main()
