#!/usr/bin/env python3
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np,pandas as pd
import plot_style as ps
OUT=Path('results/pr_lower_bound_plus_coalescence_twostep');FIG=OUT/'figures';N=['PR_lower_bound_T2_scan','energy_partition_vs_T2','highT_vs_twostep_G_rho','pore_topology_evolution_highT_vs_twostep','Chen_filled_window_map','Chen_classification_map','ablation_waterfall','fast_firing_preservation','PR_lower_bound_parameter_map']
def sv(f,n):ps.finish(f,FIG/n)
def main():
 ps.apply_style();FIG.mkdir(parents=True,exist_ok=True);d=pd.read_csv(OUT/'stageC_refined_boundary_summary.csv');win=int(d[d.complete].iloc[0].candidate_id) if d.complete.any() else int(d.iloc[0].candidate_id);best=int(d.sort_values('median_reduction',ascending=False).iloc[0].candidate_id);pts=pd.read_csv(OUT/'Chen_classification_points.csv');h=pd.read_csv(OUT/'PR_energy_partition_histories.csv');rat=pd.read_csv(OUT/'high_density_ratio_curves.csv');a=pd.read_csv(OUT/'ablation_summary.csv')
 q=pts[pts.candidate_id==win];fig,axs=plt.subplots(2,2,figsize=(8,6));axs[0,0].plot(q.T2_C,q.rho_final);axs[0,1].plot(q.T2_C,q.growth_fraction);axs[1,0].plot(q.T2_C,q.mean_w_PR,label='PR');axs[1,0].plot(q.T2_C,q.mean_w_dens,label='dens');axs[1,1].scatter(q.T2_C,q.classification.astype('category').cat.codes);[ax.set_xlabel('$T_2$ [°C]') for ax in axs.flat];axs[0,0].set_ylabel('Final density');axs[0,1].set_ylabel('Growth fraction');axs[1,0].set_ylabel('Mean partition');axs[1,1].set_ylabel('Class code');axs[1,0].legend();sv(fig,N[0])
 fig,ax=plt.subplots();ax.plot(q.T2_C,q.mean_w_PR,label='$w_{PR}$');ax.plot(q.T2_C,q.mean_w_dens,label='$w_{dens}$');ax.set(xlabel='$T_2$ [°C]',ylabel='Mean energy partition');ax.legend();sv(fig,N[1])
 z=h[h.candidate_id==best];fig,ax=plt.subplots();[ax.plot(g.rho,g.G_nm,label=k) for k,g in z.groupby('path')];ax.set(xlabel='Density',ylabel='Mean grain size [nm]');ax.legend();sv(fig,N[2])
 fig,axs=plt.subplots(2,2,figsize=(8,6));
 for path,g in z.groupby('path'):
  for ax,k,l in zip(axs.flat,['D90_nm','phi_connected_fine','phi_isolated','phi_closed'],['D90 [model nm]','Connected fine volume','Isolated volume','Closed volume']):ax.plot(g.rho,g[k],label=path);ax.set_ylabel(l)
 axs[0,0].legend();[ax.set_xlabel('Density') for ax in axs[-1]];sv(fig,N[3])
 colors={'success':'#009E73','density_exhaustion':'#6A3D9A','grain_growth':'#D55E00','mixed':'#999999'};fig,ax=plt.subplots();
 for k,g in q.groupby('classification'):ax.scatter(np.full(len(g),100),g.T2_C,color=colors[k],label=k)
 ax.set(xlabel='$G_1$ [nm]',ylabel='$T_2$ [°C]');ax.legend();sv(fig,N[4]);fig,ax=plt.subplots();
 for k,g in q.groupby('classification'):ax.scatter(g.T2_C,g.rho_final,color=colors[k],label=k)
 ax.set(xlabel='$T_2$ [°C]',ylabel='Final density');ax.legend();sv(fig,N[5])
 fig,ax=plt.subplots();
 if len(a):ax.bar(a.ablation,a.median_reduction);ax.tick_params(axis='x',rotation=45)
 ax.set(ylabel='Median reduction');sv(fig,N[6]);fig,ax=plt.subplots();ax.text(.5,.5,'Not promoted: no joint Tier A/B candidate',ha='center');ax.set(xlabel='Density',ylabel='Reference/fast grain size');sv(fig,N[7]);fig,ax=plt.subplots();ax.scatter(d.median_reduction,d.window_width_C,c=d.lower_bracketed);ax.set(xlabel='Median high-density reduction',ylabel='Chen window width [°C]');sv(fig,N[8]);ps.write_inventory(OUT/'figure_inventory.csv',[ps.inventory_row(i,n,n,'pr_lower_bound_plus_coalescence_plots.py','mechanism audit','draft') for i,n in enumerate(N,1)])
if __name__=='__main__':main()
