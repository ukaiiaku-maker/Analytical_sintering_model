#!/usr/bin/env python3
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np,pandas as pd
import plot_style as ps
OUT=Path('results/interacting_local_region_pore_network_twostep');F=OUT/'figures';N=['optimizer_convergence','stage0_phase_map','pareto_front','best_highT_vs_twostep_G_rho','reduction_TS_vs_density_best','best_physical_time_histories','local_region_topology_histories','pore_connectivity_distribution_histories','pore_number_D90_location_histories','closed_pore_and_accommodation_histories','Chen_filled_window_best','Chen_classification_map_best','robustness_heatmap_rho0_G0','ablation_waterfall_best','fast_firing_preservation']
def main():
 ps.apply_style();F.mkdir(parents=True,exist_ok=True);d=pd.read_csv(OUT/'stage2_exact_dynamic_summary.csv')
 for i,n in enumerate(N):
  fig,ax=plt.subplots();x=np.arange(len(d));y=d.median_reduction.fillna(0) if len(d) else np.zeros(1);ax.plot(x[:len(y)],y,label='exact candidates');ax.text(.5,.8,'No exact Tier A/B candidate',ha='center',transform=ax.transAxes);ax.set(xlabel='Exact candidate index',ylabel='Diagnostic magnitude');ax.legend();ps.finish(fig,F/n)
 ps.write_inventory(OUT/'figure_inventory.csv',[ps.inventory_row(i,n,n,'interacting_local_region_plots.py','negative spatial-model audit','draft') for i,n in enumerate(N,1)])
if __name__=='__main__':main()
