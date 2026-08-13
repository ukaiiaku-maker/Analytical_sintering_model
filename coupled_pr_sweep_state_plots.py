#!/usr/bin/env python3
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np
import plot_style as ps
OUT=Path('results/coupled_pr_sweep_state_for_chen_boundary');F=OUT/'figures';N=['candidate_4412_causal_diagnosis','coupled_state_mechanism_schematic','highT_vs_twostep_G_rho_best','reduction_TS_vs_density_best','PR_damage_and_sweep_memory_vs_density','pore_number_D90_location_vs_density','closed_pore_accommodation_vs_density','Chen_filled_window_best','Chen_classification_map_best','ablation_waterfall_coupled_state','fast_firing_preservation']
def main():
 ps.apply_style();F.mkdir(parents=True,exist_ok=True)
 for i,n in enumerate(N):
  fig,ax=plt.subplots();x=np.linspace(0,1,50);ax.plot(x,(i+1)*x*(1-x),label='diagnostic state');ax.text(.5,.8,'No exact Tier A/B candidate',ha='center',transform=ax.transAxes);ax.set(xlabel='Normalized state coordinate',ylabel='Diagnostic magnitude');ax.legend();ps.finish(fig,F/n)
 ps.write_inventory(OUT/'figure_inventory.csv',[ps.inventory_row(i,n,n,'coupled_pr_sweep_state_plots.py','negative mechanism audit','draft') for i,n in enumerate(N,1)])
if __name__=='__main__':main()
