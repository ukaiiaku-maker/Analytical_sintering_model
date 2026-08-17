#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
OUT=Path('results/zro2_forward_port_defined_closed_laws');FIG=OUT/'figures'
def save(fig,n):fig.tight_layout();fig.savefig(FIG/f'{n}.png',dpi=180);fig.savefig(FIG/f'{n}.pdf');plt.close(fig)
def lines(file,name,y='final_rho'):
 x=pd.read_csv(OUT/file);fig,ax=plt.subplots(figsize=(11,5));
 for mode,g in x.groupby('mode'):ax.plot(g.T2_C,g[y],label=mode)
 if y=='final_rho':ax.axhline(.976,color='k',ls='--')
 ax.set(xlabel='T2 (°C)',ylabel=y);ax.legend(fontsize=6);save(fig,name)
def main():
 FIG.mkdir(parents=True,exist_ok=True);r=pd.read_csv(OUT/'defined_closed_law_registry.csv');fig,ax=plt.subplots(figsize=(12,6));ax.axis('off');ax.table(cellText=r[['law_id','affected_process','implementation_status']].values,colLabels=['law','process','status'],loc='center',fontsize=6);save(fig,'defined_closed_law_architecture')
 f=pd.read_csv(OUT/'fixed_path_flux_integrals.csv');g=f.groupby('mode')[['Delta_rho_open','Delta_rho_closed']].max();fig,ax=plt.subplots(figsize=(11,5));g.plot.bar(ax=ax);ax.set_ylabel('integrated density channel');save(fig,'renewal_closed_vs_open_cycle')
 fig,ax=plt.subplots(figsize=(9,4));ax.text(.5,.5,'open pore → precursor/isolated → closed pore\nPR and closure transfers conserve pore volume\nonly named shrinkage changes density',ha='center',va='center',fontsize=13);ax.axis('off');save(fig,'PR_closed_preparation_flow')
 x=pd.read_csv(OUT/'natural_state_T2_scan_by_mode.csv');fig,ax=plt.subplots(figsize=(11,5));
 for mode,q in x.groupby('mode'):ax.plot(q.T2_C,q.closed_inventory_formed,label=mode)
 ax.set(xlabel='T2 (°C)',ylabel='closed inventory formed');ax.legend(fontsize=6);save(fig,'closed_inventory_accommodation_history')
 lines('candidate_state_T2_scan_by_mode.csv','candidate_state_T2_classification_by_mode');lines('natural_state_T2_scan_by_mode.csv','natural_state_T2_classification_by_mode')
 fig,ax=plt.subplots(figsize=(9,4));ax.text(.5,.5,'No naturally prepared mode passed\nthe fixed-path topology gate\nMini-map not run',ha='center',va='center');ax.axis('off');save(fig,'mini_map_window_by_mode')
 x=pd.read_csv(OUT/'phenomenological_to_mechanistic_mapping_status.csv');z=pd.crosstab(x.mapping_status,x.controls_lower_boundary);fig,ax=plt.subplots(figsize=(9,6));ax.imshow(z.values,aspect='auto',cmap='Blues');ax.set(yticks=np.arange(len(z)),yticklabels=z.index,xticks=np.arange(len(z.columns)),xticklabels=z.columns,xlabel='controls lower boundary');save(fig,'phenomenological_to_mechanistic_status_matrix')
 x=pd.read_csv(OUT/'candidate693168_defined_law_comparison.csv');fig,ax=plt.subplots(figsize=(11,5));ax.bar(x.model,x.closed_density_contribution);plt.setp(ax.get_xticklabels(),rotation=35,ha='right');ax.set_ylabel('closed density contribution');save(fig,'candidate693168_vs_forward_defined_laws')
if __name__=='__main__':main()
