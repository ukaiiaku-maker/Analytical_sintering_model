#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
OUT=Path('results/zro2_forward_mobility_envelope_audit')
def save(fig,name):fig.suptitle('Intrinsic GB mobility uncertainty audit — no validation claim',fontsize=10);fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(OUT/f'{name}.png',dpi=180);fig.savefig(OUT/f'{name}.pdf');plt.close(fig)
def main():
 f=pd.read_csv(OUT/'mobility_fast_rate_summary.csv');x=np.arange(len(f));labels=f.mobility_case_id
 fig,ax=plt.subplots(1,2,figsize=(12,4));ax[0].plot(x,f.final_rho_5,'o-',label='5 C/min');ax[0].plot(x,f.final_rho_50,'o-',label='50 C/min');ax[1].plot(x,f.final_G_um_5,'o-');ax[1].plot(x,f.final_G_um_50,'o-');ax[0].set(ylabel='final density');ax[1].set(ylabel='final G (µm)');[a.set(xticks=x,xticklabels=labels) for a in ax];[plt.setp(a.get_xticklabels(),rotation=70,ha='right',fontsize=6) for a in ax];ax[0].legend();save(fig,'mobility_envelope_CS_fast_summary')
 fig,ax=plt.subplots();ax.plot(x,f.boundary_gap_C,'o-');ax.axhline(0,c='k',ls='--');ax.axhline(-125,c='r',ls=':',label='previous required shift');ax.set(xticks=x,xticklabels=labels,ylabel='boundary gap (°C)');plt.setp(ax.get_xticklabels(),rotation=70,ha='right',fontsize=6);ax.legend();save(fig,'mobility_envelope_boundary_gap')
 c=pd.read_csv(OUT/'mobility_chen_classification_points.csv');fig,ax=plt.subplots();sc=ax.scatter(np.log10(c.M0_factor),c.T2_C,c=c.final_G_um,cmap='plasma',s=15+30*c.strict_success);fig.colorbar(sc,ax=ax,label='final G (µm)');ax.set(xlabel='log10 M0 factor',ylabel='T2 (°C)');save(fig,'mobility_envelope_chen_map')
 fig,ax=plt.subplots();ax.plot(x,f.matched_fast_smaller_G_fraction,'o-',label='smaller G fraction');ax.plot(x,f.matched_fast_smaller_D90_fraction,'o-',label='smaller D90 fraction');ax.axhline(.5,c='k',ls='--');ax.set(xticks=x,xticklabels=labels,ylabel='matched-density fraction');plt.setp(ax.get_xticklabels(),rotation=70,ha='right',fontsize=6);ax.legend();save(fig,'mobility_envelope_fast_signs')
 p=pd.read_csv(OUT/'mobility_pathway_consistency.csv');cols=['lower_boundary_present','upper_boundary_present','PR_closed_accommodation_active','pathway_consistency_flag','no_schedule_label_leakage'];fig,ax=plt.subplots(figsize=(8,6));ax.imshow(p[cols].astype(float),aspect='auto',cmap='RdYlGn',vmin=0,vmax=1);ax.set(xticks=range(len(cols)),xticklabels=cols,yticks=range(len(p)),yticklabels=p.mobility_case_id);plt.setp(ax.get_xticklabels(),rotation=35,ha='right');save(fig,'mobility_envelope_pathway_consistency')
if __name__=='__main__':main()
