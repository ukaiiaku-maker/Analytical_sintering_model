#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
OUT=Path('results/zro2_forward_resolved_rules')
def save(fig,name):fig.suptitle('Resolved-rule ZrO2 forward model — conditional, no validation claim',fontsize=10);fig.tight_layout(rect=[0,0,1,.96]);fig.savefig(OUT/f'{name}.png',dpi=180);fig.savefig(OUT/f'{name}.pdf');plt.close(fig)
def main():
 fig,ax=plt.subplots(figsize=(11,3));ax.axis('off');labels=['fixed JSON barrier\nserial cycle','open shrinkage\nconnected path','conservative PR\npreparation memory','closed transition\nbounded accommodation','intrinsic growth ×\nmigration activity'];xs=np.linspace(.08,.92,len(labels))
 for x,l in zip(xs,labels):ax.text(x,.5,l,ha='center',va='center',bbox=dict(boxstyle='round',fc='#e8f1fa'))
 for a,b in zip(xs[:-1],xs[1:]):ax.annotate('',(b-.08,.5),(a+.08,.5),arrowprops=dict(arrowstyle='->'))
 save(fig,'resolved_rule_mechanism_schematic')
 m=pd.read_csv(OUT/'resolved_rule_pdf_conditioned_matched_density.csv');fig,ax=plt.subplots();ax.plot(m.rho,m.G_um_5,label='5 C/min');ax.plot(m.rho,m.G_um_50,label='50 C/min');ax.set(xlabel='density',ylabel='G (µm)');ax.legend();save(fig,'pdf_conditioned_fast_rate_G_vs_rho')
 fig,ax=plt.subplots(1,3,figsize=(12,4));ax[0].plot(m.rho,m.pore_D90_m_5*1e9,label='5');ax[0].plot(m.rho,m.pore_D90_m_50*1e9,label='50');ax[1].plot(m.rho,m.fine_pore_fraction_5);ax[1].plot(m.rho,m.fine_pore_fraction_50);ax[2].plot(m.rho,m.R_Z_eff_m_5*1e6);ax[2].plot(m.rho,m.R_Z_eff_m_50*1e6);ax[0].set(ylabel='D90 (nm)');ax[1].set(ylabel='fine fraction');ax[2].set(ylabel='R_Z (µm)');[a.set_xlabel('density') for a in ax];ax[0].legend();save(fig,'pdf_conditioned_fast_rate_pore_state_vs_rho')
 h=pd.read_csv(OUT/'resolved_rule_smoke_histories.csv');fig,ax=plt.subplots(1,3,figsize=(12,4));
 for k,g in h.groupby('case'):ax[0].plot(g.t_s/3600,g.closed_fraction,label=k);ax[1].plot(g.t_s/3600,g.A_closed);ax[2].plot(g.t_s/3600,g.PR_memory)
 ax[0].set(ylabel='closed fraction');ax[1].set(ylabel='A_closed');ax[2].set(ylabel='PR memory');[a.set_xlabel('time (h)') for a in ax];ax[0].legend(fontsize=7);save(fig,'closed_pore_accommodation_histories')
 c=pd.read_csv(OUT/'resolved_rule_chen_classification_points.csv');codes=pd.Categorical(c.classification).codes
 fig,ax=plt.subplots();sc=ax.scatter(c.T1_C,c.T2_C,c=codes,cmap='tab10');ax.set(xlabel='T1 (°C)',ylabel='T2 (°C)');save(fig,'chen_map_T1_T2_resolved_rules')
 fig,ax=plt.subplots();ax.scatter(c.G1_um,c.T2_C,c=codes,cmap='tab10');ax.set(xlabel='G1 (µm)',ylabel='T2 (°C)');save(fig,'chen_map_G1_T2_resolved_rules')
 b=pd.read_csv(OUT/'resolved_rule_chen_window_boundaries.csv');fig,ax=plt.subplots();sc=ax.scatter(b.T1_C,b.switch_density,c=b.boundary_gap_C,cmap='coolwarm',vmin=-300,vmax=300,s=30+2*b.hold_h);fig.colorbar(sc,ax=ax,label='upper − lower boundary (°C)');ax.set(xlabel='T1 (°C)',ylabel='switch density');save(fig,'chen_boundary_gap_resolved_rules')
 a=pd.read_csv(OUT/'resolved_rule_ablation_summary.csv');fig,ax=plt.subplots(figsize=(10,5));ax.barh(a.ablation,a.success_count,color=np.where(a.controlling_expected,'#d95f02','#1b9e77'));ax.set_xlabel('strict success count in conditional 0.1× mobility context');save(fig,'ablation_summary_resolved_rules')
 v=pd.read_csv(OUT/'resolved_rules_vs_previous_forward_baseline.csv');g=v[v.metric.isin(['strict_success_count','finite_window_count','smallest_required_shift_C'])];fig,ax=plt.subplots();x=np.arange(len(g));ax.bar(x-.2,g.previous,.4,label='previous');ax.bar(x+.2,g.resolved,.4,label='resolved');ax.set(xticks=x,xticklabels=g.metric);plt.setp(ax.get_xticklabels(),rotation=25,ha='right');ax.legend();save(fig,'previous_vs_resolved_boundary_gap')
if __name__=='__main__':main()
