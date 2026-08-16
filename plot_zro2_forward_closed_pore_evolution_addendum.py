#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,matplotlib.pyplot as plt
OUT=Path('results/zro2_forward_closed_pore_evolution_addendum');FIG=OUT/'figures'
def save(fig,n):fig.tight_layout();fig.savefig(FIG/f'{n}.png',dpi=180);fig.savefig(FIG/f'{n}.pdf');plt.close(fig)
def main():
 FIG.mkdir(parents=True,exist_ok=True);x=pd.read_csv(OUT/'T2_state_scan.csv')
 for kind in ('candidate_like_injected','naturally_prepared'):
  fig,ax=plt.subplots(figsize=(10,5));
  for mode,g in x[x.state_kind.eq(kind)].groupby('mode'):ax.plot(g.T2_C,g.final_rho,label=mode)
  ax.axhline(.976,color='k',ls='--');ax.set(xlabel='T2 (°C)',ylabel='final density');ax.legend(fontsize=6);save(fig,f'{kind}_density')
 a=pd.read_csv(OUT/'boundary_topology_acceptance.csv');fig,ax=plt.subplots(figsize=(10,5));pd.crosstab(a['mode'],a['acceptance']).plot.bar(ax=ax);ax.set_ylabel('state-test count');save(fig,'boundary_topology_acceptance')
 q=pd.read_csv(OUT/'apparent_closed_rate_slopes.csv');fig,ax=plt.subplots(figsize=(9,4));ax.bar(q['mode'],q.Q_closed_app_kJ_mol);plt.setp(ax.get_xticklabels(),rotation=30,ha='right');ax.set_ylabel('Q_closed_app (kJ/mol), diagnostic only');save(fig,'apparent_closed_rate_slopes')
if __name__=='__main__':main()
