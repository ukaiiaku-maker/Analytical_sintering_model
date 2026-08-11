#!/usr/bin/env python3
"""Compact production tables and regenerate path-resolved Figures 5–6."""
from dataclasses import replace
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import production_mechanism_assessment as prod
import agentic_mechanism_model as model
import topology_constrained_sintering as aggregate


def prepared_state(row):
    p0=prod.frozen_mechanisms()[row.mechanism_id];base=replace(p0.action.location.base,G0=row.G0_nm*1e-9);p=replace(p0,action=replace(p0.action,location=replace(p0.action.location,base=base)));h=model.run(p,prod.prep.FixedBudgetRamp(row.heating_rate_C_min,row.T1_C));idx=int(np.flatnonzero(h['rho']>=row.rho_switch-1e-12)[0]);return p,h,model.final_state(h,p,idx),idx


def history_rows(label,phase,h):
    stride=max(1,len(h['rho'])//350);rows=[]
    for i in list(range(0,len(h['rho']),stride))+[len(h['rho'])-1]:
        rows.append(dict(path_label=label,phase=phase,t_s=float(h['t'][i]),T_C=float(h['T_C'][i]),rho=float(h['rho'][i]),G_nm=float(h['G'][i])*1e9,X_J=float(h['X_J'][i]) if 'X_J' in h else 0,Lambda_over_K_TJ=float(h['Lambda_over_K_TJ'][i]) if 'Lambda_over_K_TJ' in h else np.nan,P_comp_TJ=float(h['P_comp_TJ'][i]) if 'P_comp_TJ' in h else np.nan,C_GBseg=float(h['C_GBseg'][i]),C_TJ=float(h['C_TJ'][i]),f_clean_GB=float(h['f_clean_GB'][i]),P_dens=float(h['P_GBseg_dens'][i]+h['P_TJ_dens'][i]),P_clean_GB=float(h['P_clean_GB_discovery'][i]) if 'P_clean_GB_discovery' in h else float(h['P_clean_GB'][i]),P_persistent_junction_drag=float(h['P_persistent_junction_drag'][i]) if 'P_persistent_junction_drag' in h else 0,P_TJ_multihit=float(h['P_TJ_multihit'][i]) if 'P_TJ_multihit' in h else 0,sigma_base=float(h['sigma_base'][i]),sigma_GBseg=float(h['sigma_GBseg_pore'][i]),sigma_TJ=float(h['sigma_TJ_pore'][i]),sigma_act_total=float(h['sigma_act_total'][i])))
    return rows


def main():
    out=Path('results/production_mechanism_assessment');figdir=out/'figures';figdir.mkdir(exist_ok=True)
    # Preserve point-level comparisons as ignored raw evidence, then make the
    # named production table compact enough for review.
    src=out/'production_fast_firing_summary.csv';raw=out/'raw_fast_firing_comparisons.csv'
    if not raw.exists():shutil.copyfile(src,raw)
    ff=pd.read_csv(raw);compact=ff.groupby(['mechanism_id','rho_target','initial_topology','peak_T_C','heating_rate_C_min','response_class'],dropna=False).agg(n_cases=('G0_nm','size'),HR_pct_median=('HR_pct','median'),HR_pct_min=('HR_pct','min'),HR_pct_max=('HR_pct','max')).reset_index();compact.to_csv(src,index=False)
    failed=pd.read_csv(out/'failed_or_censored_cases.csv');failed.groupby(['mechanism_id','map_type','growth_tolerance','boundary_status','tier','G0_nm'],dropna=False).size().reset_index(name='n_cases').to_csv(out/'failed_or_censored_cases.csv',index=False)
    success=pd.read_csv(out/'successful_practical_windows.csv');rows=[]
    representative=success[success.tier=='Tier_A'].sort_values(['mechanism_id','window_width_C'],ascending=[True,False]).groupby('mechanism_id').head(1)
    for _,r in representative.iterrows():
        p,h1,state,idx=prepared_state(r);rows+=history_rows(f'{r.mechanism_id}_TierA','first_step', {k:(v[:idx+1] if hasattr(v,'__len__') and k!='pore_radii' else v) for k,v in h1.items()});h2=model.run(p,aggregate.Iso(r.T_first_success_C,prod.BUDGET),initial=state);rows+=history_rows(f'{r.mechanism_id}_TierA','second_step',h2)
        low=max(800,r.T_lower_density_C-10);hf=model.run(p,aggregate.Iso(low,prod.BUDGET),initial=state);rows+=history_rows(f'{r.mechanism_id}_lower_failure','second_step',hf)
    # One representative fast-firing harmful path at rho=0.90.
    rawpaths=pd.read_csv(out/'raw_fast_firing_paths.csv');harm=pd.read_csv(raw);hrow=harm[(harm.response_class=='harmful')&(harm.rho_target==.90)].sort_values('HR_pct').iloc[0];p0=prod.frozen_mechanisms()[hrow.mechanism_id];frac=prod.TOPOLOGIES[hrow.initial_topology];p=prod.fast_params(p0,hrow.G0_nm,hrow.rho0,frac);hh=model.run(p,prod.FastSchedule(hrow.heating_rate_C_min,hrow.peak_T_C,hrow.hold_time_h));rows+=history_rows(f'{hrow.mechanism_id}_fast_harmful','fast_firing',hh)
    diag=pd.DataFrame(rows);diag.to_csv(out/'representative_path_diagnostics.csv',index=False)
    fig,axs=plt.subplots(2,3,figsize=(15,8))
    for label,q in diag.groupby('path_label'):
        axs[0,0].plot(q.rho,q.X_J,label=label);axs[0,1].plot(q.rho,q.Lambda_over_K_TJ,label=label);axs[0,2].plot(q.rho,q.P_comp_TJ,label=label);axs[1,0].plot(q.rho,q.C_GBseg,label=label);axs[1,1].plot(q.rho,q.C_TJ,label=label);axs[1,2].plot(q.rho,q.f_clean_GB,label=label)
    for ax,y in zip(axs.ravel(),('X_J','Lambda/K','P_comp,TJ','C_GBseg','C_TJ','clean GB fraction')):ax.set(xlabel='relative density',ylabel=y);ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=6);fig.tight_layout();fig.savefig(figdir/'figure5_mechanism_state_diagnostics.png',dpi=200);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(13,5))
    for label,q in diag.groupby('path_label'):
        axs[0].plot(q.rho,q.P_dens,label=f'{label}: dens');axs[0].plot(q.rho,q.P_clean_GB,ls='--',label=f'{label}: clean GB');axs[0].plot(q.rho,q.P_persistent_junction_drag,ls=':',label=f'{label}: XJ drag');axs[0].plot(q.rho,q.P_TJ_multihit,ls='-.',label=f'{label}: TJ multihit');axs[1].plot(q.rho,q.sigma_act_total,label=label)
    axs[0].set(xlabel='relative density',ylabel='power channel [model units]',yscale='symlog');axs[1].set(xlabel='relative density',ylabel='activation stress [Pa]');axs[0].legend(fontsize=5,ncol=2);axs[1].legend(fontsize=6)
    for ax in axs:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(figdir/'figure6_power_stress_partition.png',dpi=200);plt.close(fig)


if __name__=='__main__':main()
