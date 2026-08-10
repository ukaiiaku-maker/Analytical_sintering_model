#!/usr/bin/env python3
"""Bounded source-grounded mechanism discovery and Chen-map escalation."""
from __future__ import annotations
import argparse,csv,math,time
from dataclasses import asdict,replace
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import agentic_mechanism_model as model
import pore_location_agentic_sensitivity as old
import topology_constrained_sintering as aggregate

TARGET=.90;BUDGET=96*3600;TOLS=(.05,.10)
G_RED=(75.,150.,300.,600.);SW_RED=(.80,.85);T2_RED=(1000.,1100.,1200.,1250.,1300.)
G_FULL=(50.,75.,100.,150.,225.,300.,450.,600.);T1_FULL=(1250.,1300.,1350.);SW_FULL=(.75,.80,.85);T2_FULL=tuple(float(x) for x in range(900,1301,25))


def write(path,rows,fields=('parameter_id','rejection_reason')):
    cols=[]
    for r in rows:
        for k in r:
            if k not in cols:cols.append(k)
    if not cols:cols=list(fields)
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)


def design():
    modes=('persistent_junction','tj_multihit_q0','tj_multihit_q1','persistent_tj_multihit_q0','persistent_tj_multihit_q1');rows=[]
    for i in range(64):
        mode=modes[i%5];A=(6.,12.,24.,48.)[(i//5)%4];tau=(5e4,1.5e5,4.5e5,1.35e6)[(i//3)%4];prod=(.5,1.,2.,4.)[(i//7)%4];lam=(.5,1.5,4.5,13.5)[(i//11)%4];K=(1.,2.,4.,8.)[(i//13)%4];Q=(340e3,380e3,420e3,460e3)[(i//17)%4];q=1 if mode.endswith('q1') else 0
        p=model.DiscoveryParams(old.base_action('action_evolving_capture'),mode,A_J=A,tau_J_ref_s=tau,XJ_prod_TJ=prod,lambda_TJ_ref=lam,K_TJ0=K,Q_TJ_event_J_mol=Q,q_TJ=q)
        rows.append((f'mech_{i:03d}',p))
    return rows


def classify(first,rho1,rho2,growth,tol):
    if not first:return 'UNATTAINABLE_FIRST_STEP'
    if rho1>=TARGET-1e-12:return 'INELIGIBLE_TARGET_ALREADY_REACHED'
    dense=rho2>=TARGET-1e-12;small=growth<=tol+1e-12
    if dense and small:return 'SUCCESS'
    if dense:return 'GRAIN_GROWTH_FAILURE'
    if small:return 'DENSIFICATION_EXHAUSTION_FAILURE'
    return 'MIXED_FAILURE'


def fractions(h,i=-1):
    a=[float(np.sum(h[k][i])) for k in ('phi_GBseg','phi_TJ','phi_iso')];z=max(sum(a),1e-300);return [x/z for x in a]


def group(pid,p,G,T1,sw,T2s):
    lp=replace(p.action.location,base=replace(p.action.location.base,G0=G*1e-9));p=replace(p,action=replace(p.action,location=lp));h1=model.run(p,aggregate.Iso(T1,BUDGET),stop_at_rho=sw);rho1=float(h1['rho'][-1]);att=rho1>=sw-1e-12;G1=float(h1['G'][-1])*1e9;state=model.final_state(h1,p) if att else None;f1=fractions(h1);rows=[]
    for T2 in T2s:
        if att:
            h=model.run(p,aggregate.Iso(T2,BUDGET),initial=state);rho2=float(h['rho'][-1]);G2=float(h['G'][-1])*1e9;growth=(G2-G1)/G1;f2=fractions(h);last=lambda k:float(h[k][-1]) if k in h else math.nan
        else:rho2=G2=growth=math.nan;f2=[math.nan]*3;last=lambda k:math.nan
        rows.append(dict(parameter_id=pid,mechanism_mode=p.mechanism_mode,G0_nm=G,T1_C=T1,rho_switch=sw,T2_C=T2,first_step_attained=att,rho1=rho1,G1_nm=G1,rho2=rho2,G2_nm=G2,growth_fraction=growth,f_GBseg_1=f1[0],f_TJ_1=f1[1],f_iso_1=f1[2],f_GBseg_2=f2[0],f_TJ_2=f2[1],f_iso_2=f2[2],C_GBseg=last('C_GBseg'),C_TJ=last('C_TJ'),f_clean_GB=last('f_clean_GB'),activity=last('activity'),growth_mobility=last('growth_mobility_discovery'),X_J=last('X_J'),X_J_production=last('X_J_production'),X_J_relaxation=last('X_J_relaxation'),Lambda_TJ=last('Lambda_TJ'),K_TJ=last('K_TJ'),Lambda_over_K_TJ=last('Lambda_over_K_TJ'),P_comp_TJ=last('P_comp_TJ'),P_persistent_junction_drag=last('P_persistent_junction_drag'),P_TJ_multihit=last('P_TJ_multihit'),P_clean_GB=last('P_clean_GB_discovery'),sigma_act_total=last('sigma_act_total'),first_budget_h=96,second_budget_h=96))
    return rows


def evaluate(item):
    pid,p=item;tra=[]
    for G in G_RED:
        for sw in SW_RED:tra.extend(group(pid,p,G,1300.,sw,T2_RED))
    counts={};cats=set()
    for tol in TOLS:
        cc=[classify(r['first_step_attained'],r['rho1'],r['rho2'],r['growth_fraction'],tol) for r in tra];counts[tol]=cc.count('SUCCESS');cats.update(cc)
    reason=''
    if counts[.05]+counts[.10]==0:reason='no_success_window'
    elif 'DENSIFICATION_EXHAUSTION_FAILURE' not in cats:reason='missing_lower_boundary'
    elif 'GRAIN_GROWTH_FAILURE' not in cats:reason='missing_upper_boundary'
    elif all(x=='SUCCESS' for x in cats):reason='universal_success'
    return pid,p,tra,counts,cats,reason


def histories(cases):
    rows=[]
    for pid,p in cases.items():
        h=model.run(p,aggregate.RampHold(.2,1350,8*3600));stride=max(1,len(h['rho'])//250)
        for i in list(range(0,len(h['rho']),stride))+[len(h['rho'])-1]:
            row=dict(parameter_id=pid,mechanism_mode=p.mechanism_mode,t_s=float(h['t'][i]),T_C=float(h['T_C'][i]),rho=float(h['rho'][i]),G_nm=float(h['G'][i])*1e9)
            for k in ('X_J','X_J_production','X_J_relaxation','Lambda_TJ','K_TJ','Lambda_over_K_TJ','P_comp_TJ','C_GBseg','C_TJ','f_clean_GB','f_iso','P_persistent_junction_drag','P_TJ_multihit','P_clean_GB_discovery','sigma_base','sigma_GBseg_pore','sigma_TJ_pore','sigma_act_total'):row[k]=float(h[k][i]) if k in h else math.nan
            f=fractions(h,i);row.update(f_GBseg=f[0],f_TJ=f[1],f_iso_location=f[2]);rows.append(row)
    return rows


def boundaries(rows,cases):
    out=[]
    for pid in cases:
      for G in G_FULL:
       for T1 in T1_FULL:
        for sw in SW_FULL:
         q=[r for r in rows if r['parameter_id']==pid and r['G0_nm']==G and r['T1_C']==T1 and r['rho_switch']==sw]
         for tol in TOLS:
          ok=[r for r in q if classify(r['first_step_attained'],r['rho1'],r['rho2'],r['growth_fraction'],tol)=='SUCCESS'];lo=min((r['T2_C'] for r in ok),default=math.nan);hi=max((r['T2_C'] for r in ok),default=math.nan)
          out.append(dict(parameter_id=pid,G0_nm=G,T1_C=T1,rho_switch=sw,growth_tolerance=tol,T_success_lower_C=lo,T_success_upper_C=hi,window_width_C=hi-lo if ok else math.nan,n_success=len(ok)))
    return out


def plots(out,full,bounds,hist):
    cases=sorted(set(r['parameter_id'] for r in full));colors={'SUCCESS':'tab:green','GRAIN_GROWTH_FAILURE':'tab:red','DENSIFICATION_EXHAUSTION_FAILURE':'tab:blue','MIXED_FAILURE':'tab:purple','UNATTAINABLE_FIRST_STEP':'.5','INELIGIBLE_TARGET_ALREADY_REACHED':'orange'}
    if cases:
      fig,axs=plt.subplots(1,len(cases),figsize=(5*len(cases),5),squeeze=False);axs=axs.ravel()
      for ax,pid in zip(axs,cases):
       q=[r for r in full if r['parameter_id']==pid and r['growth_tolerance']==.10]
       for c,col in colors.items():
        z=[r for r in q if r['classification']==c];ax.scatter([r['G0_nm'] for r in z],[r['T2_C'] for r in z],c=col,s=12,alpha=.4,label=c)
       ax.set(xscale='log',title=pid,xlabel='G0 [nm]',ylabel='T2 [C]');ax.grid(alpha=.2)
      axs[-1].legend(fontsize=6);fig.tight_layout();fig.savefig(out/'chen_maps_by_mechanism.png',dpi=150);fig.savefig(out/'failure_mode_maps.png',dpi=150);plt.close(fig)
      fig,ax=plt.subplots();
      for pid in cases:
       q=[r for r in bounds if r['parameter_id']==pid and r['growth_tolerance']==.10 and math.isfinite(r['window_width_C'])];ax.scatter([r['G0_nm'] for r in q],[r['window_width_C'] for r in q],label=pid)
      ax.set(xlabel='G0 [nm]',ylabel='window width [C]');ax.legend();fig.tight_layout();fig.savefig(out/'window_width_vs_G0.png',dpi=150);plt.close(fig)
    if hist:
      fig,axs=plt.subplots(2,3,figsize=(15,8))
      for pid in sorted(set(r['parameter_id'] for r in hist)):
       q=[r for r in hist if r['parameter_id']==pid];x=[r['rho'] for r in q];axs[0,0].plot(x,[r['X_J'] for r in q],label=pid);axs[0,1].plot(x,[r['Lambda_over_K_TJ'] for r in q],label=pid);axs[0,2].plot(x,[r['P_comp_TJ'] for r in q],label=pid);axs[1,0].plot(x,[r['f_GBseg'] for r in q],label=pid);axs[1,1].plot(x,[r['P_persistent_junction_drag'] for r in q],label=pid);axs[1,2].plot(x,[r['sigma_act_total'] for r in q],label=pid)
      for ax,t in zip(axs.ravel(),('persistent population','Lambda/K','TJ completion','pore locations','power partition','stress partition')):ax.set(title=t,xlabel='rho');ax.grid(alpha=.2)
      axs[0,0].legend(fontsize=6);fig.tight_layout();fig.savefig(out/'mechanism_state_diagnostics.png',dpi=150);plt.close(fig)
      specs=(('Lambda_over_K_TJ','lambda_over_K_class_B.png','Lambda_TJ / K_TJ'),('X_J','persistent_population_vs_density.png','X_J'),('f_GBseg','pore_location_fractions_vs_density.png','GB-segment pore fraction'),('P_persistent_junction_drag','power_channel_partition_vs_density.png','persistent drag power'),('sigma_act_total','stress_channel_partition_vs_density.png','activation stress [Pa]'))
      for key,name,ylabel in specs:
       fig,ax=plt.subplots(figsize=(7,5))
       for pid in sorted(set(r['parameter_id'] for r in hist)):
        q=[r for r in hist if r['parameter_id']==pid];ax.plot([r['rho'] for r in q],[r[key] for r in q],label=pid)
       ax.set(xlabel='rho',ylabel=ylabel);ax.grid(alpha=.2);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(out/name,dpi=150);plt.close(fig)
      # The selected iteration does not add a stress-storage state; preserve a
      # calibration-facing stress history without implying otherwise.
      fig,ax=plt.subplots(figsize=(7,5))
      for pid in sorted(set(r['parameter_id'] for r in hist)):
       q=[r for r in hist if r['parameter_id']==pid];ax.plot([r['rho'] for r in q],[r['sigma_act_total'] for r in q],label=pid)
      ax.set(xlabel='rho',ylabel='instantaneous activation stress [Pa]',title='Stress accumulation/release candidate deferred');ax.grid(alpha=.2);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(out/'stress_accumulation_release_vs_density.png',dpi=150);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',default='results/agentic_mechanism_search');ap.add_argument('--workers',type=int,default=1);ap.add_argument('--plots-only',action='store_true');a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    if a.plots_only:
      def load(name):
       with (out/name).open(newline='') as f:rows=list(csv.DictReader(f))
       for row in rows:
        for k,v in list(row.items()):
         if k not in ('parameter_id','mechanism_mode','classification','first_step_attained'):
          try:row[k]=float(v)
          except (ValueError,TypeError):pass
       return rows
      plots(out,load('full_map_classifications.csv'),load('window_boundaries.csv'),load('activity_factor_summary.csv'));print('PLOTS REFRESHED');return
    start=time.perf_counter();items=design();registry=[]
    for pid,p in items:
        row=asdict(p);row['parameter_id']=pid;registry.append(row)
    # Flatten nested action metadata to keep the registry directly tabular.
    registry=[{k:v for k,v in r.items() if k!='action'} for r in registry];write(out/'parameter_registry.csv',registry)
    if a.workers>1:
      with ProcessPoolExecutor(max_workers=a.workers) as ex:results=list(ex.map(evaluate,items))
    else:results=[evaluate(x) for x in items]
    summary=[];rejected=[];lookup={pid:p for pid,p in items}
    for pid,p,tra,counts,cats,reason in results:
      row=dict(parameter_id=pid,mechanism_mode=p.mechanism_mode,success_5pct=counts[.05],success_10pct=counts[.10],has_lower_failure='DENSIFICATION_EXHAUSTION_FAILURE' in cats,has_upper_failure='GRAIN_GROWTH_FAILURE' in cats,rejected=bool(reason),rejection_reason=reason);summary.append(row)
      if reason:rejected.append(row)
    write(out/'reduced_map_summary.csv',summary);write(out/'rejected_parameter_sets.csv',rejected)
    surviving=sorted((r for r in summary if not r['rejected']),key=lambda r:(r['success_5pct'],r['success_10pct']),reverse=True)[:2];cases={r['parameter_id']:lookup[r['parameter_id']] for r in surviving}
    # Required regression baselines are always mapped if no mechanism survives.
    if not cases:cases={'action_baseline':model.DiscoveryParams(old.base_action('action_evolving_capture'),'action_baseline')}
    fullraw=[]
    for pid,p in cases.items():
      print('full',pid,flush=True)
      for G in G_FULL:
       for T1 in T1_FULL:
        for sw in SW_FULL:fullraw.extend(group(pid,p,G,T1,sw,T2_FULL))
    full=[]
    for r in fullraw:
      for tol in TOLS:full.append({**r,'rho_target':TARGET,'growth_tolerance':tol,'classification':classify(r['first_step_attained'],r['rho1'],r['rho2'],r['growth_fraction'],tol)})
    bounds=boundaries(fullraw,cases);histcases={'nominal_persistent':model.DiscoveryParams(old.base_action('action_evolving_capture'),'persistent_junction'),'nominal_tj_q0':model.DiscoveryParams(old.base_action('action_evolving_capture'),'tj_multihit_q0'),'nominal_tj_q1':model.DiscoveryParams(old.base_action('action_evolving_capture'),'tj_multihit_q1')};hist=histories(histcases)
    write(out/'full_map_classifications.csv',full);write(out/'window_boundaries.csv',bounds);write(out/'mechanism_decision_table.csv',summary);write(out/'activity_factor_summary.csv',hist);write(out/'power_channel_summary.csv',[{k:v for k,v in r.items() if k in ('parameter_id','rho','P_persistent_junction_drag','P_TJ_multihit','P_clean_GB_discovery')} for r in hist]);write(out/'stress_channel_summary.csv',[{k:v for k,v in r.items() if k in ('parameter_id','rho','sigma_base','sigma_GBseg_pore','sigma_TJ_pore','sigma_act_total')} for r in hist]);plots(out,full,bounds,hist);write(out/'runtime_summary.csv',[{'wall_s':time.perf_counter()-start,'screen_sets':len(items),'survivors':len(surviving),'full_points':len(full)}]);print('DONE',len(surviving),'survivors',len(full),'full classifications')


if __name__=='__main__':main()
