#!/usr/bin/env python3
"""Bounded topology-action ablations, screen, histories, and Chen maps."""
from __future__ import annotations
import argparse,csv,math,time
from dataclasses import asdict,replace
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import pore_location_agentic_model as action
import pore_location_topology_model as fixed
import pore_location_topology_sensitivity as placement
import topology_constrained_sintering as aggregate

TARGET=.90;TOLS=(.05,.10);BUDGET=96*3600
MODES=('fixed_flux_baseline','action_static','action_evolving_no_capture','action_evolving_capture','action_evolving_capture_high_TJ_drag')
G_RED=(75.,150.,300.,600.);SW_RED=(.80,.85);T2_RED=(1000.,1100.,1200.,1250.,1300.)
G_FULL=(50.,75.,100.,150.,225.,300.,450.,600.);T1_FULL=(1250.,1300.,1350.);SW_FULL=(.75,.80,.85);T2_FULL=tuple(float(x) for x in range(900,1301,25))


def write_table(path,rows,empty=('parameter_id','rejection_reason')):
    fields=[]
    for row in rows:
        for key in row:
            if key not in fields:fields.append(key)
    if not fields:fields=list(empty)
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def base_action(mode='action_evolving_capture',**kw):
    lp=kw.pop('location',fixed.LocationParams(placement.base_params(),'evolving'));return action.ActionParams(lp,mode,**kw)


def classify(attained,rho2,growth,tol):
    if not attained:return 'UNATTAINABLE_FIRST_STEP'
    dense=rho2>=TARGET-1e-12;bounded=growth<=tol+1e-12
    if dense and bounded:return 'SUCCESS'
    if dense:return 'GRAIN_GROWTH_FAILURE'
    if bounded:return 'DENSIFICATION_EXHAUSTION_FAILURE'
    return 'MIXED_FAILURE'


def fractions(h,i=-1):
    v=[float(np.sum(h[k][i])) for k in ('phi_GBseg','phi_TJ','phi_iso')];z=max(sum(v),1e-300);return tuple(x/z for x in v)


ACTION_DIAG=tuple(f'action_{kind}_{n}' for kind in ('propensity','weight') for n in action.ACTION_NAMES)
FLUX_DIAG=tuple(f'action_flux_{n}' for n in ('GBseg_remove','TJ_remove','GB_smooth','GB_to_TJ','TJ_to_GBseg_capture','TJ_to_iso'))
POWERS=('P_GBseg_dens','P_TJ_dens','P_clean_GB','P_GBseg_drag','P_TJ_drag','P_GBseg_to_TJ','P_TJ_to_GBseg_capture','P_TJ_iso')
STATE_DIAG=('C_GBseg','C_TJ','f_clean_GB','sigma_GBseg_pore','sigma_TJ_pore','sigma_clean_GB','sigma_iso','sigma_act_total','activity','growth_mobility_factor','E_G')


def run_group(case_id,p,G,T1,sw,T2s):
    lp=replace(p.location,base=replace(p.location.base,G0=G*1e-9));p=replace(p,location=lp);first=action.run(p,aggregate.Iso(T1,BUDGET),stop_at_rho=sw);attained=float(first['rho'][-1])>=sw-1e-12;rho1=float(first['rho'][-1]);G1=float(first['G'][-1])*1e9;f1=fractions(first) if attained else (math.nan,)*3;state=fixed.final_state(first,action.effective_location_params(p)) if attained else None;rows=[]
    for T2 in T2s:
        if attained:
            h=action.run(p,aggregate.Iso(T2,BUDGET),initial=state);rho2=float(h['rho'][-1]);G2=float(h['G'][-1])*1e9;growth=(G2-G1)/G1;f2=fractions(h);extra={k:float(h[k][-1]) if k in h else math.nan for k in ACTION_DIAG+FLUX_DIAG+POWERS+STATE_DIAG}
        else:rho2=G2=growth=math.nan;f2=(math.nan,)*3;extra={k:math.nan for k in ACTION_DIAG+FLUX_DIAG+POWERS+STATE_DIAG}
        rows.append(dict(case_id=case_id,action_mode=p.action_mode,G0_nm=G,T1_C=T1,rho_switch=sw,T2_C=T2,first_step_attained=attained,rho1=rho1,G1_nm=G1,f_GBseg_1=f1[0],f_TJ_1=f1[1],f_iso_1=f1[2],rho2=rho2,G2_nm=G2,growth_fraction=growth,f_GBseg_2=f2[0],f_TJ_2=f2[1],f_iso_2=f2[2],first_budget_h=BUDGET/3600,second_budget_h=BUDGET/3600,**extra))
    return rows


def parameter_design():
    levels=(.5,1.,2.,3.);eta=(0.,.1,.25,.5);Qcap=(285e3,325e3,365e3,405e3);drag=(1.,1.5,2.,3.);rows=[]
    for i in range(64):
        wc=levels[i%4];wi=levels[(i//4)%4];wg=levels[(i//16)%4];e=eta[(i*3)%4];q=Qcap[(i//2)%4];td=drag[(i//8)%4]
        p=base_action('action_evolving_capture',w_TJ_capture=.8*wc,w_TJ_iso=.7*wi,w_clean_GB_growth=wg,location=replace(fixed.LocationParams(placement.base_params(),'evolving'),eta_TJ_dens=e),Q_capture_J_mol=q,high_TJ_drag_factor=td)
        rows.append((f'action_{i:02d}',p,dict(w_capture_scale=wc,w_iso_scale=wi,w_clean_growth=wg,eta_TJ_dens=e,Q_capture_J_mol=q,TJ_drag_factor=td)))
    return rows


def screen():
    registry=[];summary=[];rejected=[];lookup={}
    for i,(pid,p,meta) in enumerate(parameter_design(),1):
        print(f'screen {i}/64 {pid}',flush=True);lookup[pid]=p;tra=[]
        for G in G_RED:
            for sw in SW_RED:tra.extend(run_group(pid,p,G,1300.,sw,T2_RED))
        classes=[classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],tol) for r in tra for tol in TOLS];cats=set(classes);reason=''
        n5=sum(classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],.05)=='SUCCESS' for r in tra);n10=sum(classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],.10)=='SUCCESS' for r in tra)
        if n5+n10==0:reason='no_success_window'
        elif classes and all(c=='SUCCESS' for c in classes):reason='universal_success'
        elif 'DENSIFICATION_EXHAUSTION_FAILURE' not in cats:reason='missing_lower_boundary'
        elif 'GRAIN_GROWTH_FAILURE' not in cats:reason='missing_upper_boundary'
        row={'parameter_id':pid,**meta};registry.append(row);s={'parameter_id':pid,'success_5pct':n5,'success_10pct':n10,'has_lower_failure':'DENSIFICATION_EXHAUSTION_FAILURE' in cats,'has_upper_failure':'GRAIN_GROWTH_FAILURE' in cats,'rejected':bool(reason),'rejection_reason':reason};summary.append(s)
        if reason:rejected.append({**row,**s})
    valid=[s for s in summary if not s['rejected']];valid.sort(key=lambda x:(x['success_5pct'],x['success_10pct']),reverse=True);winner=valid[0]['parameter_id'] if valid else None
    return registry,summary,rejected,lookup,winner


def ablations_and_histories():
    summary=[];hist=[]
    for mode in MODES:
        p=base_action(mode);h=action.run(p,aggregate.Iso(1300.,BUDGET),stop_at_rho=TARGET);f=fractions(h);summary.append(dict(action_mode=mode,reached_target=float(h['rho'][-1])>=TARGET-1e-12,rho_final=float(h['rho'][-1]),G_final_nm=float(h['G'][-1])*1e9,time_h=float(h['t'][-1])/3600,f_GBseg=f[0],f_TJ=f[1],f_iso=f[2],E_G_median=float(np.median(h['E_G']))))
    for label,protocol in {'slow_0p2':aggregate.RampHold(.2,1350.,8*3600),'fast_20':aggregate.RampHold(20.,1350.,8*3600)}.items():
        for mode in ('action_evolving_no_capture','action_evolving_capture','action_evolving_capture_high_TJ_drag'):
            h=action.run(base_action(mode),protocol);stride=max(1,len(h['rho'])//300)
            for j in list(range(0,len(h['rho']),stride))+[len(h['rho'])-1]:
                f=fractions(h,j);hist.append(dict(protocol=label,action_mode=mode,t_s=float(h['t'][j]),T_C=float(h['T_C'][j]),rho=float(h['rho'][j]),G_nm=float(h['G'][j])*1e9,f_GBseg=f[0],f_TJ=f[1],f_iso=f[2],**{k:float(h[k][j]) for k in ACTION_DIAG+FLUX_DIAG+POWERS+STATE_DIAG if k in h}))
    return summary,hist


def full_maps(winner,lookup):
    cases={mode:base_action(mode) for mode in MODES}
    if winner:cases[f'screen_winner_{winner}']=lookup[winner]
    trajectories=[]
    for cid,p in cases.items():
        print('full',cid,flush=True)
        for G in G_FULL:
            for T1 in T1_FULL:
                for sw in SW_FULL:trajectories.extend(run_group(cid,p,G,T1,sw,T2_FULL))
    classes=[]
    for r in trajectories:
        for tol in TOLS:classes.append({**r,'growth_tolerance':tol,'rho_target':TARGET,'classification':classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],tol)})
    bounds=[]
    for cid in cases:
        for G in G_FULL:
            for T1 in T1_FULL:
                for sw in SW_FULL:
                    q=[r for r in trajectories if r['case_id']==cid and r['G0_nm']==G and r['T1_C']==T1 and r['rho_switch']==sw]
                    for tol in TOLS:
                        ok=[r for r in q if classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],tol)=='SUCCESS'];dense=[r for r in q if r['first_step_attained'] and r['rho2']>=TARGET-1e-12];ng=[r for r in q if r['first_step_attained'] and r['growth_fraction']<=tol+1e-12]
                        bounds.append(dict(case_id=cid,G0_nm=G,T1_C=T1,rho_switch=sw,growth_tolerance=tol,T_lower_density_C=min((r['T2_C'] for r in dense),default=math.nan),T_upper_no_growth_C=max((r['T2_C'] for r in ng),default=math.nan),T_success_lower_C=min((r['T2_C'] for r in ok),default=math.nan),T_success_upper_C=max((r['T2_C'] for r in ok),default=math.nan),window_width_C=max((r['T2_C'] for r in ok),default=math.nan)-min((r['T2_C'] for r in ok),default=math.nan) if ok else math.nan,n_success=len(ok)))
    return classes,bounds


def plots(out,hist,maps,bounds):
    cap=[r for r in hist if r['protocol']=='slow_0p2' and r['action_mode']=='action_evolving_capture'];fig,axes=plt.subplots(2,3,figsize=(16,9))
    for n in action.ACTION_NAMES:axes[0,0].plot([r['rho'] for r in cap],[r.get(f'action_weight_{n}',math.nan) for r in cap],label=n)
    for n in ('GBseg_remove','TJ_remove','GB_smooth','GB_to_TJ','TJ_to_GBseg_capture','TJ_to_iso'):axes[0,1].plot([r['rho'] for r in cap],[r.get(f'action_flux_{n}',math.nan) for r in cap],label=n)
    for k in ('f_GBseg','f_TJ','f_iso'):axes[0,2].plot([r['rho'] for r in cap],[r[k] for r in cap],label=k)
    axes[1,0].plot([r['rho'] for r in cap],[r.get('action_flux_TJ_to_GBseg_capture',0) for r in cap],label='capture');axes[1,0].plot([r['rho'] for r in cap],[r.get('action_flux_TJ_to_iso',0) for r in cap],label='isolation')
    for k in POWERS:axes[1,1].plot([r['rho'] for r in cap],[r.get(k,math.nan) for r in cap],label=k)
    for k in ('sigma_GBseg_pore','sigma_TJ_pore','sigma_clean_GB','sigma_iso','sigma_act_total'):axes[1,2].plot([r['rho'] for r in cap],[r.get(k,math.nan) for r in cap],label=k)
    for ax in axes.flat:ax.grid(alpha=.2);ax.legend(fontsize=5);ax.set_xlabel('density')
    fig.tight_layout();fig.savefig(out/'action_history_diagnostics.png',dpi=150);plt.close(fig)
    colors={'SUCCESS':'tab:green','GRAIN_GROWTH_FAILURE':'tab:red','DENSIFICATION_EXHAUSTION_FAILURE':'tab:blue','MIXED_FAILURE':'tab:purple','UNATTAINABLE_FIRST_STEP':'.5'};cases=sorted(set(r['case_id'] for r in maps));fig,axs=plt.subplots(2,math.ceil(len(cases)/2),figsize=(5*math.ceil(len(cases)/2),10),sharey=True);axs=np.asarray(axs).ravel()
    for ax,cid in zip(axs,cases):
        q=[r for r in maps if r['case_id']==cid and r['growth_tolerance']==.10]
        for cat,col in colors.items():
            x=[r for r in q if r['classification']==cat];ax.scatter([r['G0_nm'] for r in x],[r['T2_C'] for r in x],c=col,s=16,alpha=.45,label=cat)
        ax.set(xscale='log',title=cid,xlabel='G0 [nm]');ax.grid(alpha=.2)
    axs[0].set_ylabel('T2 [C]');axs[min(len(cases)-1,len(axs)-1)].legend(fontsize=6);fig.tight_layout();fig.savefig(out/'chen_style_action_modes.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));
    for cid in cases:
        q=[r for r in bounds if r['case_id']==cid and r['growth_tolerance']==.10 and np.isfinite(r['window_width_C'])];ax.scatter([r['G0_nm'] for r in q],[r['window_width_C'] for r in q],label=cid,alpha=.5)
    ax.set(xlabel='G0 [nm]',ylabel='10% window width [C]');ax.grid(alpha=.2);ax.legend(fontsize=6);fig.tight_layout();fig.savefig(out/'window_width_vs_G0_by_action_mode.png',dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--outdir',default='results/pore_location_agentic');parser.add_argument('--refresh-screen-rejections',action='store_true');args=parser.parse_args();out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    if args.refresh_screen_rejections:
        with (out/'action_parameter_screen.csv').open(newline='') as f:registry=list(csv.DictReader(f))
        with (out/'screen_summary.csv').open(newline='') as f:summary=list(csv.DictReader(f))
        rejected=[]
        for row,score in zip(registry,summary):
            if int(score['success_5pct'])+int(score['success_10pct'])==0:
                score.update(rejected=True,rejection_reason='no_success_window');rejected.append({**row,**score})
        write_table(out/'screen_summary.csv',summary);write_table(out/'rejected_action_parameter_sets.csv',rejected);print('REFRESHED',len(rejected),'rejections');return
    start=time.perf_counter();abl,hist=ablations_and_histories();registry,screen_summary,rejected,lookup,winner=screen();maps,bounds=full_maps(winner,lookup)
    write_table(out/'action_ablation_summary.csv',abl);write_table(out/'action_flux_histories.csv',hist);write_table(out/'action_parameter_screen.csv',registry);write_table(out/'selected_case_maps.csv',[r for r in maps if r['case_id'].startswith('screen_winner')]);write_table(out/'chen_style_action_map.csv',maps);write_table(out/'window_boundaries_by_action_mode.csv',bounds);write_table(out/'rejected_action_parameter_sets.csv',rejected);write_table(out/'screen_summary.csv',screen_summary);write_table(out/'runtime_summary.csv',[{'wall_s':time.perf_counter()-start,'screen_cases':64,'map_classifications':len(maps),'winner':winner}]);plots(out,hist,maps,bounds);print('DONE',len(maps),'classifications')


if __name__=='__main__':main()
