#!/usr/bin/env python3
"""Staged pore-placement ablations, bounded screen, and selected Chen maps."""
from __future__ import annotations
import argparse, math, time
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import topology_constrained_sintering as agg
import pore_location_topology_model as loc
from density_window_processing_map import write_csv
from two_step_window_map import classify

TARGET=.90;TOLS=(.05,.10);BUDGET=96*3600
G_FULL=(50.,75.,100.,150.,225.,300.,450.,600.)
T1_FULL=(1250.,1300.,1350.);SW_FULL=(.75,.80,.85);T2_FULL=tuple(float(x) for x in range(900,1301,25))
LOCATION_CLASSES={
    "GBseg_rich":(.80,.15,.05,"separated"),"TJ_rich":(.25,.70,.05,"TJ_large"),
    "isolated_rich":(.25,.15,.60,"separated"),"mixed_GBseg_TJ":(.50,.45,.05,"neutral"),
    "clean_GB_rich":(.20,.20,.60,"TJ_large"),
}


def base_params()->agg.Params:
    return agg.Params(rho0=.70,G0=150e-9,pore_radius0=25e-9,pore_ln_sigma=.65,memory_model="pore_bin_redistribution",smoothing_gate_mode="density",growth_mode="baseline")


def make_params(mode="evolving",fractions=(.60,.30,.10),bias="neutral",**kw)->loc.LocationParams:
    return loc.LocationParams(base_params(),mode,*fractions,location_size_bias=bias,**kw)


def fractions_at(h,i=-1):
    gb=float(np.sum(h['phi_GBseg'][i]));tj=float(np.sum(h['phi_TJ'][i]));iso=float(np.sum(h['phi_iso'][i]));total=max(gb+tj+iso,1e-300)
    return gb/total,tj/total,iso/total


def static_ablation()->list[dict]:
    rows=[]
    for name,(fg,ft,fi,bias) in LOCATION_CLASSES.items():
        p=make_params("static",(fg,ft,fi),bias);h=loc.run(p,agg.Iso(1300.,BUDGET),stop_at_rho=TARGET);i=-1;f=fractions_at(h,i)
        rows.append(dict(location_class=name,mode="static",f_GBseg_init=fg,f_TJ_init=ft,f_iso_init=fi,size_bias=bias,
            reached_target=float(h['rho'][i])>=TARGET-1e-12,rho_final=float(h['rho'][i]),G_final_nm=float(h['G'][i])*1e9,
            time_h=float(h['t'][i])/3600,C_GBseg=float(h['C_GBseg'][i]),C_TJ=float(h['C_TJ'][i]),f_clean_GB=float(h['f_clean_GB'][i]),
            f_GBseg_final=f[0],f_TJ_final=f[1],f_iso_final=f[2],E_G_median=float(np.median(h['E_G'])),sigma_act_MPa=float(h['sigma_act_total'][i])/1e6,
            P_GBseg_drag=float(h['P_GBseg_drag'][i]),P_TJ_drag=float(h['P_TJ_drag'][i])))
    return rows


def history_rows()->list[dict]:
    rows=[]
    protocols={"slow_0p2":agg.RampHold(.2,1350.,8*3600),"fast_20":agg.RampHold(20.,1350.,8*3600)}
    for mode in ("static","evolving"):
        p=make_params(mode)
        for label,protocol in protocols.items():
            h=loc.run(p,protocol)
            stride=max(1,len(h['t'])//350)
            for i in list(range(0,len(h['t']),stride))+[len(h['t'])-1]:
                fg,ft,fi=fractions_at(h,i)
                rows.append(dict(mode=mode,protocol=label,t_s=float(h['t'][i]),T_C=float(h['T_C'][i]),rho=float(h['rho'][i]),G_nm=float(h['G'][i])*1e9,
                    f_GBseg=fg,f_TJ=ft,f_iso=fi,C_GBseg=float(h['C_GBseg'][i]),C_TJ=float(h['C_TJ'][i]),f_clean_GB=float(h['f_clean_GB'][i]),
                    GBseg_mean_nm=float(h['GBseg_pore_mean_radius'][i])*1e9,TJ_mean_nm=float(h['TJ_pore_mean_radius'][i])*1e9,iso_mean_nm=float(h['iso_pore_mean_radius'][i])*1e9,
                    E_G=float(h['E_G'][i]),**{k:float(h[k][i]) for k in h if k.startswith('P_') or k.startswith('sigma_')}))
    return rows


def screen_design()->list[tuple[str,loc.LocationParams]]:
    fractions=((.8,.1,.1),(.5,.3,.2),(.2,.6,.2),(.5,.1,.4));biases=("neutral","separated","TJ_large")
    scales=(.35,1.,2.5);etas=(0.,.25,.5);qs=(1.,2.,3.);pin_scales=(.4,1.,2.);stress_scales=(.5,1.,1.8)
    out=[]
    for i in range(64):
        f=fractions[i%4];bias=biases[(i//4)%3];smooth=scales[(i//12)%3];reloc=scales[(i*2+1)%3];isol=scales[(i*3+2)%3]
        eta=etas[(i//3)%3];q=qs[(i//7)%3];pin=pin_scales[(i//5)%3];stress=stress_scales[(i//11)%3]
        p=make_params("evolving",f,bias,eta_TJ_dens=eta,q_GBseg=q,q_TJ=q,
            k_GB_smooth_0_s=1e5*smooth,k_GB_to_TJ_0_s=4e4*reloc,k_TJ_iso_0_s=2e5*isol,k_TJ_smooth_0_s=5e4*smooth,
            A_GBseg_pin=12*pin,A_TJ_pin=28*pin,K_GBseg=.45*stress,K_TJ=1.1*stress)
        out.append((f"screen_{i:02d}",p))
    return out


def run_second_steps(case_id,p,G0_nm,T1,switch,T2_values):
    p=replace(p,base=replace(p.base,G0=G0_nm*1e-9));first=loc.run(p,agg.Iso(T1,BUDGET),stop_at_rho=switch);attained=float(first['rho'][-1])>=switch-1e-12
    common=dict(case_id=case_id,G0_nm=G0_nm,T1_C=T1,rho_switch=switch,first_step_attained=attained,rho1=float(first['rho'][-1]),G1_nm=float(first['G'][-1])*1e9)
    if attained:state=loc.final_state(first,p);fg,ft,fi=fractions_at(first)
    else:state=None;fg=ft=fi=math.nan
    rows=[]
    for T2 in T2_values:
        if attained:
            h=loc.run(p,agg.Iso(T2,BUDGET),initial=state);rho2=float(h['rho'][-1]);G2=float(h['G'][-1])*1e9;growth=(G2-common['G1_nm'])/common['G1_nm'];f2=fractions_at(h)
            extra=dict(rho2=rho2,G2_nm=G2,growth_fraction=growth,f_GBseg_1=fg,f_TJ_1=ft,f_iso_1=fi,f_GBseg_2=f2[0],f_TJ_2=f2[1],f_iso_2=f2[2],
                C_GBseg_1=float(first['C_GBseg'][-1]),C_TJ_1=float(first['C_TJ'][-1]),C_GBseg_2=float(h['C_GBseg'][-1]),C_TJ_2=float(h['C_TJ'][-1]),
                E_G_median=float(np.median(h['E_G'])),mobility_final=float(h['growth_mobility_factor'][-1]))
        else:extra={k:math.nan for k in ('rho2','G2_nm','growth_fraction','f_GBseg_1','f_TJ_1','f_iso_1','f_GBseg_2','f_TJ_2','f_iso_2','C_GBseg_1','C_TJ_1','C_GBseg_2','C_TJ_2','E_G_median','mobility_final')}
        rows.append({**common,'T2_C':T2,**extra})
    return rows


def reduced_screen()->tuple[list[dict],list[dict],dict[str,loc.LocationParams]]:
    summaries=[];rejected=[];lookup={}
    for case_id,p in screen_design():
        lookup[case_id]=p;tra=[]
        for G in (100.,225.):tra.extend(run_second_steps(case_id,p,G,1300.,.80,(1050.,1150.,1250.,1300.)))
        classes=[]
        for r in tra:
            for tol in TOLS:classes.append(classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],TARGET,tol))
        cats=set(classes);success5=classes.count('SUCCESS') if False else sum(classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],TARGET,.05)=='SUCCESS' for r in tra)
        success10=sum(classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],TARGET,.10)=='SUCCESS' for r in tra)
        lower='DENSIFICATION_EXHAUSTION_FAILURE' in cats;upper='GRAIN_GROWTH_FAILURE' in cats
        universal=all(c=='SUCCESS' for c in classes);reason=''
        if universal:reason='universal_success'
        elif not lower:reason='missing_lower_boundary'
        elif not upper:reason='missing_upper_boundary'
        summaries.append(dict(case_id=case_id,success_5pct=success5,success_10pct=success10,has_lower_failure=lower,has_upper_failure=upper,universal_success=universal,rejected=bool(reason),rejection_reason=reason,
            f_GBseg_init=p.f_GBseg_init,f_TJ_init=p.f_TJ_init,f_iso_init=p.f_iso_init,size_bias=p.location_size_bias,eta_TJ_dens=p.eta_TJ_dens,q=p.q_GBseg,
            k_GB_smooth_scale=p.k_GB_smooth_0_s/1e5,k_GB_to_TJ_scale=p.k_GB_to_TJ_0_s/4e4,k_TJ_iso_scale=p.k_TJ_iso_0_s/2e5,
            pin_scale=p.A_GBseg_pin/12,stress_scale=p.K_GBseg/.45))
        if reason:rejected.append(summaries[-1])
    return summaries,rejected,lookup


def selected_cases(summary,lookup):
    valid=[r for r in summary if not r['rejected']];valid.sort(key=lambda r:(r['success_5pct'],r['success_10pct']),reverse=True)
    chosen=[]
    for r in valid:
        if r['case_id'] not in chosen:chosen.append(r['case_id'])
        if len(chosen)==2:break
    cases={"default_evolving":make_params("evolving"),"static_GBseg_rich":make_params("static",(.8,.15,.05),"separated")}
    cases.update({c:lookup[c] for c in chosen});return cases


def full_maps(cases):
    trajectories=[]
    for case_id,p in cases.items():
        print(f"full map {case_id}",flush=True)
        for G in G_FULL:
            for T1 in T1_FULL:
                for sw in SW_FULL:trajectories.extend(run_second_steps(case_id,p,G,T1,sw,T2_FULL))
    mapped=[]
    for r in trajectories:
        for tol in TOLS:mapped.append({**r,'rho_target':TARGET,'growth_tolerance':tol,'classification':classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],TARGET,tol)})
    boundaries=[]
    for case_id in cases:
        for G in G_FULL:
            for T1 in T1_FULL:
                for sw in SW_FULL:
                    group=[r for r in trajectories if r['case_id']==case_id and r['G0_nm']==G and r['T1_C']==T1 and r['rho_switch']==sw]
                    for tol in TOLS:
                        success=[r for r in group if classify(r['first_step_attained'],r['rho2'],r['growth_fraction'],TARGET,tol)=='SUCCESS']
                        dense=[r for r in group if r['first_step_attained'] and r['rho2']>=TARGET-1e-12];nogrow=[r for r in group if r['first_step_attained'] and r['growth_fraction']<=tol+1e-12]
                        boundaries.append(dict(case_id=case_id,G0_nm=G,T1_C=T1,rho_switch=sw,growth_tolerance=tol,
                            T_lower_density_C=min((r['T2_C'] for r in dense),default=math.nan),T_upper_no_growth_C=max((r['T2_C'] for r in nogrow),default=math.nan),
                            T_success_lower_C=min((r['T2_C'] for r in success),default=math.nan),T_success_upper_C=max((r['T2_C'] for r in success),default=math.nan),
                            window_width_C=max((r['T2_C'] for r in success),default=math.nan)-min((r['T2_C'] for r in success),default=math.nan) if success else math.nan,n_success=len(success)))
    return mapped,boundaries


def plots(out,ablation,hist,maprows,bounds):
    evo=[r for r in hist if r['mode']=='evolving'];fig,axes=plt.subplots(1,2,figsize=(12,5))
    for protocol in ('slow_0p2','fast_20'):
        x=[r for r in evo if r['protocol']==protocol]
        for key,label in (('f_GBseg','GBseg'),('f_TJ','TJ'),('f_iso','isolated')):axes[0].plot([r['rho'] for r in x],[r[key] for r in x],label=f'{protocol} {label}')
        axes[1].plot([r['T_C'] for r in x],[r['f_GBseg'] for r in x],label=f'{protocol} GBseg');axes[1].plot([r['T_C'] for r in x],[r['f_TJ'] for r in x],'--',label=f'{protocol} TJ')
    axes[0].set(xlabel='density',ylabel='pore-location volume fraction');axes[1].set(xlabel='temperature [C]',ylabel='location fraction');[a.grid(alpha=.2) for a in axes];axes[1].legend(fontsize=7);fig.tight_layout();fig.savefig(out/'pore_location_fractions.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for name,(fg,ft,fi,bias) in LOCATION_CLASSES.items():
        h=loc.run(make_params('static',(fg,ft,fi),bias),agg.Iso(1300.,BUDGET),stop_at_rho=TARGET)
        axes[0].plot(h['rho'],h['G']*1e9,label=name);axes[1].plot(h['rho'],h['E_G'],label=name)
    axes[0].set(xlabel='density',ylabel='grain size [nm]');axes[1].set(xlabel='density',ylabel='E_G');axes[0].legend(fontsize=7);[a.grid(alpha=.2) for a in axes];fig.tight_layout();fig.savefig(out/'grain_size_and_EG_by_location.png',dpi=150);plt.close(fig)
    for tol in TOLS:
        cases=sorted(set(r['case_id'] for r in maprows));fig,axes=plt.subplots(1,len(cases),figsize=(5*len(cases),5),sharey=True,squeeze=False)
        color={'SUCCESS':'tab:green','GRAIN_GROWTH_FAILURE':'tab:red','DENSIFICATION_EXHAUSTION_FAILURE':'tab:blue','MIXED_FAILURE':'tab:purple','UNATTAINABLE_FIRST_STEP':'.5'}
        for ax,case in zip(axes[0],cases):
            group=[r for r in maprows if r['case_id']==case and r['growth_tolerance']==tol]
            for cat,c in color.items():
                q=[r for r in group if r['classification']==cat];ax.scatter([r['G0_nm'] for r in q],[r['T2_C'] for r in q],s=18,alpha=.45,c=c,label=cat)
            ax.set(xscale='log',xlabel='G0 [nm]',title=case);ax.grid(alpha=.2)
        axes[0,0].set_ylabel('T2 [C]');axes[0,-1].legend(fontsize=6);fig.tight_layout();fig.savefig(out/f'chen_style_pore_location_{int(tol*100)}pct.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5.5))
    for case in sorted(set(r['case_id'] for r in bounds)):
        q=[r for r in bounds if r['case_id']==case and r['growth_tolerance']==.10 and np.isfinite(r['window_width_C'])];ax.scatter([r['G0_nm'] for r in q],[r['window_width_C'] for r in q],label=case,alpha=.55)
    ax.set(xscale='log',xlabel='G0 [nm]',ylabel='10% window width [C]');ax.grid(alpha=.2);ax.legend(fontsize=7);fig.tight_layout();fig.savefig(out/'window_width_vs_G0_by_location.png',dpi=150);plt.close(fig)
    power=[r for r in evo if r['protocol']=='slow_0p2'];fig,axes=plt.subplots(1,2,figsize=(13,5))
    pkeys=('P_GBseg_dens','P_TJ_dens','P_clean_GB','P_GBseg_drag','P_TJ_drag','P_GB_to_TJ_relocation','P_TJ_iso_conversion')
    ptotal=np.maximum(np.sum([[max(r[k],0) for k in pkeys] for r in power],axis=1),1e-300)
    for k in pkeys:axes[0].plot([r['rho'] for r in power],[max(r[k],0)/total for r,total in zip(power,ptotal)],label=k)
    for k in ('sigma_GBseg_pore','sigma_TJ_pore','sigma_clean_GB','sigma_iso','sigma_act_total'):axes[1].plot([r['rho'] for r in power],[r[k]/1e6 for r in power],label=k)
    axes[0].set(xlabel='density',ylabel='fraction of named power');axes[1].set(xlabel='density',ylabel='stress [MPa]');[a.grid(alpha=.2) for a in axes];axes[0].legend(fontsize=6);axes[1].legend(fontsize=7);fig.tight_layout();fig.savefig(out/'power_and_stress_partitions.png',dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--outdir',default='results/pore_location_topology');args=parser.parse_args();out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    ablation=static_ablation();hist=history_rows();summary,rejected,lookup=reduced_screen();cases=selected_cases(summary,lookup);mapped,bounds=full_maps(cases)
    write_csv(out/'pore_location_ablation.csv',ablation);write_csv(out/'location_evolution_histories.csv',hist)
    write_csv(out/'pore_location_power_channels.csv',[{k:r[k] for k in ('mode','protocol','t_s','T_C','rho')+tuple(x for x in r if x.startswith('P_'))} for r in hist])
    write_csv(out/'pore_location_stress_channels.csv',[{k:r[k] for k in ('mode','protocol','t_s','T_C','rho')+tuple(x for x in r if x.startswith('sigma_'))} for r in hist])
    write_csv(out/'chen_style_pore_location_map.csv',mapped);write_csv(out/'window_boundaries_by_location.csv',bounds);write_csv(out/'parameter_screen_summary.csv',summary);write_csv(out/'rejected_parameter_sets.csv',rejected)
    write_csv(out/'selected_cases.csv',[{'case_id':k,**{x:y for x,y in asdict(v).items() if x!='base'}} for k,v in cases.items()]);write_csv(out/'runtime_summary.csv',[{'wall_s':time.perf_counter()-start,'screen_cases':len(summary),'selected_cases':len(cases),'map_classifications':len(mapped)}])
    plots(out,ablation,hist,mapped,bounds);print('DONE',len(mapped),'classifications',flush=True)


if __name__=='__main__':main()
