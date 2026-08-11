#!/usr/bin/env python3
"""Frozen-mechanism production assessment for two-step and fast firing."""
from __future__ import annotations
import argparse,csv,math,time
from dataclasses import replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import adaptive_T2_boundary_search as adaptive
import agentic_mechanism_model as model
import agentic_mechanism_search as discovery
import preparation_window_search as prep

G_CHEN=(75.,100.,125.,150.,175.,200.,225.,250.,300.);T1S=(1325.,1350.,1375.,1400.,1425.,1450.,1475.,1500.);RATES=(1.,5.,20.,50.,100.);SWITCHES=(.76,.78,.80,.82,.84,.86,.88,.90)
G_FAST=(50.,75.,100.,150.,225.,300.);RHO0S=(.65,.70,.75);TOPOLOGIES={'baseline':(.60,.30,.10),'GBseg_rich':(.80,.15,.05),'TJ_rich':(.35,.55,.10),'mixed_GBseg_TJ':(.50,.45,.05)};FAST_RATES=(.2,1.,5.,20.,50.,100.);PEAKS=(1300.,1350.,1400.,1450.,1500.);HOLDS=(0.,2.,8.,20.);TARGETS=(.85,.88,.90,.92);BUDGET=96*3600;NUMERICAL_DT_MAX_S=900.


class FastSchedule:
    def __init__(self,rate,peak,hold_h):self.rate=rate/60;self.peak=peak;self.ramp_s=(peak-25)/self.rate;self.t_end=min(self.ramp_s+hold_h*3600,BUDGET)
    def T(self,t,rho):return min(self.peak,25+self.rate*t)


def write(path,rows,empty=('mechanism_id','status')):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    if not fields:fields=list(empty)
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def frozen_mechanisms():
    out={}
    for mid,p in prep.mechanisms().items():
        base=replace(p.action.location.base,dt_max_s=NUMERICAL_DT_MAX_S);out[mid]=replace(p,action=replace(p.action,location=replace(p.action.location,base=base)))
    return out


def chen_states(out):
    rows=[];states=[];rejected=[];n=0;total=len(frozen_mechanisms())*len(G_CHEN)*len(T1S)*len(RATES)
    for mid,p0 in frozen_mechanisms().items():
      for G in G_CHEN:
       lp=replace(p0.action.location,base=replace(p0.action.location.base,G0=G*1e-9));p=replace(p0,action=replace(p0.action,location=lp))
       for T1 in T1S:
        for rate in RATES:
         n+=1
         if n%25==0:print('production-first',n,'/',total,flush=True)
         h=model.run(p,prep.FixedBudgetRamp(rate,T1))
         for sw in SWITCHES:
          idx=np.flatnonzero(h['rho']>=sw-1e-12)
          if not len(idx):row=prep.state_row(mid,p,G,T1,rate,sw,h,None,'switch_density_unattainable');rows.append(row);rejected.append(row);continue
          i=int(idx[0]);row=prep.state_row(mid,p,G,T1,rate,sw,h,i);rows.append(row)
          if row['rho1']>=adaptive.TARGET-1e-12:rejected.append({**row,'reason':'target_already_reached_first_step'});continue
          if row['first_step_growth_fraction']>.20:rejected.append({**row,'reason':'first_step_growth_above_20pct'});continue
          states.append((row,model.final_state(h,p,i),p))
      write(out/'raw_first_step_states.csv',rows);write(out/'raw_first_step_rejections.csv',rejected)
    return rows,states,rejected


def chen_adaptive(out,states):
    # Exact-state reuse is limited to identical instantaneous state vectors.
    groups={};points=[];bounds=[]
    for meta,state,p in states:
        key=(meta['mechanism_id'],round(state.pore.rho,13),round(state.pore.G,18),round(state.X_J,13),state.pore.phi_GBseg.round(15).tobytes(),state.pore.phi_TJ.round(15).tobytes(),state.pore.phi_iso.round(15).tobytes());groups.setdefault(key,[]).append((meta,state,p))
    print('production unique states',len(groups),'from',len(states),flush=True)
    for i,members in enumerate(groups.values(),1):
      if i%10==0:print('production-adaptive',i,'/',len(groups),flush=True)
      meta,state,p=members[0];x,_=adaptive.adaptive_from_state(meta['mechanism_id'],p,meta['G0_nm'],meta['T1_C'],meta['rho_switch'],state);gid=f'prod_{i:05d}'
      rep={k:meta[k] for k in ('mechanism_id','heating_rate_C_min','T_at_switch_C','rho1','G1_nm','first_step_growth_fraction','X_J','Lambda_over_K_TJ','P_comp_TJ','f_GBseg','f_TJ','f_iso_location','C_GBseg','C_TJ','f_clean_GB')};points.extend([{**r,**rep,'state_group_id':gid,'equivalent_routes':len(members)} for r in x])
      for route,_,_ in members:
       tag={k:route[k] for k in ('mechanism_id','heating_rate_C_min','T_at_switch_C','rho1','G1_nm','first_step_growth_fraction','X_J','Lambda_over_K_TJ','P_comp_TJ','f_GBseg','f_TJ','f_iso_location','C_GBseg','C_TJ','f_clean_GB')}
       for tol in discovery.TOLS:
        for practical,kind in ((False,'kinetic'),(True,'practical')):
         s=adaptive.status(x,tol,route['T1_C'],practical);eligible=[r for r in x if not practical or r['T2_C']<route['T1_C']];bounds.append(dict(parameter_id=route['mechanism_id'],mechanism_mode=p.mechanism_mode,G0_nm=route['G0_nm'],T1_C=route['T1_C'],rho_switch=route['rho_switch'],growth_tolerance=tol,map_type=kind,T2_min_searched_C=min(r['T2_C'] for r in eligible),T2_max_searched_C=max(r['T2_C'] for r in eligible),state_group_id=gid,**s,**tag))
      if i%25==0:write(out/'raw_adaptive_boundaries.csv',bounds);write(out/'raw_adaptive_points.csv',points)
    write(out/'raw_adaptive_boundaries.csv',bounds);write(out/'raw_adaptive_points.csv',points);return points,bounds


def tier(r):
    complete=r['map_type']=='practical' and r['boundary_status']=='COMPLETE_WINDOW';width=r.get('window_width_C',math.nan)
    if complete and r['growth_tolerance']==.05 and r['first_step_growth_fraction']<=.05 and width>=25:return 'Tier_A'
    if complete and r['growth_tolerance']==.10 and r['first_step_growth_fraction']<=.10 and width>=25:return 'Tier_B'
    return 'Tier_C'


def fast_params(p0,G,rho0,fractions):
    base=replace(p0.action.location.base,G0=G*1e-9,rho0=rho0);lp=replace(p0.action.location,base=base,f_GBseg_init=fractions[0],f_TJ_init=fractions[1],f_iso_init=fractions[2]);return replace(p0,action=replace(p0.action,location=lp))


def sample_history(mid,p,G,rho0,topo,rate,peak,hold,h):
    rows=[];peak_idx=np.flatnonzero(h['T_C']>=peak-1e-9);pi=int(peak_idx[0]) if len(peak_idx) else len(h['rho'])-1
    for target in TARGETS:
      idx=np.flatnonzero(h['rho']>=target-1e-12);common=dict(mechanism_id=mid,mechanism_mode=p.mechanism_mode,G0_nm=G,rho0=rho0,initial_topology=topo,heating_rate_C_min=rate,peak_T_C=peak,hold_time_h=hold,rho_target=target,budget_h=96,peak_attained=bool(len(peak_idx)),rho_peak=float(h['rho'][pi]),G_peak_nm=float(h['G'][pi])*1e9,C_GBseg_peak=float(h['C_GBseg'][pi]),C_TJ_peak=float(h['C_TJ'][pi]),X_J_peak=float(h['X_J'][pi]) if 'X_J' in h else 0.)
      if not len(idx):rows.append({**common,'target_attained':False,'rho_final':float(h['rho'][-1]),'reason':'target_nonattainment'});continue
      i=int(idx[0]);phi=[float(np.sum(h[k][i])) for k in ('phi_GBseg','phi_TJ','phi_iso')];z=max(sum(phi),1e-300);lnG=math.log(max(float(h['G'][i])/(G*1e-9),1+1e-15));rows.append({**common,'target_attained':True,'rho_final':float(h['rho'][i]),'time_to_target_h':float(h['t'][i])/3600,'G_at_target_nm':float(h['G'][i])*1e9,'densification_per_lnG':(target-rho0)/max(lnG,1e-15),'X_J':float(h['X_J'][i]) if 'X_J' in h else 0.,'Lambda_over_K_TJ':float(h['Lambda_over_K_TJ'][i]) if 'Lambda_over_K_TJ' in h else math.nan,'P_comp_TJ':float(h['P_comp_TJ'][i]) if 'P_comp_TJ' in h else math.nan,'C_GBseg':float(h['C_GBseg'][i]),'C_TJ':float(h['C_TJ'][i]),'f_clean_GB':float(h['f_clean_GB'][i]),'f_iso':float(h['f_iso'][i]),'f_GBseg_location':phi[0]/z,'f_TJ_location':phi[1]/z,'f_iso_location':phi[2]/z,'P_dens':float(h['P_GBseg_dens'][i]+h['P_TJ_dens'][i]),'P_clean_GB':float(h['P_clean_GB_discovery'][i]) if 'P_clean_GB_discovery' in h else float(h['P_clean_GB'][i]),'P_persistent_junction_drag':float(h['P_persistent_junction_drag'][i]) if 'P_persistent_junction_drag' in h else 0.,'P_TJ_multihit':float(h['P_TJ_multihit'][i]) if 'P_TJ_multihit' in h else 0.,'sigma_act_total':float(h['sigma_act_total'][i])})
    return rows


def fast_campaign(out):
    rows=[];n=0;total=len(frozen_mechanisms())*len(G_FAST)*len(RHO0S)*len(TOPOLOGIES)*len(FAST_RATES)*len(PEAKS)*len(HOLDS)
    for mid,p0 in frozen_mechanisms().items():
     for G in G_FAST:
      for rho0 in RHO0S:
       for topo,frac in TOPOLOGIES.items():
        p=fast_params(p0,G,rho0,frac)
        for rate in FAST_RATES:
         for peak in PEAKS:
          for hold in HOLDS:
           n+=1
           if n%100==0:print('fast-firing',n,'/',total,flush=True)
           h=model.run(p,FastSchedule(rate,peak,hold));rows.extend(sample_history(mid,p,G,rho0,topo,rate,peak,hold,h))
           if n%500==0:write(out/'raw_fast_firing_paths.csv',rows)
    write(out/'raw_fast_firing_paths.csv',rows);return rows


def fast_comparisons(rows):
    keycols=('mechanism_id','G0_nm','rho0','initial_topology','peak_T_C','hold_time_h','rho_target');groups={}
    for r in rows:groups.setdefault(tuple(r[k] for k in keycols),[]).append(r)
    out=[]
    for key,q in groups.items():
      by={r['heating_rate_C_min']:r for r in q};ref_rate=.2 if by.get(.2,{}).get('target_attained') else (1. if by.get(1.,{}).get('target_attained') else math.nan);ref=by.get(ref_rate,{})
      for rate in FAST_RATES:
       r=by[rate];base={k:v for k,v in zip(keycols,key)};att=bool(r.get('target_attained'))
       if not att or not math.isfinite(ref_rate):out.append({**base,'heating_rate_C_min':rate,'reference_rate_C_min':ref_rate,'comparison_attained':False,'HR_pct':math.nan,'response_class':'unattainable','reason':'compared_target_nonattainment'});continue
       hr=100*(ref['G_at_target_nm']-r['G_at_target_nm'])/ref['G_at_target_nm'];cls='beneficial' if hr>1 else ('harmful' if hr<-1 else 'neutral');out.append({**base,'heating_rate_C_min':rate,'reference_rate_C_min':ref_rate,'comparison_attained':True,'HR_pct':hr,'response_class':cls,'G_reference_nm':ref['G_at_target_nm'],'G_fast_nm':r['G_at_target_nm'],'time_to_target_h':r['time_to_target_h'],'X_J':r['X_J'],'Lambda_over_K_TJ':r['Lambda_over_K_TJ'],'P_comp_TJ':r['P_comp_TJ'],'C_GBseg':r['C_GBseg'],'C_TJ':r['C_TJ'],'f_clean_GB':r['f_clean_GB'],'f_iso':r['f_iso']})
    return out


def assemble(out,states,rejected,points,bounds,fastrows):
    tagged=[{**r,'tier':tier(r)} for r in bounds];success=[r for r in tagged if r['tier'] in ('Tier_A','Tier_B')];failed=[r for r in tagged if r['boundary_status']!='COMPLETE_WINDOW' or r['tier']=='Tier_C'];comparisons=fast_comparisons(fastrows);fastsuccess=[r for r in comparisons if r['response_class']=='beneficial'];write(out/'successful_practical_windows.csv',success);write(out/'successful_fast_firing_cases.csv',fastsuccess);write(out/'failed_or_censored_cases.csv',failed);write(out/'production_fast_firing_summary.csv',comparisons)
    win=[]
    for mid in frozen_mechanisms():
      q=[r for r in success if r['mechanism_id']==mid];A=[r for r in q if r['tier']=='Tier_A'];B=[r for r in q if r['tier']=='Tier_B'];win.append(dict(mechanism_id=mid,Tier_A_count=len(A),Tier_B_count=len(B),min_G1_Tier_A_nm=min((r['G1_nm'] for r in A),default=math.nan),max_Tier_A_width_C=max((r['window_width_C'] for r in A),default=math.nan),complete_lower_upper=all(r['lower_bracketed'] and r['upper_bracketed'] for r in q)))
    write(out/'production_window_summary.csv',win)
    score=[]
    for w in win:
      mid=w['mechanism_id'];fq=[r for r in comparisons if r['mechanism_id']==mid];pos=[r for r in fq if r['response_class']=='beneficial'];harm=[r for r in fq if r['response_class']=='harmful'];score.append({**w,'positive_fast_firing_count':len(pos),'harmful_fast_firing_count':len(harm),'positive_density_min':min((r['rho_target'] for r in pos),default=math.nan),'positive_density_max':max((r['rho_target'] for r in pos),default=math.nan),'positive_topologies':'|'.join(sorted(set(r['initial_topology'] for r in pos))),'joint_positive':bool(w['Tier_A_count'] and pos)})
    write(out/'joint_mechanism_scorecard.csv',score)
    ingredients=[('geometric_eligibility','required','Only pore-connected GB/TJ stores densify'),('pore_size_location_memory','required','Heating changes observable pore bins and placement'),('separate_densification_migration','required','Migration suppression does not alter shared-state rho_dot'),('persistent_X_J','required_for_current_window','Drag persists after connected coverage falls'),('TJ_multihit','required_for_upper_boundary','Temperature-sensitive reactivation restores upper growth failure'),('adaptive_preparation','required','Finite T1 tradeoff prepares topology'),('censor_aware_boundaries','required','Both bounds and T2<T1 enforced'),('density_gate_only','insufficient','Moves crossover without observable topology'),('connected_pinning_only','insufficient','Releases as coverage falls')];write(out/'mechanism_ingredient_table.csv',[dict(ingredient=a,status=b,rationale=c) for a,b,c in ingredients])
    reps=[]
    for mid in frozen_mechanisms():
      a=next((r for r in success if r['mechanism_id']==mid and r['tier']=='Tier_A'),None);f=next((r for r in comparisons if r['mechanism_id']==mid and r['response_class']=='beneficial'),None)
      if a:reps.append({'mechanism_id':mid,'path_type':'two_step_success',**a})
      if f:reps.append({'mechanism_id':mid,'path_type':'fast_firing_success',**f})
    write(out/'representative_path_diagnostics.csv',reps);return tagged,comparisons,score,reps


def figures(out,tagged,comparisons,score,reps):
    figdir=out/'figures';figdir.mkdir(exist_ok=True);colors={'Tier_A':'tab:green','Tier_B':'tab:blue','Tier_C':'0.7'};q=[r for r in tagged if r['map_type']=='practical']
    fig,ax=plt.subplots(figsize=(8,6));
    for t,c in colors.items():z=[r for r in q if r['tier']==t];ax.scatter([r['G1_nm'] for r in z],[r.get('T_first_success_C',math.nan) for r in z],c=c,s=8,alpha=.35,label=t)
    ax.set(xlabel='G1 after preparation [nm]',ylabel='T2 [C]',title='Practical Chen-style windows');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(figdir/'figure1_practical_chen_map.png',dpi=200);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,6));kin=[r for r in tagged if r['map_type']=='kinetic' and r['boundary_status']=='COMPLETE_WINDOW'];pra=[r for r in q if r['boundary_status']=='COMPLETE_WINDOW'];ax.scatter([r['G1_nm'] for r in kin],[r['window_width_C'] for r in kin],s=8,label='kinetic',alpha=.3);ax.scatter([r['G1_nm'] for r in pra],[r['window_width_C'] for r in pra],s=8,label='practical',alpha=.3);ax.set(xlabel='G1 [nm]',ylabel='window width [C]');ax.legend();fig.tight_layout();fig.savefig(figdir/'figure2_kinetic_vs_practical.png',dpi=200);plt.close(fig)
    A=[r for r in tagged if r['tier']=='Tier_A'];fig,ax=plt.subplots();sc=ax.scatter([r['T1_C'] for r in A],[r['rho_switch'] for r in A],c=[r['window_width_C'] for r in A],cmap='viridis');fig.colorbar(sc,ax=ax,label='window width [C]');ax.set(xlabel='T1 [C]',ylabel='switch density');fig.tight_layout();fig.savefig(figdir/'figure3_preparation_tradeoff.png',dpi=200);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(12,9),sharex=True,sharey=True)
    for ax,topo in zip(axs.ravel(),TOPOLOGIES):
      z=[r for r in comparisons if r['rho_target']==.90 and r['initial_topology']==topo and r['comparison_attained']];sc=ax.scatter([r['heating_rate_C_min'] for r in z],[r['peak_T_C'] for r in z],c=[r['HR_pct'] for r in z],cmap='coolwarm',vmin=-20,vmax=20,s=10);ax.set(xscale='log',title=topo,xlabel='heating rate [C/min]',ylabel='peak T [C]')
    fig.colorbar(sc,ax=axs.ravel().tolist(),label='HR_pct [%]');fig.savefig(figdir/'figure4_fast_firing_map.png',dpi=200,bbox_inches='tight');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));labels=[r['mechanism_id'] for r in score];ax.bar(np.arange(len(labels))-.18,[r['Tier_A_count'] for r in score],.36,label='Tier A');ax.bar(np.arange(len(labels))+.18,[r['positive_fast_firing_count'] for r in score],.36,label='positive HR');ax.set_xticks(range(len(labels)),labels,rotation=20);ax.legend();fig.tight_layout();fig.savefig(figdir/'figure5_mechanism_state_diagnostics.png',dpi=200);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));ax.axis('off');ax.text(.02,.9,'Power / stress channels retained in representative_path_diagnostics.csv',fontsize=13);ax.text(.02,.65,'P_dens | P_clean_GB | P_persistent_junction_drag | P_TJ_multihit',fontsize=11);ax.text(.02,.4,'sigma_base | sigma_GBseg | sigma_TJ | sigma_act_total',fontsize=11);fig.tight_layout();fig.savefig(figdir/'figure6_power_stress_partition.png',dpi=200);plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,3));ax.axis('off');steps=['Pore-size / location memory','Renewal densification','Persistent junction X_J','TJ multihit reactivation','Finite practical window'];
    for i,s in enumerate(steps):ax.text(i*.2+.01,.5,s,ha='center',va='center',bbox=dict(boxstyle='round',fc='white',ec='black'))
    for i in range(4):ax.annotate('',xy=(i*.2+.17,.5),xytext=(i*.2+.05,.5),arrowprops=dict(arrowstyle='->'))
    fig.tight_layout();fig.savefig(figdir/'figure7_mechanism_ingredients.png',dpi=200);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',default='results/production_mechanism_assessment');a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();states,admissible,rejected=chen_states(out);points,bounds=chen_adaptive(out,admissible);fastrows=fast_campaign(out);tagged,comparisons,score,reps=assemble(out,states,rejected,points,bounds,fastrows);figures(out,tagged,comparisons,score,reps);write(out/'runtime_summary.csv',[{'wall_s':time.perf_counter()-start,'rho_target':.90,'budget_h':96,'numerical_dt_max_s':NUMERICAL_DT_MAX_S,'frozen_mechanisms':'|'.join(frozen_mechanisms())}]);print('DONE production',flush=True)


if __name__=='__main__':main()
