#!/usr/bin/env python3
"""First-step preparation search using the adaptive Chen-window engine."""
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
import pore_location_topology_model as location

T1S=(1250.,1300.,1350.,1400.,1450.,1500.);RATES=(.2,1.,5.,20.,100.);SWITCHES=(.72,.75,.80,.85,.88,.90);G0S=(50.,75.,100.,150.,225.,300.);PREP_TOLS=(.05,.10,.20)


class FixedBudgetRamp:
    def __init__(self,rate,target,start=25.,budget=adaptive.BUDGET):self.rate=rate/60;self.target=target;self.start=start;self.t_end=budget
    def T(self,t,rho):return min(self.target,self.start+self.rate*t)


def write(path,rows,empty=('mechanism_id','reason')):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    if not fields:fields=list(empty)
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def mechanisms():
    base=adaptive.survivor_params();return {'mech_009':base['mech_009'],'mech_019':base['mech_019'],'mech_009_q0':replace(base['mech_009'],mechanism_mode='persistent_tj_multihit_q0',q_TJ=0),'mech_019_q0':replace(base['mech_019'],mechanism_mode='persistent_tj_multihit_q0',q_TJ=0)}


def state_row(mid,p,G,T1,rate,sw,h,index=None,reason=''):
    common=dict(mechanism_id=mid,mechanism_mode=p.mechanism_mode,G0_nm=G,T1_C=T1,heating_rate_C_min=rate,rho_switch=sw,first_step_budget_h=96)
    if index is None:return {**common,'first_step_attained':False,'reason':reason}
    phi=[float(np.sum(h[k][index])) for k in ('phi_GBseg','phi_TJ','phi_iso')];z=max(sum(phi),1e-300);rho=float(h['rho'][index]);G1=float(h['G'][index])*1e9
    return {**common,'first_step_attained':True,'reason':reason,'t1_s':float(h['t'][index]),'T_at_switch_C':float(h['T_C'][index]),'reached_nominal_T1':float(h['T_C'][index])>=T1-1e-9,'rho1':rho,'G1_nm':G1,'first_step_growth_fraction':(G1-G)/G,'X_J':float(h['X_J'][index]),'X_J_production':float(h['X_J_production'][index]),'X_J_relaxation':float(h['X_J_relaxation'][index]),'Lambda_TJ':float(h['Lambda_TJ'][index]),'K_TJ':float(h['K_TJ'][index]),'Lambda_over_K_TJ':float(h['Lambda_over_K_TJ'][index]),'P_comp_TJ':float(h['P_comp_TJ'][index]),'f_GBseg':phi[0]/z,'f_TJ':phi[1]/z,'f_iso_location':phi[2]/z,'C_GBseg':float(h['C_GBseg'][index]),'C_TJ':float(h['C_TJ'][index]),'f_clean_GB':float(h['f_clean_GB'][index]),'f_iso':float(h['f_iso'][index])}


def first_states(out):
    rows=[];states=[];rejected=[];total=len(mechanisms())*len(G0S)*len(T1S)*len(RATES);n=0
    for mid,p0 in mechanisms().items():
      for G in G0S:
       lp=replace(p0.action.location,base=replace(p0.action.location.base,G0=G*1e-9));p=replace(p0,action=replace(p0.action,location=lp))
       for T1 in T1S:
        for rate in RATES:
         n+=1
         if n%20==0:print('first-step',n,'/',total,flush=True)
         h=model.run(p,FixedBudgetRamp(rate,T1))
         for sw in SWITCHES:
          idx=np.flatnonzero(h['rho']>=sw-1e-12)
          if not len(idx):
           row=state_row(mid,p,G,T1,rate,sw,h,None,'switch_density_unattainable');rows.append(row);rejected.append(row);continue
          i=int(idx[0]);row=state_row(mid,p,G,T1,rate,sw,h,i);rows.append(row)
          if row['rho1']>=adaptive.TARGET-1e-12:rejected.append({**row,'reason':'target_already_reached_first_step'});continue
          if row['first_step_growth_fraction']>.20:rejected.append({**row,'reason':'first_step_growth_above_20pct'});continue
          states.append((row,model.final_state(h,p,i),p))
      write(out/'first_step_preparation_states.csv',rows);write(out/'preparation_rejected_cases.csv',rejected)
    return rows,states,rejected


def second_steps(out,states):
    points=[];bounds=[];groups={}
    for meta,state,p in states:
        arrays=(state.pore.phi_GBseg,state.pore.phi_TJ,state.pore.phi_iso)
        key=(meta['mechanism_id'],round(state.pore.rho,13),round(state.pore.G,18),round(state.X_J,13),*(a.round(15).tobytes() for a in arrays))
        groups.setdefault(key,[]).append((meta,state,p))
    total=len(groups);print('unique adaptive states',total,'from',len(states),'routes',flush=True)
    for i,members in enumerate(groups.values(),1):
        if i%10==0:print('adaptive-state',i,'/',total,flush=True)
        meta,state,p=members[0];x,_=adaptive.adaptive_from_state(meta['mechanism_id'],p,meta['G0_nm'],meta['T1_C'],meta['rho_switch'],state);group_id=f"state_{i:04d}"
        rep={k:meta[k] for k in ('mechanism_id','heating_rate_C_min','T_at_switch_C','rho1','G1_nm','first_step_growth_fraction','X_J','Lambda_over_K_TJ','P_comp_TJ','f_GBseg','f_TJ','f_iso_location','C_GBseg','C_TJ','f_clean_GB')};points.extend([{**r,**rep,'state_group_id':group_id,'equivalent_route_count':len(members)} for r in x])
        for route,_,_ in members:
          tag={k:route[k] for k in ('mechanism_id','heating_rate_C_min','T_at_switch_C','rho1','G1_nm','first_step_growth_fraction','X_J','Lambda_over_K_TJ','P_comp_TJ','f_GBseg','f_TJ','f_iso_location','C_GBseg','C_TJ','f_clean_GB')}
          for tol in discovery.TOLS:
           for practical,kind in ((False,'kinetic'),(True,'practical')):
            s=adaptive.status(x,tol,route['T1_C'],practical);eligible=[r for r in x if not practical or r['T2_C']<route['T1_C']]
            bounds.append(dict(parameter_id=route['mechanism_id'],mechanism_mode=p.mechanism_mode,G0_nm=route['G0_nm'],T1_C=route['T1_C'],rho_switch=route['rho_switch'],growth_tolerance=tol,map_type=kind,T2_min_searched_C=min(r['T2_C'] for r in eligible),T2_max_searched_C=max(r['T2_C'] for r in eligible),T2_cap_C=adaptive.CAP,extension_model_extrapolative=True,state_group_id=group_id,equivalent_route_count=len(members),**s,**tag))
        if i%25==0:write(out/'adaptive_second_step_boundaries.csv',bounds)
    write(out/'adaptive_second_step_boundaries.csv',bounds);return points,bounds


def score(bounds):
    practical=[];kinetic_only=[];censored=[]
    key=lambda r:(r['parameter_id'],r['G0_nm'],r['T1_C'],r['rho_switch'],r['heating_rate_C_min'],r['growth_tolerance'])
    practical_lookup={key(r):r for r in bounds if r['map_type']=='practical'}
    for r in bounds:
        if 'CENSORED' in r['boundary_status']:censored.append(r)
        if r['map_type']!='kinetic' or r['boundary_status']!='COMPLETE_WINDOW':continue
        pr=practical_lookup[key(r)]
        if pr['boundary_status']!='COMPLETE_WINDOW':kinetic_only.append({**r,'practical_boundary_status':pr['boundary_status']})
    for r in bounds:
        if r['map_type']=='practical' and r['boundary_status']=='COMPLETE_WINDOW':
            for ptol in PREP_TOLS:
                practical.append({**r,'preparation_growth_tolerance':ptol,'preparation_admissible':r['first_step_growth_fraction']<=ptol+1e-12,'practical_success':r['first_step_growth_fraction']<=ptol+1e-12})
    return practical,kinetic_only,censored


def summary(states,bounds,practical):
    rows=[]
    for mid in mechanisms():
      for pt in PREP_TOLS:
       for st in discovery.TOLS:
        q=[r for r in practical if r['mechanism_id']==mid and r['preparation_growth_tolerance']==pt and r['growth_tolerance']==st and r['practical_success']]
        rows.append(dict(mechanism_id=mid,preparation_growth_tolerance=pt,second_step_growth_tolerance=st,n_complete_practical=len(q),n_below_150nm=sum(r['G0_nm']<150 for r in q),min_G0_nm=min((r['G0_nm'] for r in q),default=math.nan),max_window_width_C=max((r.get('window_width_C',math.nan) for r in q),default=math.nan)))
    return rows


def plots(out,states,points,bounds,practical,summary_rows):
    ok=[r for r in practical if r['practical_success'] and r['second_step_growth_tolerance' if 'second_step_growth_tolerance' in r else 'growth_tolerance']==.05 and r['preparation_growth_tolerance']==.05]
    fig,ax=plt.subplots(figsize=(8,6));ax.scatter([r['G1_nm'] for r in ok],[r['T_first_success_C'] for r in ok],c=[r['T1_C'] for r in ok],cmap='viridis',alpha=.6);ax.set(xlabel='G1 [nm]',ylabel='first practical success T2 [C]',title='Practical Chen-style map');ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'T2_vs_G1_practical_map.png',dpi=150);plt.close(fig)
    kin=[r for r in bounds if r['map_type']=='kinetic' and r['growth_tolerance']==.05];fig,ax=plt.subplots(figsize=(8,6));ax.scatter([r.get('G1_nm',math.nan) for r in kin],[r.get('T_first_success_C',math.nan) for r in kin],c=['orange' if 'CENSORED' in r['boundary_status'] else 'tab:green' for r in kin],alpha=.5);ax.set(xlabel='G1 [nm]',ylabel='T2 [C]',title='Kinetic map (orange=censored)');fig.tight_layout();fig.savefig(out/'T2_vs_G1_kinetic_censored.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots();ax.scatter([r['rho1'] for r in ok],[r['T_first_success_C'] for r in ok],c=[r['G1_nm'] for r in ok],cmap='plasma');ax.set(xlabel='rho1',ylabel='T2 [C]');fig.tight_layout();fig.savefig(out/'T2_vs_rho1_practical_map.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots();ax.scatter([r['G1_nm'] for r in ok],[r['window_width_C'] for r in ok],alpha=.5);ax.set(xlabel='G1 [nm]',ylabel='practical window width [C]');fig.tight_layout();fig.savefig(out/'practical_window_width_vs_G1.png',dpi=150);plt.close(fig)
    attained=[r for r in states if r.get('first_step_attained')];
    for key,name,ylabel in (('first_step_growth_fraction','first_step_growth_vs_T1_rate.png','first-step growth fraction'),('X_J','XJ_vs_first_step_history.png','X_J'),('Lambda_over_K_TJ','Lambda_over_K_vs_first_step_history.png','Lambda/K')):
      fig,ax=plt.subplots();sc=ax.scatter([r['T1_C'] for r in attained],[r[key] for r in attained],c=[math.log10(r['heating_rate_C_min']) for r in attained],alpha=.35);ax.set(xlabel='T1 [C]',ylabel=ylabel);fig.colorbar(sc,ax=ax,label='log10 heating rate');fig.tight_layout();fig.savefig(out/name,dpi=150);plt.close(fig)
    fig,ax=plt.subplots();sample=points[::max(1,len(points)//3000)];ax.scatter([r['T2_C'] for r in sample],[r['P_comp_TJ'] for r in sample],s=7,alpha=.3);ax.set(xlabel='T2 [C]',ylabel='P_comp_TJ');fig.tight_layout();fig.savefig(out/'Pcomp_TJ_vs_T2.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots();z=[r for r in ok if r['mechanism_id']=='mech_019'];ax.scatter([r['T1_C'] for r in z],[r['rho_switch'] for r in z],c=[r['heating_rate_C_min'] for r in z],norm=matplotlib.colors.LogNorm(),alpha=.6);ax.set(xlabel='T1 [C]',ylabel='rho switch');fig.tight_layout();fig.savefig(out/'practical_window_T1_rho_switch.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots();
    for mid in mechanisms():
      q=[r for r in summary_rows if r['mechanism_id']==mid and r['preparation_growth_tolerance']==.05];ax.plot([r['second_step_growth_tolerance'] for r in q],[r['n_complete_practical'] for r in q],marker='o',label=mid)
    ax.set(xlabel='second-step tolerance',ylabel='complete practical cases');ax.legend();fig.tight_layout();fig.savefig(out/'mechanism_preparation_comparison.png',dpi=150);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',default='results/preparation_window_search');a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();states,admissible,rejected=first_states(out);points,bounds=second_steps(out,admissible);practical,kinetic,censored=score(bounds);summ=summary(states,bounds,practical);write(out/'practical_two_step_windows.csv',practical);write(out/'kinetic_only_windows.csv',kinetic);write(out/'censored_preparation_cases.csv',censored);write(out/'mechanism_preparation_summary.csv',summ);write(out/'adaptive_second_step_points.csv',points);plots(out,states,points,bounds,practical,summ);write(out/'runtime_summary.csv',[{'wall_s':time.perf_counter()-start,'first_step_states':len(states),'adaptive_states':len(admissible),'rho_target':adaptive.TARGET,'budget_h':96}]);print('DONE',len(states),len(admissible),len(practical),flush=True)


if __name__=='__main__':main()
