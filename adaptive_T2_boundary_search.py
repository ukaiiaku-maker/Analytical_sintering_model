#!/usr/bin/env python3
"""Adaptive, censor-aware Chen-style T2 boundary search."""
from __future__ import annotations
import argparse,csv,math,time
from dataclasses import replace
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import agentic_mechanism_model as model
import agentic_mechanism_search as search
import topology_constrained_sintering as aggregate

COARSE=tuple(float(x) for x in range(900,1301,25));DOWN=tuple(float(x) for x in range(800,900,25));UP=tuple(float(x) for x in range(1325,1551,25));CAP=1550.;TARGET=search.TARGET;BUDGET=search.BUDGET


def write(path,rows,empty=('parameter_id','boundary_status')):
    fields=[]
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    if not fields:fields=list(empty)
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def survivor_params():
    lookup={pid:p for pid,p in search.design()};return {k:lookup[k] for k in ('mech_009','mech_019')}


def point(state,p,pid,G,T1,sw,T2,stage):
    h=model.run(p,aggregate.Iso(T2,BUDGET),initial=state);rho2=float(h['rho'][-1]);G2=float(h['G'][-1])*1e9;G1=state.pore.G*1e9;growth=(G2-G1)/G1
    return dict(parameter_id=pid,mechanism_mode=p.mechanism_mode,G0_nm=G,T1_C=T1,rho_switch=sw,T2_C=T2,scan_stage=stage,rho_target=TARGET,first_budget_h=96,second_budget_h=96,rho1=state.pore.rho,G1_nm=G1,rho2=rho2,G2_nm=G2,growth_fraction=growth,X_J=float(h['X_J'][-1]),Lambda_TJ=float(h['Lambda_TJ'][-1]),K_TJ=float(h['K_TJ'][-1]),P_comp_TJ=float(h['P_comp_TJ'][-1]))


def status(points,tol,T1,practical=False):
    q=sorted((r for r in points if not practical or r['T2_C']<T1),key=lambda r:r['T2_C'])
    if not q:return dict(boundary_status='NO_OVERLAP',n_success=0)
    for r in q:r['_class']=search.classify(True,r['rho1'],r['rho2'],r['growth_fraction'],tol)
    if q[0]['rho1']>=TARGET-1e-12:return dict(boundary_status='INELIGIBLE_TARGET_ALREADY_REACHED',n_success=0)
    dense=[r for r in q if r['rho2']>=TARGET-1e-12];success=[r for r in q if r['_class']=='SUCCESS']
    if not dense:return dict(boundary_status='LOWER_BOUND_RIGHT_CENSORED',n_success=0)
    lower=min(dense,key=lambda r:r['T2_C']);lower_fail=[r for r in q if r['T2_C']<lower['T2_C'] and r['rho2']<TARGET-1e-12]
    lower_bracket=bool(lower_fail)
    if not success:
        return dict(boundary_status='NO_OVERLAP',n_success=0,T_lower_density_C=lower['T2_C'],lower_bracketed=lower_bracket)
    first=min(success,key=lambda r:r['T2_C']);last=max(success,key=lambda r:r['T2_C']);upper=[r for r in q if r['T2_C']>last['T2_C'] and r['_class']=='GRAIN_GROWTH_FAILURE']
    upper_bracket=bool(upper);upperT=min((r['T2_C'] for r in upper),default=math.nan)
    if not lower_bracket:label='LOWER_BOUND_LEFT_CENSORED'
    elif not upper_bracket:label='UPPER_BOUND_RIGHT_CENSORED'
    else:label='COMPLETE_WINDOW'
    return dict(boundary_status=label,n_success=len(success),T_lower_density_C=lower['T2_C'],T_first_success_C=first['T2_C'],T_last_success_C=last['T2_C'],T_upper_growth_C=upperT,window_width_C=last['T2_C']-first['T2_C'],lower_bracketed=lower_bracket,upper_bracketed=upper_bracket)


def intervals_for_refinement(points,tol):
    q=sorted(points,key=lambda r:r['T2_C']);intervals=[]
    for a,b in zip(q,q[1:]):
        ca=search.classify(True,a['rho1'],a['rho2'],a['growth_fraction'],tol);cb=search.classify(True,b['rho1'],b['rho2'],b['growth_fraction'],tol)
        da=a['rho2']>=TARGET-1e-12;db=b['rho2']>=TARGET-1e-12;ga=a['growth_fraction']<=tol;gb=b['growth_fraction']<=tol
        if ca!=cb or da!=db or ga!=gb:intervals.append((a['T2_C'],b['T2_C']))
    return intervals


def adaptive_from_state(pid,p,G,T1,sw,state):
    """Run the adaptive second step from an externally prepared local state."""
    pts=[point(state,p,pid,G,T1,sw,T,'coarse') for T in COARSE]
    # Downward extension is needed only when density is already attained at the left edge.
    if pts[0]['rho2']>=TARGET-1e-12:
        for T in reversed(DOWN):
            r=point(state,p,pid,G,T1,sw,T,'downward_extension');pts.append(r)
            if r['rho2']<TARGET-1e-12:break
    # Extend upward while any tolerance has a success without its upper failure.
    for T in UP:
        needs=False
        for tol in search.TOLS:
            s=status(pts,tol,T1,False)
            if s.get('n_success',0)>0 and not s.get('upper_bracketed',False):needs=True
        if not needs:break
        pts.append(point(state,p,pid,G,T1,sw,T,'upward_extension'))
    # Local 10 C refinement only inside changing 25 C brackets.
    intervals=set()
    for tol in search.TOLS:intervals.update(intervals_for_refinement(pts,tol))
    existing={r['T2_C'] for r in pts}
    for lo,hi in sorted(intervals):
        for T in np.arange(lo+10,hi,10):
            if float(T) not in existing:pts.append(point(state,p,pid,G,T1,sw,float(T),'local_refinement'));existing.add(float(T))
    rows=[]
    for tol in search.TOLS:
        for practical,kind in ((False,'kinetic'),(True,'practical')):
            s=status(pts,tol,T1,practical);rows.append(dict(parameter_id=pid,mechanism_mode=p.mechanism_mode,G0_nm=G,T1_C=T1,rho_switch=sw,growth_tolerance=tol,map_type=kind,T2_min_searched_C=min(r['T2_C'] for r in pts if not practical or r['T2_C']<T1),T2_max_searched_C=max(r['T2_C'] for r in pts if not practical or r['T2_C']<T1),T2_cap_C=CAP,extension_model_extrapolative=True,**s))
    for r in pts:
        for tol in search.TOLS:r[f'classification_{int(tol*100)}pct']=search.classify(True,r['rho1'],r['rho2'],r['growth_fraction'],tol)
    return pts,rows


def adaptive_group(pid,p,G,T1,sw):
    lp=replace(p.action.location,base=replace(p.action.location.base,G0=G*1e-9));p=replace(p,action=replace(p.action,location=lp));h1=model.run(p,aggregate.Iso(T1,BUDGET),stop_at_rho=sw);rho1=float(h1['rho'][-1])
    if rho1<sw-1e-12:
        base=dict(parameter_id=pid,G0_nm=G,T1_C=T1,rho_switch=sw,boundary_status='UNATTAINABLE_FIRST_STEP',n_success=0);return [],[{**base,'growth_tolerance':tol,'map_type':kind} for tol in search.TOLS for kind in ('kinetic','practical')]
    return adaptive_from_state(pid,p,G,T1,sw,model.final_state(h1,p))


def refinement_design():
    bases=survivor_params();out=[]
    for bid,p in bases.items():
        variants=(('A_half',{'A_J':p.A_J*.5}),('A_double',{'A_J':p.A_J*2}),('tau_low',{'tau_J_ref_s':p.tau_J_ref_s*.3}),('tau_high',{'tau_J_ref_s':p.tau_J_ref_s*3}),('prod_half',{'XJ_prod_TJ':p.XJ_prod_TJ*.5}),('prod_double',{'XJ_prod_TJ':p.XJ_prod_TJ*2}),('lambda_half',{'lambda_TJ_ref':p.lambda_TJ_ref*.5}),('lambda_double',{'lambda_TJ_ref':p.lambda_TJ_ref*2}),('Q_minus40',{'Q_TJ_event_J_mol':p.Q_TJ_event_J_mol-40e3}),('Q_plus40',{'Q_TJ_event_J_mol':p.Q_TJ_event_J_mol+40e3}),('q0_visible',{'mechanism_mode':'persistent_tj_multihit_q0','q_TJ':0}))
        for name,kw in variants:out.append((f'{bid}_{name}',replace(p,**kw),bid,name))
    return out


def run_refinement():
    rows=[]
    for pid,p,bid,name in refinement_design():
        _,_,tra,counts,cats,reason=search.evaluate((pid,p));rows.append(dict(parameter_id=pid,parent=bid,variation=name,mechanism_mode=p.mechanism_mode,A_J=p.A_J,tau_J_ref_s=p.tau_J_ref_s,XJ_prod_TJ=p.XJ_prod_TJ,lambda_TJ_ref=p.lambda_TJ_ref,K_TJ0=p.K_TJ0,Q_TJ_event_J_mol=p.Q_TJ_event_J_mol,q_TJ=p.q_TJ,success_5pct=counts[.05],success_10pct=counts[.10],has_lower_failure='DENSIFICATION_EXHAUSTION_FAILURE' in cats,has_upper_failure='GRAIN_GROWTH_FAILURE' in cats,decision='reject' if reason else 'promising_reduced_only',reason=reason))
    return rows


def plots(out,points,bounds,refine):
    complete=[r for r in bounds if r['map_type']=='kinetic' and r['boundary_status']=='COMPLETE_WINDOW'];censored=[r for r in bounds if r['map_type']=='kinetic' and 'CENSORED' in r['boundary_status']]
    for kind,name in (('kinetic','chen_kinetic_T2_vs_G_censored.png'),('practical','chen_practical_T2_less_T1.png')):
        fig,ax=plt.subplots(figsize=(8,6));q=[r for r in bounds if r['map_type']==kind and r['growth_tolerance']==.05]
        colors={'COMPLETE_WINDOW':'tab:green','UPPER_BOUND_RIGHT_CENSORED':'orange','LOWER_BOUND_LEFT_CENSORED':'gold','LOWER_BOUND_RIGHT_CENSORED':'tab:blue','NO_OVERLAP':'tab:red','UNATTAINABLE_FIRST_STEP':'.5'}
        for c,col in colors.items():
            z=[r for r in q if r['boundary_status']==c];ax.scatter([r['G0_nm'] for r in z],[r.get('T_first_success_C',math.nan) for r in z],c=col,label=c,marker='>' if 'CENSORED' in c else 'o',alpha=.65)
        ax.set(xscale='log',xlabel='G0 [nm]',ylabel='first success T2 [C]',title=kind);ax.grid(alpha=.2);ax.legend(fontsize=6);fig.tight_layout();fig.savefig(out/name,dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5))
    for pid in sorted(set(r['parameter_id'] for r in complete)):
        q=[r for r in complete if r['parameter_id']==pid and r['growth_tolerance']==.05];ax.scatter([r['G0_nm'] for r in q],[r['T_lower_density_C'] for r in q],label=f'{pid} lower');ax.scatter([r['G0_nm'] for r in q],[r['T_upper_growth_C'] for r in q],marker='x',label=f'{pid} upper')
    ax.set(xscale='log',xlabel='G0 [nm]',ylabel='boundary T2 [C]');ax.legend(fontsize=7);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'T_lower_T_upper_vs_G0.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5))
    for pid in sorted(set(r['parameter_id'] for r in bounds)):
        q=[r for r in bounds if r['parameter_id']==pid and r['map_type']=='kinetic' and r['growth_tolerance']==.05 and r.get('n_success',0)>0];ax.scatter([r['G0_nm'] for r in q],[r.get('window_width_C',math.nan) for r in q],marker='>' if any('CENSORED' in r['boundary_status'] for r in q) else 'o',label=pid)
    ax.set(xscale='log',xlabel='G0 [nm]',ylabel='sampled success width [C]');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'window_width_vs_G0_censored.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));reps=[r for r in points if r['G0_nm'] in (150.,300.) and r['T1_C']==1300 and r['rho_switch']==.85]
    for pid in sorted(set(r['parameter_id'] for r in reps)):
        q=sorted([r for r in reps if r['parameter_id']==pid],key=lambda r:r['T2_C']);ax.plot([r['T2_C'] for r in q],[r['growth_fraction'] for r in q],marker='.',label=pid)
    ax.axhline(.05,color='k',ls='--');ax.set(xlabel='T2 [C]',ylabel='growth fraction',title='Representative adaptive classifications');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'classification_vs_T2_representative.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));colors={'promising_reduced_only':'tab:green','reject':'tab:red'}
    for dec,col in colors.items():
        q=[r for r in refine if r['decision']==dec];ax.scatter([r['A_J'] for r in q],[r['lambda_TJ_ref'] for r in q],c=col,label=dec,alpha=.7)
    ax.set(xlabel='A_J',ylabel='lambda_TJ_ref',title='OAT refinement decisions');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'parameter_refinement_decision_map.png',dpi=150);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',default='results/agentic_mechanism_search');a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();points=[];bounds=[]
    for pid,p in survivor_params().items():
        print('adaptive',pid,flush=True)
        for G in search.G_FULL:
            for T1 in search.T1_FULL:
                for sw in search.SW_FULL:
                    x,y=adaptive_group(pid,p,G,T1,sw);points.extend(x);bounds.extend(y)
    refine=run_refinement();kin=[];practical=[]
    for r in points:
        for tol in search.TOLS:
            kin.append({**r,'growth_tolerance':tol,'classification':r[f'classification_{int(tol*100)}pct']})
            if r['T2_C']<r['T1_C']:practical.append({**r,'growth_tolerance':tol,'classification':r[f'classification_{int(tol*100)}pct']})
    cens=[r for r in bounds if 'CENSORED' in r['boundary_status']];complete_practical=[r for r in bounds if r['map_type']=='practical' and r['boundary_status']=='COMPLETE_WINDOW']
    write(out/'adaptive_T2_boundary_points.csv',points);write(out/'adaptive_window_boundaries.csv',bounds);write(out/'censored_window_cases.csv',cens);write(out/'practical_T2_less_than_T1_windows.csv',complete_practical);write(out/'kinetic_window_map.csv',kin);write(out/'practical_two_step_map.csv',practical);write(out/'extended_parameter_refinement.csv',refine);plots(out,points,bounds,refine);write(out/'adaptive_runtime_summary.csv',[{'wall_s':time.perf_counter()-start,'T2_cap_C':CAP,'rho_target':TARGET,'step_budget_h':96,'parameter_sets_refined':len(refine)}]);print('DONE',len(points),len(bounds),len(cens),len(complete_practical))


if __name__=='__main__':main()
