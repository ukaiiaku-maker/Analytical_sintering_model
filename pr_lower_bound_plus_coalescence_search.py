#!/usr/bin/env python3
"""Focused staged search for a PR-controlled Chen lower boundary."""
from pathlib import Path
import argparse,csv,gzip,json,time
import numpy as np,pandas as pd
import pr_lower_bound_coalescence_model as m
import grain_growth_pore_coalescence_search as prior
import massive_latent_topology_optimizers as opt
import massive_latent_topology_objectives as obj
OUT=Path('results/pr_lower_bound_plus_coalescence_twostep');BASE_ID=540979
SPECS=(('k_PR_ref',1e-8,1e-2,'log'),('Q_PR',100e3,400e3,'linear'),('T_PR_ref_C',900,1200,'linear'),('low_activity_gate_mid',.05,.5,'linear'),('low_activity_gate_width',.03,.2,'linear'),('activity_power',1,4,'discrete'),('PR_to_smoothing_fraction',0,1,'linear'),('PR_to_large_tail_fraction',0,1,'linear'),('PR_to_TJ_fraction',0,1,'linear'),('PR_to_isolated_fraction',0,1,'linear'),('PR_to_closed_fraction',0,1,'linear'),('drive_loss_coupling',0,10,'linear'),('connected_sink_loss_coupling',0,10,'linear'),('large_tail_penalty_power',.5,4,'linear'),('PR_damage_persistence_tau',1e3,1e7,'log'),('Q_PR_damage_relax',100e3,500e3,'linear'),('coal_factor',.3,10,'log'),('detach_factor',.3,10,'log'),('recapture_factor',.1,3,'log'),('closed_factor',.3,3,'log'))
def write(path,rows,gz=False):
 rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status'];op=gzip.open if gz else open
 with op(path,'wt' if gz else 'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def base_params():
 x=opt.latin_hypercube(1_000_000,len(prior.SPECS));return prior.decode(x[BASE_ID])
def decode(x,b):
 p={**m.defaults(),**b}
 for a,(n,lo,hi,k) in zip(x,SPECS):p[n]=10**(np.log10(lo)+a*(np.log10(hi)-np.log10(lo))) if k=='log' else (round(lo+a*(hi-lo)) if k=='discrete' else lo+a*(hi-lo))
 p['k_sweep_coalesce']*=p.pop('coal_factor');p['k_drag_detach']*=p.pop('detach_factor');p['k_recapture']*=p.pop('recapture_factor');p['k_closed_transition']*=p.pop('closed_factor');return p
def pair(cid,p,dt=7200,T2=1200):
 h=m.simulate(p,T1=1400,dt=dt);t=m.simulate(p,T1=1400,T2=T2,dt=dt);z=obj.trajectory_score({'rho':h['rho'],'G_nm':h['G_nm']},{'rho':t['rho'],'G_nm':t['G_nm']});return {**z,'candidate_id':cid,'rho_high_final':h['rho'][-1],'rho_two_final':t['rho'][-1],'G_high_final_nm':h['G_nm'][-1],'G_two_final_nm':t['G_nm'][-1]},h,t
def screen(x):
 k=10**(-8+6*x[:,0]);Q=100+300*x[:,1];gate=.05+.45*x[:,3];drive=10*x[:,11];iso=x[:,9]+x[:,10];damage=np.clip(np.log10(k)+8,0,6)/6*np.exp(-(Q-250)**2/80000)*gate;lower=damage*drive*(.3+iso);preserve=np.exp(-damage*.4);score=2*lower+preserve-np.abs(lower-.7);return [dict(candidate_id=i,lower_bound_likelihood=float(lower[i]),trajectory_preservation=float(preserve[i]),score=float(score[i]),projection_only=True,tier='unscored') for i in np.argpartition(score,-min(50000,len(score)))[-min(50000,len(score)):]]
def run(a):
 start=time.time();OUT.mkdir(parents=True,exist_ok=True);n=max(a.stageA_min,a.stageA_count);x=opt.latin_hypercube(n,len(SPECS),20260813);b=base_params();sa=sorted(screen(x),key=lambda r:r['score'],reverse=True);write(OUT/'stageA_PR_feasibility_screen.csv.gz',sa,True);write(OUT/'parameter_registry.csv',[dict(parameter=n,minimum=l,maximum=h,scale=k) for n,l,h,k in SPECS]);prom=sa[:min(a.stageB_count,len(sa))];params={};sb=[]
 for r in prom:cid=int(r['candidate_id']);params[cid]=decode(x[cid],b);q,_,_=pair(cid,params[cid]);sb.append(q)
 write(OUT/'stageB_reduced_T2_scan_summary.csv',sb);rank=sorted(sb,key=lambda r:(np.nan_to_num(r['median_reduction'],nan=-9),r['rho_two_final']),reverse=True)[:a.stageC_count];sc=[];points=[];hist=[];rat=[]
 for r in rank:
  cid=int(r['candidate_id']);p=params[cid];q,hi,ts=pair(cid,p,1800);classes=[]
  for T2 in range(800,1551,25):
   z=m.simulate(p,T1=1400,T2=T2,dt=3600);att=z['rho'][-1]>=.98;growth=z['G_nm'][-1]/z['G_nm'][0]-1;cl='success' if att and growth<=.2 else ('grain_growth' if att else ('mixed' if growth>.2 else 'density_exhaustion'));row=dict(candidate_id=cid,T2_C=T2,rho_final=z['rho'][-1],growth_fraction=growth,classification=cl,mean_w_PR=float(np.mean(z['w_PR'])),mean_w_dens=float(np.mean(z['w_dens'])));points.append(row);classes.append(row)
  w=obj.chen_window(classes);sc.append({**q,**w,'tier':'Tier_B' if q['span20']>=.02 and w['complete'] else ('Tier_C' if q['max_reduction']>0 else 'reject')})
  for path,z in [('highT',hi),('two_step',ts)]:
   for j in range(0,len(z['rho']),max(1,len(z['rho'])//250)):hist.append({'candidate_id':cid,'path':path,**{k:float(v[j]) for k,v in z.items() if isinstance(v,np.ndarray)}})
  top=min(hi['rho'].max(),ts['rho'].max(),.98)
  if top>=.95:
   rr=np.arange(.95,top+5e-4,.001);gh=np.interp(rr,hi['rho'],hi['G_nm']);gt=np.interp(rr,ts['rho'],ts['G_nm']);rat += [dict(candidate_id=cid,rho=x,G_highT_nm=y,G_two_step_nm=z,reduction=1-z/y) for x,y,z in zip(rr,gh,gt)]
 write(OUT/'stageC_refined_boundary_summary.csv',sc);write(OUT/'Chen_classification_points.csv',points);write(OUT/'Chen_window_boundaries.csv',[{k:r[k] for k in ('candidate_id','window_width_C','lower_bracketed','upper_bracketed','complete')} for r in sc]);write(OUT/'PR_energy_partition_histories.csv',hist);write(OUT/'pore_topology_histories.csv',hist);write(OUT/'high_density_ratio_curves.csv',rat);prod=sorted(sc,key=lambda r:(r['tier']=='Tier_B',r['median_reduction']),reverse=True)[:a.production_count];write(OUT/'production_candidate_summary.csv',prod);write(OUT/'accepted_tier_candidates.csv',[r for r in prod if r['tier']=='Tier_B']);write(OUT/'rejected_candidates.csv',[{**r,'rejection_reason':'incomplete_window' if not r['complete'] else 'reduction_below_threshold'} for r in prod if r['tier']!='Tier_B']);write(OUT/'ablation_summary.csv',[]);write(OUT/'fast_firing_preservation.csv',[{'status':'pending_for_accepted_candidate' if any(r['tier']=='Tier_B' for r in prod) else 'not_promoted','material_regime':'frozen nucleation-limited'}]);(OUT/'run_state.json').write_text(json.dumps(dict(status='complete',runtime_s=time.time()-start,stageA=n,stageB=len(sb),stageC=len(sc),production=len(prod)),indent=2)+'\n');print((OUT/'run_state.json').read_text())
def main():
 p=argparse.ArgumentParser();p.add_argument('--max-hours',type=float,default=10);p.add_argument('--stageA-min',type=int,default=100000);p.add_argument('--stageA-count',type=int,default=250000);p.add_argument('--stageB-count',type=int,default=5000);p.add_argument('--stageC-count',type=int,default=500);p.add_argument('--production-count',type=int,default=20);run(p.parse_args())
if __name__=='__main__':main()
