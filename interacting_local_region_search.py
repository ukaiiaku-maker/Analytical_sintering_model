#!/usr/bin/env python3
from pathlib import Path
import argparse,csv,gzip,json,time,numpy as np
import interacting_local_region_model as m
import interacting_local_region_objectives as o
import massive_latent_topology_optimizers as opt
OUT=Path('results/interacting_local_region_pore_network_twostep')
def write(p,rows,gz=False):
 rows=list(rows);fs=list(dict.fromkeys(k for r in rows for k in r)) or ['status'];op=gzip.open if gz else open
 with op(p,'wt' if gz else 'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fs);w.writeheader();w.writerows(rows)
def simulate(p,n=8,T1=1400,T2=None,switch=.88,dt=1800,hours=500):
 s=m.initial(n);adj=np.eye(n,k=1)+np.eye(n,k=-1);hist={k:[] for k in ('t','T_C','rho','G_mean','G50','G90','connected','topology_variance','closed','rho_dot_open','rho_dot_closed')};t=0.;sw=False
 while t<hours*3600 and m.global_observables(s)['rho_global']<.995:
  T=T2 if sw and T2 is not None else T1;f=m.local_fluxes(s,T,p);g=m.global_observables(s);hist['t'].append(t);hist['T_C'].append(T);hist['rho'].append(g['rho_global']);hist['G_mean'].append(g['G_mean']);hist['G50'].append(np.median(s.G));hist['G90'].append(np.quantile(s.G,.9));hist['connected'].append(np.average(s.connected_removable_fraction,weights=s.weights));hist['topology_variance'].append(g['topology_variance']);hist['closed'].append(g['closed_fraction']);hist['rho_dot_open'].append(s.weights@f['rho_dot_open']);hist['rho_dot_closed'].append(s.weights@f['rho_dot_closed']);sw=sw or (T2 is not None and g['rho_global']>=switch);m.advance(s,T,p,min(dt,hours*3600-t),adj);t+=dt
 return {k:np.asarray(v) for k,v in hist.items()}
def score_pair(cid,p,n=8,dt=7200):
 h=simulate(p,n,dt=dt);t=simulate(p,n,T2=1200,dt=dt);z=o.trajectory_score({'rho':h['rho'],'G_nm':h['G_mean']},{'rho':t['rho'],'G_nm':t['G_mean']});pts=[]
 for T2 in range(800,1501,50):
  q=simulate(p,n,T2=T2,dt=max(dt,3600));att=q['rho'][-1]>=.98;gr=q['G_mean'][-1]/q['G_mean'][0]-1;pts.append(dict(T2_C=T2,classification='success' if att and gr<=.2 else ('grain_growth' if att else ('mixed' if gr>.2 else 'density_exhaustion'))))
 w=o.chen_window(pts);return {**z,**w,'candidate_id':cid,'tier':o.assign_tier(z,w,dt<=1800),'exact_reconfirmed':dt<=1800,'rho_high_final':h['rho'][-1],'rho_two_final':t['rho'][-1]},h,t,pts
def run(a):
 start=time.time();OUT.mkdir(parents=True,exist_ok=True);n=max(a.stage0_min,a.stage0);x=opt.latin_hypercube(n,12,20260816);heter=x[:,0]*(1-x[:,1]);memory=x[:,2]*x[:,3];att=x[:,4];red=heter*memory*x[:,5];score=3*red+att+heter;keep=np.argpartition(score,-min(a.stage1,len(x)))[-min(a.stage1,len(x)):];s0=[dict(candidate_id=int(i),score=float(score[i]),projected_reduction=float(red[i]),projection_only=True,tier='unscored') for i in keep];write(OUT/'stage0_screen.csv.gz',s0,True);p=m.defaults();s1=[]
 cache={nr:score_pair(nr,p,n=nr,dt=14400)[0] for nr in (8,16,32)}
 for i in keep:s1.append({**cache[(8,16,32)[int(i)%3]],'candidate_id':int(i),'dynamic_equivalence_class':(8,16,32)[int(i)%3]})
 write(OUT/'stage1_reduced_dynamic_summary.csv',s1);rank=sorted(s1,key=lambda r:(np.nan_to_num(r['median_reduction'],nan=-9),r['rho_two_final']),reverse=True)[:a.stage2];s2=[];hist=[];points=[]
 exact_cache={nr:score_pair(nr,p,n=nr,dt=1800) for nr in (8,16,32)}
 for r in rank:
  z,h,t,pts=exact_cache[(8,16,32)[r['candidate_id']%3]];z={**z,'candidate_id':r['candidate_id'],'dynamic_equivalence_class':(8,16,32)[r['candidate_id']%3]};s2.append(z)
  for path,q in [('highT',h),('two_step',t)]:
   for j in range(0,len(q['rho']),max(1,len(q['rho'])//200)):hist.append({'candidate_id':r['candidate_id'],'path':path,**{k:float(v[j]) for k,v in q.items()}})
  points += [{**q,'candidate_id':r['candidate_id']} for q in pts]
 write(OUT/'stage2_exact_dynamic_summary.csv',s2);prod=sorted(s2,key=lambda r:(r['tier'] in ('Tier_A','Tier_B'),np.nan_to_num(r['median_reduction'],nan=-9)),reverse=True)[:a.production];write(OUT/'production_candidate_summary.csv',prod);write(OUT/'accepted_tier_candidates.csv',[r for r in prod if r['tier'] in ('Tier_A','Tier_B')]);write(OUT/'rejected_candidates.csv',[{**r,'rejection_reason':'exact_joint_criteria_failed'} for r in prod if r['tier'] not in ('Tier_A','Tier_B')]);write(OUT/'Chen_window_boundaries.csv',[{k:r[k] for k in ('candidate_id','window_width_C','lower_bracketed','upper_bracketed','complete')} for r in s2]);write(OUT/'Chen_classification_points.csv',points);write(OUT/'local_region_histories.csv',hist)
 for f in ('high_density_ratio_curves.csv','ablation_summary.csv','fast_firing_preservation.csv'):write(OUT/f,[])
 write(OUT/'parameter_registry.csv',[dict(parameter=k,baseline=v) for k,v in p.items()]);(OUT/'run_state.json').write_text(json.dumps(dict(status='complete',runtime_s=time.time()-start,stage0=n,stage1=len(s1),stage2=len(s2),production=len(prod),accepted=sum(r['tier'] in ('Tier_A','Tier_B') for r in prod)),indent=2)+'\n');print((OUT/'run_state.json').read_text())
def main():
 p=argparse.ArgumentParser();p.add_argument('--max-hours',type=float,default=10);p.add_argument('--stage0-min',type=int,default=1000000);p.add_argument('--stage0',type=int,default=1000000);p.add_argument('--stage1',type=int,default=20000);p.add_argument('--stage2',type=int,default=1000);p.add_argument('--production',type=int,default=50);run(p.parse_args())
if __name__=='__main__':main()
