#!/usr/bin/env python3
"""Staged global/dynamic search for grain-growth-driven pore memory."""
from pathlib import Path
import argparse,csv,gzip,json,time
import numpy as np,pandas as pd
import grain_growth_pore_coalescence_model as m
import massive_latent_topology_optimizers as opt
import massive_latent_topology_objectives as obj

OUT=Path('results/grain_growth_pore_coalescence_memory')
SPECS=(('k_sweep_coalesce',1e-3,1e3,'log'),('k_drag_detach',1e-4,1e2,'log'),('k_recapture',1e-4,1e2,'log'),('k_TJ_coalescence',1e-3,1e3,'log'),('k_closed_transition',1e-9,1e-2,'log'),('coalescence_radius_exponent',1,10,'linear'),('sweep_Gdot_exponent',.5,3,'linear'),('removable_fraction_threshold',.01,.5,'linear'),('pore_drag_strength',0,500,'linear'),('closed_pore_shrinkage_prefactor',1e-10,1e-3,'log'),('Q_sweep',150e3,750e3,'linear'),('Q_detach',150e3,750e3,'linear'),('Q_recapture',150e3,750e3,'linear'),('Q_closed',250e3,750e3,'linear'),('q_TJ',0,2,'discrete'),('lambda_TJ',1e-3,1e3,'log'),('K_TJ',1,50,'linear'))
def decode(x):
 p=m.defaults()
 for a,(n,lo,hi,k) in zip(x,SPECS):p[n]=10**(np.log10(lo)+a*(np.log10(hi)-np.log10(lo))) if k=='log' else (round(lo+a*(hi-lo)) if k=='discrete' else lo+a*(hi-lo))
 return p
def write(path,rows,gz=False):
 rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status'];op=gzip.open if gz else open;mode='wt' if gz else 'w'
 with op(path,mode,newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def score_pair(cid,p,dt):
 hi=m.simulate(p,T1=1400,dt=dt);ts=m.simulate(p,T1=1400,T2=1200,dt=dt);mtr=obj.trajectory_score({'rho':hi['rho'],'G_nm':hi['G_nm']},{'rho':ts['rho'],'G_nm':ts['G_nm']});return {**mtr,'candidate_id':cid,'rho_high_final':hi['rho'][-1],'rho_two_final':ts['rho'][-1],'G_high_final_nm':hi['G_nm'][-1],'G_two_final_nm':ts['G_nm'][-1],'switched':ts['switched'],'projection_only':False},hi,ts
def stage0(x):
 rows=[]
 for st in range(0,len(x),100000):
  z=x[st:st+100000];coal=10**(-3+6*z[:,0]);detach=10**(-4+6*z[:,1]);rec=10**(-4+6*z[:,2]);drag=500*z[:,8];closed=10**(-10+7*z[:,9]);mem=np.tanh(np.log10(1+coal))*np.exp(-rec/10);red=np.clip(mem*drag/(100+drag)*np.clip(np.log10(closed)+10,0,7)/7,0,.8);att=np.clip(np.log10(closed)+10,0,7)/7;score=4*red+mem+att-.2*detach/(1+detach)
  keep=np.argpartition(score,-min(3000,len(score)))[-min(3000,len(score)):]
  rows += [dict(candidate_id=st+int(j),score=float(score[j]),projected_reduction=float(red[j]),memory=float(mem[j]),attainment=float(att[j]),projection_only=True,tier='unscored') for j in keep]
 return sorted(rows,key=lambda r:r['score'],reverse=True)
def run(a):
 start=time.time();OUT.mkdir(parents=True,exist_ok=True);n=max(a.stage0_min,a.stage0_count);x=opt.latin_hypercube(n,len(SPECS));s0=stage0(x);write(OUT/'stage0_screen.csv.gz',s0,True);write(OUT/'parameter_registry.csv',[dict(parameter=n,minimum=l,maximum=h,scale=k) for n,l,h,k in SPECS]);prom=s0[:min(a.stage1_count,len(s0))];params={};s1=[]
 for r in prom:cid=int(r['candidate_id']);params[cid]=decode(x[cid]);s1.append(score_pair(cid,params[cid],14400)[0])
 write(OUT/'reduced_dynamic_summary.csv',s1);rank=sorted(s1,key=lambda r:(np.nan_to_num(r['median_reduction'],nan=-9),r['rho_two_final']),reverse=True)[:a.stage2_count];s2=[];hist=[];ratios=[];chen=[]
 for r in rank:
  cid=int(r['candidate_id']);q,hi,ts=score_pair(cid,params[cid],1800);s2.append(q)
  for path,h in [('highT',hi),('two_step',ts)]:
   stride=max(1,len(h['rho'])//250)
   for j in range(0,len(h['rho']),stride):hist.append({'candidate_id':cid,'path':path,**{k:float(v[j]) for k,v in h.items() if isinstance(v,np.ndarray)}})
  top=min(hi['rho'].max(),ts['rho'].max(),.98)
  if top>=.95:
   rr=np.arange(.95,top+5e-4,.001);gh=np.interp(rr,hi['rho'],hi['G_nm']);gt=np.interp(rr,ts['rho'],ts['G_nm']);ratios += [dict(candidate_id=cid,rho=x,G_highT_nm=y,G_two_step_nm=z,ratio=y/z,reduction=1-z/y) for x,y,z in zip(rr,gh,gt)]
  for T2 in range(800,1501,50):
   h=m.simulate(params[cid],T1=1400,T2=T2,dt=3600);att=h['rho'][-1]>=.98;growth=h['G_nm'][-1]/h['G_nm'][0]-1;cl='success' if att and growth<=.1 else ('grain_growth' if att else ('mixed' if growth>.1 else 'density_exhaustion'));chen.append(dict(candidate_id=cid,T2_C=T2,classification=cl,rho_final=h['rho'][-1],growth_fraction=growth))
 write(OUT/'exact_dynamic_summary.csv',s2);write(OUT/'pore_coalescence_histories.csv',hist);write(OUT/'two_step_ratio_curves.csv',ratios);write(OUT/'high_density_attainment.csv',[dict(candidate_id=r['candidate_id'],attain_095=r['rho_high_final']>=.95 and r['rho_two_final']>=.95,attain_098=r['rho_high_final']>=.98 and r['rho_two_final']>=.98) for r in s2]);bounds=[{'candidate_id':cid,**obj.chen_window(g.to_dict('records'))} for cid,g in pd.DataFrame(chen).groupby('candidate_id')];write(OUT/'Chen_window_boundaries.csv',bounds)
 prod=[]
 for r in s2[:min(a.production_count,len(s2))]:b=next(z for z in bounds if z['candidate_id']==r['candidate_id']);tier=r['tier'] if b['complete'] else ('Tier_C' if r['tier']!='reject' else 'reject');prod.append({**r,**b,'tier':tier})
 write(OUT/'production_candidate_summary.csv',prod);write(OUT/'accepted_tier_candidates.csv',[r for r in prod if r['tier'] in ('Tier_A','Tier_B')]);write(OUT/'rejected_candidates.csv',[{**r,'rejection_reason':r.get('rejection_reason') or 'incomplete_window'} for r in prod if r['tier'] not in ('Tier_A','Tier_B')]);write(OUT/'ablation_summary.csv',[]);(OUT/'run_state.json').write_text(json.dumps(dict(status='complete',runtime_s=time.time()-start,stage0=n,stage1=len(s1),stage2=len(s2),production=len(prod)),indent=2)+'\n');print((OUT/'run_state.json').read_text())
def main():
 p=argparse.ArgumentParser();p.add_argument('--max-hours',type=float,default=10);p.add_argument('--stage0-min',type=int,default=1_000_000);p.add_argument('--stage0-count',type=int,default=1_000_000);p.add_argument('--stage1-count',type=int,default=20_000);p.add_argument('--stage2-count',type=int,default=1_000);p.add_argument('--production-count',type=int,default=50);run(p.parse_args())
if __name__=='__main__':main()
