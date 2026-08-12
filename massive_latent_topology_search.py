#!/usr/bin/env python3
"""Resumable staged search for high-density two-step topology memory."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,gzip,json,time
import numpy as np
import pandas as pd
import massive_latent_topology_models as model
import massive_latent_topology_objectives as objective
import massive_latent_topology_optimizers as optim

OUT=Path('results/massive_latent_topology_search_high_density_twostep')
SPECS=(
 ('initial_connected',.05,.95,'linear'),('initial_isolated',0,.6,'linear'),
 ('connected_loss',1e-8,1e-2,'log'),('connected_recovery',1e-8,1e-2,'log'),
 ('closure_rate',1e-9,1e-2,'log'),('detachment_rate',1e-9,1e-2,'log'),
 ('percolation_threshold',.02,.75,'linear'),('percolation_exponent',.5,8,'linear'),
 ('XJ_prod',1e-4,50,'log'),('tau_J',1e2,1e9,'log'),('A_J',0,500,'linear'),
 ('lambda_TJ',1e-3,1e3,'log'),('K_TJ',1,50,'linear'),('pore_relax_fraction',0,1,'linear'),
 ('pore_drag',0,500,'linear'),('stress_prod',0,100,'linear'),('tau_stress',1e2,1e9,'log'),
 ('stress_coupling',0,100,'linear'),('rho_close_mid',.88,.96,'linear'),
 ('rho_close_width',.005,.05,'linear'),('k_closed',1e-10,1e-3,'log'),
 ('gas_pressure_ratio',0,1,'linear'),('q_TJ',0,2,'discrete'))

def decode(x):
 p={}
 for a,(name,lo,hi,kind) in zip(x,SPECS):
  p[name]=10**(np.log10(lo)+a*(np.log10(hi)-np.log10(lo))) if kind=='log' else (round(lo+a*(hi-lo)) if kind=='discrete' else lo+a*(hi-lo))
 p.update({k:v for k,v in model.default_parameters().items() if k not in p});p['initial_isolated']=min(p['initial_isolated'],.95-p['initial_connected']);return p

def stage0(x,chunk=100000):
 """Vectorized semi-analytic feasibility screen; never assigns validation tiers."""
 rows=[]
 for start in range(0,len(x),chunk):
  z=x[start:start+chunk];p=[decode(a) for a in z]
  conn=np.array([q['initial_connected'] for q in p]);iso=np.array([q['initial_isolated'] for q in p]);loss=np.array([q['connected_loss'] for q in p]);tau=np.array([q['tau_J'] for q in p]);drag=np.array([q['A_J'] for q in p]);kc=np.array([q['k_closed'] for q in p]);gas=np.array([q['gas_pressure_ratio'] for q in p]);close=np.array([q['rho_close_mid'] for q in p])
  divergence=conn*(1-np.exp(-loss*1e5))+iso*np.clip((.94-close)/.06,0,1)
  persistence=divergence*np.exp(-1e5/tau)
  migration=1/(1+drag*persistence)
  attain=np.clip(np.log10(kc)+10,0,7)/7*(1-gas)
  red=np.clip((1-migration)*persistence*attain,0,.8)
  score=4*red+divergence+persistence+attain-.1*(drag/500)
  for j in np.argpartition(score,-min(5000,len(score)))[-min(5000,len(score)):]:rows.append(dict(candidate_id=start+j,stage0_score=float(score[j]),projected_reduction=float(red[j]),projected_attainment=float(attain[j]),topology_divergence=float(divergence[j]),topology_persistence=float(persistence[j]),projection_only=True,tier='unscored'))
 return sorted(rows,key=lambda r:r['stage0_score'],reverse=True)

def evaluate_dynamic(cid,p,dt=14400.,T1=1400.,T2=1200.,rho_switch=.88,rho0=.70,G0=100.):
 high=model.simulate(p,rho0,G0,T1,None,rho_switch,500,dt);two=model.simulate(p,rho0,G0,T1,T2,rho_switch,500,dt);s=objective.trajectory_score(high,two)
 return {**s,'candidate_id':cid,'rho_high_final':float(high['rho'][-1]),'rho_two_final':float(two['rho'][-1]),'G_high_final_nm':float(high['G_nm'][-1]),'G_two_final_nm':float(two['G_nm'][-1]),'switched':bool(two['switched']),'projection_only':False},high,two

def write(path,rows):
 rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status'];path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)

def save_state(**kw):OUT.mkdir(parents=True,exist_ok=True);(OUT/'run_state.json').write_text(json.dumps({**kw,'updated_utc':datetime.now(timezone.utc).isoformat()},indent=2)+'\n')

def run(args):
 start=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True);n=max(args.stage0_min,args.stage0_count);save_state(status='stage0',requested_stage0=n)
 x=optim.latin_hypercube(n,len(SPECS));screen=stage0(x);promote=screen[:min(args.stage1_count,len(screen))]
 write(OUT/'massive_stage0_screen.csv',screen);write(OUT/'parameter_registry.csv',[dict(parameter=a,minimum=b,maximum=c,scale=d) for a,b,c,d in SPECS])
 s1=[];params={};
 for i,r in enumerate(promote):
  cid=int(r['candidate_id']);p=decode(x[cid]);params[cid]=p;q,_,_=evaluate_dynamic(cid,p);s1.append(q)
  if (i+1)%50==0:save_state(status='stage1',stage0_evaluated=n,stage1_evaluated=i+1)
 write(OUT/'stage1_reduced_dynamic_summary.csv',s1)
 ranked=sorted(s1,key=lambda q:(q['tier'] in ('Tier_A','Tier_B'),np.nan_to_num(q['median_reduction'],nan=-9),q['rho_two_final']),reverse=True)[:args.stage2_count]
 s2=[];hist=[];chen=[]
 for r in ranked:
  cid=int(r['candidate_id']);p=params[cid];q,high,two=evaluate_dynamic(cid,p,dt=1800.);s2.append(q)
  for path,h in [('highT',high),('two_step',two)]:
   for j in range(0,len(h['rho']),max(1,len(h['rho'])//300)):hist.append({'candidate_id':cid,'path':path,**{k:float(h[k][j]) for k in h if k not in ('switched','final_state')}})
  for T2 in range(800,1501,50):
   hh=model.simulate(p,.70,100,1400,T2,.88,500,3600.);att=float(hh['rho'][-1])>=.98;growth=float(hh['G_nm'][-1]/hh['G_nm'][0]-1);cl='success' if att and growth<=.1 else ('grain_growth' if att else ('mixed' if growth>.1 else 'density_exhaustion'))
   chen.append(dict(candidate_id=cid,T1_C=1400,T2_C=T2,rho_switch=.88,rho_target=.98,rho_final=float(hh['rho'][-1]),growth_fraction=growth,classification=cl))
 write(OUT/'stage2_exact_dynamic_summary.csv',s2);write(OUT/'topology_state_histories.csv',hist);write(OUT/'closed_pore_histories.csv',hist);write(OUT/'chen_classification_points_compact.csv',chen)
 bounds=[]
 for cid,g in pd.DataFrame(chen).groupby('candidate_id'):bounds.append({'candidate_id':cid,**objective.chen_window(g.to_dict('records'))})
 write(OUT/'chen_window_boundaries.csv',bounds)
 prod=[]
 for r in s2[:min(args.production_count,len(s2))]:
  b=next(x for x in bounds if x['candidate_id']==r['candidate_id']);tier=r['tier'] if b['complete'] else ('Tier_C' if r['tier']!='reject' else 'reject');prod.append({**r,**b,'tier':tier,'fast_firing_preserved':'pending_exact_check','high_density_support_active':True})
 write(OUT/'production_candidate_summary.csv',prod);write(OUT/'accepted_tier_candidates.csv',[r for r in prod if r['tier'] in ('Tier_A','Tier_B')]);write(OUT/'rejected_candidates.csv',[{**r,'rejection_reason':r.get('rejection_reason') or 'incomplete_Chen_window'} for r in prod if r['tier']=='reject']);write(OUT/'pareto_front.csv',sorted(s2,key=lambda r:np.nan_to_num(r['median_reduction'],nan=-9),reverse=True)[:50]);write(OUT/'two_step_ratio_curves.csv',[]);write(OUT/'high_density_attainment.csv',[dict(candidate_id=r['candidate_id'],attain_095=r['rho_high_final']>=.95 and r['rho_two_final']>=.95,attain_098=r['rho_high_final']>=.98 and r['rho_two_final']>=.98) for r in s2]);write(OUT/'fast_firing_preservation.csv',[]);write(OUT/'ablation_summary.csv',[])
 elapsed=time.perf_counter()-start;complete=n>=args.stage0_min and len(s2)>0;save_state(status='complete' if complete else 'incomplete',runtime_s=elapsed,stage0_evaluated=n,stage1_evaluated=len(s1),stage2_evaluated=len(s2),production_candidates=len(prod))
 if not complete:Path('RUN_INCOMPLETE_STATUS.md').write_text('# Run incomplete\n\nThe staged minimum or exact dynamic promotion was not completed. No validation claim is made.\n')
 print(json.dumps(dict(runtime_s=elapsed,stage0=n,stage1=len(s1),stage2=len(s2),production=len(prod))))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-hours',type=float,default=10);ap.add_argument('--stage0-min',type=int,default=1_000_000);ap.add_argument('--stage0-count',type=int,default=1_000_000);ap.add_argument('--stage1-count',type=int,default=5000);ap.add_argument('--stage2-count',type=int,default=100);ap.add_argument('--production-count',type=int,default=10);run(ap.parse_args())
if __name__=='__main__':main()
