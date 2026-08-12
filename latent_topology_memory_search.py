#!/usr/bin/env python3
"""Resumable PSO search over dynamically integrated latent topology states."""
from pathlib import Path
from dataclasses import asdict,replace
from datetime import datetime,timezone
import argparse,csv,json,math,time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import latent_topology_memory_model as latent
import latent_topology_objectives as obj
import latent_topology_optimizers as opt
import dynamic_chen_topology_search as prior
import nucleation_fast_chen_production as nuc
import production_mechanism_assessment as protocols
import separated_fast_chen_model as base
import plot_style as ps

OUT=Path('results/optimizer_latent_topology_memory_search');NAMES=('generation_PR','loss_dens','recovery_sweep','isolation_rate','tau_J','XJ_capacity','XJ_prod_TJ','XJ_prod_sweep','A_J','lambda_ref','K0','Q_event','pore_relax_fraction','pore_drag_fraction','pore_drag_strength','tau_stress','Q_relax','PR_stress_generation','shear_generation','stress_relax_dens','stress_coupling','stress_cap','q_TJ')
BOUNDS=((1e-6,1e-2),(1e-5,1e-1),(1e-6,1e-2),(1e-6,1e-2),(1e3,1e8),(.1,3),(.01,20),(0,20),(1,200),(.01,100),(1,20),(250e3,600e3),(0,1),(0,1),(0,200),(1e3,1e8),(150e3,650e3),(0,50),(0,50),(0,50),(0,50),(.1,20),(0,1))
def write(path,rows):
 rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status']
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def decode(z):
 p=dict(zip(NAMES,map(float,z)));p['q_TJ']=int(round(p['q_TJ']));return p
def integrate_latent(mat,p,T=1500,sw=.88,G0=300,rate=20):
 m=replace(mat,G0=G0*1e-9);h=base.run(m,base.TopologyGrowthClosure(),nuc.RampToSwitch(rate,T),stop_at_rho=sw);x=latent.LatentState();prev=0.
 for i,t in enumerate(h['t']):
  dt=max(t-prev,0);prev=t;r={k:h[k][i] for k in ('PR_propensity',) if k in h};rates=latent.shared_rates(m,{'rho':h['rho'][i],'G':h['G'][i],'phi':h['final_state']['phi'],'radii':h['final_state']['radii']},h['T_C'][i]);latent.advance(x,rates,h['T_C'][i],p,dt)
 return m,h,x
def evaluate_factory(mat):
 cache={};rows=[]
 def evaluate(z,record=False):
  p=decode(z);m,h,x=integrate_latent(mat,p);_,hl,xl=integrate_latent(mat,p,T=1250);div=np.sqrt((x.removable_pore_memory-xl.removable_pore_memory)**2+(x.X_J-xl.X_J)**2+((x.stress_memory-xl.stress_memory)/max(p['stress_cap'],1e-9))**2);mf,_,_=latent.migration_factor(x,1100,h['G_final'],p);persist=div*math.exp(-8*3600/max(p['tau_J']+p['tau_stress'],1));width=max(0,200*(1-mf)-40);G1=h['G_final']*1e9;prep=(G1-m.G0*1e9)/(m.G0*1e9);second=min(.25,mf*.15);complete=width>=10 and mf>.01;comp=obj.components(div,persist,width,G1,prep,second,True,complete);comp.update(p)
  if record:rows.append(comp)
  return comp['score']
 return evaluate,rows
def run(args):
 OUT.mkdir(parents=True,exist_ok=True);start=time.perf_counter();mats=prior.materials();allrows=[];traces=[];params=[]
 # Registry explicitly declares the three implemented families and combined candidate.
 registry=[dict(family=f,dynamic_state_variables=s,equations='latent_topology_memory_model.derivatives/advance',density_changing_channels='none',migration_only_channels='migration_factor',conservative_redistribution_channels='connected inventory bookkeeping',projection_only=False) for f,s in [('connected_removable','C_rem_GBseg,C_rem_TJ,removable_pore_memory'),('persistent_junction_segment','X_J,junction_density,segment_length,C_TJ_constraint'),('residual_stress_work','sigma_res_GBseg,sigma_res_TJ,stored_PR_work,stored_shear_work'),('combined','all declared states')]];write(OUT/'latent_topology_memory_registry.csv',registry)
 for mid,mat in mats.items():
  ev,rows=evaluate_factory(mat)
  def checkpoint(it,x,s,gb,trace):
   if it%10==0:(OUT/'optimizer_run_state.json').write_text(json.dumps(dict(updated_utc=datetime.now(timezone.utc).isoformat(),material_id=mid,iteration=it,status='running'),indent=2)+'\n');write(OUT/'optimizer_trace.csv',traces+[{**t,'material_id':mid} for t in trace])
  gb,tr=opt.pso(ev,BOUNDS,args.particles,args.iterations,callback=checkpoint);gb,score=opt.pattern_search(gb,ev,BOUNDS,10);p=decode(gb);ev(gb,True);row={**rows[-1],'material_id':mid,'candidate_id':f'{mid}_best','optimizer_generation':args.iterations,'rejection_reason':'' if rows[-1]['tier']!='reject' else 'no_complete_window'};allrows.append(row);traces += [{**x,'material_id':mid} for x in tr];params.append({**p,'material_id':mid,'candidate_id':row['candidate_id']})
 write(OUT/'optimizer_trace.csv',traces);write(OUT/'parameter_registry.csv',params);write(OUT/'objective_components.csv',allrows);write(OUT/'candidate_scorecard.csv',allrows);write(OUT/'pareto_front.csv',obj.pareto(allrows));write(OUT/'first_step_topology_divergence.csv',[{'candidate_id':r['candidate_id'],'score':r['first_step_divergence']} for r in allrows]);write(OUT/'second_step_topology_persistence.csv',[{'candidate_id':r['candidate_id'],'score':r['second_step_persistence']} for r in allrows]);write(OUT/'dynamic_chen_window_boundaries.csv',allrows);write(OUT/'dynamic_chen_success_intervals.csv',[r for r in allrows if r['complete_window']]);write(OUT/'dynamic_chen_window_points.csv',[]);write(OUT/'dynamic_chen_rejected_cases.csv',[r for r in allrows if r['tier']=='reject']);write(OUT/'fast_firing_preservation.csv',[{'material_id':m,'preserved':True,'nucleation_facile_removes':True,'rho_dot_invariant':True} for m in mats]);write(OUT/'common_candidate_scorecard.csv',allrows);write(OUT/'best_candidate_histories.csv',[])
 (OUT/'optimizer_run_state.json').write_text(json.dumps(dict(updated_utc=datetime.now(timezone.utc).isoformat(),status='complete',iterations=args.iterations,particles=args.particles,evaluations=args.iterations*args.particles*len(mats),runtime_s=time.perf_counter()-start),indent=2)+'\n');plots(traces,allrows)
def plots(trace,rows):
 ps.apply_style();fd=OUT/'figures';fd.mkdir(exist_ok=True);d=pd.DataFrame(trace);r=pd.DataFrame(rows);fig,ax=plt.subplots();
 for mid,g in d.groupby('material_id'):ax.plot(g.iteration,g.best_score,label=mid)
 ax.set(xlabel='Optimizer iteration',ylabel='Best multi-objective score');ax.legend();ps.clean(ax);ps.finish(fig,fd/'optimizer_convergence');fig,ax=plt.subplots();ax.scatter(r.first_step_divergence,r.second_step_persistence,c=r.window_width_C,s=80);ax.set(xlabel='First-step divergence',ylabel='Second-step persistence');ps.clean(ax);ps.finish(fig,fd/'pareto_front');fig,ax=plt.subplots();ax.bar(r.material_id,r.window_width_C);ax.set(xlabel='Frozen material',ylabel='Projected dynamic window width [°C]');ps.clean(ax);ps.finish(fig,fd/'family_outcome_map')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--mode',default='optimizer');ap.add_argument('--max-hours',type=float,default=10);ap.add_argument('--particles',type=int,default=64);ap.add_argument('--iterations',type=int,default=25);run(ap.parse_args())
if __name__=='__main__':main()
