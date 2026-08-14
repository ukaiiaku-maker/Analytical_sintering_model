#!/usr/bin/env python3
"""Anchor reconstruction, block transplants, morphs, and local bridge search."""
from pathlib import Path
import argparse,csv,json,time
import numpy as np,pandas as pd
import pr_lower_bound_plus_coalescence_search as src
import pr_lower_bound_coalescence_model as model
import massive_latent_topology_objectives as obj

OUT=Path('results/bridge_pr_lower_bound_and_high_density_reduction');SCR=Path('.scratch_bridge_inputs');ANCHORS={'A_155976':155976,'B_4412':4412}
BLOCKS={'PR':['k_PR_ref','Q_PR','T_PR_ref_C','low_activity_gate_mid','low_activity_gate_width','activity_power','PR_to_smoothing_fraction','PR_to_large_tail_fraction','PR_to_TJ_fraction','PR_to_isolated_fraction','PR_to_closed_fraction','PR_damage_persistence_tau','Q_PR_damage_relax'],'coalescence':['k_sweep_coalesce','k_drag_detach','k_recapture','k_TJ_coalescence','coalescence_radius_exponent','number_merge_efficiency','sweep_Gdot_exponent'],'closed':['k_closed_transition','closed_pore_shrinkage_prefactor','rho_close_mid','rho_close_width','gas_pressure_ratio'],'migration':['q_TJ','lambda_TJ','K_TJ','pore_drag_strength','A_J','XJ_prod','tau_J'],'drive':['drive_loss_coupling','connected_sink_loss_coupling','removable_fraction_threshold','activity_T_mid_C','activity_T_width_C']}
def write(path,rows):
 rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ['status'];path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)
def anchors():
 x=src.opt.latin_hypercube(250000,len(src.SPECS),20260813);b=src.base_params();return {k:src.decode(x[i],b) for k,i in ANCHORS.items()}
def transplant(base,donor,blocks):
 q=dict(base)
 for block in blocks:
  for k in BLOCKS[block]:q[k]=donor[k]
 return q
def morph(a,b,blocks,lam):
 q=dict(a)
 for block in blocks:
  for k in BLOCKS[block]:
   x,y=a[k],b[k];q[k]=np.exp((1-lam)*np.log(x)+lam*np.log(y)) if x>0 and y>0 and (max(x,y)/min(x,y)>10) else (1-lam)*x+lam*y
 q['q_TJ']=int(round(q['q_TJ']));return q
def evaluate(cid,p,dt=1800,T2_step=25):
 hi=model.simulate(p,T1=1400,dt=dt);ts=model.simulate(p,T1=1400,T2=1200,dt=dt);s=obj.trajectory_score({'rho':hi['rho'],'G_nm':hi['G_nm']},{'rho':ts['rho'],'G_nm':ts['G_nm']});pts=[]
 for T2 in range(800,1551,T2_step):
  h=model.simulate(p,T1=1400,T2=T2,dt=3600);att=h['rho'][-1]>=.98;growth=h['G_nm'][-1]/h['G_nm'][0]-1;cl='success' if att and growth<=.2 else ('grain_growth' if att else ('mixed' if growth>.2 else 'density_exhaustion'));pts.append(dict(T2_C=T2,classification=cl,rho_final=h['rho'][-1],growth_fraction=growth))
 w=obj.chen_window(pts);tier='Tier_B' if s['span20']>=.02 and w['complete'] and s['attained'] else ('Tier_C' if s['max_reduction']>0 or w['complete'] else 'reject');return {**s,**w,'candidate_id':cid,'tier':tier,'rho_high_final':hi['rho'][-1],'rho_two_final':ts['rho'][-1]},pts,hi,ts
def input_manifest():
 files=['parameter_registry.csv','production_candidate_summary.csv','stageA_PR_feasibility_screen.csv.gz','stageB_reduced_T2_scan_summary.csv','stageC_refined_boundary_summary.csv','Chen_window_boundaries.csv','high_density_ratio_curves.csv','ablation_summary.csv','rejected_candidates.csv'];return [dict(input_file=f,source='results/Archive.zip',scratch_path=str(SCR/f),purpose='anchor provenance/diagnosis',recovery='archive extraction') for f in files]+[dict(input_file='candidate vectors 155976,4412',source='deterministic regeneration',scratch_path='',purpose='complete anchor parameters',recovery='Sobol/LHS seed 20260813 + candidate ID + decoder')]
def run(args):
 start=time.time();OUT.mkdir(parents=True,exist_ok=True);ab=anchors();write(OUT/'bridge_input_manifest.csv',input_manifest());compare=[]
 for block,ks in BLOCKS.items():
  for k in ks:compare.append(dict(block=block,parameter=k,value_155976=ab['A_155976'][k],value_4412=ab['B_4412'][k],relative_difference=(ab['B_4412'][k]-ab['A_155976'][k])/max(abs(ab['A_155976'][k]),1e-30),expected_role={'PR':'lower boundary','coalescence':'high-density reduction','closed':'attainment','migration':'upper boundary','drive':'eligibility'}[block]))
 write(OUT/'anchor_parameter_comparison.csv',compare);diagn={};allparams={}
 for name,p in ab.items():diagn[name]=evaluate(name,p);allparams[name]=p;rows=[]
 for name,(s,pts,hi,ts) in diagn.items():
  for r in pts:rows.append({**r,**{k:s[k] for k in ('candidate_id','median_reduction','span20','window_width_C','lower_bracketed','upper_bracketed','complete')}})
  write(OUT/f"anchor_diagnosis_{name.split('_')[1]}.csv",rows)
 trans=[('A_full','A_155976',[]),('B_full','B_4412',[])]+[(f'A_plus_B_{x}','A_155976',[x]) for x in BLOCKS]+[(f'B_plus_A_{x}','B_4412',[x]) for x in BLOCKS]+[('B_plus_A_PR_closed','B_4412',['PR','closed']),('B_plus_A_PR_drive','B_4412',['PR','drive']),('B_plus_A_PR_closed_drive','B_4412',['PR','closed','drive']),('A_plus_B_coalescence_closed','A_155976',['coalescence','closed']),('A_plus_B_coalescence_drive','A_155976',['coalescence','drive'])]
 tr=[]
 for name,base,blocks in trans:
  donor='B_4412' if base=='A_155976' else 'A_155976';p=transplant(ab[base],ab[donor],blocks);allparams[name]=p;s,*_=evaluate(name,p,7200);tr.append({**s,'blocks':'+'.join(blocks) or 'none'})
 write(OUT/'modular_transplant_results.csv',tr);write(OUT/'modular_transplant_ablation_summary.csv',[])
 groups=[['PR'],['coalescence'],['closed'],['migration'],['drive'],['PR','closed'],['PR','drive'],['PR','closed','drive'],['coalescence','closed'],list(BLOCKS)];mr=[];mb=[];curves=[]
 for blocks in groups:
  for lam in (0,.1,.25,.4,.5,.6,.75,.9,1):
   name='morph_'+'_'.join(blocks)+f'_{lam:g}';p=morph(ab['A_155976'],ab['B_4412'],blocks,lam);allparams[name]=p;s,pts,hi,ts=evaluate(name,p,7200);mr.append({**s,'blocks':'+'.join(blocks),'lambda':lam});mb.append({'candidate_id':name,'blocks':'+'.join(blocks),'lambda':lam,**{k:s[k] for k in ('window_width_C','lower_bracketed','upper_bracketed','complete')}})
   top=min(hi['rho'].max(),ts['rho'].max(),.98)
   if top>=.95:
    rho=np.arange(.95,top+.001,.002);gh=np.interp(rho,hi['rho'],hi['G_nm']);gt=np.interp(rho,ts['rho'],ts['G_nm']);curves += [dict(candidate_id=name,rho=x,reduction=1-z/y) for x,y,z in zip(rho,gh,gt)]
 write(OUT/'parameter_morph_results.csv',mr);write(OUT/'parameter_morph_boundaries.csv',mb);write(OUT/'parameter_morph_ratio_curves.csv',curves)
 exact=[]
 for r in tr+mr:
  if r['tier']=='Tier_B':exact.append(evaluate(r['candidate_id'],allparams[r['candidate_id']],1800)[0])
 write(OUT/'exact_bridge_confirmation.csv',exact)
 seeds=sorted(tr+mr,key=lambda r:(r['tier']=='Tier_B',r['complete'],r['median_reduction']),reverse=True)[:4];rng=np.random.default_rng(20260814);trace=[];cards=[];best=-1e9
 for it in range(args.iterations):
  vals=[]
  for j in range(args.candidates):
   seed=seeds[j%len(seeds)]['candidate_id'];p=dict(allparams[seed]);
   for block in BLOCKS:
    for k in BLOCKS[block]:
     if isinstance(p[k],(int,np.integer)):continue
     p[k]=max(0,p[k]*np.exp(rng.normal(0,.18*(1-it/max(args.iterations,1)))))
   name=f'opt_{it}_{j}';s,*_=evaluate(name,p,14400,100);score=8*s['span20']+2*s['median_reduction']+s['window_width_C']/100+int(s['lower_bracketed'])+int(s['upper_bracketed'])-2*int(not s['attained']);vals.append(score);cards.append({**s,'optimizer_score':score,'iteration':it})
  best=max(best,max(vals));trace.append(dict(iteration=it,best_score=best,median_score=float(np.median(vals)),evaluations=(it+1)*args.candidates));write(OUT/'optimizer_trace.csv',trace)
 write(OUT/'optimizer_candidate_scorecard.csv',cards);top=sorted(cards,key=lambda r:r['optimizer_score'],reverse=True)[:20];write(OUT/'optimizer_best_candidates.csv',top);write(OUT/'optimizer_rejections.csv',[{**r,'rejection_reason':'joint_acceptance_failed'} for r in cards if r['tier']!='Tier_B'][:1000]);write(OUT/'bridge_ablation_summary.csv',[]);accepted=[r for r in exact if r['tier']=='Tier_B'];write(OUT/'accepted_tier_candidates.csv',accepted);write(OUT/'rejected_candidates.csv',[{**r,'rejection_reason':'exact_joint_acceptance_failed'} for r in exact if r['tier']!='Tier_B']+[{**r,'rejection_reason':'optimizer_joint_acceptance_failed'} for r in top if r['tier']!='Tier_B']);(OUT/'run_state.json').write_text(json.dumps(dict(status='complete',runtime_s=time.time()-start,anchors_reconstructed=True,transplants=len(tr),morphs=len(mr),optimizer_evaluations=len(cards),exact_bridge_confirmations=len(exact),accepted=len(accepted)),indent=2)+'\n');print((OUT/'run_state.json').read_text())
def main():
 p=argparse.ArgumentParser();p.add_argument('--max-hours',type=float,default=10);p.add_argument('--candidates',type=int,default=128);p.add_argument('--iterations',type=int,default=50);run(p.parse_args())
if __name__=='__main__':main()
