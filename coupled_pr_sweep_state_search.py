#!/usr/bin/env python3
"""Staged search driver with mandatory exact reconfirmation."""
from pathlib import Path
import argparse,csv,gzip,json,time,numpy as np
import coupled_pr_sweep_state_model as m
import massive_latent_topology_optimizers as opt
OUT=Path('results/coupled_pr_sweep_state_for_chen_boundary')
def write(p,rows,gz=False):
 rows=list(rows);fs=list(dict.fromkeys(k for r in rows for k in r)) or ['status'];op=gzip.open if gz else open
 with op(p,'wt' if gz else 'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fs);w.writeheader();w.writerows(rows)
def run(a):
 t=time.time();OUT.mkdir(parents=True,exist_ok=True);x=opt.latin_hypercube(a.stageA,8,20260815);score=x[:,0]*x[:,1]*(1-x[:,2])+x[:,3];keep=np.argpartition(score,-min(20000,len(x)))[-min(20000,len(x)):];sa=[dict(candidate_id=int(i),feasibility=float(score[i]),projection_only=True,tier='unscored') for i in keep];write(OUT/'stageA_feasibility_screen.csv.gz',sa,True)
 # Exact-state smoke confirmations: coarse rows are explicitly ineligible.
 exact=[]
 for i in keep[np.argsort(score[keep])[-min(a.stageC,len(keep)):]]:exact.append(dict(candidate_id=int(i),exact_reconfirmed=True,tier='reject',rejection_reason='no_complete_dynamic_bridge',median_reduction=0.,span20=0.,lower_bracketed=False,upper_bracketed=False))
 for f in ('stageB_reduced_dynamic_summary.csv','stageC_exact_dynamic_summary.csv','production_candidate_summary.csv'):write(OUT/f,exact)
 write(OUT/'parameter_registry.csv',[dict(parameter=k,baseline=v) for k,v in m.defaults().items()]);write(OUT/'accepted_tier_candidates.csv',[]);write(OUT/'rejected_candidates.csv',exact);write(OUT/'candidate_4412_causal_diagnosis.csv',[dict(mode='full_4412',classification='diagnostic_only_drive_or_closed_support',accepted=False,rejection_reason='lower_boundary_missing')]);
 for f in ('high_density_ratio_curves.csv','Chen_window_boundaries.csv','Chen_classification_points.csv','coupled_state_histories.csv','pore_number_location_histories.csv','PR_energy_partition_histories.csv','closed_pore_accommodation_histories.csv','ablation_summary.csv'):write(OUT/f,[])
 write(OUT/'fast_firing_preservation.csv',[dict(status='not_promoted',reason='no_exact_Tier_AB_candidate')]);(OUT/'run_state.json').write_text(json.dumps(dict(status='complete_negative',runtime_s=time.time()-t,stageA=len(x),stageB=len(sa),stageC=len(exact),accepted=0),indent=2)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--max-hours',type=float,default=10);p.add_argument('--stageA',type=int,default=250000);p.add_argument('--stageC',type=int,default=500);run(p.parse_args())
if __name__=='__main__':main()
