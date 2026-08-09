#!/usr/bin/env python3
import argparse,csv,random
from dataclasses import replace
from pathlib import Path
import numpy as np
import topology_constrained_sintering as m
def main():
    a=argparse.ArgumentParser();a.add_argument('--n',type=int,default=12);a.add_argument('--seed',type=int,default=17);a.add_argument('--target',type=float,default=.90);a.add_argument('--out',default='results/initial/search_summary.csv');x=a.parse_args();rng=random.Random(x.seed);rows=[]
    for i in range(x.n):
        p=m.Params() if i==0 else m.Params(Q_growth=rng.uniform(385e3,455e3),Q_nucleation=rng.uniform(410e3,470e3),pore_drag_resistance=10**rng.uniform(-.2,1),tl_drag_resistance=10**rng.uniform(-.2,1),connectivity_rho_mid=rng.uniform(.85,.90)); ps=[p,p,replace(p,G0=75e-9),replace(p,G0=75e-9)];protos=[m.RampHold(.2),m.RampHold(20),m.Iso(1350),m.TwoStep()];rr=[m.run(pp,pr,x.target) for pp,pr in zip(ps,protos)];v=[m.value_at_density(r,x.target) for r in rr];reached=all(ok for _,ok in v);g=[z*1e9 for z,_ in v];hr=m.percent_gain(g[0],g[1]);ts=m.percent_gain(g[2],g[3]);eg=float(np.median(np.concatenate([r['E_G'] for r in rr])));rows.append({'sample':i,'all_reached':reached,'HR_pct':hr,'TS_pct':ts,'median_E_G':eg,'combined_score':min(hr,ts)+eg if reached else -1e9,'Q_growth_kJ':p.Q_growth/1e3,'Q_nucleation_kJ':p.Q_nucleation/1e3,'pore_drag_resistance':p.pore_drag_resistance,'tl_drag_resistance':p.tl_drag_resistance,'connectivity_rho_mid':p.connectivity_rho_mid})
    rows.sort(key=lambda z:z['combined_score'],reverse=True);path=Path(x.out);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    print('best',rows[0])
if __name__=='__main__':main()
