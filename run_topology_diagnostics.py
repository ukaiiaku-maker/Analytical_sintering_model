#!/usr/bin/env python3
import argparse,csv
from dataclasses import replace
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import topology_constrained_sintering as m

def main():
    a=argparse.ArgumentParser();a.add_argument('--outdir',default='results/initial');a.add_argument('--target',type=float,default=.90);args=a.parse_args()
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);p=m.Params()
    runs={'slow_0p2':m.run(p,m.RampHold(.2),args.target),'fast_20':m.run(p,m.RampHold(20),args.target),'high_1350':m.run(replace(p,G0=75e-9),m.Iso(1350),args.target),'two_1350_1250':m.run(replace(p,G0=75e-9),m.TwoStep(),args.target)}
    rows=[]
    for name,r in runs.items():
        G,reached=m.value_at_density(r,args.target);rows.append({'protocol':name,'reached_target':reached,'rho_final':r['rho'][-1],'G_at_target_nm':G*1e9,'time_final_h':r['t'][-1]/3600,'median_E_G':float(np.median(r['E_G']))})
    q={x['protocol']:x for x in rows};metrics={'HR_pct':m.percent_gain(q['slow_0p2']['G_at_target_nm'],q['fast_20']['G_at_target_nm']),'TS_pct':m.percent_gain(q['high_1350']['G_at_target_nm'],q['two_1350_1250']['G_at_target_nm'])}
    for name,data in [('protocol_summary.csv',rows),('percent_metrics.csv',[metrics])]:
        with (out/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
    fig,ax=plt.subplots(2,2,figsize=(11,8))
    for name,r in runs.items():ax[0,0].plot(r['t']/3600,r['rho'],label=name);ax[0,1].plot(r['G']*1e9,r['rho'],label=name);ax[1,0].plot(r['rho'],r['E_G'],label=name);ax[1,1].plot(r['rho'],r['f_pore'],label=name)
    for x,title in zip(ax.flat,['Density history','Density-grain trajectory','Trajectory efficiency','Pore-boundary coverage']):x.set_title(title);x.grid(alpha=.25)
    ax[0,0].legend(fontsize=8);fig.tight_layout();fig.savefig(out/'protocol_diagnostics.png',dpi=140);plt.close(fig)
    window=[]
    for G in (50,100,250,500):
        for T in (1150,1250,1350):
            r=m.run(replace(p,rho0=.83,G0=G*1e-9),m.Iso(T,12*3600),.90);window.append({'G0_nm':G,'T_C':T,'rho_final':r['rho'][-1],'activity_median':float(np.median(r['activity'])),'E_G_median':float(np.median(r['E_G']))})
    with (out/'second_step_window.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=window[0].keys());w.writeheader();w.writerows(window)
    print({**metrics,'all_reached':all(x['reached_target'] for x in rows)})
if __name__=='__main__':main()
