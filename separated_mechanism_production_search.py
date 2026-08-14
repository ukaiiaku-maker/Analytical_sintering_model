#!/usr/bin/env python3
"""Checkpointed causal production audit for the separated reduced model."""
from __future__ import annotations
from dataclasses import asdict,replace
from datetime import datetime,timezone
from pathlib import Path
import argparse,csv,json,math
import numpy as np
import pandas as pd

import separated_fast_chen_model as model
import separated_fast_chen_search as discovery
import production_mechanism_assessment as protocols
import observable_trajectory_effect_audit as observable

OUT=Path("results/separated_mechanism_production_search")
MODES=("full_material_model","no_PR_redistribution","no_nucleation_limitation","transport_only","no_surface_PR_low_activity_gate","PR_only_no_densification_activation","exchange_limited_variant")

def write(path,rows):
    rows=list(rows);fields=list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def state(phase,completed,best=None,failures=0,note=""):
    data=dict(updated_utc=datetime.now(timezone.utc).isoformat(),current_phase=phase,completed_material_sets=completed,completed_topology_sets=0,checkpoints_written=completed//1,best_candidates_so_far=best or [],failures_rejections_so_far=failures,estimated_remaining_work=note)
    (OUT/"run_state.json").write_text(json.dumps(data,indent=2)+"\n")

def seed_materials():
    old=pd.read_csv("results/separate_fast_firing_and_chen_mechanisms/fast_firing_successful_material_sets.csv")
    ids=sorted(old.material_id.unique());registry=dict(discovery.materials(256));return [(mid,registry[mid]) for mid in ids]

def metrics(ref,fast):
    if ref["numerical_censored"] or fast["numerical_censored"]:return dict(attained=False,rejection_reason="numerical_censor")
    c=discovery.curve(ref,fast)
    attained=len(c)>1 and c.rho.min()<=.751 and c.rho.max()>=.919
    if not attained:return dict(attained=False,rejection_reason="unattained_interval")
    rho=c.rho.to_numpy();ratio=c.ratio.to_numpy();span15=observable.longest_span(rho,ratio,1.5)
    positive=np.maximum(np.log(np.maximum(ratio,1e-300)),0)
    return dict(attained=True,max_ratio=float(ratio.max()),median_ratio=float(np.median(ratio)),span_ge_1p2=observable.longest_span(rho,ratio,1.2),span_ge_1p5=span15,span_ge_2p0=observable.longest_span(rho,ratio,2),integrated_separation=float(np.trapezoid(positive,rho)),effect_rho_min=float(rho[ratio>=1.5].min()) if np.any(ratio>=1.5) else math.nan,effect_rho_max=float(rho[ratio>=1.5].max()) if np.any(ratio>=1.5) else math.nan,meaningful=span15>=.03,extreme_warning=bool(ratio.max()>10),promotion_blocked=bool(ratio.max()>100),rejection_reason="" if span15>=.03 else "ratio_or_span_below_threshold",curve=c)

def ablations():
    OUT.mkdir(parents=True,exist_ok=True);summary=[];curves=[];histories=[];reject=[];seeds=seed_materials();state("phase_A",0,note=f"{len(seeds)} seed materials")
    for j,(mid,p0) in enumerate(seeds,1):
      # Re-run the schedule that produced the seed's strongest non-extreme hit.
      old=pd.read_csv("results/separate_fast_firing_and_chen_mechanisms/fast_firing_successful_material_sets.csv");q=old[old.material_id==mid].sort_values("max_ratio");q=q[q.max_ratio<=10] if np.any(q.max_ratio<=10) else q;best=q.iloc[-1]
      for mode in MODES:
        p=replace(p0,ablation_mode=mode);ref=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(1,float(best.peak_T_C),float(best.hold_h)));fast=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(float(best.fast_rate_C_min),float(best.peak_T_C),float(best.hold_h)));m=metrics(ref,fast);curve=m.pop("curve",None);row=dict(material_id=mid,ablation_mode=mode,peak_T_C=best.peak_T_C,hold_h=best.hold_h,fast_rate_C_min=best.fast_rate_C_min,**m);summary.append(row)
        if row.get("rejection_reason"):reject.append(row)
        if curve is not None:
          for x in curve.iloc[::max(1,len(curve)//200)].to_dict("records"):curves.append(dict(material_id=mid,ablation_mode=mode,**x))
        for path,h in (("reference",ref),("fast",fast)):
          stride=max(1,len(h["rho"])//120)
          for i in range(0,len(h["rho"]),stride):histories.append(dict(material_id=mid,ablation_mode=mode,path=path,t_h=h["t"][i]/3600,T_C=h["T_C"][i],rho=h["rho"][i],G_nm=h["G"][i]*1e9,tau_nuc=h["tau_nuc"][i],tau_exchange=h["tau_exchange"][i],tau_transport=h["tau_transport"][i],activity=h["activity"][i],PR_exposure=h["PR_exposure"][i],pore_D50_nm=h["pore_D50"][i]*1e9,pore_D90_nm=h["pore_D90"][i]*1e9,connected_fine=h["connected_fine"][i]))
      write(OUT/"fast_ablation_summary.csv",summary);write(OUT/"fast_ablation_ratio_curves.csv",curves);write(OUT/"fast_ablation_state_histories.csv",histories);write(OUT/"fast_ablation_rejections.csv",reject);survive=[r["material_id"] for r in summary if r["ablation_mode"]=="full_material_model" and r.get("meaningful") and not r.get("promotion_blocked")];state("phase_A",j,survive,len(reject),"Phase A complete" if j==len(seeds) else f"{len(seeds)-j} seeds remain")
    return summary

def compact_summary(rows):
    df=pd.DataFrame(rows);out=[]
    for mid,g in df.groupby("material_id"):
        q={r.ablation_mode:r for _,r in g.iterrows()};full=q["full_material_model"]
        causal=bool(full.get("meaningful",False) and not q["no_PR_redistribution"].get("meaningful",False) and not q["no_nucleation_limitation"].get("meaningful",False))
        out.append(dict(material_id=mid,full_meaningful=full.get("meaningful",False),full_max_ratio=full.get("max_ratio",math.nan),full_span_ge_1p5=full.get("span_ge_1p5",0),no_PR_meaningful=q["no_PR_redistribution"].get("meaningful",False),no_nucleation_meaningful=q["no_nucleation_limitation"].get("meaningful",False),transport_only_meaningful=q["transport_only"].get("meaningful",False),exchange_limited_meaningful=q["exchange_limited_variant"].get("meaningful",False),survives_causal_gate=causal,plausibility_flag="stiff" if full.get("extreme_warning",False) else "prototype_scale"))
    write(OUT/"causal_survivor_scorecard.csv",out);return out

def expanded_materials(n):
    levels=dict(Q_GB_diffusion=(350,425,500,575,650),Q_surface_diffusion=(250,325,400,475,550,625,700),Q_disconnection_nucleation=(350,400,450,500,550,600,650),Q_exchange=(180,245,325,425),Q_transport=(280,360,450,550),v_star=(3e-29,8e-29,2e-28,5e-28),stress_concentration=(1,3,10,30,60),PR_prefactor=(2e-6,2e-5,2e-4,2e-3),PR_partition=("smoothing","GB_to_TJ","isolation","balanced","large_pore_tail"),zeta_eta_ratio=(.25,.5,1,2),pore_ln_sigma=(.35,.65,.95),pore_radius0=(15,22,35,50),G0=(50,75,100,150,225,300),rho0=(.65,.70,.75));primes=(1,3,5,7,11,13,17,19,23,29,31,37,41,43,47);out=[]
    for i in range(n):
        kw={};
        for (name,vals),prime in zip(levels.items(),primes):
            value=vals[(i*prime+i//max(len(vals),1))%len(vals)]
            if name.startswith("Q_"):value*=1e3
            if name=="pore_radius0" or name=="G0":value*=1e-9
            kw[name]=value
        out.append((f"E{i:04d}",model.MaterialKinetics(**kw)))
    return out

def expanded_screen(n=512):
    registry=[];screen=[];near=[];reject=[];ratios=[];schedules=((1450,8),(1450,20),(1550,8),(1550,20));state("phase_B",0,note=f"broad material screen n={n}")
    for j,(mid,p) in enumerate(expanded_materials(n),1):
        registry.append(dict(material_id=mid,**asdict(p)));best=None
        for peak,hold in schedules:
            ref=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(1,peak,hold))
            for rate in (20,50,100):
                fast=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(rate,peak,hold));m=metrics(ref,fast);c=m.pop("curve",None);row=dict(material_id=mid,peak_T_C=peak,hold_h=hold,fast_rate_C_min=rate,**m);screen.append(row)
                if c is not None and (best is None or row.get("span_ge_1p5",0)>best.get("span_ge_1p5",0)):best=row
        hit=best is not None and best.get("meaningful",False) and not best.get("promotion_blocked",False)
        causal=False
        if hit:
            checks={}
            for mode in ("no_PR_redistribution","no_nucleation_limitation"):
                pa=replace(p,ablation_mode=mode);r=model.run(pa,model.TopologyGrowthClosure(),protocols.FastSchedule(1,best["peak_T_C"],best["hold_h"]));f=model.run(pa,model.TopologyGrowthClosure(),protocols.FastSchedule(best["fast_rate_C_min"],best["peak_T_C"],best["hold_h"]));am=metrics(r,f);am.pop("curve",None);checks[mode]=am
            causal=not checks["no_PR_redistribution"].get("meaningful",False) and not checks["no_nucleation_limitation"].get("meaningful",False)
            best.update(no_PR_meaningful=checks["no_PR_redistribution"].get("meaningful",False),no_nucleation_meaningful=checks["no_nucleation_limitation"].get("meaningful",False),survives_causal_gate=causal)
        if best:
            target=near if (hit and not causal) or best.get("max_ratio",0)>=1.2 else reject;target.append({**best,**asdict(p),"rejection_reason":"" if causal else ("PR_ablation_did_not_collapse" if hit else best.get("rejection_reason","below_threshold"))})
        if j%25==0 or j==n:
            write(OUT/"material_parameter_registry.csv",registry);write(OUT/"fast_firing_material_screen.csv",screen);write(OUT/"fast_firing_successful_material_sets.csv",[r for r in near if r.get("survives_causal_gate")]);write(OUT/"fast_firing_near_hit_material_sets.csv",near);write(OUT/"fast_firing_rejected_material_sets.csv",reject);state("phase_B",j,[r["material_id"] for r in near if r.get("survives_causal_gate")],len(reject),f"{n-j} material sets remain")
    return [r for r in near if r.get("survives_causal_gate")]

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--phase",choices=("ablations","expanded"),default="ablations");ap.add_argument("--materials",type=int,default=512);a=ap.parse_args()
    if a.phase=="ablations":
        rows=ablations();survive=compact_summary(rows);state("phase_A_complete",len(seed_materials()),[r["material_id"] for r in survive if r["survives_causal_gate"]],sum(not r["survives_causal_gate"] for r in survive),"Expanded material and dynamic Chen phases pending causal survivors");print(pd.DataFrame(survive).to_string(index=False))
    else:
        hits=expanded_screen(a.materials);state("phase_B_complete",a.materials,[r["material_id"] for r in hits],0,"Dynamic Chen blocked if no causal survivors");print("causal survivors",len(hits))

if __name__=="__main__":main()
