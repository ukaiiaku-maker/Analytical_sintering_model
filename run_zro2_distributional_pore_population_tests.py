#!/usr/bin/env python3
"""Run bounded synthetic, heating-rate, topology, and ablation forward tests."""
from __future__ import annotations
from dataclasses import replace
import json, math
from pathlib import Path
import runpy, sys, traceback
import numpy as np
import pandas as pd

import build_zro2_distributional_pore_population_model as model

OUT=model.OUT


def prepared_state(name,representation):
    specs={
      "narrow_fine_connected":dict(family="lognormal",D50_nm=15,sigma_ln=.30),
      "broad_fine":dict(family="lognormal",D50_nm=15,sigma_ln=.85),
      "bimodal_agglomerate_tail":dict(family="bimodal",D50_nm=15,sigma_ln=.45,tail_weight=.25,D50_2_nm=150),
      "high_precursor_low_closed":dict(family="bimodal",D50_nm=24.5,sigma_ln=.65,tail_weight=.10),
      "high_closed_low_accommodation":dict(family="lognormal",D50_nm=24.5,sigma_ln=.45),
      "high_useful_closed":dict(family="bimodal",D50_nm=15,sigma_ln=.45,tail_weight=.05),
      "gas_limited_closed":dict(family="bimodal",D50_nm=15,sigma_ln=.45,tail_weight=.05),
      "candidate_response_target_like":dict(family="bimodal",D50_nm=15,sigma_ln=.65,tail_weight=.10),
    }
    s=model.initial_state(rho0=.82,G0_nm=65,**specs[name]);s.representation=representation
    total=s.total_pore;shape=s.phi[0]/max(s.phi[0].sum(),model.EPS)
    target={
      "high_precursor_low_closed":(.055,.110,.010,.005,.20),
      "high_closed_low_accommodation":(.045,.010,.005,.120,.02),
      "high_useful_closed":(.040,.010,.010,.120,.70),
      "gas_limited_closed":(.040,.010,.010,.120,.70),
      "candidate_response_target_like":(.015,.020,.015,.130,.85),
    }
    if name in target:
        vals=target[name];s.phi[:]=np.array(vals[:4])[:,None]*shape[None,:]
        # Keep exact rho identity; any unassigned void is connected open.
        s.phi[0]+=(total-s.phi.sum())*shape;s.closed_A[:]=vals[4]
    if name=="candidate_response_target_like":s.PR_memory=.8;s.candidate_response_target_only=True;s.provenance="candidate_response_target_only"
    if name=="gas_limited_closed":s.provenance="synthetic_forward_state_gas_limited"
    s.rho=1-s.phi.sum();s.N=s.phi/((4/3)*math.pi*s.radii_m**3)
    return s


def window_summary(points):
    rows=[]
    for keys,g in points.groupby(["state_id","representation","closed_kernel","m"]):
        g=g.sort_values("T2_C");succ=g[g.classification.eq("SUCCESS")]
        lower=bool(len(succ) and ((g.T2_C<succ.T2_C.min())&g.classification.eq("DENSIFICATION_EXHAUSTION_FAILURE")).any())
        upper=bool(len(succ) and ((g.T2_C>succ.T2_C.max())&g.classification.isin(["GRAIN_GROWTH_FAILURE","MIXED_FAILURE"])).any())
        finite=bool(len(succ)>1 and lower and upper and succ.T2_C.max()>succ.T2_C.min())
        rows.append({"state_id":keys[0],"representation":keys[1],"closed_kernel":keys[2],"m":keys[3],
                     "success_points":len(succ),"lower_boundary":lower,"upper_boundary":upper,"strict_finite_window":finite,
                     "lower_success_T2_C":succ.T2_C.min() if len(succ) else np.nan,"upper_success_T2_C":succ.T2_C.max() if len(succ) else np.nan,
                     "window_width_C":succ.T2_C.max()-succ.T2_C.min() if len(succ) else 0,"cloned_first_step_state":True})
    return pd.DataFrame(rows)


def synthetic_boundary_test():
    names=["narrow_fine_connected","broad_fine","bimodal_agglomerate_tail","high_precursor_low_closed","high_closed_low_accommodation","high_useful_closed","gas_limited_closed","candidate_response_target_like"]
    rows=[]
    for name in names:
      for rep in ("lognormal","bimodal","discrete_bin"):
       initial=prepared_state(name,rep);fingerprint=json.dumps({"rho":initial.rho,"G":initial.G_m,"phi":initial.phi.round(14).tolist(),"A":initial.closed_A.round(12).tolist()},sort_keys=True)
       for kernel in ("renewal","GB_diffusion"):
        for m in (3,4):
         gas=.9 if name=="gas_limited_closed" else .25
         p=model.PopulationParameters(closed_kernel=kernel,radius_exponent=m,gas_fraction=gas)
         for T2 in range(850,1351,25):
          final,_=model.evolve(initial.clone(),T2,96*3600,p,dt_max=7200,record_every=1e99)
          cls,growth=model.classify(final,initial);met=model.metrics(final)
          rows.append({"state_id":name,"representation":rep,"closed_kernel":kernel,"m":m,"T2_C":T2,"hold_h":96,
                       "initial_state_fingerprint":fingerprint,"candidate_response_target_only":initial.candidate_response_target_only,
                       "initial_rho":initial.rho,"initial_G_nm":initial.G_m*1e9,"final_rho":final.rho,"final_G_nm":final.G_m*1e9,
                       "growth_fraction":growth,"classification":cls,"strict_success":cls=="SUCCESS",**met})
    q=pd.DataFrame(rows);q.to_csv(OUT/"synthetic_distribution_boundary_test.csv",index=False)
    w=window_summary(q);w.to_csv(OUT/"synthetic_distribution_window_boundaries.csv",index=False);return q,w


def heating_rate_tests():
    histories=[];summary=[]
    for family in ("lognormal","bimodal","discrete_bin"):
      for peak in (1300,1350,1400,1450,1500):
       for rate in (.2,1,5,10,20,50,100):
        for hold in (0,1,2,8,20):
         s=model.initial_state(rho0=.66,G0_nm=24.5,family=family,D50_nm=24.5,sigma_ln=.45,tail_weight=.10)
         s.representation=family;p=model.PopulationParameters(closed_kernel="GB_diffusion",radius_exponent=3)
         run=f"{family}:{peak}:{rate}:{hold}";t_global=0.
         for TC in np.arange(950,peak+1e-9,10):
          dt=10/rate*60;s,h=model.evolve(s,float(TC),dt,p,dt_max=dt,record_every=0)
          if len(h):
           row=h.iloc[-1].to_dict();row.update(run_id=run,family=family,peak_C=peak,rate_C_min=rate,hold_h=hold,thermal_phase="heating");histories.append(row)
         if hold>0:
          s,h=model.evolve(s,peak,hold*3600,p,dt_max=1800,record_every=3600)
          for _,x in h.iterrows():
           row=x.to_dict();row.update(run_id=run,family=family,peak_C=peak,rate_C_min=rate,hold_h=hold,thermal_phase="hold");histories.append(row)
         met=model.metrics(s);summary.append({"run_id":run,"family":family,"peak_C":peak,"rate_C_min":rate,"hold_h":hold,**met})
    h=pd.DataFrame(histories);s=pd.DataFrame(summary);h.to_csv(OUT/"distribution_heating_rate_histories.csv",index=False);s.to_csv(OUT/"distribution_heating_rate_state_summary.csv",index=False)
    ratios=[]
    targets=(.75,.80,.85,.88,.90,.92,.95)
    for (fam,peak,hold),g in h.groupby(["family","peak_C","hold_h"]):
      ref=g[g.rate_C_min.eq(.2)];fast=g[g.rate_C_min.eq(100)]
      for target in targets:
       def at(x):
        if len(x)==0 or x.rho.max()<target:return None
        return x.iloc[(x.rho-target).abs().argmin()]
       a,b=at(ref),at(fast)
       ratios.append({"family":fam,"peak_C":peak,"hold_h":hold,"rho_target":target,"both_attain":a is not None and b is not None,
                      "G_reference_over_G_fast":a.G_nm/b.G_nm if a is not None and b is not None else np.nan,
                      "D50_difference_nm":a.D50_nm-b.D50_nm if a is not None and b is not None else np.nan,
                      "D90_difference_nm":a.D90_nm-b.D90_nm if a is not None and b is not None else np.nan,
                      "large_tail_difference":a.large_pore_tail_fraction-b.large_pore_tail_fraction if a is not None and b is not None else np.nan,
                      "connected_fine_difference":a.connected_fine_pore_fraction-b.connected_fine_pore_fraction if a is not None and b is not None else np.nan,
                      "useful_closed_difference":a.useful_closed_inventory-b.useful_closed_inventory if a is not None and b is not None else np.nan,
                      "Zener_difference":a.Zener_pinning_metric_minv-b.Zener_pinning_metric_minv if a is not None and b is not None else np.nan})
    r=pd.DataFrame(ratios)
    for threshold in (1.2,1.5,2.0):r[f"density_span_ratio_above_{str(threshold).replace('.','p')}"]=r.groupby(["family","peak_C","hold_h"])["G_reference_over_G_fast"].transform(lambda x: float((x>=threshold).sum())/len(targets))
    r.to_csv(OUT/"distribution_heating_rate_ratios.csv",index=False);return h,r,s


def ablations(points):
    gap=(.90-points.final_rho).clip(lower=0)+(points.growth_fraction-.10).clip(lower=0)
    best=points.iloc[gap.argmin()];base=prepared_state(best.state_id,best.representation)
    definitions={
      "baseline":{},"no_regularization":{"regularization_enabled":False},"no_damaging_coarsening":{"damage_enabled":False},
      "no_PR_pinch_off":{"pinch_enabled":False},"no_precursor_to_closed":{"transition_enabled":False},"no_closed_shrinkage":{"closed_shrinkage_enabled":False},
      "no_accommodation_recovery":{"accommodation_recovery_enabled":False},"infinite_accommodation":{"infinite_accommodation":True},
      "no_gas_pressure":{"gas_enabled":False},"high_gas_pressure":{"gas_fraction":1.1},"no_distributional_Zener":{"distributional_zener":False},
      "mean_radius_Zener_only":{"mean_radius_zener_only":True},"no_energy_ledger_coupling":{"energy_ledger_coupling":False},
      "strict_GB_area_only_power_balance_diagnostic":{"energy_ledger_coupling":False},
    }
    histories=[];rows=[]
    for name,kw in definitions.items():
      p=model.PopulationParameters(closed_kernel=best.closed_kernel,radius_exponent=int(best.m),**kw)
      final,h=model.evolve(base.clone(),best.T2_C,96*3600,p,dt_max=7200,record_every=4*3600)
      cls,growth=model.classify(final,base)
      rows.append({"ablation":name,"source_state":best.state_id,"T2_C":best.T2_C,"final_rho":final.rho,"final_G_nm":final.G_m*1e9,"growth_fraction":growth,"classification":cls,
                   "density_gain":final.rho-base.rho,"closed_density_gain":final.cumulative_closed_shrink,"PR_memory":final.PR_memory,**model.metrics(final)})
      h["ablation"]=name;histories.append(h)
    q=pd.DataFrame(rows);hh=pd.concat(histories,ignore_index=True);q.to_csv(OUT/"distribution_ablation_matrix.csv",index=False);hh.to_csv(OUT/"distribution_ablation_histories.csv",index=False);return q,hh


def write_process_outputs(points,windows):
    # A broad process grid is gated on a strict synthetic finite window.
    justified=bool(windows.strict_finite_window.any())
    if justified:
        p=points.copy();p["map_scope"]="bounded_synthetic_state_map_after_gate"
    else:
        p=points.copy();p["map_scope"]="synthetic_failure_modes_only_process_grid_not_justified"
    p.to_csv(OUT/"distribution_process_map_points.csv",index=False)
    windows.to_csv(OUT/"distribution_window_boundaries.csv",index=False)
    p.sort_values(["strict_success","final_rho","growth_fraction"],ascending=[False,False,True]).head(50).to_csv(OUT/"distribution_best_cases.csv",index=False)
    windows.assign(process_map_run=justified).to_csv(OUT/"distribution_failure_modes.csv",index=False)
    return justified


def main():
    points,windows=synthetic_boundary_test();justified=write_process_outputs(points,windows)
    heat,ratios,states=heating_rate_tests();ab,abh=ablations(points)
    state=json.loads((OUT/"run_state.json").read_text());state.update({"synthetic_points":len(points),"synthetic_successes":int(points.strict_success.sum()),
      "synthetic_finite_windows":int(windows.strict_finite_window.sum()),"broad_process_grid_run":justified,"heating_paths":int(states.shape[0]),
      "ablation_cases":len(ab),"final_interpretation":"bounded_semi_phenomenological_candidate" if windows.strict_finite_window.any() else "diagnostic_negative_result",
      "validation":False});(OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n");print(state)
    pd.DataFrame([{"final_interpretation":state["final_interpretation"],"limiting_tests_pass":state["all_conservation_tests_pass"],
      "synthetic_successes":state["synthetic_successes"],"finite_Chen_windows":state["synthetic_finite_windows"],
      "broad_process_grid_run":state["broad_process_grid_run"],"accepted_model_physics_changed":False,"physical_Q_closed_introduced":False,
      "dominant_failure":"successes begin at minimum T2; no lower boundary; heating-rate sign mostly reversed; distributional Zener/regularization/damage noncausal in best ablation",
      "promoted":False,"validation_claim":False}]).to_csv(OUT/"distribution_promotion_decision.csv",index=False)

def direct_guardrails():
    ns=runpy.run_path(str(model.ROOT/"tests/test_zro2_forward_distributional_pore_population_model.py"));rows=[]
    for name in sorted(x for x in ns if x.startswith("test_")):
        try:ns[name]();rows.append({"test":name,"passed":True,"error":""})
        except Exception as exc:rows.append({"test":name,"passed":False,"error":f"{type(exc).__name__}: {exc}"})
    q=pd.DataFrame(rows);q.to_csv(OUT/"direct_guardrail_results.csv",index=False);print(q.to_string(index=False));print({"direct_guardrails_pass":bool(q.passed.all()),"count":len(q)})
    if not q.passed.all():raise SystemExit(1)

if __name__=="__main__":
    direct_guardrails() if "--direct-guardrails" in sys.argv else main()
