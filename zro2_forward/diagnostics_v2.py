from __future__ import annotations
from dataclasses import asdict, replace
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .barrier_json import BarrierModel
from .campaign import BARRIER
from .integrator import ForwardModel, ModelParameters
from .material_zro2 import MaterialParameters
from .schedules import RampNoHold, TwoStep

OUT=Path("results/zro2_forward_diagnostic_calibration_v2")
TARGET=Path("data/targets/mazaheri_8ysz_2008")
MODES=["nearest_slice_clamp","pchip_extrapolate","fixed_lowT_slope","generic_anchor_barrier"]

def build(mode="nearest_slice_clamp",params=None,M0=5.8e-3):
    barrier=BarrierModel.load(BARRIER).with_mode(mode)
    return ForwardModel(barrier,MaterialParameters(M0_m4_J_s=M0),params or ModelParameters())

def run(model,path,dt=180,record=180,max_steps=20000):
    state,rows=model.run(path,dt_s=dt,record_every_s=record,max_steps=max_steps)
    f=pd.DataFrame(rows); f["G_um"]=f.G_m*1e6; f["T_C"]=f.T_K-273.15
    return state,f

def density_metrics(frame,target=None):
    target=target if target is not None else pd.read_csv(TARGET/"density_vs_temperature_digitized.csv").query("method=='CS'")
    g=frame.sort_values("T_C").drop_duplicates("T_C")
    prediction=np.interp(target.T_C,g.T_C,g.rho)
    onset_target=float(target.loc[target.fractional_density>=.65,"T_C"].min())
    eligible=g[g.rho>=.65]; onset_model=float(eligible.T_C.min()) if len(eligible) else np.nan
    return {"density_curve_rmse":float(np.sqrt(np.mean((prediction-target.fractional_density)**2))),
            "density_curve_mae":float(np.mean(abs(prediction-target.fractional_density))),
            "onset_target_C":onset_target,"onset_model_C":onset_model,"onset_error_C":onset_model-onset_target}

def coarse_map(model):
    rows=[]
    for T2 in [1050,1150,1250,1300]:
      for switch in [.75,.90]:
        path=TwoStep(1450,T2,switch,20); state,_=run(model,path,600,1e9)
        complete=state.t_s>=path.t_end_s-1
        density_ok=abs(state.rho-.976)<=.01; growth_ok=abs(state.G_m*1e6-.29)<=.10
        rows.append({"T1_C":1450,"T2_C":T2,"switch_density":switch,"hold_h":20,
                     "final_rho":state.rho,"final_G_um":state.G_m*1e6,"density_ok":density_ok,
                     "growth_ok":growth_ok,"integration_complete":complete,"chen_success":complete and density_ok and growth_ok})
    return pd.DataFrame(rows)

def baseline_summary():
    root=Path("results/zro2_forward_natural_pore_evolution_target_8ysz")
    h=pd.read_csv(root/"dense_histories.csv"); cs=h[h.case.eq("CS_thermal")]; rapid=h[h.case.eq("rate50_thermal")]
    chen=pd.read_csv(root/"chen_classification_points.csv"); bounds=pd.read_csv(root/"chen_window_boundaries.csv")
    metric=density_metrics(cs); counts=chen.apply(lambda x:"success" if x.chen_success else "density_only" if x.density_ok else "growth_only" if x.growth_ok else "neither",axis=1).value_counts()
    rows=[("CS_final_density",cs.iloc[-1].rho,"dimensionless"),("CS_final_grain_size",cs.iloc[-1].G_um,"um"),
          ("CS_density_onset_error",metric["onset_error_C"],"C"),("rate50_final_density",rapid.iloc[-1].rho,"dimensionless"),
          ("rate50_final_grain_size",rapid.iloc[-1].G_um,"um"),("HMS_density_miss",rapid.iloc[-1].rho-.98,"dimensionless"),
          ("chen_success_count",int(chen.chen_success.sum()),"cases"),("chen_density_only_count",int(counts.get("density_only",0)),"cases"),
          ("chen_growth_only_count",int(counts.get("growth_only",0)),"cases"),("chen_neither_count",int(counts.get("neither",0)),"cases"),
          ("lower_boundary_min_C",bounds.lower_density_exhaustion_max_T2_C.min(),"C"),("lower_boundary_max_C",bounds.lower_density_exhaustion_max_T2_C.max(),"C"),
          ("upper_boundary_min_C",bounds.upper_growth_min_T2_C.min(),"C"),("upper_boundary_max_C",bounds.upper_growth_min_T2_C.max(),"C"),
          ("finite_success_band",bool(chen.chen_success.any()),"boolean"),("M0",5.8e-3,"m4_J_s"),("site_density_multiplier",39.5,"dimensionless"),
          ("barrier_mode","nearest_slice_clamp","label")]
    out=pd.DataFrame(rows,columns=["metric","value","units"]); out.to_csv(OUT/"baseline_failure_summary.csv",index=False)

def barrier_audit():
    comparisons=[]; curves=[]; fast=[]; counts=[]
    for mode in MODES:
        m=build(mode); cs,h=run(m,RampNoHold(5,1500)); rf,hf=run(m,RampNoHold(50,1500),60,60); two,_=run(m,TwoStep(1450,1150,.85,20),600,1e9)
        met=density_metrics(h); comparisons.append({"barrier_mode":mode,"CS_final_rho":cs.rho,"CS_final_G_um":cs.G_m*1e6,
            "rate50_final_rho":rf.rho,"rate50_final_G_um":rf.G_m*1e6,"representative_two_level_rho":two.rho,"representative_two_level_G_um":two.G_m*1e6,
            "JSON_final_evidence":mode!="generic_anchor_barrier","lowT_policy_flagged":True})
        curves.append({"barrier_mode":mode,**met}); fast.append({"barrier_mode":mode,"final_rho":rf.rho,"final_G_um":rf.G_m*1e6,
            "density_error_vs_HMS":rf.rho-.98,"smaller_than_CS":rf.G_m<cs.G_m})
        cm=coarse_map(m); counts.append({"barrier_mode":mode,"success_count":int(cm.chen_success.sum()),"density_ok_count":int(cm.density_ok.sum()),"growth_ok_count":int(cm.growth_ok.sum()),"cases":len(cm)})
        h.assign(barrier_mode=mode).to_csv(OUT/f"history_CS_{mode}.csv",index=False)
    pd.DataFrame(comparisons).to_csv(OUT/"barrier_mode_comparison.csv",index=False)
    pd.DataFrame(curves).to_csv(OUT/"barrier_mode_CS_curve_metrics.csv",index=False)
    pd.DataFrame(fast).to_csv(OUT/"barrier_mode_fast_rate_metrics.csv",index=False)
    pd.DataFrame(counts).to_csv(OUT/"barrier_mode_chen_counts.csv",index=False)

def calibration_audit():
    target=pd.read_csv(TARGET/"density_vs_temperature_digitized.csv").query("method=='CS'"); rows=[]; residuals=[]
    candidates=[("endpoint_only",39.5,1.,5.8e-3)]
    trials=[]
    for site in [10.,40.,160.]:
      for work in [.1,1.,10.]: trials.append((site,work,5.8e-3))
    scored=[]
    for site,work,M0 in trials:
        q=replace(ModelParameters(),site_density_multiplier=site,stress_work_factor=work); s,h=run(build(params=q,M0=M0),RampNoHold(5,1500))
        met=density_metrics(h,target); scored.append((met["density_curve_rmse"],site,work,M0,s,h,met))
    best=min(scored,key=lambda x:x[0]); candidates.append(("density_curve_only",best[1],best[2],best[3]))
    joint=[]
    for site,work,_ in [(x[1],x[2],x[3]) for x in sorted(scored)[:4]]:
      for M0 in [2e-3,5.8e-3,2e-2]:
        q=replace(ModelParameters(),site_density_multiplier=site,stress_work_factor=work); s,h=run(build(params=q,M0=M0),RampNoHold(5,1500)); met=density_metrics(h,target)
        objective=met["density_curve_rmse"]**2+(.05*(s.G_m*1e6-2.14))**2+.002*((np.log10(site/39.5))**2+(np.log10(work))**2+(np.log10(M0/5.8e-3))**2)
        joint.append((objective,site,work,M0,s,h,met))
    x=min(joint,key=lambda y:y[0]); candidates.append(("density_plus_grain_endpoint",x[1],x[2],x[3]))
    params=[]; objectives=[]; grain=[]
    for label,site,work,M0 in candidates:
        q=replace(ModelParameters(),site_density_multiplier=site,stress_work_factor=work); s,h=run(build(params=q,M0=M0),RampNoHold(5,1500)); met=density_metrics(h,target)
        params.append({"objective_mode":label,"site_density_multiplier":site,"stress_work_factor":work,"M0_m4_J_s":M0})
        objectives.append({"objective_mode":label,**met,"final_density":s.rho,"final_G_um":s.G_m*1e6,"density_endpoint_error":s.rho-.975,"G_endpoint_error_um":s.G_m*1e6-2.14})
        grain.append({"objective_mode":label,"model_G_um":s.G_m*1e6,"target_G_um":2.14,"residual_um":s.G_m*1e6-2.14})
        pred=np.interp(target.T_C,h.T_C,h.rho)
        for (_,t),p in zip(target.iterrows(),pred): residuals.append({"objective_mode":label,"T_C":t.T_C,"target_rho":t.fractional_density,"model_rho":p,"residual":p-t.fractional_density})
        h.assign(objective_mode=label).to_csv(OUT/f"history_calibration_{label}.csv",index=False)
    pd.DataFrame(objectives).to_csv(OUT/"CS_calibration_objective_comparison.csv",index=False)
    pd.DataFrame(params).to_csv(OUT/"CS_calibrated_parameter_sets.csv",index=False)
    pd.DataFrame(residuals).to_csv(OUT/"CS_density_curve_residuals.csv",index=False)
    pd.DataFrame(grain).to_csv(OUT/"CS_grain_endpoint_residuals.csv",index=False)

def failure_class(row,first,final):
    if not first: return "first_step_unattainable"
    if first.get("G_um",0)>.39: return "first_step_overgrown"
    if first.get("rho",0)>=.966: return "target_reached_in_first_step"
    if final["final_rho"]<.966:
        if final.get("activity",1)<.05: return "low_T2_nucleation_or_sink_limited"
        if final.get("closed_fraction",0)>.03: return "low_T2_closed_shrinkage_limited"
        if final.get("P_excess_W_m3",0)>final.get("P_dens_W_m3",0): return "low_T2_excess_power_PR_loss"
        return "low_T2_open_pore_removal_limited"
    if final["final_G_um"]>.39:
        return "high_T2_pore_coarsening_pin_loss" if final.get("S_Z",1)>.5 else "high_T2_grain_growth_activation"
    return "mixed"

def chen_decomposition():
    source=pd.read_csv("results/zro2_forward_natural_pore_evolution_target_8ysz/chen_classification_points.csv"); rows=[]; reps=[]
    m=build()
    for _,x in source.iterrows():
        state,h=run(m,TwoStep(x.T1_C,x.T2_C,x.switch_density,x.hold_h),600,600)
        second=h[h.T_C.eq(x.T2_C)]; firstrow=second.iloc[0].to_dict() if len(second) else None
        last=h.iloc[-1].to_dict(); start2=second.iloc[0] if len(second) else h.iloc[-1]
        last.update({"final_rho":state.rho,"final_G_um":state.G_m*1e6})
        tau=[]
        try: tau=[v for v in json.loads(last["tau_remove_s_json"]) if np.isfinite(v)]
        except Exception: pass
        row={"T1_C":x.T1_C,"T2_C":x.T2_C,"switch_density":x.switch_density,"hold_h":x.hold_h,
             "first_step_attained":bool(len(second)),"first_step_rho":firstrow.get("rho",np.nan) if firstrow else np.nan,
             "first_step_G1_um":firstrow.get("G_um",np.nan) if firstrow else np.nan,"first_step_fine_pore_fraction":firstrow.get("fine_pore_fraction",np.nan) if firstrow else np.nan,
             "first_step_D50_nm":firstrow.get("pore_D50_m",np.nan)*1e9 if firstrow else np.nan,"first_step_D90_nm":firstrow.get("pore_D90_m",np.nan)*1e9 if firstrow else np.nan,
             "first_step_closed_fraction":firstrow.get("closed_fraction",np.nan) if firstrow else np.nan,"first_step_R_Z_eff_um":firstrow.get("R_Z_eff_m",np.nan)*1e6 if firstrow else np.nan,
             "first_step_Gamma_growth":firstrow.get("Gamma_growth",np.nan) if firstrow else np.nan,"final_rho":state.rho,"final_G_um":state.G_m*1e6,
             "second_step_density_gain":state.rho-start2.rho,"second_step_growth_fraction":(state.G_m-start2.G_m)/max(state.G_m,1e-30),
             "open_shrinkage_rate":last.get("rho_dot_open_sinv",np.nan),"closed_shrinkage_rate":last.get("rho_dot_closed_sinv",np.nan),
             "PR_flux":last.get("bin_crossing_rate",np.nan),"isolation_flux":last.get("isolation_rate",np.nan),"closure_flux":last.get("closure_rate",np.nan),
             "mean_tau_remove_s":np.mean(tau) if tau else np.inf,"D90_tau_remove_s":np.max(tau) if tau else np.inf,
             **{k:last.get(k,np.nan) for k in ["P_surf_W_m3","P_dens_W_m3","P_excess_W_m3","sigma_eff_Pa","activity","S_Z","Gamma_mobile","Gamma_growth","closed_fraction"]}}
        row["failure_class"]=failure_class(row["first_step_attained"],firstrow,last); rows.append(row)
        if x.T1_C==1450 and x.switch_density==.85 and x.hold_h==20 and x.T2_C in [1050,1150,1300]: reps.append(h.assign(path_label=f"T2_{int(x.T2_C)}C"))
    out=pd.DataFrame(rows); out.to_csv(OUT/"chen_failure_decomposition.csv",index=False)
    out.groupby(["T1_C","switch_density","hold_h","failure_class"],dropna=False).size().reset_index(name="count").to_csv(OUT/"chen_boundary_decomposition.csv",index=False)
    pd.concat(reps,ignore_index=True).to_csv(OUT/"representative_chen_path_histories.csv",index=False)

def oat_audit():
    base=ModelParameters(); cases=[]
    def add(parameter,values,fn):
        for label,value in values: cases.append((parameter,label,value,fn(value)))
    factors=[("0.1x",.1),("0.3x",.3),("1x",1.),("3x",3.),("10x",10.)]
    add("C_PR",factors,lambda f:replace(base,C_PR_m2=base.C_PR_m2*f)); add("tau_remove0",factors,lambda f:replace(base,open_removal_multiplier=1/f))
    add("rho_close_mid",[(str(x),x) for x in [.84,.87,.90,.93]],lambda x:replace(base,rho_close_mid=x))
    add("rho_close_width",[(str(x),x) for x in [.015,.03,.06]],lambda x:replace(base,rho_close_width=x))
    add("closed_shrinkage",factors,lambda f:replace(base,closed_tau0_s=base.closed_tau0_s/f)); add("A_closed_capacity",[("0.3x",.3),("1x",1.),("3x",3.)],lambda x:replace(base,A_closed_capacity=x))
    add("C_Z",factors,lambda f:replace(base,zener_length_factor=f)); add("C_pd",factors,lambda f:replace(base,mobile_drag_coefficient=base.mobile_drag_coefficient*f))
    summary=[]; counts=[]; pathways=[]
    for parameter,label,value,q in cases:
        print(f"OAT {parameter} {label}",flush=True)
        m=build(params=q); cs,h=run(m,RampNoHold(5,1500),240,600); rapid,hf=run(m,RampNoHold(50,1500),60,120); two,ht=run(m,TwoStep(1450,1150,.85,20),600,1200); cm=coarse_map(m)
        complete=bool(cs.t_s>=RampNoHold(5,1500).t_end_s-1 and rapid.t_s>=RampNoHold(50,1500).t_end_s-1 and two.t_s>=TwoStep(1450,1150,.85,20).t_end_s-1)
        summary.append({"parameter":parameter,"level":label,"value":value,"integration_complete":complete,"CS_rho":cs.rho,"CS_G_um":cs.G_m*1e6,"rate50_rho":rapid.rho,"rate50_G_um":rapid.G_m*1e6,"two_rho":two.rho,"two_G_um":two.G_m*1e6})
        counts.append({"parameter":parameter,"level":label,"chen_success_count":int(cm.chen_success.sum()),"lower_failure_count":int((~cm.density_ok).sum()),"upper_failure_count":int((~cm.growth_ok).sum()),"cases":len(cm)})
        common=max(h.rho.min(),hf.rho.min()),min(h.rho.max(),hf.rho.max()); rho=common[1]
        slowD=np.interp(rho,h.rho,h.pore_D90_m); fastD=np.interp(rho,hf.rho,hf.pore_D90_m); slowG=np.interp(rho,h.rho,h.G_um); fastG=np.interp(rho,hf.rho,hf.G_um)
        pathways.append({"parameter":parameter,"level":label,"matched_rho":rho,"both_paths_attain_interval":common[1]>common[0],"fast_smaller_D90":fastD<slowD,"fast_smaller_G":fastG<slowG,
                         "integration_complete":complete,"two_density_target":two.rho>=.966,"two_grain_target":two.G_m*1e6<=.39,"two_finer_than_highT":two.G_m<cs.G_m,"chen_success_count":int(cm.chen_success.sum())})
        pd.DataFrame(summary).to_csv(OUT/"pore_evolution_OAT_summary.csv",index=False); pd.DataFrame(counts).to_csv(OUT/"pore_evolution_OAT_chen_counts.csv",index=False); pd.DataFrame(pathways).to_csv(OUT/"pore_evolution_OAT_pathway_metrics.csv",index=False)

def microwave_audit():
    rows=[]
    for multiplier in [1,3,10,30,100,300,1000]:
        q=replace(ModelParameters(),effective_nucleation_activity_multiplier=multiplier,microwave_mode="effective_nucleation_activity_multiplier")
        s,_=run(build(params=q),RampNoHold(50,1500),60,1e9)
        rows.append({"mode":q.microwave_mode,"multiplier":multiplier,"final_rho":s.rho,"final_G_um":s.G_m*1e6,"target_reached":s.rho>=.98,"thermal_validation_evidence":False})
    frame=pd.DataFrame(rows)
    below=frame[frame.final_rho<.98].tail(1); above=frame[frame.final_rho>=.98].head(1)
    if len(below) and len(above):
        a,b=below.iloc[0],above.iloc[0]; frac=(.98-a.final_rho)/(b.final_rho-a.final_rho)
        estimate=float(np.exp(np.log(a.multiplier)+frac*(np.log(b.multiplier)-np.log(a.multiplier))))
        frame=pd.concat([frame,pd.DataFrame([{"mode":frame.iloc[0]["mode"],"multiplier":estimate,"final_rho":.98,"final_G_um":np.nan,"target_reached":True,"thermal_validation_evidence":False}])],ignore_index=True)
    frame.to_csv(OUT/"microwave_multiplier_required.csv",index=False)

def run_all():
    OUT.mkdir(parents=True,exist_ok=True)
    if not (OUT/"baseline_failure_summary.csv").is_file(): baseline_summary()
    if not (OUT/"barrier_mode_comparison.csv").is_file(): barrier_audit()
    if not (OUT/"CS_calibration_objective_comparison.csv").is_file(): calibration_audit()
    if not (OUT/"chen_failure_decomposition.csv").is_file(): chen_decomposition()
    oat_audit(); microwave_audit()
    state={"status":"diagnostic_complete_not_validated","model_physics_changed":False,"barrier_modes":MODES,"validation_claim":False}
    (OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n"); return state
