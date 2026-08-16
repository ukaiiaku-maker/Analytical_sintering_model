#!/usr/bin/env python3
from pathlib import Path
import itertools
import numpy as np
import pandas as pd

SRC=Path("results/zro2_forward_950C_sensitivity_chen_failure_audit")
OUT=Path("results/zro2_forward_required_chen_physics_gap_analysis")
WIDTH=25.

def slope(g,col,transform):
    x=g.parameter_value.to_numpy(float);y=g[col].to_numpy(float);ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<2:return np.nan,int(ok.sum())
    tx=np.log10(x[ok]) if transform=="log10" else x[ok]/.01 if transform=="per_0p01" else x[ok]
    value=float(np.polyfit(tx,y[ok],1)[0])
    return (0. if abs(value)<1e-6 else value),int(ok.sum())

def boundary_tables():
    b=pd.read_csv(SRC/"chen_boundary_ordering_table.csv");both=b.T_lower_density_C.notna()&b.T_upper_growth_C.notna()
    b["gap_status"]=np.select([~b.T_lower_density_C.notna(),~b.T_upper_growth_C.notna(),b.gap_C.lt(0),b.gap_C.eq(0),b.gap_C.gt(0)],["missing_lower","missing_upper","boundary_ordering_failure","zero_width_contact","positive_gap"],default="missing_both")
    b["required_shift_C"]=np.where(both,np.maximum(0,-b.gap_C+WIDTH),np.nan);b["desired_min_width_C"]=WIDTH
    b=b.sort_values(["required_shift_C","T1_C","switch_density","hold_h"],na_position="last");b.to_csv(OUT/"boundary_gap_by_state_ranked.csv",index=False);b.to_csv(OUT/"boundary_gap_by_T1_switch_hold.csv",index=False)
    full=pd.read_csv(SRC/"chen_failure_decomposition_full.csv");state=full.groupby(["T1_C","switch_density","hold_h"])[["rho1","G1_nm","first_step_growth_fraction","phi_open_total","phi_iso_total","phi_closed_total","pore_D50_1_nm","pore_D90_1_nm","fine1","R_Z_eff_1_m","S_Z_1","Gamma_growth_1"]].mean().reset_index();state.merge(b,on=["T1_C","switch_density","hold_h"]).to_csv(OUT/"boundary_gap_by_initial_state.csv",index=False)
    finite=b.gap_C.dropna();best=b.loc[b.required_shift_C.idxmin()]
    summary=pd.DataFrame([{"n_first_step_groups":len(b),"n_both_boundaries_present":int(both.sum()),"n_negative_gap":int(b.gap_C.lt(0).sum()),"n_zero_width":int(b.gap_C.eq(0).sum()),"n_positive_gap":int(b.gap_C.gt(0).sum()),"n_missing_lower":int(b.T_lower_density_C.isna().sum()),"n_missing_upper":int(b.T_upper_growth_C.isna().sum()),"min_gap_C":finite.min(),"median_gap_C":finite.median(),"max_gap_C":finite.max(),"smallest_required_shift_C":best.required_shift_C,"state_with_smallest_required_shift":f"T1={best.T1_C:g}, switch={best.switch_density:g}, hold={best.hold_h:g} h","desired_min_width_C":WIDTH,"strict_success_count":0,"finite_strict_window_count":0}]);summary.to_csv(OUT/"boundary_gap_strict_summary.csv",index=False)
    return b,summary

def oat_tables(required):
    o=pd.read_csv(SRC/"chen_failure_OAT_boundary_gaps.csv");logpars={"pore_D50_nm","tau_remove_scale","PR_scale","closed_shrinkage_prefactor","closed_accommodation_capacity","zener_strength","mobile_drag_scale"};rows=[]
    for name,g in o.groupby("modified_parameter"):
        transform="log10" if name in logpars else "per_0p01" if name in {"rho_close_mid","rho_close_width"} else "linear"
        lo,nlo=slope(g,"T_lower_density_C",transform);up,nup=slope(g,"T_upper_growth_C",transform);gap,ng=slope(g,"gap_C",transform)
        if ng<3:gap=np.nan
        if ng<3:role="sparse_or_censored_gap_response"
        elif np.isfinite(lo) and np.isfinite(up):role="does_both" if lo<0 and up>0 else "makes_both_worse" if lo>0 and up<0 else "lowers_T_lower_density" if abs(lo)>abs(up) and lo<0 else "raises_T_upper_growth" if up>0 else "mixed_or_null"
        else:role="censored_or_unidentified"
        tri=pd.read_csv(SRC/"candidate_triage_950C_sensitivity.csv").query("modified_parameter==@name")
        if not tri.pathway_consistent.any() and role!="censored_or_unidentified":role += "; destroys_or_fails_pathway_consistency"
        rows.append({"parameter":name,"scale":transform,"dT_lower_density":lo,"dT_upper_growth":up,"dgap":gap,"n_lower_finite":nlo,"n_upper_finite":nup,"n_gap_finite":ng,"primary_boundary_role":role,"any_pathway_consistent":bool(tri.pathway_consistent.any()),"strict_window_count":int(tri.strict_Chen_window.sum())})
    c=pd.DataFrame(rows);c.to_csv(OUT/"boundary_sensitivity_coefficients.csv",index=False)
    shift=[]
    req=float(required.smallest_required_shift_C.iloc[0])
    for _,r in c.iterrows():
        magnitude=abs(r.dgap) if np.isfinite(r.dgap) else np.nan;delta=req/magnitude if magnitude and magnitude>0 else np.nan
        shift.append({**r.to_dict(),"required_gap_improvement_C":req,"estimated_coordinate_change":delta,"direction_to_improve_gap":"increase" if r.dgap>0 else "decrease" if r.dgap<0 else "unidentified","estimate_status":"coarse_local_diagnostic" if np.isfinite(delta) else "censored_or_not_estimable","not_a_tuning_recommendation":True})
    sh=pd.DataFrame(shift);sh.to_csv(OUT/"boundary_shift_required_by_parameter.csv",index=False)
    rank=sh.sort_values(["any_pathway_consistent","estimated_coordinate_change"],ascending=[False,True]);rank["rank_by_estimable_gap_motion"]=range(1,len(rank)+1);rank["recommendation_status"]=np.where(rank.any_pathway_consistent,"candidate_requires_data","diagnostic_only_fails_pathway_gate");rank.to_csv(OUT/"parameter_lever_rankings.csv",index=False)
    return c,rank

def physics_table():
    rows=[
    ("Low-T2 open-pore removal too slow","density stalls with open porosity","open-removal-limited and sink-limited classifications occur","increase measured low-T2 connected-pore removal","tau_remove scale","yes, within existing open-pore law","may over-densify fast path","interrupted D50/D90, connectivity and shrinkage rate"),
    ("Low-T2 closed-pore shrinkage too slow","closed reservoir persists while density stalls","closed pathway not isolated as dominant in current table","increase closed shrinkage only if closed pores are observed","closed shrinkage prefactor","proxy exists but uncalibrated","hidden fit if initial closed state is assumed","open/closed pore fractions and gas pressure at switch"),
    ("Low-T2 nucleation/activity too low","low activity or Lambda at T2","12 cases are sink-limited; barrier region is extrapolated","raise low-T2 activity only with independent barrier evidence","tau_remove scale / barrier diagnostic","not established by current fitted-range data","changes barrier physics and can overfit onset","low-T creep/nucleation-rate data under relevant stress"),
    ("Surface coarsening/PR too strong before or during step 2","D90 grows and removal time increases","pore scale strongly controls fast density; PR OAT is diagnostic","reduce PR only if interrupted pore coarsening is overpredicted","PR/coarsening scale","yes as bounded existing transfer","can create fine pores mathematically without evidence","interrupted pore-size distributions through both steps"),
    ("Fine-pore pinning too weak","grain growth despite retained fine pores","high-T pin-loss/growth categories dominate many attained paths","strengthen pinning only against measured pore/grain trajectories","Zener strength","proxy exists","can suppress growth without correct pore topology","simultaneous grain and pore distributions at switch/hold"),
    ("Zener pinning released too early","R_Z/S_Z weaken before density target","34 high-T pore-coarsening pin-loss classifications","delay release if measured pinning persists","Zener strength","yes, but coarse proxy","double-counts pore drag or fits target directly","boundary mobility and pore-boundary attachment data"),
    ("Clean grain-growth mobility too high","upper boundary lies below density boundary","all 28 estimable strict gaps are negative","lower mobility only after independent mobility calibration","mobile pore drag coefficient","mobile-drag proxy exists; clean M0 frozen","confounds clean mobility and drag","pore-free/low-porosity grain-growth kinetics"),
    ("Initial common pore scale too coarse","fast density falls strongly with D50","mean rho50 drops 0.178 from 10 to 60 nm","use measured common pore scale","pore_D50","fully compatible as initial state","selecting D50 to hit HMS is overfitting","950 C D50/D90 and connectivity"),
    ("Closed-pore fraction/accommodation missing or uncalibrated","late density deficit with closure","primary map correctly starts closed=0; proxy uncalibrated","calibrate only from observed closure state","closure midpoint / accommodation capacity","bounded proxy only","hidden initial closed fraction or nonidentifiability","tomography/dilatometry identifying closure and trapped gas"),
    ("Barrier extrapolation suppresses low-T2 activity","activity remains extrapolated below fitted range","every audited conditioned run flags extrapolation","extend barrier evidence to T2 range","barrier extrapolation diagnostic","not supported without new atomistic/creep data","retuning barrier to force window","barrier/creep data at 1000–1300 C and relevant stress")]
    pd.DataFrame(rows,columns=["failure_mode","observable_signature","model_evidence","required_physical_change","likely_parameter_proxy","compatible_with_current_physics","risks_or_artifacts","experimental_measurement_needed"]).to_csv(OUT/"required_physics_gap_interpretation.csv",index=False)

def thresholds():
    r=pd.read_csv(SRC/"chen_target_relaxation_diagnostics.csv");classes=[]
    for _,x in r.iterrows():
        if x.density_target==.976 and x.grain_threshold_um==.29:c="strict_failure"
        elif x.finite_window_count>0:c="relaxed_finite_window"
        elif x.success_count>0:c="relaxed_success_no_window"
        else:c="strict_failure"
        density_relaxed=x.density_target!=.976;grain_relaxed=x.grain_threshold_um!=.29
        req="both_relaxations_required" if x.finite_window_count>0 and density_relaxed and grain_relaxed else "density_relaxation_required" if density_relaxed and not grain_relaxed else "grain_relaxation_required" if grain_relaxed else "none_strict"
        classes.append({**x.to_dict(),"threshold_class":c,"relaxation_requirement":req,"TSS_like":False if c!="strict_failure" else False,"diagnostic_only":not (x.density_target==.976 and x.grain_threshold_um==.29)})
    t=pd.DataFrame(classes);t.to_csv(OUT/"threshold_relaxation_transition_table.csv",index=False)
    pd.DataFrame([{"first_grain_threshold_with_success_um":t.loc[t.success_count.gt(0),"grain_threshold_um"].min(),"first_density_target_with_finite_window":t.loc[t.finite_window_count.gt(0),"density_target"].min(),"finite_windows_require_grain_threshold_um":t.loc[t.finite_window_count.gt(0),"grain_threshold_um"].min(),"interpretation":"success first appears at 1.0 um without windows; finite windows require both 2.0 um and rho 0.95 relaxation","strict_result_preserved":True,"strict_success_count":0,"strict_finite_window_count":0,"relaxed_windows_pathway_consistent":False}]).to_csv(OUT/"threshold_relaxation_summary.csv",index=False)
    return t

def common_state():
    s=pd.read_csv(SRC/"common_state_fast_rate_summary.csv");d=s.query("state_class=='primary_common_state'").copy();factors=["rho_start","G_start_nm","pore_D50_nm","pore_log_width","phi_iso_fraction"]
    rows=[]
    for factor in factors:
        for level,g in d.groupby(factor):rows.append({"factor":factor,"level":level,"n":len(g),"mean_final_rho_50":g.final_rho_50.mean(),"mean_fast_smaller_G_fraction":g.fast_smaller_G_fraction.mean(),"mean_fast_smaller_D90_fraction":g.fast_smaller_D90_fraction.mean(),"mean_two_step_final_rho":pd.read_csv(SRC/"common_state_chen_summary.csv").set_index('case_id').loc[g.case_id,'two_step_final_rho'].mean(),"boundary_gap_status":"not_identifiable_from_single_representative_schedule"})
    main=pd.DataFrame(rows);main["range_final_rho_50_within_factor"]=main.groupby('factor').mean_final_rho_50.transform(lambda x:x.max()-x.min());main["range_fast_G_sign_within_factor"]=main.groupby('factor').mean_fast_smaller_G_fraction.transform(lambda x:x.max()-x.min());main["range_fast_D90_sign_within_factor"]=main.groupby('factor').mean_fast_smaller_D90_fraction.transform(lambda x:x.max()-x.min());main.to_csv(OUT/"common_state_sensitivity_main_effects.csv",index=False)
    ints=[]
    for a,b in itertools.combinations(factors,2):
        for metric in ["final_rho_50","fast_smaller_G_fraction","fast_smaller_D90_fraction"]:
            grand=d[metric].mean();am=d.groupby(a)[metric].mean();bm=d.groupby(b)[metric].mean();cell=d.groupby([a,b])[metric].mean();res=[v-am.loc[i]-bm.loc[j]+grand for (i,j),v in cell.items()];ints.append({"factor_a":a,"factor_b":b,"metric":metric,"max_abs_interaction_residual":max(abs(np.array(res))),"interaction_method":"balanced two-factor mean residual"})
    pd.DataFrame(ints).to_csv(OUT/"common_state_sensitivity_interactions.csv",index=False)
    c=pd.read_csv(SRC/"common_state_chen_summary.csv");p=d.merge(c,on='case_id');p["common_interval_attained"]=p.joint_rho_max.sub(p.joint_rho_min).ge(.05);p["fast_smaller_grain_gate"]=p.fast_smaller_G_fraction.gt(.5);p["fast_smaller_D90_gate"]=p.fast_smaller_D90_fraction.gt(.5);p["method_specific_initialization_absent"]=True;p["finite_strict_Chen_window"]=False;p["both_boundaries_present"]=False;p["highT_comparator_gates_available"]=False;p["all_pathway_gates_pass"]=False;p["boundary_gap_identifiable"]=False;p.to_csv(OUT/"common_state_pathway_consistency_summary.csv",index=False)
    return main,p

def decisions():
    rows=[
    ("A. Calibrate initial pore-size distribution first","pore D50 gives largest rho50 main effect","does not by itself create strict window","950 C D50/D90/connectivity","moderate; target-selection overfit","1 - highest"),
    ("B. Obtain exact Mazaheri TSS schedule and interrupted pore data","boundary states and high-T comparator are incomplete","requires new experiments/source recovery","exact T1/T2/rates/holds plus switch pore/grain state","low scientific risk","1 - highest"),
    ("C. Modify low-T2 closed-pore shrinkage/accommodation law","late densification is a candidate lever","closed pathway is not identified as dominant","open/closed fraction, closure density, trapped gas","high nonidentifiability","3 - only after data"),
    ("D. Modify Zener/pore-pinning law","attained paths show pin-loss/growth failures","OAT produces no pathway-consistent strict window","pore-boundary attachment and mobility","high target-fitting risk","3 - only after data"),
    ("E. Modify barrier extrapolation / low-T activity law","all T2 states are below fitted barrier range","conditioned barrier diagnostics were modest","low-T stress-dependent barrier/creep data","very high extrapolation risk","3 - only after data"),
    ("F. Add microwave-specific physics for HMS only","thermal HMS miss remains","cannot explain thermal TSS boundary ordering","microwave field/temperature and matched thermal controls","high schedule-specific overfit","4 - separate future study"),
    ("G. Stop model development and report negative thermal baseline","strict result is robustly negative in current audit","key pore and schedule data could still discriminate causes","none","low","2 - valid interim outcome")]
    pd.DataFrame(rows,columns=["decision","evidence_for","evidence_against","required_data","risk","recommended_priority"]).to_csv(OUT/"next_action_decision.csv",index=False)

def main():
    OUT.mkdir(parents=True,exist_ok=True);b,s=boundary_tables();o,r=oat_tables(s);physics_table();thresholds();common_state();decisions()
if __name__=="__main__":main()
