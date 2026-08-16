#!/usr/bin/env python3
from dataclasses import replace
from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
from zro2_forward.conditioned_950c import ConditionedTwoStep,make_pdf_conditioned_initial_state,run_path
from zro2_forward.sensitivity_audit import calibrated_model,chen_success,finite_window,tau_d90,matched_metrics
from zro2_forward.schedules import RampNoHold,Iso

OUT=Path("results/zro2_forward_950C_sensitivity_chen_failure_audit")
T2S=[1000,1100,1200,1300]

def integrate(h,col):
    return float(np.trapezoid(h[col],h.t_s)) if len(h)>1 else 0.

def classify(first,final,attained):
    density=final.rho>=.966; grain=final.G_um<=.29
    if not attained:return "first_step_unattainable"
    if first.rho>=.976:return "first_step_already_target"
    if not density and not grain:return "mixed_density_and_growth_failure"
    if not density:
        if first.activity<1e-3:return "low_T2_nucleation_limited"
        if first.Lambda<1e-2:return "low_T2_sink_limited"
        if first.closed_fraction>first.open_fraction:return "low_T2_closed_shrinkage_limited"
        if first.bin_crossing_rate>abs(first.rho_dot_open_sinv):return "low_T2_excess_PR_coarsening"
        return "low_T2_open_pore_removal_limited"
    if not grain:
        if first.Gamma_growth>.8:return "high_T2_growth_activation"
        if first.S_Z>.5:return "high_T2_pore_coarsening_pin_loss"
        return "high_T2_mobile_drag_loss"
    return "success_point"

def boundaries(x,density_target,grain_threshold):
    rows=[]
    for key,g in x.groupby(["T1_C","switch_density","hold_h"]):
        g=g.sort_values("T2_C");dok=g.final_rho.ge(density_target-.01);gok=g.final_G_um.le(grain_threshold);success=dok&gok
        lower=float(g.loc[dok,"T2_C"].min()) if dok.any() else np.nan;upper=float(g.loc[gok,"T2_C"].max()) if gok.any() else np.nan
        low_present=bool(np.isfinite(lower) and (g.loc[g.T2_C<lower,"final_rho"]<density_target-.01).any())
        up_present=bool(np.isfinite(upper) and (g.loc[g.T2_C>upper,"final_G_um"]>grain_threshold).any())
        vals=g.loc[success,"T2_C"].tolist();isolated=len(vals)==1;gap=upper-lower if np.isfinite(lower) and np.isfinite(upper) else np.nan
        rows.append(dict(T1_C=key[0],switch_density=key[1],hold_h=key[2],density_target=density_target,grain_threshold_um=grain_threshold,T_lower_density_C=lower,T_upper_growth_C=upper,gap_C=gap,lower_boundary_present=low_present,upper_boundary_present=up_present,success_band_present=len(vals)>=2,finite_window_present=finite_window(vals,low_present,up_present),zero_width_or_isolated=isolated,lower_above_upper=bool(np.isfinite(gap) and gap<0),lower_censored=not np.isfinite(lower),upper_censored=not np.isfinite(upper),success_count=len(vals)))
    return pd.DataFrame(rows)

def run_full_decomposition():
    base=pd.read_csv("results/zro2_forward_pdf_conditioned_950C_comparison/pdf_conditioned_chen_classification_points.csv");m,_=calibrated_model();rows=[];hist=[]
    for i,r in base.iterrows():
        final,h=run_path(m,ConditionedTwoStep(r.T1_C,r.T2_C,r.switch_density,r.hold_h),make_pdf_conditioned_initial_state(),900,900,"chen_audit")
        second=h[h.T_C.eq(r.T2_C)];attained=len(second)>0;first=second.iloc[0] if attained else h.iloc[-1];post=second if attained else h.iloc[0:0]
        pores=[json.loads(first[c]) for c in ["phi_open_json","phi_iso_json","phi_closed_json"]];radii=np.array(json.loads(first.pore_radii_m_json));openp=np.array(pores[0]);large=openp[radii>40e-9].sum()/max(openp.sum(),1e-300)
        rows.append({**r.to_dict(),"T_start_C":950,"G1_nm":first.G_um*1000,"first_step_growth_fraction":first.G_um/.05-1,"phi_open_total":sum(pores[0]),"phi_iso_total":sum(pores[1]),"phi_closed_total":sum(pores[2]),"pore_D50_1_nm":first.pore_D50_m*1e9,"pore_D90_1_nm":first.pore_D90_m*1e9,"large_pore_fraction":large,"tau_remove_mean_s":np.nanmean([v for v in json.loads(first.tau_remove_s_json) if np.isfinite(v)]),"tau_remove_D90_s":tau_d90(first),"R_Z_eff_1_m":first.R_Z_eff_m,"S_Z_1":first.S_Z,"Gamma_growth_1":first.Gamma_growth,"sigma_eff_1_Pa":first.sigma_eff_Pa,"activity_1":first.activity,"Lambda_1":first.Lambda,"P_surf_1":first.P_surf_W_m3,"P_dens_1":first.P_dens_W_m3,"P_excess_1":first.P_excess_W_m3,"target_density":.976,"density_gain":final.rho-first.rho,"growth_fraction_second":final.G_m/(first.G_um*1e-6)-1,"open_shrinkage_contribution":integrate(post,"rho_dot_open_sinv"),"closed_shrinkage_contribution":integrate(post,"rho_dot_closed_sinv"),"PR_coarsening_flux_integral":integrate(post,"bin_crossing_rate"),"isolation_flux_integral":integrate(post,"isolation_rate"),"closure_flux_integral":integrate(post,"closure_rate"),"final_pore_D50_nm":h.iloc[-1].pore_D50_m*1e9,"final_pore_D90_nm":h.iloc[-1].pore_D90_m*1e9,"final_fine_pore_fraction":h.iloc[-1].fine_pore_fraction,"final_R_Z_eff_m":h.iloc[-1].R_Z_eff_m,"final_Gamma_growth":h.iloc[-1].Gamma_growth,"failure_mechanism":classify(first,h.iloc[-1],attained),"barrier_mode":"nearest_slice_clamp","barrier_extrapolated":True})
        if i in [0,base.final_rho.sub(.976).abs().add(base.final_G_um.sub(.29).abs()).idxmin(),len(base)-1]:hist.append(h.assign(representative_index=i))
    x=pd.DataFrame(rows);x.to_csv(OUT/"chen_failure_decomposition_full.csv",index=False);pd.concat(hist).to_csv(OUT/"representative_chen_path_histories.csv",index=False)
    b=boundaries(x,.976,.29);b.to_csv(OUT/"chen_boundary_ordering_table.csv",index=False);b.to_csv(OUT/"chen_boundary_gap_by_state.csv",index=False)
    score=abs(x.final_rho-.976)/.01+abs(x.final_G_um-.29)/.1;x.nsmallest(25,score.name if score.name else "final_rho") if False else x.loc[score.nsmallest(25).index].to_csv(OUT/"chen_near_miss_table.csv",index=False)
    relax=[];maps=[]
    for density,grain in itertools.product([.95,.976,.98],[.29,.5,1.,2.]):
        bb=boundaries(x,density,grain);bb["evidence_type"]="strict" if density==.976 and grain==.29 else "relaxed_diagnostic";maps.append(bb)
        ok=x.final_rho.ge(density-.01)&x.final_G_um.le(grain);best=x.loc[((x.final_rho-density)/.01).pow(2)+((x.final_G_um-grain)/max(grain,.1)).pow(2).idxmin()] if False else x.loc[(((x.final_rho-density)/.01)**2+((x.final_G_um-grain)/max(grain,.1))**2).idxmin()]
        relax.append({"density_target":density,"grain_threshold_um":grain,"evidence_type":"strict" if density==.976 and grain==.29 else "relaxed_diagnostic","success_count":int(ok.sum()),"finite_window_count":int(bb.finite_window_present.sum()),"lower_boundary_count":int(bb.lower_boundary_present.sum()),"upper_boundary_count":int(bb.upper_boundary_present.sum()),"best_T1_C":best.T1_C,"best_T2_C":best.T2_C,"best_rho":best.final_rho,"best_G_um":best.final_G_um,"overwrites_strict":False})
    pd.DataFrame(relax).to_csv(OUT/"chen_target_relaxation_diagnostics.csv",index=False);pd.concat(maps).to_csv(OUT/"chen_relaxed_threshold_boundary_maps.csv",index=False)

def oat():
    specs={"pore_D50_nm":[10,15,25,40,60],"pore_log_width":[.45,.65,.85],"tau_remove_scale":[.1,.3,1,3,10],"PR_scale":[.1,.3,1,3,10],"rho_close_mid":[.84,.87,.90,.93],"rho_close_width":[.015,.03,.06],"closed_shrinkage_prefactor":[.1,.3,1,3,10],"closed_accommodation_capacity":[.3,1,3],"zener_strength":[.1,.3,1,3,10],"mobile_drag_scale":[.1,.3,1,3,10]}
    design=[]
    for name,vals in specs.items():
        for v in vals:design.append({"case_id":f"{name}_{v}","modified_parameter":name,"parameter_value":v})
    pd.DataFrame(design).to_csv(OUT/"chen_failure_OAT_design.csv",index=False);summary=[];gaps=[];effects=[]
    for r in design:
        state_kw={};over={}
        if r["modified_parameter"] in ("pore_D50_nm","pore_log_width"):state_kw[r["modified_parameter"]]=r["parameter_value"]
        maps={"tau_remove_scale":("sink_time_factor",r["parameter_value"]),"PR_scale":("C_PR_m2",1e-23*r["parameter_value"]),"rho_close_mid":("rho_close_mid",r["parameter_value"]),"rho_close_width":("rho_close_width",r["parameter_value"]),"closed_shrinkage_prefactor":("closed_tau0_s",1e5/r["parameter_value"]),"closed_accommodation_capacity":("closed_tau0_s",1e5/r["parameter_value"]),"zener_strength":("zener_strength",r["parameter_value"]),"mobile_drag_scale":("mobile_drag_scale",r["parameter_value"])}
        if r["modified_parameter"] in maps:over[maps[r["modified_parameter"]][0]]=maps[r["modified_parameter"]][1]
        m,_=calibrated_model(over);make=lambda:make_pdf_conditioned_initial_state(**state_kw)
        f5,h5=run_path(m,RampNoHold(5,1500,start_C=950),make(),180,300,r["case_id"]+"_5");f50,h50=run_path(m,RampNoHold(50,1500,start_C=950),make(),60,120,r["case_id"]+"_50");mm=matched_metrics(r["case_id"],h5,h50)
        pts=[]
        for T2 in T2S:
            f,_=run_path(m,ConditionedTwoStep(1400,T2,.9,30),make(),900,999999,r["case_id"]);pts.append({"T1_C":1400,"T2_C":T2,"switch_density":.9,"hold_h":30,"final_rho":f.rho,"final_G_um":f.G_m*1e6})
        px=pd.DataFrame(pts);bb=boundaries(px,.976,.29).iloc[0];summary.append({**r,"final_rho_5":f5.rho,"final_rho_50":f50.rho,"fast_smaller_grain_sign":(mm.G_5_over_G_50>=1).mean()>.5,"fast_smaller_D90_sign":(mm.pore_D90_m_50<=mm.pore_D90_m_5).mean()>.5,"strict_Chen_window":bb.finite_window_present,"boundary_gap_C":bb.gap_C,"lower_boundary_present":bb.lower_boundary_present,"upper_boundary_present":bb.upper_boundary_present,"barrier_mode":"nearest_slice_clamp","barrier_extrapolated":True});gaps.append({**r,**bb.to_dict()});effects.append({**r,"HMS_density_miss":f50.rho-.98,"fast_G_sign_fraction":(mm.G_5_over_G_50>=1).mean(),"fast_D90_sign_fraction":(mm.pore_D90_m_50<=mm.pore_D90_m_5).mean()})
    s=pd.DataFrame(summary);s.to_csv(OUT/"chen_failure_OAT_summary.csv",index=False);pd.DataFrame(gaps).to_csv(OUT/"chen_failure_OAT_boundary_gaps.csv",index=False);pd.DataFrame(effects).to_csv(OUT/"chen_failure_OAT_fast_rate_effects.csv",index=False)
    tri=s.copy();tri["fast_density_attained"]=tri.final_rho_50>=.97;tri["relaxed_Chen_window"]=False;tri["failure_category"]=np.where(tri.lower_boundary_present,"growth_or_ordering","density_or_censored");tri["pathway_consistent"]=tri.fast_density_attained&tri.fast_smaller_grain_sign&tri.fast_smaller_D90_sign&tri.strict_Chen_window;tri["accepted_for_future_calibration"]=tri.pathway_consistent;tri["reason"]=np.where(tri.pathway_consistent,"bounded case worth future calibration","fails one or more physical pathway gates");tri.to_csv(OUT/"candidate_triage_950C_sensitivity.csv",index=False)

def main():run_full_decomposition();oat()
if __name__=="__main__":main()
