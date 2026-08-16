#!/usr/bin/env python3
from dataclasses import replace,asdict
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from zro2_forward.conditioned_950c import make_pdf_conditioned_initial_state,ConditionedTwoStep,run_path,matched
from zro2_forward.resolved_rules import ResolvedRuleModel,ResolvedRuleParameters,resolved_initial_state,ABLATIONS
from zro2_forward.schedules import RampNoHold

OUT=Path("results/zro2_forward_resolved_rules");TARGET_RHO=.976;TARGET_G=.29
def initial():return resolved_initial_state(make_pdf_conditioned_initial_state())
def params(**kw):
    cal=pd.read_csv("results/zro2_forward_pdf_conditioned_950C_comparison/pdf_conditioned_calibrated_parameters.csv").query("calibration_mode=='density_plus_grain_trajectory_conditioned'").iloc[0]
    return replace(ResolvedRuleParameters(),site_density_multiplier=cal.site,surface_power_length2_m2=1e-19*cal.work,**kw)
def model(q=None):return ResolvedRuleModel(parameters=q or params())
def classify(final,first,attained):
    density=final.rho>=TARGET_RHO;growth=final.G_m*1e6<=TARGET_G
    if not attained:return "FIRST_STEP_UNATTAINABLE"
    if first.rho>=TARGET_RHO:return "FIRST_STEP_ALREADY_TARGET"
    if density and growth:return "SUCCESS"
    if not density and growth:return "DENSIFICATION_EXHAUSTION_FAILURE"
    if density and not growth:return "GRAIN_GROWTH_FAILURE"
    return "MIXED_FAILURE"
def window_rows(x):
    rows=[]
    for key,g in x.groupby(["T1_C","switch_density","hold_h"]):
        g=g.sort_values('T2_C');s=g[g.classification.eq('SUCCESS')];density=g[g.density_ok];grain=g[g.grain_ok]
        lower=density.T2_C.min() if len(density) else np.nan;upper=grain.T2_C.max() if len(grain) else np.nan
        low_present=bool(np.isfinite(lower) and (~g.loc[g.T2_C<lower,'density_ok']).any());up_present=bool(np.isfinite(upper) and (~g.loc[g.T2_C>upper,'grain_ok']).any())
        finite=len(s)>=2 and low_present and up_present and s.T2_C.max()>s.T2_C.min() and s.T2_C.max()<key[0]
        rows.append({"T1_C":key[0],"switch_density":key[1],"hold_h":key[2],"lower_boundary_present":low_present,"upper_boundary_present":up_present,"T_lower_density_C":lower,"T_upper_growth_C":upper,"success_count":len(s),"success_min_T2_C":s.T2_C.min() if len(s) else np.nan,"success_max_T2_C":s.T2_C.max() if len(s) else np.nan,"finite_window":finite,"window_width_C":s.T2_C.max()-s.T2_C.min() if finite else np.nan,"boundary_gap_C":upper-lower if np.isfinite(lower) and np.isfinite(upper) else np.nan})
    return pd.DataFrame(rows)
def run_chen(m,T2_values=None,label="resolved_rules"):
    rows=[]
    for T1 in [1200,1300,1400,1500]:
     for T2 in (T2_values or [900,1000,1100,1200,1300]):
      if T2>=T1:continue
      for switch in [.7,.8,.9]:
       for hold in [5,10,20,30,40]:
        f,h=run_path(m,ConditionedTwoStep(T1,T2,switch,hold),initial(),900,900,label);second=h[h.T_C.eq(T2)];attained=len(second)>0;first=second.iloc[0] if attained else h.iloc[-1];last=h.iloc[-1];c=classify(f,first,attained)
        rows.append({"T1_C":T1,"T2_C":T2,"switch_density":switch,"hold_h":hold,"classification":c,"first_step_attained":attained,"rho1":first.rho,"G1_um":first.G_um,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"density_ok":f.rho>=TARGET_RHO,"grain_ok":f.G_m*1e6<=TARGET_G,"strict_success":c=="SUCCESS","closed_fraction_at_switch":first.closed_fraction,"A_closed_at_switch":first.A_closed,"PR_memory_at_switch":first.PR_memory,"final_closed_fraction":last.closed_fraction,"final_A_closed":last.A_closed,"final_PR_memory":last.PR_memory,"barrier_extrapolated":True,"mechanism_mode":"resolved_rules","gb_mobility_mode":m.parameters.gb_mobility_mode})
    return pd.DataFrame(rows)
def main():
    OUT.mkdir(parents=True,exist_ok=True);q=params();m=model(q)
    pd.DataFrame([{"parameter":k,"value":v,"role":"global resolved-rule parameter"} for k,v in asdict(q).items()]).to_csv(OUT/"resolved_rule_parameters.csv",index=False)
    eq=[("serial densification","tau_cycle=tau_nuc+tau_exchange+tau_transport","open shrinkage only"),("PR preparation","adjacent-bin conservative transfer plus named precursor transition","no direct density"),("closed shrinkage","phi_closed/tau_closed times bounded accommodation","closed density flux"),("growth","M_GB gamma/G times Gamma_migration","growth only"),("density identity","1-sum(open+isolated+closed)","all states")]
    pd.DataFrame(eq,columns=["rule","equation","density_effect"]).to_csv(OUT/"resolved_rule_equation_map.csv",index=False)
    smoke=[]
    for T2 in [950,1100,1300]:
        _,h=run_path(m,ConditionedTwoStep(1400,T2,.85,20),initial(),600,600,f"smoke_T2_{T2}");smoke.append(h)
    pd.concat(smoke).to_csv(OUT/"resolved_rule_smoke_histories.csv",index=False)
    f5,h5=run_path(m,RampNoHold(5,1500,start_C=950),initial(),120,120,"resolved_5C");f50,h50=run_path(m,RampNoHold(50,1500,start_C=950),initial(),60,60,"resolved_50C");mm=matched(h5,h50);mm.to_csv(OUT/"resolved_rule_pdf_conditioned_matched_density.csv",index=False)
    pd.DataFrame([{"case":"5C","final_rho":f5.rho,"final_G_um":f5.G_m*1e6,"final_closed_fraction":h5.iloc[-1].closed_fraction,"final_A_closed":h5.iloc[-1].A_closed,"final_PR_memory":h5.iloc[-1].PR_memory,"fast_smaller_G_fraction":np.nan,"fast_smaller_D90_fraction":np.nan},{"case":"50C","final_rho":f50.rho,"final_G_um":f50.G_m*1e6,"final_closed_fraction":h50.iloc[-1].closed_fraction,"final_A_closed":h50.iloc[-1].A_closed,"final_PR_memory":h50.iloc[-1].PR_memory,"fast_smaller_G_fraction":float((mm.G_5_over_G_50>1).mean()),"fast_smaller_D90_fraction":float((mm.pore_D90_m_50<mm.pore_D90_m_5).mean())}]).to_csv(OUT/"resolved_rule_pdf_conditioned_fast_rate_summary.csv",index=False)
    chen=run_chen(m);bounds=window_rows(chen)
    if chen.strict_success.any():
        lo=max(900,int(chen.loc[chen.strict_success,'T2_C'].min()-50));hi=min(1300,int(chen.loc[chen.strict_success,'T2_C'].max()+50));ref=run_chen(m,list(range(lo,hi+1,10)),"resolved_refined");chen=pd.concat([chen.assign(map_resolution_C=100),ref.assign(map_resolution_C=10)],ignore_index=True);bounds=window_rows(chen)
    else:chen["map_resolution_C"]=100
    chen.to_csv(OUT/"resolved_rule_chen_classification_points.csv",index=False);bounds.to_csv(OUT/"resolved_rule_chen_window_boundaries.csv",index=False)
    pd.DataFrame([{"groups":len(bounds),"strict_success_count":int(chen.strict_success.sum()),"finite_window_count":int(bounds.finite_window.sum()),"both_boundaries_count":int((bounds.lower_boundary_present&bounds.upper_boundary_present).sum()),"min_gap_C":bounds.boundary_gap_C.min(),"median_gap_C":bounds.boundary_gap_C.median(),"max_gap_C":bounds.boundary_gap_C.max()}]).to_csv(OUT/"resolved_rule_boundary_gap_summary.csv",index=False)
    ab=[]
    for name in ABLATIONS:
        aq=params(M0_factor=.1,gb_mobility_mode="bounded_uncertainty_factor",**{name:True});am=model(aq);pts=[]
        for T2 in range(1050,1300,20):
            f,h=run_path(am,ConditionedTwoStep(1400,T2,.8,40),initial(),900,900,name);second=h[h.T_C.eq(T2)];first=second.iloc[0] if len(second) else h.iloc[-1];c=classify(f,first,len(second)>0);pts.append({"T1_C":1400,"T2_C":T2,"switch_density":.8,"hold_h":40,"classification":c,"density_ok":f.rho>=TARGET_RHO,"grain_ok":f.G_m*1e6<=TARGET_G})
        px=pd.DataFrame(pts);wb=window_rows(px).iloc[0]
        ab.append({"ablation":name,"mobility_context":"bounded_uncertainty_factor_0p1","parent_window_present":False,"controlling_expected":name in ABLATIONS[:4],"success_count":int(px.classification.eq("SUCCESS").sum()),"finite_window":bool(wb.finite_window),"lower_boundary_present":bool(wb.lower_boundary_present),"upper_boundary_present":bool(wb.upper_boundary_present),"boundary_gap_C":wb.boundary_gap_C,"destructive_as_expected":np.nan,"noncontrolling_preserves_window":np.nan,"ablation_interpretation":"not testable because parent context has no finite window","classifications":"|".join(px.classification),"resolved_pathway_active":not (name in ABLATIONS[:4])})
    pd.DataFrame(ab).to_csv(OUT/"resolved_rule_ablation_summary.csv",index=False)
    prev=pd.read_csv("results/zro2_forward_required_chen_physics_gap_analysis/boundary_gap_strict_summary.csv").iloc[0];summary=pd.read_csv(OUT/"resolved_rule_boundary_gap_summary.csv").iloc[0]
    pd.DataFrame([{"metric":"strict_success_count","previous":0,"resolved":summary.strict_success_count},{"metric":"finite_window_count","previous":0,"resolved":summary.finite_window_count},{"metric":"smallest_required_shift_C","previous":125,"resolved":max(0,-summary.max_gap_C+25) if np.isfinite(summary.max_gap_C) else np.nan},{"metric":"fast_50_final_rho","previous":.959215,"resolved":f50.rho},{"metric":"fast_50_final_G_um","previous":.687985,"resolved":f50.G_m*1e6},{"metric":"closed_fraction_at_switch","previous":0,"resolved":chen.closed_fraction_at_switch.max()},{"metric":"accommodation_at_switch","previous":0,"resolved":chen.A_closed_at_switch.max()}]).to_csv(OUT/"resolved_rules_vs_previous_forward_baseline.csv",index=False)
    pd.DataFrame([{"limitation":"barrier extrapolation below fitted range","status":"unchanged"},{"limitation":"resolved parameters not independently calibrated","status":"conditional implementation"},{"limitation":"exact Mazaheri TSS schedule unavailable","status":"map diagnostic"},{"limitation":"GB mobility uncertain","status":"separate bounded audit"},{"limitation":"no validation","status":"explicit"}]).to_csv(OUT/"unresolved_limitations.csv",index=False)
    print(pd.read_csv(OUT/"resolved_rule_pdf_conditioned_fast_rate_summary.csv").to_string(index=False));print(pd.read_csv(OUT/"resolved_rule_boundary_gap_summary.csv").to_string(index=False))
if __name__=="__main__":main()
