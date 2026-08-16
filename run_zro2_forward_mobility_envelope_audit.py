#!/usr/bin/env python3
from dataclasses import replace
from pathlib import Path
import numpy as np
import pandas as pd
from zro2_forward.conditioned_950c import ConditionedTwoStep,run_path,matched
from zro2_forward.material_zro2 import MaterialParameters,R
from zro2_forward.resolved_rules import ResolvedRuleModel
from zro2_forward.schedules import RampNoHold,Iso
from run_zro2_forward_resolved_rules import params,initial,classify,window_rows,TARGET_RHO,TARGET_G

OUT=Path("results/zro2_forward_mobility_envelope_audit");TANCHOR=1773.15
def design():
    base=MaterialParameters();rows=[("fixed_highT_literature",1.,base.Q_M_J_mol/1000,False,"current high-T estimate"),("CS_endpoint_calibrated",1.,base.Q_M_J_mol/1000,True,"inherited global CS endpoint calibration"),("CS_curve_regularized",1.,base.Q_M_J_mol/1000,True,"inherited CS curve/final-G regularized calibration")]
    rows += [("bounded_uncertainty_factor",x,base.Q_M_J_mol/1000,False,"bounded intrinsic M0 uncertainty") for x in [.03,.1,.3,1,3,10,30]]
    for q in [300,350,400,450,500,550]:
        factor=np.exp((q*1000-base.Q_M_J_mol)/(R*TANCHOR));rows.append(("activation_energy_envelope",factor,q,True,"M0 refit to preserve M_GB at 1500 C anchor"))
    out=[]
    for i,(mode,factor,q,refit,target) in enumerate(rows):out.append({"mobility_case_id":f"mobility_{i:02d}","gb_mobility_mode":mode,"M0_factor":factor,"Q_M_kJ_mol":q,"M0_or_K0_value":base.M0_m4_J_s*factor,"M0_refit":refit,"calibration_target_used":target,"non_validation_flag":True})
    return pd.DataFrame(out)
def chen_mini(m,drow):
    rows=[]
    for T1 in [1300,1400]:
     for T2 in [900,1000,1100,1200,1300]:
      if T2>=T1:continue
      for switch in [.8,.9]:
       for hold in [20,40]:
        f,h=run_path(m,ConditionedTwoStep(T1,T2,switch,hold),initial(),900,900,drow.mobility_case_id);second=h[h.T_C.eq(T2)];att=len(second)>0;first=second.iloc[0] if att else h.iloc[-1];c=classify(f,first,att)
        rows.append({**drow.to_dict(),"T1_C":T1,"T2_C":T2,"switch_density":switch,"hold_h":hold,"classification":c,"first_step_attained":att,"rho1":first.rho,"G1_um":first.G_um,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"density_ok":f.rho>=TARGET_RHO,"grain_ok":f.G_m*1e6<=TARGET_G,"strict_success":c=="SUCCESS","closed_fraction_at_switch":first.closed_fraction,"A_closed_at_switch":first.A_closed,"PR_memory_at_switch":first.PR_memory})
    return pd.DataFrame(rows)
def refine_success_groups(m,drow,coarse):
    rows=[]
    for key,g in coarse[coarse.strict_success].groupby(["T1_C","switch_density","hold_h"]):
        center=int(g.T2_C.median())
        for T2 in range(max(900,center-100),min(key[0]-10,center+100)+1,10):
            f,h=run_path(m,ConditionedTwoStep(key[0],T2,key[1],key[2]),initial(),900,900,drow.mobility_case_id+"_refined");second=h[h.T_C.eq(T2)];att=len(second)>0;first=second.iloc[0] if att else h.iloc[-1];c=classify(f,first,att)
            rows.append({**drow.to_dict(),"T1_C":key[0],"T2_C":T2,"switch_density":key[1],"hold_h":key[2],"classification":c,"first_step_attained":att,"rho1":first.rho,"G1_um":first.G_um,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"density_ok":f.rho>=TARGET_RHO,"grain_ok":f.G_m*1e6<=TARGET_G,"strict_success":c=="SUCCESS","closed_fraction_at_switch":first.closed_fraction,"A_closed_at_switch":first.A_closed,"PR_memory_at_switch":first.PR_memory,"map_resolution_C":10})
    return pd.DataFrame(rows)
def main():
    OUT.mkdir(parents=True,exist_ok=True);d=design();d.to_csv(OUT/"mobility_envelope_design.csv",index=False);d[["mobility_case_id","gb_mobility_mode","M0_factor","Q_M_kJ_mol","M0_or_K0_value","M0_refit","calibration_target_used"]].to_csv(OUT/"mobility_parameter_values.csv",index=False)
    cs=[];fast=[];points=[];bounds=[];path=[]
    for _,r in d.iterrows():
        q=params(gb_mobility_mode=r.gb_mobility_mode,M0_factor=r.M0_factor,Q_M_J_mol_override=r.Q_M_kJ_mol*1000);m=ResolvedRuleModel(parameters=q)
        f5,h5=run_path(m,RampNoHold(5,1500,start_C=950),initial(),120,120,r.mobility_case_id+"_5");f50,h50=run_path(m,RampNoHold(50,1500,start_C=950),initial(),60,60,r.mobility_case_id+"_50");mm=matched(h5,h50)
        iso,hi=run_path(m,Iso(1400,10),initial(),600,600,r.mobility_case_id+"_iso");two,ht=run_path(m,ConditionedTwoStep(1400,1100,.85,20),initial(),900,900,r.mobility_case_id+"_two")
        x=chen_mini(m,r);x["map_resolution_C"]=100;refined=refine_success_groups(m,r,x)
        if len(refined):x=pd.concat([x,refined],ignore_index=True).drop_duplicates(["T1_C","T2_C","switch_density","hold_h"],keep="last")
        b=window_rows(x);b.insert(0,"mobility_case_id",r.mobility_case_id);b["gb_mobility_mode"]=r.gb_mobility_mode;b["M0_factor"]=r.M0_factor;b["Q_M_kJ_mol"]=r.Q_M_kJ_mol;b["M0_or_K0_value"]=r.M0_or_K0_value;b["M0_refit"]=r.M0_refit;b["calibration_target_used"]=r.calibration_target_used;b["non_validation_flag"]=True
        points.append(x);bounds.append(b);succ=int(x.strict_success.sum());wins=int(b.finite_window.sum());lower=bool(b.lower_boundary_present.any());upper=bool(b.upper_boundary_present.any());active=bool(x.PR_memory_at_switch.max()>1e-3 and x.A_closed_at_switch.max()>1e-3 and x.closed_fraction_at_switch.max()>1e-6);consistent=bool(wins>0 and lower and upper and active)
        common={**r.to_dict(),"final_rho_5":f5.rho,"final_G_um_5":f5.G_m*1e6,"final_rho_50":f50.rho,"final_G_um_50":f50.G_m*1e6,"matched_fast_smaller_G_fraction":float((mm.G_5_over_G_50>1).mean()),"matched_fast_smaller_D90_fraction":float((mm.pore_D90_m_50<mm.pore_D90_m_5).mean()),"strict_Chen_success_count":succ,"finite_strict_Chen_window_count":wins,"lower_boundary_present":lower,"upper_boundary_present":upper,"T_lower_density_C":b.T_lower_density_C.max(),"T_upper_growth_C":b.T_upper_growth_C.min(),"boundary_gap_C":b.boundary_gap_C.max(),"PR_closed_accommodation_active":active,"pathway_consistency_flag":consistent,"non_validation_flag":True}
        cs.append({**common,"CS_endpoint_G_error_um":f5.G_m*1e6-2.14,"CS_endpoint_density_error":f5.rho-.975,"CS_trajectory_status":"not destroyed" if abs(f5.G_m*1e6-2.14)<=.3 else "outside endpoint tolerance"})
        fast.append(common);path.append({**common,"highT_iso_final_rho":iso.rho,"highT_iso_final_G_um":iso.G_m*1e6,"representative_two_step_rho":two.rho,"representative_two_step_G_um":two.G_m*1e6,"two_step_finer_than_highT":two.G_m<iso.G_m,"no_schedule_label_leakage":True,"mobility_direct_density_effect":False})
        print(r.mobility_case_id,r.gb_mobility_mode,r.M0_factor,succ,wins,flush=True)
    pd.DataFrame(cs).to_csv(OUT/"mobility_CS_conditioned_summary.csv",index=False);pd.DataFrame(fast).to_csv(OUT/"mobility_fast_rate_summary.csv",index=False);pd.concat(points).to_csv(OUT/"mobility_chen_classification_points.csv",index=False);pd.concat(bounds).to_csv(OUT/"mobility_chen_window_boundaries.csv",index=False);pd.DataFrame(path).to_csv(OUT/"mobility_pathway_consistency.csv",index=False)
    f=pd.DataFrame(fast);f[["mobility_case_id","gb_mobility_mode","M0_factor","Q_M_kJ_mol","strict_Chen_success_count","finite_strict_Chen_window_count","T_lower_density_C","T_upper_growth_C","boundary_gap_C","pathway_consistency_flag","non_validation_flag"]].to_csv(OUT/"mobility_boundary_gap_summary.csv",index=False)
    prior=pd.read_csv("results/zro2_forward_resolved_rules/resolved_rule_boundary_gap_summary.csv").iloc[0];f.assign(previous_required_shift_C=125,previous_strict_success_count=0,previous_finite_window_count=0,resolved_default_success_count=prior.strict_success_count,resolved_default_window_count=prior.finite_window_count).to_csv(OUT/"mobility_vs_previous_baseline.csv",index=False)
if __name__=="__main__":main()
