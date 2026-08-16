#!/usr/bin/env python3
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
import pandas as pd
from zro2_forward.conditioned_950c import *
from zro2_forward.schedules import RampNoHold

def interp(frame,xcol,ycol,x):
    g=frame.sort_values(xcol).drop_duplicates(xcol); return np.interp(x,g[xcol],g[ycol])
def main():
    OUT.mkdir(parents=True,exist_ok=True); states=initial_state_table(); targets=combine_targets()
    old=Path("results/zro2_forward_natural_pore_evolution_target_8ysz"); base=pd.read_csv(old/"benchmark_summary.csv"); chen=pd.read_csv(old/"chen_classification_points.csv")
    pd.DataFrame([
      {"metric":"original_CS_final_density","value":base.query("case=='CS_thermal'").final_rho.iloc[0],"corrected_interpretation":"full-process endpoint"},
      {"metric":"original_CS_final_G_um","value":base.query("case=='CS_thermal'").final_G_um.iloc[0],"corrected_interpretation":"full-process endpoint"},
      {"metric":"original_CS_density_curve_delay_C","value":330.17,"corrected_interpretation":"initialization/pre-950 C prediction mismatch plus possible barrier extrapolation effect"},
      {"metric":"original_rate50_density","value":base.query("case=='rate50_thermal'").final_rho.iloc[0],"corrected_interpretation":"full-process thermal prediction"},
      {"metric":"original_rate50_G_um","value":base.query("case=='rate50_thermal'").final_G_um.iloc[0],"corrected_interpretation":"full-process thermal prediction"},
      {"metric":"original_HMS_density_miss","value":base.query("case=='rate50_thermal'").final_rho.iloc[0]-.98,"corrected_interpretation":"thermal model versus microwave experiment"},
      {"metric":"original_Chen_success_count","value":chen.chen_success.sum(),"corrected_interpretation":"full-process map; no finite band"},
      {"metric":"original_onset_error_label","value":np.nan,"corrected_interpretation":"replaced: initialization/pre-950 C prediction mismatch plus possible barrier extrapolation effect"}
    ]).to_csv(OUT/"baseline_reinterpretation.csv",index=False)
    m=model(); full=[]; arrivals=[]
    for name,rate in [("CS_full_process",5),("rate50_full_process",50)]:
        s,h=run_path(m,RampNoHold(rate,1500),m.initial_state(),60,60,name); full.append(h)
        vals={col:interp(h,"T_C",col,950) for col in ["rho","G_um","pore_D50_m","pore_D90_m"]}
        arrivals.append({"case":name,"T_nearest_C":950,"rho_at_950":vals["rho"],"G_at_950_nm":vals["G_um"]*1000,"density_mismatch_vs_common":vals["rho"]-.66,"G_mismatch_nm":vals["G_um"]*1000-50,
                         "pore_D50_nm":vals["pore_D50_m"]*1e9,"pore_D90_nm":vals["pore_D90_m"]*1e9,"barrier_status":"nearest-slice clamp; below fitted 1557 C"})
    pd.concat(full,ignore_index=True).to_csv(OUT/"full_process_dense_histories.csv",index=False); pd.DataFrame(arrivals).to_csv(OUT/"full_process_arrival_to_950C_state.csv",index=False)
    pd.DataFrame([{**x,"final_rho":f.rho,"final_G_um":f.G_m*1e6} for x,(f,_) in zip([{"case":"CS_full_process"},{"case":"rate50_full_process"}],[run_path(m,RampNoHold(5,1500),m.initial_state(),120,99999),run_path(m,RampNoHold(50,1500),m.initial_state(),60,99999)])]).to_csv(OUT/"full_process_prediction_summary.csv",index=False)

    td=pd.read_csv(TARGET/"density_vs_temperature_digitized.csv").query("method=='CS' and T_C>=950"); tg=pd.read_csv(TARGET/"grain_size_vs_temperature_digitized.csv").query("method=='CS' and T_C>=950")
    trials=[]
    for site in [15.,30.,39.5,60.,100.]:
      for work in [.3,1.,3.]:
       for M0 in [3e-3,5.8e-3,1e-2]:
        q=replace(ModelParameters(),site_density_multiplier=site, surface_power_length2_m2=1e-19*work)
        s,h=run_path(model(q,M0),RampNoHold(5,1500,start_C=950),make_pdf_conditioned_initial_state(),120,120,"trial")
        dr=np.sqrt(np.mean((np.interp(td.T_C,h.T_C,h.rho)-td.fractional_density)**2)); gr=np.sqrt(np.mean((np.interp(tg.T_C,h.T_C,h.G_um)-tg.G_um)**2))
        endpoint=((s.rho-.975)/.01)**2+((s.G_m*1e6-2.14)/.15)**2
        trials.append({"site":site,"work":work,"M0":M0,"density_rmse":dr,"grain_rmse_um":gr,"endpoint_objective":endpoint,"curve_objective":dr/.03+abs(s.G_m*1e6-2.14)/.3,"joint_objective":dr/.03+gr/.3,"final_rho":s.rho,"final_G_um":s.G_m*1e6})
    tr=pd.DataFrame(trials); modes=[]; histories=[]; residuals=[]
    for label,obj in [("endpoint_only_conditioned","endpoint_objective"),("observed_curve_conditioned","curve_objective"),("density_plus_grain_trajectory_conditioned","joint_objective")]:
        x=tr.loc[tr[obj].idxmin()]; q=replace(ModelParameters(),site_density_multiplier=x.site,surface_power_length2_m2=1e-19*x.work)
        s,h=run_path(model(q,x.M0),RampNoHold(5,1500,start_C=950),make_pdf_conditioned_initial_state(),60,60,label); histories.append(h)
        modes.append({"calibration_mode":label,**x.to_dict()})
        for _,t in td.iterrows():residuals.append({"calibration_mode":label,"observable":"rho_T","x":t.T_C,"target":t.fractional_density,"model":interp(h,"T_C","rho",t.T_C),"residual":interp(h,"T_C","rho",t.T_C)-t.fractional_density})
        for _,t in tg.iterrows():residuals.append({"calibration_mode":label,"observable":"G_T","x":t.T_C,"target":t.G_um,"model":interp(h,"T_C","G_um",t.T_C),"residual":interp(h,"T_C","G_um",t.T_C)-t.G_um})
    pd.DataFrame(modes).to_csv(OUT/"pdf_conditioned_CS_calibration_modes.csv",index=False); pd.DataFrame(modes)[["calibration_mode","site","work","M0"]].to_csv(OUT/"pdf_conditioned_calibrated_parameters.csv",index=False); pd.DataFrame(residuals).to_csv(OUT/"pdf_conditioned_CS_curve_residuals.csv",index=False); pd.concat(histories).to_csv(OUT/"pdf_conditioned_CS_calibration_histories.csv",index=False)
    chosen=pd.DataFrame(modes).query("calibration_mode=='density_plus_grain_trajectory_conditioned'").iloc[0]; q=replace(ModelParameters(),site_density_multiplier=chosen.site,surface_power_length2_m2=1e-19*chosen.work); cm=model(q,chosen.M0)
    s5,h5=run_path(cm,RampNoHold(5,1500,start_C=950),make_pdf_conditioned_initial_state(),60,60,"conditioned_5C_min"); s50,h50=run_path(cm,RampNoHold(50,1500,start_C=950),make_pdf_conditioned_initial_state(),20,20,"conditioned_50C_min")
    both=pd.concat([h5,h50]); both.to_csv(OUT/"pdf_conditioned_dense_histories.csv",index=False); match=matched(h5,h50); match.to_csv(OUT/"pdf_conditioned_matched_density_curves.csv",index=False)
    pd.DataFrame([{"case":"conditioned_5C_min","final_rho":s5.rho,"final_G_um":s5.G_m*1e6},{"case":"conditioned_50C_min","final_rho":s50.rho,"final_G_um":s50.G_m*1e6}]).to_csv(OUT/"pdf_conditioned_fast_rate_summary.csv",index=False)
    both[["case","T_C","rho","G_um","pore_D50_m","pore_D90_m","fine_pore_fraction","R_Z_eff_m","Gamma_growth","tau_remove_s_json"]].to_csv(OUT/"pdf_conditioned_pore_trajectory_comparison.csv",index=False)
    both[["case","T_C","rho","activity","Lambda","sigma_eff_Pa","P_surf_W_m3","P_dens_W_m3","P_excess_W_m3"]].to_csv(OUT/"pdf_conditioned_energy_balance_comparison.csv",index=False)
    sensitivity=[]
    for _,x in states.iterrows():
        st=make_pdf_conditioned_initial_state(pore_D50_nm={"nominal":25,"fine_pore":15,"coarse_pore":40,"partial_isolation":25}[x.state_id],phi_iso_fraction=.10 if x.state_id=="partial_isolation" else 0)
        for rate in [5,50]:
            sf,_=run_path(cm,RampNoHold(rate,1500,start_C=950),st,60,99999,f"{x.state_id}_{rate}")
            sensitivity.append({"state_id":x.state_id,"rate_C_min":rate,"final_rho":sf.rho,"final_G_um":sf.G_m*1e6})
    pd.DataFrame(sensitivity).to_csv(OUT/"pdf_conditioned_initial_pore_sensitivity.csv",index=False)
    barrier_rows=[]; activity_rows=[]; g_rows=[]
    for mode in ["nearest_slice_clamp","pchip_extrapolate","fixed_lowT_slope","generic_anchor_barrier"]:
        bm=model(q,chosen.M0,mode); a,ha=run_path(bm,RampNoHold(5,1500,start_C=950),make_pdf_conditioned_initial_state(),120,300,mode); b,hb=run_path(bm,RampNoHold(50,1500,start_C=950),make_pdf_conditioned_initial_state(),60,300,mode)
        t,_=run_path(bm,ConditionedTwoStep(1400,1100,.8,20),make_pdf_conditioned_initial_state(),600,99999,mode)
        barrier_rows.append({"barrier_mode":mode,"CS_rho":a.rho,"CS_G_um":a.G_m*1e6,"rate50_rho":b.rho,"rate50_G_um":b.G_m*1e6,"representative_two_step_rho":t.rho,"representative_two_step_G_um":t.G_m*1e6,"generic_non_JSON":mode=="generic_anchor_barrier"})
        for _,r in pd.concat([ha.assign(path="5C"),hb.assign(path="50C")]).iterrows():activity_rows.append({"barrier_mode":mode,"path":r.path,"T_C":r.T_C,"rho":r.rho,"activity":r.activity,"extrapolated":r.T_C<1557})
        for T in np.linspace(950,1600,66):
            for sigma in [1e7,2.5e8]:g_rows.append({"barrier_mode":mode,"T_C":T,"sigma_Pa":sigma,"Gstar_eV":bm.barrier.Gstar(sigma,T+273.15)/1.602176634e-19,"extrapolated":T<1557})
    pd.DataFrame(barrier_rows).to_csv(OUT/"pdf_conditioned_barrier_mode_comparison.csv",index=False);pd.DataFrame(activity_rows).to_csv(OUT/"barrier_mode_activity_vs_T_post950.csv",index=False);pd.DataFrame(g_rows).to_csv(OUT/"barrier_mode_Gstar_vs_T_relevant_sigma.csv",index=False)
    print(pd.DataFrame(modes).to_string(index=False)); print(pd.read_csv(OUT/"pdf_conditioned_fast_rate_summary.csv").to_string(index=False))
if __name__=="__main__":main()
