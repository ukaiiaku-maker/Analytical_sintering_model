#!/usr/bin/env python3
from dataclasses import replace
import numpy as np
import pandas as pd
from zro2_forward.conditioned_950c import *
from zro2_forward.integrator import ModelParameters

def chen_success(final_rho, final_G_um):
    return abs(final_rho-.976)<=.01 and abs(final_G_um-.29)<=.10

def finite_window(success_T2, lower_boundary_present, upper_boundary_present):
    success=sorted(success_T2)
    return bool(lower_boundary_present and upper_boundary_present and len(success)>=2 and
                any(b-a<=100 for a,b in zip(success,success[1:])))

def main():
    p=pd.read_csv(OUT/"pdf_conditioned_calibrated_parameters.csv").query("calibration_mode=='density_plus_grain_trajectory_conditioned'").iloc[0]
    q=replace(ModelParameters(),site_density_multiplier=p.site,surface_power_length2_m2=1e-19*p.work); m=model(q,p.M0); rows=[]
    for T1 in [1200,1300,1400,1500]:
      for T2 in [1000,1100,1200,1300]:
       if T2>=T1:continue
       for switch in [.70,.80,.90]:
        for hold in [5,10,20,30,40]:
         path=ConditionedTwoStep(T1,T2,switch,hold); final,h=run_path(m,path,make_pdf_conditioned_initial_state(),600,600,"conditioned_map")
         second=h[h.T_C.eq(T2)]; first=second.iloc[0] if len(second) else h.iloc[-1]; last=h.iloc[-1]
         density_ok=abs(final.rho-.976)<=.01; strict=abs(final.G_m*1e6-.29)<=.10
         if not len(second):fail="first_step_unattainable"
         elif first.G_um>.5:fail="first_step_overgrown"
         elif not density_ok:fail="lower_density_exhaustion"
         elif not strict:fail="upper_grain_growth"
         else:fail="success"
         rows.append({"T1_C":T1,"T2_C":T2,"switch_density":switch,"hold_h":hold,"barrier_mode":"nearest_slice_clamp",
           "first_step_attained":bool(len(second)),"G1_um":first.G_um,"rho1":first.rho,"fine1":first.fine_pore_fraction,"D90_1_nm":first.pore_D90_m*1e9,
           "final_rho":final.rho,"final_G_um":final.G_m*1e6,"density_ok":density_ok,"grain_0p3_ok":final.G_m*1e6<=.3,"grain_0p5_ok":final.G_m*1e6<=.5,"grain_1p0_ok":final.G_m*1e6<=1.,
           "strict_target_ok":strict,"chen_success":chen_success(final.rho,final.G_m*1e6),"failure_class":fail,"final_D90_nm":last.pore_D90_m*1e9,"final_fine_fraction":last.fine_pore_fraction,"final_R_Z_um":last.R_Z_eff_m*1e6,"final_Gamma_growth":last.Gamma_growth})
    x=pd.DataFrame(rows); x.to_csv(OUT/"pdf_conditioned_chen_classification_points.csv",index=False)
    bounds=[]
    for key,g in x.groupby(["T1_C","switch_density","hold_h"]):
        low=g.loc[~g.density_ok,"T2_C"]; high=g.loc[~g.strict_target_ok,"T2_C"]; success=sorted(g.loc[g.chen_success,"T2_C"].tolist())
        finite=finite_window(success,bool(len(low)),bool(len(high)))
        bounds.append({"T1_C":key[0],"switch_density":key[1],"hold_h":key[2],"lower_boundary_max_T2_C":low.max() if len(low) else np.nan,"upper_boundary_min_T2_C":high.min() if len(high) else np.nan,
                       "lower_boundary_present":bool(len(low)),"upper_boundary_present":bool(len(high)),"success_count":len(success),"finite_window":finite})
    pd.DataFrame(bounds).to_csv(OUT/"pdf_conditioned_chen_window_boundaries.csv",index=False)
    score=((x.final_rho-.976)/.01)**2+((x.final_G_um-.29)/.1)**2; best=x.loc[score.nsmallest(15).index].copy();best["normalized_error"]=score.loc[best.index];best.to_csv(OUT/"pdf_conditioned_best_TSS_like_paths.csv",index=False)
    x.groupby("failure_class").size().reset_index(name="count").to_csv(OUT/"pdf_conditioned_chen_failure_decomposition.csv",index=False)
    print(best.head(10).to_string(index=False)); print("success",int(x.chen_success.sum()),"finite windows",int(pd.DataFrame(bounds).finite_window.sum()))
if __name__=="__main__":main()
