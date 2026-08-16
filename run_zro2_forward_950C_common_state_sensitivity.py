#!/usr/bin/env python3
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from zro2_forward.sensitivity_audit import *

OUT=Path("results/zro2_forward_950C_sensitivity_chen_failure_audit")

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    prior=Path("results/zro2_forward_pdf_conditioned_950C_comparison")
    fast=pd.read_csv(prior/"pdf_conditioned_fast_rate_summary.csv");match=pd.read_csv(prior/"pdf_conditioned_matched_density_curves.csv");chen=pd.read_csv(prior/"pdf_conditioned_chen_classification_points.csv");win=pd.read_csv(prior/"pdf_conditioned_chen_window_boundaries.csv")
    pd.DataFrame([{"T_start_C":950,"rho_start":.66,"G_start_nm":50,"pore_D50_nominal_nm":25,"rate5_final_rho":fast.iloc[0].final_rho,"rate5_final_G_um":fast.iloc[0].final_G_um,"rate50_final_rho":fast.iloc[1].final_rho,"rate50_final_G_um":fast.iloc[1].final_G_um,"HMS_density_miss":fast.iloc[1].final_rho-.98,"fast_smaller_G_fraction":float((match.G_5_over_G_50>1).mean()),"fast_smaller_D90_fraction":float((match.pore_D90_m_50<match.pore_D90_m_5).mean()),"conditioned_Chen_cases":len(chen),"strict_success_count":int(chen.chen_success.sum()),"finite_window_count":int(win.finite_window.sum()),"barrier_mode":"nearest_slice_clamp","non_claims":"controlled sensitivity; microwave targets are not thermal validation; no validation claim"}]).to_csv(OUT/"baseline_pdf_conditioned_summary.csv",index=False)
    design=[]
    for i,v in enumerate(itertools.product([.62,.66,.70],[35,50,75],[10,15,25,40,60],[.45,.65,.85],[0,.05,.10])):
        design.append(dict(case_id=f"primary_{i:03d}",rho_start=v[0],G_start_nm=v[1],pore_D50_nm=v[2],pore_log_width=v[3],phi_iso_fraction=v[4],phi_closed_fraction=0.,state_class="primary_common_state"))
    for closed in [.02,.05]:design.append(dict(case_id=f"diagnostic_closed_{closed:.2f}",rho_start=.66,G_start_nm=50,pore_D50_nm=25,pore_log_width=.65,phi_iso_fraction=0.,phi_closed_fraction=closed,state_class="diagnostic_closed_initial_state"))
    d=pd.DataFrame(design);d.to_csv(OUT/"common_state_factorial_design.csv",index=False)
    m,_=calibrated_model(); summaries=[]; matches=[]; pores=[]; energy=[]; chenrows=[]; initials=[]
    for j,r in d.iterrows():
        state=state_from_row(r); p=state.pores
        initials.append({**r.to_dict(),"phi_open_total":p.phi_open.sum(),"phi_iso_total":p.phi_iso.sum(),"phi_closed_total":p.phi_closed.sum(),"pore_arrays_equal_across_paths":True})
        s5,h5=run_path(m,RampNoHold(5,1500,start_C=950),state_from_row(r),120,300,f"{r.case_id}_5")
        s50,h50=run_path(m,RampNoHold(50,1500,start_C=950),state_from_row(r),60,120,f"{r.case_id}_50")
        mm=matched_metrics(r.case_id,h5,h50);matches.append(mm)
        summaries.append({**r.to_dict(),"final_rho_5":s5.rho,"final_G_um_5":s5.G_m*1e6,"final_rho_50":s50.rho,"final_G_um_50":s50.G_m*1e6,"HMS_density_miss":s50.rho-.98,"joint_rho_min":mm.rho.min(),"joint_rho_max":mm.rho.max(),"fast_smaller_G_fraction":(mm.G_5_over_G_50>1).mean(),"fast_smaller_D90_fraction":(mm.pore_D90_m_50<mm.pore_D90_m_5).mean(),"span_ratio_ge_1p2":longest_span(mm.rho,mm.G_5_over_G_50,1.2),"span_ratio_ge_1p5":longest_span(mm.rho,mm.G_5_over_G_50,1.5),"barrier_mode":"nearest_slice_clamp","barrier_extrapolated":True})
        for h in (h5,h50):
            pores.append(h[["case","T_C","rho","G_um","pore_D50_m","pore_D90_m","fine_pore_fraction","R_Z_eff_m","S_Z","Gamma_growth","tau_remove_s_json"]])
            energy.append(h[["case","T_C","rho","activity","Lambda","sigma_eff_Pa","P_surf_W_m3","P_dens_W_m3","P_excess_W_m3"]])
        ts,ht=run_path(m,ConditionedTwoStep(1400,1100,.9,30),state_from_row(r),900,1800,f"{r.case_id}_two_step")
        chenrows.append({"case_id":r.case_id,"two_step_final_rho":ts.rho,"two_step_final_G_um":ts.G_m*1e6,"density_ok":ts.rho>=.966,"grain_ok":ts.G_m*1e6<=.29,"strict_success":chen_success(ts.rho,ts.G_m*1e6),"comparator_status":"optional high-temperature isothermal comparator not previously defined"})
        if j%25==0:print(f"completed {j+1}/{len(d)}",flush=True)
    pd.DataFrame(initials).to_csv(OUT/"common_state_initial_microstructures.csv",index=False);pd.DataFrame(summaries).to_csv(OUT/"common_state_fast_rate_summary.csv",index=False);pd.concat(matches).to_csv(OUT/"common_state_matched_density_curves.csv",index=False);pd.concat(pores).to_csv(OUT/"common_state_pore_trajectory_metrics.csv",index=False);pd.concat(energy).to_csv(OUT/"common_state_energy_balance_metrics.csv",index=False);pd.DataFrame(chenrows).to_csv(OUT/"common_state_chen_summary.csv",index=False)
    s=pd.DataFrame(summaries);score=abs(s.final_rho_50-.98)+abs(s.final_G_um_50-.70)/10;s.loc[score.nsmallest(20).index].to_csv(OUT/"common_state_best_cases.csv",index=False)
    print(s.groupby("pore_D50_nm").final_rho_50.mean())
if __name__=="__main__":main()
