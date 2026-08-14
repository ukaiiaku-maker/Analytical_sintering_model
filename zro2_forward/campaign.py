from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .barrier_json import BarrierModel
from .integrator import ForwardModel
from .schedules import RampNoHold, Iso, TwoStep

ROOT=Path("results/zro2_forward_natural_pore_evolution_target_8ysz")
BARRIER=Path("data/zro2/bicrystal_creep_barrier_export.json")


def model() -> ForwardModel:
    return ForwardModel(BarrierModel.load(BARRIER))


def run_case(name, path, dt_s=120., record_s=120.):
    final, rows=model().run(path,dt_s=dt_s,record_every_s=record_s)
    frame=pd.DataFrame(rows); frame.insert(0,"case",name)
    frame["G_um"]=frame.G_m*1e6; frame["T_C"]=frame.T_K-273.15
    frame["barrier_temperature_extrapolated"]=~frame.T_K.between(1830.15,2325.15)
    return final,frame


def _matched(histories):
    rows=[]
    for case,g in histories.groupby("case"):
        g=g.sort_values("rho").drop_duplicates("rho")
        for rho in np.arange(.65,.971,.01):
            if g.rho.min() <= rho <= g.rho.max():
                rows.append({"case":case,"rho":rho,"G_um":np.interp(rho,g.rho,g.G_um),
                             "pore_D90_nm":np.interp(rho,g.rho,g.pore_D90_m)*1e9,
                             "fine_pore_fraction":np.interp(rho,g.rho,g.fine_pore_fraction),
                             "R_Z_eff_um":np.interp(rho,g.rho,g.R_Z_eff_m)*1e6})
    return pd.DataFrame(rows)


def _explode_pores(histories):
    rows=[]
    for _,x in histories.iterrows():
        arrays={k:json.loads(x[k]) for k in ["pore_radii_m_json","phi_open_json","phi_iso_json","phi_closed_json","tau_remove_s_json"]}
        for i,r in enumerate(arrays["pore_radii_m_json"]):
            rows.append({"case":x.case,"t_s":x.t_s,"T_K":x.T_K,"rho":x.rho,"bin":i,"radius_m":r,
                         "phi_open":arrays["phi_open_json"][i],"phi_iso":arrays["phi_iso_json"][i],
                         "phi_closed":arrays["phi_closed_json"][i],"tau_remove_s":arrays["tau_remove_s_json"][i]})
    return pd.DataFrame(rows)


def run_benchmarks():
    ROOT.mkdir(parents=True,exist_ok=True)
    definitions=[("CS_thermal",RampNoHold(5,1500),120.,120.),
                 ("rate50_thermal",RampNoHold(50,1500),30.,30.),
                 ("iso_1100C_40h",Iso(1100,40),300.,600.),
                 ("iso_1500C_1h",Iso(1500,1),30.,60.),
                 ("two_level_best_map_case",TwoStep(1250,1150,.90,20),300.,600.)]
    frames=[]; summaries=[]
    for name,path,dt,record in definitions:
        state,frame=run_case(name,path,dt,record); frames.append(frame)
        summaries.append({"case":name,"final_rho":state.rho,"final_G_um":state.G_m*1e6,
                          "duration_h":state.t_s/3600})
    dense=pd.concat(frames,ignore_index=True); dense.to_csv(ROOT/"dense_histories.csv",index=False)
    _matched(dense).to_csv(ROOT/"matched_density_curves.csv",index=False)
    pores=_explode_pores(dense); pores.drop(columns="tau_remove_s").to_csv(ROOT/"pore_distribution_histories.csv",index=False)
    pores[["case","t_s","rho","bin","radius_m","tau_remove_s"]].to_csv(ROOT/"tau_remove_by_bin.csv",index=False)
    dense[["case","t_s","rho","rho_dot_open_sinv","rho_dot_closed_sinv","bin_crossing_rate","isolation_rate","closure_rate"]].to_csv(ROOT/"pore_state_fluxes.csv",index=False)
    dense[["case","t_s","rho","G_um","R_Z_eff_m","S_Z","Gamma_growth","P_Z_Pa"]].to_csv(ROOT/"zener_pinning_histories.csv",index=False)
    dense[["case","t_s","T_K","rho","sigma_eff_Pa","P_surf_W_m3","P_dens_W_m3","P_excess_W_m3","efficiency","stress_bound_hit"]].to_csv(ROOT/"energy_balance_histories.csv",index=False)
    summary=pd.DataFrame(summaries); summary.to_csv(ROOT/"benchmark_summary.csv",index=False)
    slow=summary[summary.case.eq("CS_thermal")].iloc[0]; rapid=summary[summary.case.eq("rate50_thermal")].iloc[0]
    pd.DataFrame([{"comparison":"thermal_rate_50_vs_5_C_min","rho_5":slow.final_rho,"rho_50":rapid.final_rho,
                   "G_5_um":slow.final_G_um,"G_50_um":rapid.final_G_um,
                   "smaller_grain_sign":bool(rapid.final_G_um<slow.final_G_um),
                   "density_target_met_at_50":bool(abs(rapid.final_rho-.98)<=.01),
                   "interpretation":"thermal path only; experiment is microwave-assisted"}]).to_csv(ROOT/"fast_firing_summary.csv",index=False)
    return summary


def run_map():
    ROOT.mkdir(parents=True,exist_ok=True); rows=[]
    for T1 in [1250,1350,1450,1500]:
      for T2 in [1050,1150,1250,1300]:
       if T2>=T1: continue
       for switch in [.75,.85,.90]:
        for hold in [20,40]:
         state,_=run_case("map",TwoStep(T1,T2,switch,hold),600.,1e9)
         density_ok=abs(state.rho-.976)<=.01; growth_ok=abs(state.G_m*1e6-.29)<=.10
         rows.append({"T1_C":T1,"T2_C":T2,"switch_density":switch,"hold_h":hold,
                      "final_rho":state.rho,"final_G_um":state.G_m*1e6,
                      "density_ok":density_ok,"growth_ok":growth_ok,"chen_success":density_ok and growth_ok})
    frame=pd.DataFrame(rows); frame.to_csv(ROOT/"chen_classification_points.csv",index=False)
    boundaries=[]
    for key,g in frame.groupby(["T1_C","switch_density","hold_h"]):
        lower=g.loc[g.final_rho<.966,"T2_C"]; upper=g.loc[g.final_G_um>.39,"T2_C"]
        boundaries.append({"T1_C":key[0],"switch_density":key[1],"hold_h":key[2],
                           "lower_density_exhaustion_max_T2_C":lower.max() if len(lower) else np.nan,
                           "upper_growth_min_T2_C":upper.min() if len(upper) else np.nan,
                           "lower_boundary_present":bool(len(lower)),"upper_boundary_present":bool(len(upper)),
                           "success_count":int(g.chen_success.sum())})
    pd.DataFrame(boundaries).to_csv(ROOT/"chen_window_boundaries.csv",index=False)
    score=((frame.final_rho-.976)/.01)**2+((frame.final_G_um-.29)/.10)**2
    best=frame.loc[score.nsmallest(10).index].copy(); best["normalized_target_error"]=score.loc[best.index]
    best.to_csv(ROOT/"two_step_summary.csv",index=False)
    return frame,best


def write_state():
    fast=pd.read_csv(ROOT/"fast_firing_summary.csv").iloc[0]
    chen=pd.read_csv(ROOT/"chen_classification_points.csv")
    cs=pd.read_csv(ROOT/"benchmark_summary.csv").query("case=='CS_thermal'").iloc[0]
    state={"status":"complete_not_validated","barrier_json_loaded":True,"target_pdf_loaded":True,
           "CS_final_rho":cs.final_rho,"CS_final_G_um":cs.final_G_um,
           "CS_density_within_0p01":bool(abs(cs.final_rho-.975)<=.01),"CS_G_within_0p15_um":bool(abs(cs.final_G_um-2.14)<=.15),
           "thermal_rate_smaller_grain_sign":bool(fast.smaller_grain_sign),"thermal_rate_density_target_met":bool(fast.density_target_met_at_50),
           "chen_success_count":int(chen.chen_success.sum()),"validation_claim":False,
           "barrier_temperature_policy":"PCHIP inside fit; nearest-slice clamp below 1557 C with explicit flag",
           "TSS_schedule_status":"reported final state only; schedule mapped, not reproduced"}
    (ROOT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n")
    return state
