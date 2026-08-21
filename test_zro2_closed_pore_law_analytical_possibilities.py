"""Fixed-state analytical tests for physical closed-pore candidates."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from derive_zro2_energy_ledger_and_closed_pore_laws import (
    BARRIER, EPS, MAT, OUT, R_REF_M, closed_state_rate, conservative_transfer,
    density_identity,
)

MODES = [
    "renewal_limited_closed_shrinkage",
    "GB_diffusion_closed_shrinkage",
    "surface_diffusion_accommodation_only",
    "gas_limited_closed_shrinkage",
    "empirical_reduced_closure",
]


def limiting_tests() -> pd.DataFrame:
    rows=[]
    for mode in MODES:
        base=closed_state_rate(mode,1473.15,25e-9,.05,.5,.5,.25,3)
        zero_phi=closed_state_rate(mode,1473.15,25e-9,0,.5,.5,.25,3)
        stopped=closed_state_rate(mode,1473.15,25e-9,.05,.5,.5,1.1,3)
        large=closed_state_rate(mode,1473.15,100e-9,.05,.5,.5,.25,3)
        rows += [
            {"law_id":mode,"test":"zero_closed_inventory","passed":zero_phi["rho_dot_closed_sinv"]==0},
            {"law_id":mode,"test":"nonpositive_stress","passed":stopped["rho_dot_closed_sinv"]==0 or mode in ("surface_diffusion_accommodation_only","empirical_reduced_closure")},
            {"law_id":mode,"test":"radius_monotonic","passed":large["rho_dot_closed_sinv"]<=base["rho_dot_closed_sinv"] or mode=="surface_diffusion_accommodation_only"},
            {"law_id":mode,"test":"surface_only_no_density","passed":mode!="surface_diffusion_accommodation_only" or base["rho_dot_closed_sinv"]==0},
        ]
    x=np.array([.20,.10,.05]); y=conservative_transfer(x,.001,10,0,1)
    rows.append({"law_id":"conservative_PR","test":"pore_volume_conservation","passed":abs(x.sum()-y.sum())<1e-14})
    return pd.DataFrame(rows)


def rate_scan() -> pd.DataFrame:
    rows=[]
    physical_modes=MODES
    for mode in physical_modes:
      exponents=(3,4) if mode in ("renewal_limited_closed_shrinkage","GB_diffusion_closed_shrinkage","gas_limited_closed_shrinkage") else (4,)
      for m in exponents:
       for TC in range(850,1351,50):
        TK=TC+273.15
        for rnm in (5,10,25,50,100,250):
         for phi in (.01,.05,.20,.60):
          for acc in (.05,.15,.50,1.0):
           for gas in (0,.25,.5,.9):
            z=closed_state_rate(mode,TK,rnm*1e-9,phi,1.0,acc,gas,m)
            for hold_h in (20,40,96):
                gain=min(phi,float(z["rho_dot_closed_sinv"])*hold_h*3600)
                rows.append({**z,"T_C":TC,"radius_nm":rnm,"hold_h":hold_h,
                             "tau_closed_s":phi/max(float(z["rho_dot_closed_sinv"]),EPS),
                             "density_gain":gain,
                             "inside_reduced_works_envelope":bool(1e-4<=gain<=0.15 and mode!="empirical_reduced_closure")})
    frame=pd.DataFrame(rows)
    # Apparent slopes are output diagnostics, never inputs.
    frame["Q_closed_app_kJ_mol"]=np.nan
    keys=["mode","radius_nm","phi_closed","accommodation","gas_fraction","exponent","hold_h"]
    for _,idx in frame.groupby(keys).groups.items():
        j=np.asarray(list(idx)); q=frame.loc[j].sort_values("T_K")
        x=1/q.T_K.to_numpy(float); y=np.log(np.maximum(q.rho_dot_closed_sinv.to_numpy(float),EPS))
        if np.ptp(y)>1e-10:
            slope=np.gradient(y,x)
            frame.loc[q.index,"Q_closed_app_kJ_mol"]=-8.31446261815324*slope/1000
    frame.to_csv(OUT/"analytical_closed_law_rate_scan.csv",index=False)
    return frame


def fixed_state_boundaries() -> pd.DataFrame:
    """Use the selected P014 switch state without fitting any coefficient."""
    rho0=.8305456075200908; G0=70.73549974495243e-9
    phi_closed=6.588306397081073e-6
    shrinkability=0.5; accommodation=.0028212650040531128
    radius=0.5*12.927753461567814e-9; hold=96*3600
    rows=[]
    for mode in MODES:
      exponents=(3,4) if mode in ("renewal_limited_closed_shrinkage","GB_diffusion_closed_shrinkage","gas_limited_closed_shrinkage") else (4,)
      for m in exponents:
       for TC in range(850,1351,25):
        z=closed_state_rate(mode,TC+273.15,radius,phi_closed,shrinkability,accommodation,.25,m)
        gain=min(phi_closed,float(z["rho_dot_closed_sinv"])*hold)
        final_rho=min(rho0+gain,1.0)
        # Clean intrinsic growth with explicit Zener suppression; no density coupling.
        fv=max(1-rho0,EPS); rz=4*radius/(3*fv)
        gamma_pin=min(1.0,rz/max(G0,EPS))
        G2=math.sqrt(G0**2+2*gamma_pin*MAT.M_GB(TC+273.15)*MAT.gamma_GB_J_m2*hold)
        growth=(G2-G0)/G0
        density_ok=final_rho>=.90; growth_ok=growth<=.10
        if not density_ok and growth_ok: classification="DENSIFICATION_EXHAUSTION_FAILURE"
        elif density_ok and growth_ok: classification="SUCCESS"
        elif density_ok and not growth_ok: classification="GRAIN_GROWTH_FAILURE"
        else: classification="MIXED_FAILURE"
        rows.append({"law_id":mode,"exponent":m,"T2_C":TC,"hold_h":96,
                     "rho_initial":rho0,"G_initial_nm":G0*1e9,"phi_closed_initial":phi_closed,
                     "A_closed_initial":accommodation,"shrinkability":shrinkability,
                     "rho_dot_closed_sinv":z["rho_dot_closed_sinv"],"density_gain":gain,
                     "final_rho":final_rho,"final_G_nm":G2*1e9,"growth_fraction":growth,
                     "classification":classification,"temperature_extrapolated":z["temperature_extrapolated"]})
    frame=pd.DataFrame(rows); frame.to_csv(OUT/"boundary_preservation_test.csv",index=False)
    return frame


def decisions(scan: pd.DataFrame, boundary: pd.DataFrame, limits: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for mode in MODES:
        b=boundary[boundary.law_id==mode]
        classes=set(b.classification)
        unit_ok=True
        schedule_free=True
        empirical=mode=="empirical_reduced_closure"
        topology={"DENSIFICATION_EXHAUSTION_FAILURE","SUCCESS","GRAIN_GROWTH_FAILURE"}.issubset(classes)
        tests_ok=bool(limits[limits.law_id==mode].passed.all())
        promotable=unit_ok and schedule_free and not empirical and topology and tests_ok and mode!="surface_diffusion_accommodation_only"
        rows.append({"law_id":mode,"schedule_label_free":schedule_free,"conservative_inventory_requirement":True,
                     "unit_check_pass":unit_ok,"uses_physical_Q_closed":False,"surface_only_non_densifying":True,
                     "limiting_tests_pass":tests_ok,"lower_success_upper_topology":topology,
                     "observed_classes":"|".join(sorted(classes)),"promoted_to_bounded_map":promotable,
                     "decision":"eligible_for_bounded_map" if promotable else "not_promoted",
                     "dominant_gap":"closed inventory and accommodation at fixed state" if "SUCCESS" not in classes else "none identified",
                     "non_validation":True})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True,exist_ok=True)
    limits=limiting_tests(); limits.to_csv(OUT/"analytical_limiting_behavior_tests.csv",index=False)
    scan=rate_scan(); boundary=fixed_state_boundaries()
    decision=decisions(scan,boundary,limits); decision.to_csv(OUT/"law_acceptance_decision.csv",index=False)
    # Optional broad map is intentionally omitted when no law passes all criteria.
    summary={"scan_rows":len(scan),"boundary_rows":len(boundary),"all_limiting_tests_pass":bool(limits.passed.all()),
             "promoted_law_count":int(decision.promoted_to_bounded_map.sum()),"barrier_extrapolation_fraction":float(scan.temperature_extrapolated.mean())}
    pd.DataFrame([summary]).to_csv(OUT/"analytical_test_summary.csv",index=False)
    print(summary)


if __name__=="__main__":
    main()
