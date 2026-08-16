from __future__ import annotations
from dataclasses import replace
import json
import numpy as np
import pandas as pd

from .conditioned_950c import make_pdf_conditioned_initial_state, model, run_path, ConditionedTwoStep
from .integrator import ModelParameters
from .schedules import RampNoHold, Iso


def calibrated_model(parameter_overrides=None):
    p=pd.read_csv("results/zro2_forward_pdf_conditioned_950C_comparison/pdf_conditioned_calibrated_parameters.csv").query("calibration_mode=='density_plus_grain_trajectory_conditioned'").iloc[0]
    q=replace(ModelParameters(),site_density_multiplier=p.site,surface_power_length2_m2=1e-19*p.work,**(parameter_overrides or {}))
    return model(q,p.M0),q


def state_from_row(row):
    return make_pdf_conditioned_initial_state(rho=row.rho_start,G_nm=row.G_start_nm,pore_D50_nm=row.pore_D50_nm,
        pore_log_width=row.pore_log_width,phi_iso_fraction=row.phi_iso_fraction,phi_closed_fraction=row.phi_closed_fraction)


def tau_d90(row):
    values=[x for x in json.loads(row.tau_remove_s_json) if np.isfinite(x)]
    return max(values or [np.nan])


def longest_span(x,y,threshold):
    good=np.asarray(y)>=threshold; x=np.asarray(x); best=run=None
    for i,ok in enumerate(good):
        if ok and run is None: run=i
        if run is not None and (not ok or i==len(good)-1):
            end=i if ok else i-1; span=x[end]-x[run]
            if best is None or span>best: best=span
            run=None
    return float(best or 0.)


def matched_metrics(case_id,h5,h50,n=25):
    lo=max(h5.rho.min(),h50.rho.min()); hi=min(h5.rho.max(),h50.rho.max()); rows=[]
    cols=["G_um","pore_D50_m","pore_D90_m","fine_pore_fraction","R_Z_eff_m","S_Z","Gamma_growth","activity","Lambda","sigma_eff_Pa","P_surf_W_m3","P_dens_W_m3","P_excess_W_m3"]
    for rho in np.linspace(lo,hi,n):
        z={"case_id":case_id,"rho":rho,"barrier_mode":"nearest_slice_clamp","barrier_extrapolated":True}
        for tag,h in (("5",h5),("50",h50)):
            g=h.sort_values("rho").drop_duplicates("rho")
            for col in cols:z[f"{col}_{tag}"]=np.interp(rho,g.rho,g[col])
            z[f"tau_remove_D90_s_{tag}"]=np.interp(rho,g.rho,[tau_d90(r) for _,r in g.iterrows()])
        z["G_5_over_G_50"]=z["G_um_5"]/max(z["G_um_50"],1e-30);rows.append(z)
    return pd.DataFrame(rows)


def chen_success(rho,G_um,density_target=.976,grain_threshold=.29):
    return rho>=density_target-.01 and G_um<=grain_threshold


def finite_window(success_T2,lower,upper):
    s=sorted(success_T2)
    return bool(lower and upper and len(s)>=2 and max(np.diff(s),default=np.inf)<=100)
