from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .barrier_json import BarrierModel
from .integrator import ForwardModel, ModelState, ModelParameters
from .material_zro2 import MaterialParameters
from .pore_population import initial_population, diagnostics
from .grain_growth import growth_state
from .schedules import RampNoHold, TwoStep

OUT=Path("results/zro2_forward_pdf_conditioned_950C_comparison")
BARRIER=Path("data/zro2/bicrystal_creep_barrier_export.json")
TARGET=Path("data/targets/mazaheri_8ysz_2008")

class ConditionedTwoStep:
    def __init__(self,T1_C,T2_C,switch_rho,hold_h,rate_C_min=5.,start_C=950.):
        self.T1_C=T1_C; self.T2_C=T2_C; self.switch_rho=switch_rho; self.start_C=start_C; self.rate=rate_C_min/60
        self.ramp_s=(T1_C-start_C)/self.rate; self.t_end_s=self.ramp_s+hold_h*3600
    def temperature_K(self,t_s,rho):
        if rho>=self.switch_rho:return self.T2_C+273.15
        return min(self.T1_C,self.start_C+self.rate*t_s)+273.15

def make_pdf_conditioned_initial_state(rho=.66,G_nm=50.,T_C=950.,pore_mode="common_open_lognormal",
                                       pore_D50_nm=25.,pore_log_width=.65,phi_iso_fraction=0.,phi_closed_fraction=0.):
    if pore_mode!="common_open_lognormal": raise ValueError("unsupported common pore mode")
    pop=initial_population(rho0=rho,center_m=.5*pore_D50_nm*1e-9,ln_sigma=pore_log_width)
    if phi_iso_fraction+phi_closed_fraction >= 1: raise ValueError("non-open fractions must sum below one")
    isolated=pop.phi_open*phi_iso_fraction; closed=pop.phi_open*phi_closed_fraction
    pop.phi_open-=isolated+closed; pop.phi_iso=isolated; pop.phi_closed=closed
    return ModelState(0.,T_C+273.15,rho,G_nm*1e-9,pop,1.)

def model(params=None,M0=5.8e-3,barrier_mode="nearest_slice_clamp"):
    b=BarrierModel.load(BARRIER)
    # This source branch supports the fitted-range PCHIP plus the conservative clamp.
    if barrier_mode!="nearest_slice_clamp":
        from .conditioned_barriers import barrier_with_mode
        b=barrier_with_mode(b,barrier_mode)
    return ForwardModel(b,MaterialParameters(M0_m4_J_s=M0),params or ModelParameters())

def run_path(m,path,state,dt=120,record=120,label="conditioned"):
    final,rows=m.run(path,dt_s=dt,record_every_s=record,initial_state=state)
    f=pd.DataFrame(rows); f.insert(0,"case",label); f["T_C"]=f.T_K-273.15; f["G_um"]=f.G_m*1e6
    return final,f

def initial_state_table():
    rows=[]; m=model()
    definitions=[("nominal",25,.65,0.),("fine_pore",15,.65,0.),("coarse_pore",40,.65,0.),("partial_isolation",25,.65,.10)]
    for name,D,width,iso in definitions:
        s=make_pdf_conditioned_initial_state(pore_D50_nm=D,pore_log_width=width,phi_iso_fraction=iso); d=diagnostics(s.pores)
        g=growth_state(s.G_m,s.pores.radii_m,s.pores.phi_open,s.T_K,m.material)
        rows.append({"state_id":name,"T_start_C":950,"rho_start":s.rho,"G_start_nm":s.G_m*1e9,
                     "phi_open_total":s.pores.phi_open.sum(),"phi_iso_total":s.pores.phi_iso.sum(),"phi_closed_total":s.pores.phi_closed.sum(),
                     "pore_D50_nm":d["pore_D50_m"]*1e9,"pore_D90_nm":d["pore_D90_m"]*1e9,"fine_pore_fraction":d["fine_pore_fraction"],
                     "large_pore_fraction":float(s.pores.phi_open[s.pores.radii_m>40e-9].sum()/max(s.pores.phi_open.sum(),1e-300)),
                     "R_Z_eff":g["R_Z_eff_m"],"Gamma_growth":g["Gamma_growth"],"notes":"controlled common-entry state; not fitted per method"})
    out=pd.DataFrame(rows); out.to_csv(OUT/"pdf_conditioned_initial_states.csv",index=False); return out

def combine_targets():
    frames=[]
    specs=[("density_vs_temperature_digitized.csv","Figure 2"),("density_vs_time_digitized.csv","Figure 4"),("grain_size_vs_temperature_digitized.csv","Figure 5"),("grain_size_vs_density_digitized.csv","Figure 6")]
    for file,fig in specs:
        x=pd.read_csv(TARGET/file); x["source_figure"]=fig
        x=x.rename(columns={"fractional_density":"rho","G_um":"G_um","uncertainty_density":"digitization_uncertainty"})
        for col in ["T_C","time_min","rho","G_um","digitization_uncertainty"]:
            if col not in x:x[col]=np.nan
        frames.append(x[["method","source_figure","T_C","time_min","rho","G_um","digitization_uncertainty"]])
    full=pd.concat(frames,ignore_index=True); full["included_in_pdf_conditioned_fit"]=full.T_C.isna()|full.T_C.ge(950); full["excluded_reason"]=np.where(full.included_in_pdf_conditioned_fit,"","below visible 950 C interval")
    full.to_csv(OUT/"target_curves_full_digitized.csv",index=False); full[full.included_in_pdf_conditioned_fit].to_csv(OUT/"target_curves_observed_interval_950C.csv",index=False)
    return full

def matched(frame5,frame50):
    lo=max(frame5.rho.min(),frame50.rho.min()); hi=min(frame5.rho.max(),frame50.rho.max()); rows=[]
    for rho in np.linspace(lo,hi,50):
        row={"rho":rho}
        for tag,f in [("5",frame5),("50",frame50)]:
            g=f.sort_values("rho").drop_duplicates("rho")
            for col in ["G_um","pore_D50_m","pore_D90_m","fine_pore_fraction","R_Z_eff_m","Gamma_growth","activity","Lambda","sigma_eff_Pa"]: row[f"{col}_{tag}"]=np.interp(rho,g.rho,g[col])
            tau=[max([v for v in json.loads(x) if np.isfinite(v)] or [np.nan]) for x in g.tau_remove_s_json]; row[f"tau_remove_D90_s_{tag}"]=np.interp(rho,g.rho,tau)
        row["G_5_over_G_50"]=row["G_um_5"]/max(row["G_um_50"],1e-30); rows.append(row)
    return pd.DataFrame(rows)
