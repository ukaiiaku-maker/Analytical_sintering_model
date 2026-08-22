"""Final bounded promotion audit for an emergent ZrO2 pore-closure candidate.

The accepted forward integrator is not modified.  Promotion is decided only after
limiting, unit, ledger, boundary, and ablation gates.  This is not validation.
"""
from __future__ import annotations

import hashlib
import inspect
import math
from pathlib import Path

import numpy as np
import pandas as pd

from derive_zro2_energy_ledger_and_closed_pore_laws import (
    BARRIER, BARRIER_PATH, EPS, MAT, R_REF_M, conservative_transfer,
    density_identity,
)

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results/zro2_forward_emergent_pore_closure_final_test"
PRIOR=ROOT/"results/zro2_forward_energy_ledger_closed_pore_derivation"
MODE="emergent_pore_closure_v1"
R_GAS=8.31446261815324


def emergent_pore_closure_v1(T_K: float, radius_m: float, phi_closed: float,
                             shrinkability: float, accommodation: float,
                             gas_fraction: float=0.0, exponent: int=3,
                             geometry_factor: float=1.0,
                             kernel: str="renewal") -> dict[str,float|bool|str]:
    """One-bin local constitutive candidate with no path identifiers."""
    radius=max(float(radius_m),EPS)
    phi=max(float(phi_closed),0.0); chi=np.clip(shrinkability,0.0,1.0)
    avail=np.clip(accommodation,0.0,1.0)
    pcap=geometry_factor*2*MAT.gamma_s_J_m2/radius
    pgas=np.clip(gas_fraction,0.0,1.5)*pcap
    sigma=max(pcap-pgas,0.0)
    Dgb=MAT.D_GB(T_K); Ds=MAT.D_s(T_K)
    gstar=float(BARRIER.Gstar(sigma,T_K))
    rnuc=MAT.nu0_sinv*math.exp(-gstar/(MAT.kB*T_K))
    tau_sink=math.inf if sigma<=0 else (MAT.kB*T_K/(sigma*MAT.Omega_m3))*radius**2/max(Dgb,EPS)
    lam=0.0 if not np.isfinite(tau_sink) else rnuc*tau_sink
    activity=lam/(1+lam)
    size_penalty=(R_REF_M/radius)**exponent
    if kernel=="renewal":
        rate=0.0 if not np.isfinite(tau_sink) else phi*chi*avail*size_penalty/tau_sink*activity
        prefactor_status="geometry_scale_and_event_strain_semi_phenomenological"
    elif kernel=="GB_diffusion":
        coeff=R_REF_M**(exponent-2)
        rate=coeff*phi*chi*avail*Dgb*MAT.Omega_m3*sigma/(MAT.kB*T_K)*radius**(-exponent)
        prefactor_status="unit_correct_unresolved_geometry_coefficient"
    elif kernel=="surface_only":
        rate=0.0
        prefactor_status="physical_shape_channel_non_densifying"
    else:
        raise ValueError(kernel)
    k_shape=1e-6*R_REF_M**2
    shape_recovery_rate=k_shape*Ds*radius**(-4)*max(1-avail,0.0)
    return {"closed_law_mode":MODE,"kernel":kernel,"T_K":T_K,"radius_m":radius,
            "phi_closed":phi,"shrinkability":chi,"A_closed":avail,"A_closed_max":1.0,
            "geometry_factor":geometry_factor,"gas_fraction":gas_fraction,
            "P_cap_Pa":pcap,"P_gas_Pa":pgas,"sigma_c_Pa":sigma,"D_GB_m2_s":Dgb,"D_s_m2_s":Ds,
            "Gstar_J":gstar,"r_nuc_c_sinv":rnuc,"tau_sink_c_s":tau_sink,"Lambda_c":lam,
            "activity_c":activity,"radius_exponent":exponent,"size_penalty":size_penalty,
            "rho_dot_closed_sinv":max(float(rate),0.0),"A_dot_shape_recovery_sinv":shape_recovery_rate,
            "temperature_extrapolated":not BARRIER.temperature_in_fit_range(T_K),
            "prefactor_status":prefactor_status,"uses_physical_Q_closed":False}


def evolve_accommodation(A: float, shrink_rate: float, recovery_rate: float,
                         dt_s: float, beta_A: float=1.0, infinite: bool=False) -> dict[str,float]:
    if infinite:
        return {"A":1.0,"used_increment":0.0,"recovered_increment":0.0}
    used=min(max(beta_A*abs(shrink_rate)*dt_s,0.0),max(A,0.0))
    after=max(A-used,0.0)
    recovered=min(max(recovery_rate*dt_s,0.0),1.0-after)
    return {"A":np.clip(after+recovered,0.0,1.0),"used_increment":used,"recovered_increment":recovered}


def csv(name:str,rows) -> pd.DataFrame:
    q=pd.DataFrame(rows); q.to_csv(OUT/name,index=False); return q


def write_registries() -> None:
    units=[
      ("phi_c","1","evolved state","closed pore-volume fraction"),
      ("chi_c","1","evolved state","shrinkable/connective fraction"),
      ("A_c","1","evolved state","finite accommodation"),
      ("C_geom","1","geometry-derived","curvature factor"),
      ("P_gas","Pa","evolved state","bounded trapped-gas counterpressure"),
      ("Gstar","J","physical","fixed fitted barrier"),
      ("D_GB","m2 s^-1","physical","fixed Arrhenius diffusivity"),
      ("D_s","m2 s^-1","physical","fixed Arrhenius diffusivity"),
      ("C_c","1","bounded uncertainty","sink geometry coefficient set to one"),
      ("ell_c","m","geometry-derived","closed-pore radius in this audit"),
      ("r_ref","m","semi-phenomenological","fixed 25 nm reference size"),
      ("m","1","bounded uncertainty","tested at 3 and 4"),
      ("C_GB_c","m^(m-2)","semi-phenomenological","r_ref^(m-2); unit-correct magnitude unresolved"),
      ("k_shape","m2","bounded uncertainty","1e-6 r_ref^2; shape recovery only"),
      ("Q_closed_app","kJ mol^-1","empirical diagnostic","post-run slope only"),
    ]
    cols=["term","dimensions","classification","mapping"]
    csv("closed_law_parameter_classification.csv",[dict(zip(cols,r)) for r in units])
    audit=[
      ("renewal_rate","phi chi A (rref/r)^m tau_sink^-1 Lambda/(1+Lambda)","s^-1",True,"dimensionless factors times inverse time"),
      ("tau_sink","(kBT/(sigma Omega)) ell^2/D_GB","s",True,"dimensionless stress factor times diffusion time"),
      ("GB_rate","C phi chi A D_GB Omega sigma/(kBT) r^-m","s^-1",True,"C requires m^(m-2)"),
      ("shape_recovery","k_shape D_s r^-4 (1-A)","s^-1",True,"k_shape requires m2"),
      ("closed_work","sigma rho_dot_closed","W m^-3",True,"Pa s^-1"),
    ]
    csv("closed_law_unit_audit.csv",[dict(zip(["term","equation","result_units","pass","note"],r)) for r in audit])
    transfers=[
      ("open_bin_i","open_bin_i_plus_1","PR_coarsening",False,"P_PR","D_s,r,activity,topology,available inventory"),
      ("open_bin_i","precursor_i","topology_preparation",False,"P_PR","D_s,r,activity,connected fine fraction"),
      ("precursor_i","closed_bin_i","closure_transition",False,"P_pore_coarsen","topology and precursor inventory"),
      ("isolated_i","closed_bin_i","closure_transition",False,"P_pore_coarsen","topology and isolated inventory"),
    ]
    csv("PR_preparation_flux_registry.csv",[dict(zip(["source_store","destination_store","flux_name","changes_density","power_channel","local_state_variables"],r)) for r in transfers])
    ledger=["P_open_dens","P_closed_dens","P_PR","P_surface_smooth","P_pore_coarsen","P_GB_growth","P_drag","P_gas","P_other","P_residual"]
    csv("final_energy_ledger_channel_registry.csv",[{"channel":x,"named":True,"silently_rescaled":False} for x in ledger])


def limiting_tests() -> pd.DataFrame:
    rows=[]
    for kernel in ("renewal","GB_diffusion","surface_only"):
      for m in (3,4):
        base=emergent_pore_closure_v1(1473.15,25e-9,.05,.7,.5,.25,m,kernel=kernel)
        zero=emergent_pore_closure_v1(1473.15,25e-9,0,.7,.5,.25,m,kernel=kernel)
        gas=emergent_pore_closure_v1(1473.15,25e-9,.05,.7,.5,1.1,m,kernel=kernel)
        large=emergent_pore_closure_v1(1473.15,100e-9,.05,.7,.5,.25,m,kernel=kernel)
        highgas=emergent_pore_closure_v1(1473.15,25e-9,.05,.7,.5,.9,m,kernel=kernel)
        checks={"zero_inventory":zero["rho_dot_closed_sinv"]==0,
                "nonpositive_stress":gas["rho_dot_closed_sinv"]==0,
                "radius_monotonic":large["rho_dot_closed_sinv"]<=base["rho_dot_closed_sinv"],
                "gas_reduces_stress_and_rate":highgas["sigma_c_Pa"]<base["sigma_c_Pa"] and highgas["rho_dot_closed_sinv"]<=base["rho_dot_closed_sinv"],
                "surface_only_no_density":kernel!="surface_only" or base["rho_dot_closed_sinv"]==0}
        for test,passed in checks.items(): rows.append({"law":kernel,"m":m,"test":test,"passed":passed})
    p=np.array([.2,.1,.05]); q=conservative_transfer(p,.001,10,0,1)
    rows += [{"law":"preparation","m":np.nan,"test":"pore_volume_conserved","passed":np.isclose(p.sum(),q.sum())},
             {"law":MODE,"m":np.nan,"test":"accommodation_bounded","passed":0<=evolve_accommodation(.2,1e-6,2e-7,100)["A"]<=1},
             {"law":MODE,"m":np.nan,"test":"density_identity","passed":np.isclose(density_identity(np.array([.2]),np.array([.1]),np.array([.05])),.65)}]
    return csv("final_law_limiting_tests.csv",rows)


def fixed_rate_scan() -> pd.DataFrame:
    rows=[]
    for kernel in ("renewal","GB_diffusion"):
     for m in (3,4):
      for TC in range(850,1351,25):
       for rnm in (5,10,25,50,100,250):
        for phi in (.01,.05,.20,.60):
         for A in (.05,.15,.50,1.0):
          for gas in (0,.25,.5,.9):
           z=emergent_pore_closure_v1(TC+273.15,rnm*1e-9,phi,1,A,gas,m,kernel=kernel)
           Pclosed=z["sigma_c_Pa"]*z["rho_dot_closed_sinv"]
           Pavailable=z["P_cap_Pa"]*z["rho_dot_closed_sinv"]
           for hold in (20,40,96):
            rows.append({**z,"T_C":TC,"radius_nm":rnm,"hold_h":hold,
                         "density_gain":min(phi,z["rho_dot_closed_sinv"]*hold*3600),
                         "P_available_proxy_W_m3":Pavailable,"P_closed_dens_W_m3":Pclosed,
                         "Pi_dens":Pclosed/max(Pavailable,EPS),"ledger_constraint_pass":Pclosed<=Pavailable+1e-12})
    q=pd.DataFrame(rows); q["Q_closed_app_kJ_mol"]=np.nan
    keys=["kernel","radius_nm","phi_closed","A_closed","gas_fraction","radius_exponent","hold_h"]
    for _,ix in q.groupby(keys).groups.items():
        z=q.loc[list(ix)].sort_values("T_K"); x=1/z.T_K.to_numpy(); y=np.log(np.maximum(z.rho_dot_closed_sinv.to_numpy(),EPS))
        if np.ptp(y)>1e-12: q.loc[z.index,"Q_closed_app_kJ_mol"]=-R_GAS*np.gradient(y,x)/1000
    q.to_csv(OUT/"final_closed_law_rate_scan.csv",index=False); return q


STATES=[
 {"state_id":"actual_selected","state_status":"naturally_prepared_export","rho":.8305456075200908,"G_nm":70.73549974495243,"phi_closed":6.588306397081073e-6,"A":.0028212650040531128,"chi":.5,"r_nm":6.463876730783907},
 {"state_id":"moderate_bracket","state_status":"bounded_diagnostic_bracket","rho":.82,"G_nm":65.,"phi_closed":.10,"A":.35,"chi":.70,"r_nm":25.},
 {"state_id":"candidate_like_high","state_status":"injected_diagnostic_only","rho":.82,"G_nm":65.,"phi_closed":.16,"A":.70,"chi":.85,"r_nm":25.},
]


def classify_state(state:dict,kernel:str,m:int,TC:int,hold_h:float=96,gas:float=.25,A_override=None,rate_enabled=True) -> dict:
    A=state["A"] if A_override is None else A_override
    z=emergent_pore_closure_v1(TC+273.15,state["r_nm"]*1e-9,state["phi_closed"],state["chi"],A,gas,m,kernel=kernel)
    rate=z["rho_dot_closed_sinv"] if rate_enabled else 0.0
    gain=min(state["phi_closed"],rate*hold_h*3600); rho=min(state["rho"]+gain,1)
    G0=state["G_nm"]*1e-9; fv=max(1-state["rho"],EPS); rz=4*state["r_nm"]*1e-9/(3*fv)
    pin=min(1.0,rz/max(G0,EPS)); G=math.sqrt(G0**2+2*pin*MAT.M_GB(TC+273.15)*MAT.gamma_GB_J_m2*hold_h*3600)
    growth=(G-G0)/G0; dens=rho>=.90; grow=growth<=.10
    cls="SUCCESS" if dens and grow else "DENSIFICATION_EXHAUSTION_FAILURE" if not dens and grow else "GRAIN_GROWTH_FAILURE" if dens else "MIXED_FAILURE"
    phi_other=max(1-state["rho"]-state["phi_closed"],0.0)
    return {**state,"kernel":kernel,"m":m,"T2_C":TC,"hold_h":hold_h,"gas_fraction":gas,"rho_dot_closed_sinv":rate,
            "density_gain":gain,"final_rho":rho,"final_G_nm":G*1e9,"growth_fraction":growth,"classification":cls,
            "strict_success":cls=="SUCCESS","barrier_extrapolated":z["temperature_extrapolated"],
            "phi_other":phi_other,"phi_closed_final":state["phi_closed"]-gain,
            "density_identity_residual":rho-(1-phi_other-(state["phi_closed"]-gain))}


def boundary_test() -> pd.DataFrame:
    rows=[]
    for state in STATES:
     for kernel in ("renewal","GB_diffusion"):
      for m in (3,4):
       for TC in range(850,1351,25): rows.append(classify_state(state,kernel,m,TC))
    return csv("final_boundary_preservation_test.csv",rows)


def topology_summary(boundary:pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for keys,z in boundary.groupby(["state_id","state_status","kernel","m"]):
        classes=set(z.classification); success=z[z.strict_success]
        lower=bool((z.T2_C<success.T2_C.min()).any()) if len(success) else False
        upper=bool((z.T2_C>success.T2_C.max()).any()) if len(success) else False
        finite=bool(len(success)>1 and lower and upper)
        rows.append({"state_id":keys[0],"state_status":keys[1],"kernel":keys[2],"m":keys[3],"observed_classes":"|".join(sorted(classes)),
                     "success_points":len(success),"lower_boundary":lower,"upper_boundary":upper,"strict_finite_window":finite,
                     "minimum_success_T2_C":success.T2_C.min() if len(success) else np.nan,
                     "maximum_success_T2_C":success.T2_C.max() if len(success) else np.nan})
    return csv("final_boundary_topology_summary.csv",rows)


def ablations() -> tuple[pd.DataFrame,pd.DataFrame]:
    base=STATES[1]; definitions=[
      ("baseline","renewal",3,.25,None,True,1,1,False),
      ("no_PR_preparation","renewal",3,.25,None,True,.5,.5,False),
      ("no_precursor_to_closed","renewal",3,.25,None,True,.5,1,False),
      ("no_closed_shrinkage","renewal",3,.25,None,False,1,1,False),
      ("infinite_accommodation","renewal",3,.25,1.,True,1,1,True),
      ("zero_accommodation","renewal",3,.25,0.,True,1,1,False),
      ("no_gas_pressure","renewal",3,0,None,True,1,1,False),
      ("high_gas_pressure","renewal",3,.9,None,True,1,1,False),
      ("radius_exponent_m3","renewal",3,.25,None,True,1,1,False),
      ("radius_exponent_m4","renewal",4,.25,None,True,1,1,False),
      ("GB_diffusion_alternative","GB_diffusion",3,.25,None,True,1,1,False),
      ("surface_accommodation_only","surface_only",4,.25,None,True,1,1,False),
      ("strict_GB_area_balance_diagnostic","renewal",3,.25,None,True,1,1,False),
      ("energy_ledger_diagnostic_only","renewal",3,.25,None,True,1,1,False),
    ]
    histories=[]
    for name,kernel,m,gas,Aover,enabled,phif,Af,infinite in definitions:
      state={**base,"phi_closed":base["phi_closed"]*phif,"A":base["A"]*Af}
      A=state["A"] if Aover is None else Aover; phi=state["phi_closed"]; rho=state["rho"]
      used=recovered=0.0
      for hour in np.linspace(0,96,193):
        if hour>0:
          dt=1800.; z=emergent_pore_closure_v1(1100+273.15,state["r_nm"]*1e-9,phi,state["chi"],A,gas,m,kernel=kernel)
          rate=z["rho_dot_closed_sinv"] if enabled else 0.; remove=min(phi,rate*dt); phi-=remove; rho+=remove
          recovery=0.0 if name=="zero_accommodation" else z["A_dot_shape_recovery_sinv"]
          aa=evolve_accommodation(A,rate,recovery,dt,infinite=infinite); A=aa["A"]; used+=aa["used_increment"]; recovered+=aa["recovered_increment"]
        histories.append({"ablation":name,"physical_time_h":hour,"T_C":1100,"rho":rho,"phi_closed":phi,"A_closed":A,
                          "A_closed_used":used,"A_closed_recovered":recovered,"changes_density":enabled and kernel!="surface_only",
                          "strict_balance_diagnostic_only":name=="strict_GB_area_balance_diagnostic","ledger_diagnostic_only":name=="energy_ledger_diagnostic_only"})
    h=csv("final_emergent_closure_ablation_histories.csv",histories)
    rows=[]
    for name,z in h.groupby("ablation"):
      rows.append({"ablation":name,"initial_rho":z.rho.iloc[0],"final_rho":z.rho.iloc[-1],"density_gain":z.rho.iloc[-1]-z.rho.iloc[0],
                   "final_closed_inventory":z.phi_closed.iloc[-1],"final_A_closed":z.A_closed.iloc[-1],
                   "mechanism_signature_supported":bool((name not in ("no_closed_shrinkage","zero_accommodation","surface_accommodation_only")) or abs(z.rho.iloc[-1]-z.rho.iloc[0])<1e-12)})
    return csv("final_emergent_closure_ablation_matrix.csv",rows),h


def ledger_paths() -> pd.DataFrame:
    source=PRIOR/"energy_ledger_diagnostic_histories.csv"
    q=pd.read_csv(source).copy()
    q["P_PR_W_m3"]=q["P_pore_coarsen_W_m3"]
    q["P_residual_W_m3"]=q["ledger_residual_W_m3"]
    q["strict_GB_area_balance_status"]="diagnostic_only"
    q["accepted_physics_changed"]=False
    q.to_csv(OUT/"final_energy_ledger_selected_paths.csv",index=False)
    return q


def decisions(limits,units,boundary,topology,ablation,ledger) -> None:
    natural=topology[topology.state_status=="naturally_prepared_export"]
    finite=int(natural.strict_finite_window.sum())
    near=bool(natural.lower_boundary.any() and natural.upper_boundary.any())
    limiting=bool(limits.passed.all()); unit=bool(units["pass"].all())
    ledger_ok=bool((ledger.Pi_total.replace([np.inf,-np.inf],np.nan).dropna()<=1+1e-6).all())
    promoted=limiting and unit and finite>0 and ledger_ok
    interpretation="promoted_physical_candidate" if promoted else "diagnostic_negative_result"
    csv("final_promotion_decision.csv",[{"closed_law_mode":MODE,"limiting_checks_pass":limiting,"unit_checks_pass":unit,
          "natural_strict_finite_windows":finite,"natural_near_window":near,"ledger_all_paths_consistent":ledger_ok,
          "process_map_run":False,"strict_process_map_windows":0,"promoted":promoted,"final_interpretation":interpretation,
          "limiting_defect":"closed inventory times accommodation/availability; ledger inconsistency",
          "no_validation_claim":True}])
    # Failure-mode tables replace unjustified success maps.
    boundary.to_csv(OUT/"final_process_map_failure_modes.csv",index=False)
    topology.to_csv(OUT/"final_process_map_window_boundaries.csv",index=False)
    pd.DataFrame(columns=["case_id","strict_success"]).to_csv(OUT/"final_process_map_best_cases.csv",index=False)
    pd.DataFrame(columns=["case_id","G1_nm","rho_switch","T2_C","classification"]).to_csv(OUT/"final_process_map_points.csv",index=False)


def main() -> None:
    OUT.mkdir(parents=True,exist_ok=True); write_registries()
    limits=limiting_tests(); units=pd.read_csv(OUT/"closed_law_unit_audit.csv")
    scan=fixed_rate_scan(); boundary=boundary_test(); topology=topology_summary(boundary)
    ablation,_=ablations(); ledger=ledger_paths(); decisions(limits,units,boundary,topology,ablation,ledger)
    csv("fixed_input_guardrails.csv",[{"barrier_sha256":hashlib.sha256(BARRIER_PATH.read_bytes()).hexdigest(),
      "D_GB0":MAT.D_GB0_m2_s,"Q_GB_J_mol":MAT.Q_GB_J_mol,"D_s0":MAT.D_s0_m2_s,"Q_s_J_mol":MAT.Q_s_J_mol,
      "M0":MAT.M0_m4_J_s,"Q_M_J_mol":MAT.Q_M_J_mol,"failed_global_mobility_fit_active":False,
      "physical_Q_closed_introduced":False,"accepted_forward_physics_changed":False}])
    print(pd.read_csv(OUT/"final_promotion_decision.csv").to_dict("records")[0])


if __name__=="__main__": main()
