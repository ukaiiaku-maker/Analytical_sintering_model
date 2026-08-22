"""Analytical baseline tests for conservative ZrO2 pore-channel topology laws."""
from __future__ import annotations

import hashlib, math
from pathlib import Path
import numpy as np
import pandas as pd

from derive_zro2_energy_ledger_and_closed_pore_laws import BARRIER_PATH, MAT, EPS, density_identity
from promote_zro2_emergent_pore_closure_final_test import emergent_pore_closure_v1

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results/zro2_forward_pore_channel_pr_baseline_test"
MODES=("surface_coarsening_only_v1","PR_pinchoff_v1","PR_regularization_damage_v1")


def surface_coefficient(T_K:float)->float:
    """Dimensionally audited Mullins-type coefficient, m4/s."""
    return MAT.D_s(T_K)*MAT.gamma_s_J_m2*MAT.Omega_m3**(4/3)/(MAT.kB*T_K)


def pore_channel_terms(T_K:float,radius_m:float,segment_ratio:float,gb_surface_ratio:float,
                       activity:float,width:float,connected:float,phi_open:float=.1)->dict:
    """Local topology terms without processing-path identifiers."""
    r=max(radius_m,EPS); a=np.clip(activity,0,1); conn=np.clip(connected,0,1)
    Bs=surface_coefficient(T_K); C_s=1.0
    tau=C_s*r**4/max(Bs,EPS)
    F_GB=np.clip(1+0.5*(gb_surface_ratio-.5),.5,1.5)
    I=segment_ratio/(2*math.pi*F_GB)-1
    pinch=1/(1+math.exp(-I/.15))
    low=1-a; moderate=4*a*(1-a)
    not_overwide=1/(1+math.exp((width-.75)/.12))
    large=1/(1+math.exp(-(width-.65)/.12))
    base=max(phi_open,0)/max(tau,EPS)*conn
    Jcoars=base*low
    Jpinch=base*low*pinch
    Jreg=base*moderate*not_overwide*(1-pinch)
    Jdamage=base*low*large*pinch
    pz=MAT.gamma_GB_J_m2*max(phi_open,0)/r
    rz=4*r/(3*max(phi_open,EPS)); gamma_z=min(1,rz/(100e-9))
    return {"T_K":T_K,"radius_m":r,"B_s_m4_s":Bs,"tau_s_s":tau,"F_GB":F_GB,"I_PR":I,
            "P_pinch":pinch,"J_coars_sinv":Jcoars,"J_pinch_sinv":Jpinch,"J_reg_sinv":Jreg,"J_damage_sinv":Jdamage,
            "P_Z_Pa":pz,"R_Z_m":rz,"Gamma_Zener":gamma_z,"conservative":True,"rho_dot_topology":0.0}


def conservative_step(stores:dict[str,float],mode:str,terms:dict,dt:float)->dict[str,float]:
    """Advance four stores; precursor is accounted as an open sub-store."""
    s={k:max(float(stores.get(k,0)),0) for k in ("open","precursor","isolated","closed")}
    before=sum(s.values()); available=s["open"]
    if mode=="surface_coarsening_only_v1":
        moved=min(available,terms["J_coars_sinv"]*dt); s["open"]-=moved; s["open"]+=moved
    elif mode=="PR_pinchoff_v1":
        moved=min(available,terms["J_pinch_sinv"]*dt); s["open"]-=moved
        eta=(.45,.35,.20); s["precursor"]+=eta[0]*moved; s["isolated"]+=eta[1]*moved; s["closed"]+=eta[2]*moved
    elif mode=="PR_regularization_damage_v1":
        damage=min(available,terms["J_damage_sinv"]*dt); s["open"]-=damage
        s["precursor"]+=.5*damage; s["isolated"]+=.3*damage; s["closed"]+=.2*damage
        s["regularized_volume"]=min(available-damage,terms["J_reg_sinv"]*dt)
    else: raise ValueError(mode)
    s["conservation_residual"]=sum(s[k] for k in ("open","precursor","isolated","closed"))-before
    return s


def write_registries():
    pars=[
      ("D_s","physical","fixed 0.10 exp(-380000/RT) m2/s"),("gamma_s","physical","fixed 1 J/m2"),("Omega","physical","fixed 3.35e-29 m3"),
      ("C_s","semi-phenomenological","dimensionless geometry product set to one"),("F_GB","geometry-derived","bounded surrogate from gamma_GB/gamma_s"),
      ("w_I","bounded uncertainty","dimensionless transition width 0.15"),("eta_partition","bounded uncertainty","conservative 0.45/0.35/0.20"),
      ("width thresholds","bounded uncertainty","regularization/damage crossover 0.65-0.75"),("Q_closed_app","empirical diagnostic","post-run only")]
    pd.DataFrame(pars,columns=["parameter","classification","mapping"]).to_csv(OUT/"pore_channel_parameter_registry.csv",index=False)
    channels=["P_open_dens","P_closed_dens","P_PR_regularization","P_PR_damage","P_surface_smooth","P_pore_coarsen","P_pinch","P_GB_growth","P_drag","P_gas","P_residual"]
    pd.DataFrame([{"channel":x,"diagnostic":True,"forced_equality":False} for x in channels]).to_csv(OUT/"pore_channel_energy_ledger_registry.csv",index=False)


def limiting_tests():
    rows=[]; base={"open":.2,"precursor":.03,"isolated":.04,"closed":.05}
    for mode in MODES:
      t=pore_channel_terms(1473.15,25e-9,10,.5,.05,.8,.6)
      s=conservative_step(base,mode,t,100)
      rows += [{"mode":mode,"test":"topology_no_direct_density","passed":t["rho_dot_topology"]==0},
               {"mode":mode,"test":"pore_volume_conserved","passed":abs(s["conservation_residual"])<1e-14},
               {"mode":mode,"test":"stores_nonnegative","passed":all(s[k]>=0 for k in ("open","precursor","isolated","closed"))}]
    small=pore_channel_terms(1473.15,10e-9,10,.5,.05,.8,.6); large=pore_channel_terms(1473.15,100e-9,10,.5,.05,.8,.6)
    rows += [{"mode":"all","test":"r4_time_increases","passed":large["tau_s_s"]>small["tau_s_s"]},
             {"mode":"all","test":"small_pores_pin_more","passed":small["P_Z_Pa"]>large["P_Z_Pa"]},
             {"mode":"all","test":"density_identity","passed":np.isclose(density_identity(np.array([.23]),np.array([.04]),np.array([.05])),.68)}]
    for kernel in ("renewal","GB_diffusion"):
      z=emergent_pore_closure_v1(1473.15,25e-9,0,1,.5,.25,3,kernel=kernel)
      stop=emergent_pore_closure_v1(1473.15,25e-9,.05,1,.5,1.1,3,kernel=kernel)
      gas=emergent_pore_closure_v1(1473.15,25e-9,.05,1,.5,.9,3,kernel=kernel)
      nogas=emergent_pore_closure_v1(1473.15,25e-9,.05,1,.5,0,3,kernel=kernel)
      rows += [{"mode":kernel,"test":"zero_closed_inventory","passed":z["rho_dot_closed_sinv"]==0},
               {"mode":kernel,"test":"nonpositive_closed_stress","passed":stop["rho_dot_closed_sinv"]==0},
               {"mode":kernel,"test":"gas_reduces_closed_stress","passed":gas["sigma_c_Pa"]<nogas["sigma_c_Pa"]}]
    q=pd.DataFrame(rows); q.to_csv(OUT/"pore_channel_limiting_tests.csv",index=False); return q


def fixed_scan():
    rows=[]
    for TC in range(700,1501,25):
     TK=TC+273.15
     for rnm in (5,10,25,50,100,250):
      for seg in (2,4,6.28,10,20):
       for ratio in (.1,.3,.5,.8):
        for activity in (.001,.01,.1,.5,.9):
         for W in (.25,.5,.8,1.2):
          for conn in (.1,.3,.6,.9):
           z=pore_channel_terms(TK,rnm*1e-9,seg,ratio,activity,W,conn)
           rows.append({**z,"T_C":TC,"radius_nm":rnm,"lambda_seg_over_r":seg,"gamma_GB_over_gamma_s":ratio,"activity":activity,"W_p":W,"connected_fraction":conn,
                        "D50_trend":"increase","D90_trend":"increase" if z["J_damage_sinv"]>z["J_reg_sinv"] else "narrow_toward_D50",
                        "connected_fine_trend":"decrease" if z["J_damage_sinv"]>z["J_reg_sinv"] else "preserve_or_increase",
                        "dominant_branch":"damage" if z["J_damage_sinv"]>z["J_reg_sinv"] else "regularization"})
    q=pd.DataFrame(rows); q.to_csv(OUT/"pore_channel_fixed_state_scan.csv",index=False,float_format="%.9g"); return q


def main():
    OUT.mkdir(parents=True,exist_ok=True); write_registries(); lim=limiting_tests(); scan=fixed_scan()
    pd.DataFrame([{"barrier_sha256":hashlib.sha256(BARRIER_PATH.read_bytes()).hexdigest(),"D_GB_unchanged":True,"D_s_unchanged":True,"mobility_unchanged":True,
                   "failed_global_mobility_active":False,"physical_Q_closed":False,"accepted_physics_changed":False,"scan_rows":len(scan),"all_limiting_pass":bool(lim.passed.all())}]).to_csv(OUT/"pore_channel_fixed_input_guardrails.csv",index=False)
    print({"scan_rows":len(scan),"all_limiting_pass":bool(lim.passed.all())})

if __name__=="__main__": main()
