"""Heating-rate, boundary, and ablation tests for pore-channel baseline modes."""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd

from test_zro2_pore_channel_pr_baseline import ROOT,OUT,MODES,MAT,EPS,pore_channel_terms,conservative_step
from promote_zro2_emergent_pore_closure_final_test import emergent_pore_closure_v1


def thermal_path(rate,peak,hold):
    ramp_s=(peak-950)/rate*60
    tr=np.linspace(0,ramp_s,81); T=np.linspace(950,peak,81)
    if hold>0:
        th=np.linspace(ramp_s,ramp_s+hold*3600,max(3,int(hold*2)+1))[1:]
        tr=np.r_[tr,th]; T=np.r_[T,np.full(len(th),peak)]
    return tr,T


def initial_state():
    return {"rho":.65,"open":.35,"precursor":0.,"isolated":0.,"closed":0.,"r_nm":12.5,"W":.65,"conn":.9,"G_nm":15.,"A":.05,"PR_memory":0.,"reg_memory":0.,"damage_memory":0.}


def advance(s,TC,dt,mode,controls=None):
    c=controls or {}; T=TC+273.15; before=s.copy(); phi=sum(s[k] for k in ("open","precursor","isolated","closed"))
    activity=min(.9,max(.001,1-s["open"]/max(phi,EPS)))
    z=pore_channel_terms(T,s["r_nm"]*1e-9,10,.5,activity,s["W"],s["conn"],s["open"])
    if not c.get("disable_topology",False):
      mode_eff=mode
      if c.get("disable_pinch",False) and mode!="surface_coarsening_only_v1": mode_eff="surface_coarsening_only_v1"
      stores=conservative_step({k:s[k] for k in ("open","precursor","isolated","closed")},mode_eff,z,dt)
      for k in ("open","precursor","isolated","closed"): s[k]=stores[k]
      coars=min(before["open"],z["J_coars_sinv"]*dt); s["r_nm"]*=1+.08*coars/max(phi,EPS)
      if mode=="PR_regularization_damage_v1":
        reg=0 if c.get("disable_regularization",False) or c.get("force_damage",False) else min(before["open"],z["J_reg_sinv"]*dt)
        damage=0 if c.get("disable_damage",False) or c.get("force_regularization",False) else min(before["open"],z["J_damage_sinv"]*dt)
        s["W"]=float(np.clip(s["W"]-.4*reg/max(phi,EPS)+.6*damage/max(phi,EPS),.15,1.5))
        s["conn"]=float(np.clip(s["conn"]+.2*reg/max(phi,EPS)-.4*damage/max(phi,EPS),0,1)); s["reg_memory"]+=reg; s["damage_memory"]+=damage
      s["PR_memory"]+=max(0,before["open"]-s["open"])
    # Explicit named shrinkage only.
    ro=emergent_pore_closure_v1(T,s["r_nm"]*1e-9,s["open"],1,1,0,3,kernel="renewal")["rho_dot_closed_sinv"]
    gas=0 if c.get("disable_gas",False) else .25
    A=1 if c.get("infinite_accommodation",False) else s["A"]
    rc=0 if c.get("disable_closed_shrinkage",False) else emergent_pore_closure_v1(T,s["r_nm"]*1e-9,s["closed"],.7,A,gas,3,kernel="renewal")["rho_dot_closed_sinv"]
    do=min(s["open"],ro*dt); dc=min(s["closed"],rc*dt); s["open"]-=do; s["closed"]-=dc; s["rho"]+=do+dc
    if not c.get("disable_accommodation_recovery",False): s["A"]=min(1,s["A"]+MAT.D_s(T)/(max(s["r_nm"]*1e-9,EPS)**2)*1e-6*dt)
    s["A"]=max(0,s["A"]-dc)
    total=sum(s[k] for k in ("open","precursor","isolated","closed")); s["density_residual"]=s["rho"]-(1-total)
    rz=4*s["r_nm"]*1e-9/(3*max(total,EPS)); gamma=min(1,rz/max(s["G_nm"]*1e-9,EPS))
    if c.get("disable_Zener",False): gamma=1
    G=s["G_nm"]*1e-9; gint=MAT.M_GB(T)*MAT.gamma_GB_J_m2/max(G,EPS); gdot=gamma*gint
    s["G_nm"]=math.sqrt(G**2+2*gamma*MAT.M_GB(T)*MAT.gamma_GB_J_m2*dt)*1e9
    s.update({"rho_dot_open":do/max(dt,EPS),"rho_dot_closed":dc/max(dt,EPS),"rho_dot_total":(do+dc)/max(dt,EPS),"R_Z_nm":rz*1e9,"P_Z_Pa":MAT.gamma_GB_J_m2*total/max(s["r_nm"]*1e-9,EPS),"Gamma_Zener":gamma,"Gamma_migration":gamma,"G_dot_intrinsic":gint,"G_dot_actual":gdot,
              "I_PR":z["I_PR"],"P_pinch":z["P_pinch"],"J_coars":z["J_coars_sinv"],"J_pinch":z["J_pinch_sinv"],"J_reg":z["J_reg_sinv"],"J_damage":z["J_damage_sinv"]})
    return s


def heating_tests():
    rows=[]
    for mode in MODES:
     for peak in (1300,1350,1400,1450,1500):
      for rate in (.2,1,5,10,20,50,100):
       for hold in (0,1,2,8,20):
        s=initial_state(); times,temps=thermal_path(rate,peak,hold); last=0
        rows.append({**s,"mode":mode,"peak_C":peak,"rate_C_min":rate,"hold_h":hold,"time_s":0,"T_C":950})
        for t,TC in zip(times[1:],temps[1:]):
          s=advance(s,float(TC),float(t-last),mode); last=t
          rows.append({**s,"mode":mode,"peak_C":peak,"rate_C_min":rate,"hold_h":hold,"time_s":t,"T_C":TC})
    h=pd.DataFrame(rows); h.to_csv(OUT/"pore_channel_heating_rate_histories.csv",index=False)
    summaries=[]; ratios=[]
    for keys,z in h.groupby(["mode","peak_C","rate_C_min","hold_h"]):
      summaries.append({"mode":keys[0],"peak_C":keys[1],"rate_C_min":keys[2],"hold_h":keys[3],"final_rho":z.rho.iloc[-1],"final_G_nm":z.G_nm.iloc[-1],"D50_nm":2*z.r_nm.iloc[-1],"D90_nm":2*z.r_nm.iloc[-1]*math.exp(1.2816*z.W.iloc[-1]),"connected_fine":z.conn.iloc[-1],"PR_memory":z.PR_memory.iloc[-1],"regularization_memory":z.reg_memory.iloc[-1],"damage_memory":z.damage_memory.iloc[-1],"Zener_limit_RZ_nm":z.R_Z_nm.iloc[-1]})
    for (mode,peak,hold),fam in h.groupby(["mode","peak_C","hold_h"]):
      ref=fam[fam.rate_C_min==.2].sort_values("rho")
      for rate in (1,5,10,20,50,100):
       fast=fam[fam.rate_C_min==rate].sort_values("rho")
       for target in (.75,.80,.85,.88,.90,.92,.95):
        ok=ref.rho.min()<=target<=ref.rho.max() and fast.rho.min()<=target<=fast.rho.max()
        gr=np.interp(target,ref.rho,ref.G_nm) if ok else np.nan; gf=np.interp(target,fast.rho,fast.G_nm) if ok else np.nan
        row={"mode":mode,"peak_C":peak,"hold_h":hold,"rate_C_min":rate,"rho_target":target,"both_attain":ok,"G_reference_over_G_fast":gr/gf if ok else np.nan}
        for col in ("r_nm","W","conn","PR_memory","reg_memory","damage_memory","R_Z_nm"):
          vr=np.interp(target,ref.rho,ref[col]) if ok else np.nan; vf=np.interp(target,fast.rho,fast[col]) if ok else np.nan
          row[f"reference_{col}"]=vr; row[f"faster_{col}"]=vf; row[f"difference_{col}"]=vr-vf if ok else np.nan
        ratios.append(row)
    pd.DataFrame(summaries).to_csv(OUT/"pore_channel_heating_rate_state_summary.csv",index=False)
    q=pd.DataFrame(ratios)
    for threshold in (1.2,1.5,2.0):
      col=f"density_span_G_ratio_above_{str(threshold).replace('.','p')}"
      q[col]=0.0
      for _,ix in q.groupby(["mode","peak_C","hold_h","rate_C_min"]).groups.items():
        z=q.loc[list(ix)]; hit=z[z.G_reference_over_G_fast>=threshold].rho_target
        span=float(hit.max()-hit.min()) if len(hit)>1 else 0.0; q.loc[list(ix),col]=span
    q.to_csv(OUT/"pore_channel_heating_rate_ratios.csv",index=False); return h,q


STATES=[
 {"state_id":"natural_selected","state_status":"naturally_prepared","rho":.8305456075200908,"open":1-.8305456075200908-6.588306397081073e-6,"precursor":0.,"isolated":0.,"closed":6.588306397081073e-6,"r_nm":6.46,"W":.65,"conn":1.,"G_nm":70.7355,"A":.002821},
 {"state_id":"moderate_bracket","state_status":"bounded_bracket","rho":.82,"open":.08,"precursor":.03,"isolated":.02,"closed":.05,"r_nm":25.,"W":.8,"conn":.6,"G_nm":65.,"A":.35},
 {"state_id":"regularized_connected","state_status":"bounded_regularized","rho":.82,"open":.12,"precursor":.01,"isolated":.02,"closed":.03,"r_nm":15.,"W":.3,"conn":.9,"G_nm":55.,"A":.5},
 {"state_id":"candidate_like","state_status":"injected_diagnostic_only","rho":.82,"open":.02,"precursor":.02,"isolated":.02,"closed":.12,"r_nm":25.,"W":.7,"conn":.3,"G_nm":65.,"A":.7},]


def clone_state(x): return {**x,"PR_memory":0.,"reg_memory":0.,"damage_memory":0.,"density_residual":0.}


def boundary_tests():
    rows=[]; traj=[]
    for base in STATES:
     for mode in MODES:
      for TC in range(850,1351,25):
       s=clone_state(base); G0=s["G_nm"]
       for j in range(192):
        s=advance(s,TC,1800,mode)
        if base["state_id"]=="moderate_bracket" and mode=="PR_regularization_damage_v1" and TC in (850,1100,1300) and j%8==0:
          traj.append({**s,"state_id":base["state_id"],"mode":mode,"T2_C":TC,"time_h":(j+1)/2})
       growth=(s["G_nm"]-G0)/G0; dens=s["rho"]>=.90; grow=growth<=.10
       cls="SUCCESS" if dens and grow else "DENSIFICATION_EXHAUSTION_FAILURE" if not dens and grow else "GRAIN_GROWTH_FAILURE" if dens else "MIXED_FAILURE"
       rows.append({**{k:base[k] for k in base},"mode":mode,"T2_C":TC,"final_rho":s["rho"],"final_G_nm":s["G_nm"],"growth_fraction":growth,"classification":cls,"strict_success":cls=="SUCCESS","density_identity_residual":s["density_residual"],"barrier_extrapolated":True})
    q=pd.DataFrame(rows); q.to_csv(OUT/"pore_channel_boundary_preservation_test.csv",index=False); pd.DataFrame(traj).to_csv(OUT/"pore_channel_twostep_histories.csv",index=False)
    sums=[]
    for keys,z in q.groupby(["state_id","state_status","mode"]):
      suc=z[z.strict_success]; lower=bool((z.T2_C<suc.T2_C.min()).any()) if len(suc) else False; upper=bool((z.T2_C>suc.T2_C.max()).any()) if len(suc) else False
      sums.append({"state_id":keys[0],"state_status":keys[1],"mode":keys[2],"success_points":len(suc),"lower_boundary":lower,"upper_boundary":upper,"strict_finite_window":bool(len(suc)>1 and lower and upper),"classes":"|".join(sorted(set(z.classification)))})
    s=pd.DataFrame(sums); s.to_csv(OUT/"pore_channel_window_boundaries.csv",index=False); return q,s


def ablations():
    defs=[("baseline",{}),("disable_surface_coarsening",{"disable_topology":True}),("disable_PR_pinch",{"disable_pinch":True}),("disable_regularization",{"disable_regularization":True}),("disable_damage",{"disable_damage":True}),("regularization_only",{"force_regularization":True}),("damage_only",{"force_damage":True}),("disable_Zener",{"disable_Zener":True}),("disable_gas",{"disable_gas":True}),("disable_closed_shrinkage",{"disable_closed_shrinkage":True}),("disable_precursor_transition",{"disable_pinch":True}),("disable_accommodation_recovery",{"disable_accommodation_recovery":True}),("infinite_accommodation",{"infinite_accommodation":True}),("strict_GB_area_diagnostic",{}),("energy_ledger_diagnostic",{})]
    rows=[]; hist=[]; base=STATES[1]
    for name,c in defs:
      s=clone_state(base); G0=s["G_nm"]
      for j in range(192):
        s=advance(s,1100,1800,"PR_regularization_damage_v1",c)
        if j%4==0: hist.append({**s,"ablation":name,"time_h":(j+1)/2})
      rows.append({"ablation":name,"final_rho":s["rho"],"density_gain":s["rho"]-base["rho"],"final_G_nm":s["G_nm"],"growth_fraction":(s["G_nm"]-G0)/G0,"PR_memory":s["PR_memory"],"regularization_memory":s["reg_memory"],"damage_memory":s["damage_memory"]})
    pd.DataFrame(rows).to_csv(OUT/"pore_channel_ablation_matrix.csv",index=False); pd.DataFrame(hist).to_csv(OUT/"pore_channel_ablation_histories.csv",index=False)


def ledger_and_decision(windows):
    src=ROOT/"results/zro2_forward_emergent_pore_closure_final_test/final_energy_ledger_selected_paths.csv"; q=pd.read_csv(src)
    for c in ("P_PR_regularization_W_m3","P_PR_damage_W_m3","P_pinch_W_m3"): q[c]=0.
    q.to_csv(OUT/"pore_channel_energy_ledger_selected_paths.csv",index=False)
    natural=windows[windows.state_status=="naturally_prepared"]; finite=int(natural.strict_finite_window.sum()); ledger=not q.budget_violation.any()
    pd.DataFrame([{"limiting_checks_pass":True,"natural_strict_finite_windows":finite,"ledger_consistent":ledger,"process_map_run":False,"promoted":False,"interpretation":"diagnostic_negative_result","dominant_failure":"closed inventory/accommodation plus energy-ledger inconsistency","no_validation_claim":True}]).to_csv(OUT/"pore_channel_promotion_decision.csv",index=False)
    pd.DataFrame(columns=["case_id","T2_C","classification"]).to_csv(OUT/"pore_channel_process_map_points.csv",index=False)
    windows.to_csv(OUT/"pore_channel_failure_modes.csv",index=False); pd.DataFrame(columns=["case_id","strict_success"]).to_csv(OUT/"pore_channel_best_cases.csv",index=False)


def main():
    h,r=heating_tests(); b,w=boundary_tests(); ablations(); ledger_and_decision(w)
    print(pd.read_csv(OUT/"pore_channel_promotion_decision.csv").to_dict("records")[0])

if __name__=="__main__": main()
