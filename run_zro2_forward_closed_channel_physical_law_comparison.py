#!/usr/bin/env python3
from __future__ import annotations
from dataclasses import asdict,replace
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import pandas as pd

from zro2_forward.closed_channel_laws import ALLOWED_CLOSED_CHANNEL_LAWS
from zro2_forward.conditioned_950c import make_pdf_conditioned_initial_state,ConditionedTwoStep,run_path,matched,BARRIER
from zro2_forward.pore_population import diagnostics
from zro2_forward.resolved_rules import ResolvedRuleModel,ResolvedRuleParameters,resolved_initial_state
from zro2_forward.schedules import RampNoHold,Iso

OUT=Path("results/zro2_forward_closed_channel_physical_law_comparison")
OLD=Path("results/zro2_forward_closed_channel_apparent_law_correction")
TARGET_RHO=.976;TARGET_G_UM=.29;FACTORS=(.03,.1,.3,1.,3.,10.,30.,100.)

def params(**kw):
    cal=pd.read_csv("results/zro2_forward_pdf_conditioned_950C_comparison/pdf_conditioned_calibrated_parameters.csv").query("calibration_mode=='density_plus_grain_trajectory_conditioned'").iloc[0]
    return replace(ResolvedRuleParameters(),site_density_multiplier=cal.site,surface_power_length2_m2=1e-19*cal.work,**kw)
def model(**kw):return ResolvedRuleModel(parameters=params(**kw))
def initial():return resolved_initial_state(make_pdf_conditioned_initial_state())
def integrate(h,col):return float(np.trapezoid(np.maximum(h[col],0),h.t_s)) if len(h)>1 else 0.
def classify(final,first,attained=True):
    if not attained:return "UNATTAINABLE_FIRST_STEP"
    if first.rho>=TARGET_RHO:return "INELIGIBLE_TARGET_ALREADY_REACHED"
    d=final.rho>=TARGET_RHO;g=final.G_m*1e6<=TARGET_G_UM
    return "SUCCESS" if d and g else "DENSIFICATION_EXHAUSTION_FAILURE" if (not d and g) else "GRAIN_GROWTH_FAILURE" if (d and not g) else "MIXED_FAILURE"
def state_from_row(row,rho=None,G_nm=None,closed_fraction=None,A=None,PR=None):
    s=initial();p=s.pores.copy()
    if "phi_open_json" in row:
        p.phi_open=np.array(json.loads(row.phi_open_json));p.phi_iso=np.array(json.loads(row.phi_iso_json));p.phi_closed=np.array(json.loads(row.phi_closed_json))
    total=1-(rho if rho is not None else float(row.rho))
    if closed_fraction is not None:
        shape=p.phi_open/max(p.phi_open.sum(),1e-300);p.phi_closed=shape*total*closed_fraction;p.phi_iso=np.zeros_like(shape);p.phi_open=shape*total*(1-closed_fraction)
    return replace(s,t_s=0.,rho=1-p.total,G_m=(G_nm*1e-9 if G_nm else float(row.G_m)),pores=p,A_closed=float(A if A is not None else row.A_closed),PR_memory=float(PR if PR is not None else row.PR_memory))
def prepared_states():
    m=model();_,h=run_path(m,ConditionedTwoStep(1400,1100,.88,1),initial(),600,600,"natural_prepare")
    pre=h[h.rho<.88];row=(pre.iloc[-1] if len(pre) else h.iloc[-1]);natural=state_from_row(row)
    candidate=state_from_row(row,rho=.88,G_nm=117,closed_fraction=.649,A=.152,PR=1.)
    return candidate,natural
def configurations():
    rows=[("resolved_proxy_current",3,1.,"inherited baseline"),("renewal_limited_closed_shrinkage",3,1.,"physical candidate"),("surface_diffusion_accommodation_only",4,1.,"conservative accommodation only")]
    for m in (3,4):rows += [("GB_diffusion_closed_shrinkage",m,f,"prefactor_uncertainty") for f in FACTORS]
    rows += [("gas_accommodation_limited",3,f,"prefactor_uncertainty") for f in FACTORS]
    rows += [("empirical_closed_rate_scale",4,f,"diagnostic empirical scale") for f in (1.,30.,100.,300.)]
    return rows
def tag(law,m,f):return f"{law}__m{m}__x{f:g}".replace(".","p")
def run_fixed(law,mexp,factor):
    ident=tag(law,mexp,factor);mod=model(closed_channel_law=law,closed_radius_exponent=mexp,closed_prefactor_factor=factor)
    rows=[];flux=[];histories={}
    for rate in (5,50):
        final,h=run_path(mod,RampNoHold(rate,1500,start_C=950),initial(),300 if rate==5 else 60,300,ident)
        histories[rate]=h;d=diagnostics(final.pores);closed_history=h.closed_fraction*(1-h.rho);rows.append({"config_id":ident,"law":law,"m":mexp,"prefactor_factor":factor,"path":f"PDF_conditioned_{rate}C_min","final_rho":final.rho,"final_G_um":final.G_m*1e6,"pore_D50_nm":d["pore_D50_m"]*1e9,"pore_D90_nm":d["pore_D90_m"]*1e9,"closed_inventory_formed":closed_history.max(),"final_A_closed":final.A_closed,"candidate_state_injected":False})
        closure=integrate(h,"closure_rate");closed=min(integrate(h,"rho_dot_closed_sinv"),closed_history.iloc[0]+closure)
        flux.append({"config_id":ident,"path":f"PDF_conditioned_{rate}C_min","Delta_rho_open":integrate(h,"rho_dot_open_sinv"),"Delta_rho_closed":closed,"state_density_gain":final.rho-.66,"integral_closure_transfer":closure})
    ratio=float(matched(histories[5],histories[50]).G_5_over_G_50.median())
    for row in rows:row["matched_density_grain_ratio_median"]=ratio
    return rows,flux
def scan_state(state,state_kind,configs):
    out=[]
    for law,mexp,factor,_ in configs:
        ident=tag(law,mexp,factor);mod=model(closed_channel_law=law,closed_radius_exponent=mexp,closed_prefactor_factor=factor)
        for T2 in range(900,1301,25):
            final,h=run_path(mod,Iso(T2,40),replace(state,t_s=0.),1800,1800,ident);first=replace(state)
            closed=min(integrate(h,"rho_dot_closed_sinv"),1-state.rho)
            out.append({"config_id":ident,"law":law,"m":mexp,"prefactor_factor":factor,"state_kind":state_kind,"T2_C":T2,"hold_h":40,"initial_rho":state.rho,"initial_G_nm":state.G_m*1e9,"initial_closed_fraction":state.pores.phi_closed.sum()/max(state.pores.total,1e-300),"initial_A_closed":state.A_closed,"initial_PR_memory":state.PR_memory,"final_rho":final.rho,"final_G_um":final.G_m*1e6,"Delta_rho_open":integrate(h,"rho_dot_open_sinv"),"Delta_rho_closed":closed,"classification":classify(final,first),"strict_success":classify(final,first)=="SUCCESS","candidate_state_injected":state_kind=="candidate_like_injected","diagnostic_only":state_kind=="candidate_like_injected"})
    return pd.DataFrame(out)
def window_rows(x):
    rows=[]
    for key,g in x.groupby(["config_id","law","T1_C","switch_density","hold_h"]):
        g=g.sort_values("T2_C");s=g[g.classification.eq("SUCCESS")];dens=g[g.final_rho.ge(TARGET_RHO)];grain=g[g.final_G_um.le(TARGET_G_UM)]
        lo=dens.T2_C.min() if len(dens) else np.nan;hi=grain.T2_C.max() if len(grain) else np.nan
        lp=bool(np.isfinite(lo) and (g.loc[g.T2_C<lo,"classification"].isin(["DENSIFICATION_EXHAUSTION_FAILURE","MIXED_FAILURE"])).any());up=bool(np.isfinite(hi) and (g.loc[g.T2_C>hi,"classification"].isin(["GRAIN_GROWTH_FAILURE","MIXED_FAILURE"])).any())
        finite=bool(len(s)>=2 and lp and up and s.T2_C.max()>s.T2_C.min() and s.T2_C.max()<key[2])
        rows.append({"config_id":key[0],"law":key[1],"T1_C":key[2],"switch_density":key[3],"hold_h":key[4],"lower_boundary_present":lp,"upper_boundary_present":up,"T_lower_density_C":lo,"T_upper_growth_C":hi,"success_count":len(s),"success_min_T2_C":s.T2_C.min() if len(s) else np.nan,"success_max_T2_C":s.T2_C.max() if len(s) else np.nan,"finite_window":finite,"window_width_C":s.T2_C.max()-s.T2_C.min() if finite else np.nan,"boundary_gap_C":hi-lo if np.isfinite(lo) and np.isfinite(hi) else np.nan})
    return pd.DataFrame(rows)
def mini_map(configs):
    rows=[]
    for law,mexp,factor,_ in configs:
      ident=tag(law,mexp,factor);mod=model(closed_channel_law=law,closed_radius_exponent=mexp,closed_prefactor_factor=factor)
      for T1 in (1250,1300,1350,1400,1450,1500):
       for switch in (.75,.80,.85,.88,.90):
        # Prepare a switch state once; T2 is causally downstream and must not
        # alter this first-step state.  This also avoids redundant ramp solves.
        _,prep=run_path(mod,ConditionedTwoStep(T1,900,switch,.25),initial(),900,900,ident)
        second=prep[prep.T_C.eq(900)];att=len(second)>0
        first=second.iloc[0] if att else prep.iloc[-1]
        switch_state=state_from_row(first) if att else None
        for hold in (20,40):
         for T2 in range(900,1301,25):
          if T2>=T1:continue
          if att:final,h=run_path(mod,Iso(T2,hold),replace(switch_state,t_s=0.,T_K=T2+273.15),1800,1800,ident)
          else:final,h=initial(),prep
          c=classify(final,first,att)
          rows.append({"config_id":ident,"law":law,"m":mexp,"prefactor_factor":factor,"T1_C":T1,"rho_switch":switch,"switch_density":switch,"T2_C":T2,"hold_h":hold,"classification":c,"first_step_attained":att,"rho1":first.rho,"G1_um":first.G_um,"final_rho":final.rho,"final_G_um":final.G_m*1e6,"strict_success":c=="SUCCESS","closed_fraction_at_switch":first.closed_fraction,"A_closed_at_switch":first.A_closed,"PR_memory_at_switch":first.PR_memory,"map_resolution_C":25})
    return pd.DataFrame(rows)
def main():
    OUT.mkdir(parents=True,exist_ok=True);configs=configurations()
    inherited=pd.DataFrame([{"item":"correction_branch","value":"model physics was not changed; barrier/diffusivities/mobility unchanged"},{"item":"Q_closed_app","value":"apparent diagnostic only, not a material property"},{"item":"closed_channel","value":"current rate too small"},{"item":"candidate_693168","value":"conditional comparator only"},{"item":"validation","value":"no validation claim"}]);inherited.to_csv(OUT/"inherited_closed_channel_correction_summary.csv",index=False)
    registry=[]
    descriptions={"resolved_proxy_current":"inherited reduced proxy","GB_diffusion_closed_shrinkage":"dimensional GB diffusion and capillary pressure","renewal_limited_closed_shrinkage":"serial nucleation/exchange/GB transport","surface_diffusion_accommodation_only":"bounded shape relaxation; no volume removal","gas_accommodation_limited":"bounded gas-pressure factor times GB transport","empirical_closed_rate_scale":"diagnostic Arrhenius reduced closure"}
    for law in ALLOWED_CLOSED_CHANNEL_LAWS:registry.append({"closed_channel_law":law,"description":descriptions[law],"direct_density_change":law!="surface_diffusion_accommodation_only","schedule_blind":True,"status":"default inherited" if law=="resolved_proxy_current" else "candidate; not calibrated; not validation","empirical":law=="empirical_closed_rate_scale"})
    pd.DataFrame(registry).to_csv(OUT/"closed_law_registry.csv",index=False)
    prow=[]
    for law,m,f,note in configs:
        for k,v in asdict(params(closed_channel_law=law,closed_radius_exponent=m,closed_prefactor_factor=f)).items():
            if k.startswith("closed_") or k in ("C_closed_GB","C_sigma_closed","C_transport_closed","k0_closed_emp_sinv","Q_closed_emp_J_mol","C_surface_accommodation","accommodation_max"):prow.append({"config_id":tag(law,m,f),"law":law,"parameter":k,"value":v,"audit_label":note,"validation_status":"not validated"})
    pd.DataFrame(prow).to_csv(OUT/"closed_law_parameter_values.csv",index=False)
    fixed_file=OUT/"fixed_path_closed_law_summary.csv";flux_file=OUT/"fixed_path_closed_law_flux_integrals.csv"
    if fixed_file.exists() and flux_file.exists() and pd.read_csv(fixed_file).get("matched_density_grain_ratio_median",pd.Series(dtype=float)).notna().all():fixed=pd.read_csv(fixed_file);flux=pd.read_csv(flux_file)
    else:
        fixed=[];flux=[]
        for c in configs:
            a,b=run_fixed(*c[:3]);fixed+=a;flux+=b
        fixed=pd.DataFrame(fixed);flux=pd.DataFrame(flux)
        fixed.to_csv(fixed_file,index=False);flux.to_csv(flux_file,index=False)
    cs_file=OUT/"candidate_state_closed_law_T2_scan.csv";ns_file=OUT/"natural_state_closed_law_T2_scan.csv"
    if cs_file.exists() and ns_file.exists():
        cs=pd.read_csv(cs_file);ns=pd.read_csv(ns_file)
        for frame,path in ((cs,cs_file),(ns,ns_file)):
            frame["Delta_rho_closed"]=np.minimum(frame.Delta_rho_closed,1-frame.initial_rho)
            frame["inventory_bounded_integral"]=True;frame.to_csv(path,index=False)
    else:
        cand,natural=prepared_states();cs=scan_state(cand,"candidate_like_injected",configs);ns=scan_state(natural,"naturally_prepared_resolved",configs)
        cs.to_csv(cs_file,index=False);ns.to_csv(ns_file,index=False)
    # A mode is plausible enough for the map only if a physical, naturally prepared path gains closed density without universal low-T success.
    phys=[c for c in configs if c[0] in ("resolved_proxy_current","GB_diffusion_closed_shrinkage","renewal_limited_closed_shrinkage","gas_accommodation_limited") and c[2]==1.]
    plausible=[]
    for c in phys:
        g=ns[ns.config_id.eq(tag(*c[:3]))];
        if c[0]=="resolved_proxy_current" or g.Delta_rho_closed.max()>1e-6:plausible.append(c)
    chen=mini_map(plausible);chen.to_csv(OUT/"closed_law_chen_classification_points.csv",index=False);bounds=window_rows(chen);bounds.to_csv(OUT/"closed_law_chen_window_boundaries.csv",index=False)
    gap=bounds.groupby(["config_id","law"]).agg(groups=("config_id","size"),strict_success_groups=("success_count",lambda x:int((x>0).sum())),finite_window_count=("finite_window","sum"),lower_boundary_count=("lower_boundary_present","sum"),upper_boundary_count=("upper_boundary_present","sum"),min_boundary_gap_C=("boundary_gap_C","min"),max_boundary_gap_C=("boundary_gap_C","max")).reset_index();gap.to_csv(OUT/"closed_law_boundary_gap_summary.csv",index=False)
    lower=[]
    for c in configs:
        g=cs[cs.config_id.eq(tag(*c[:3]))];low=g[g.T2_C.le(950)];mid=g[g.T2_C.between(975,1175)];high=g[g.T2_C.ge(1200)];loses=bool(low.final_rho.ge(TARGET_RHO).all());success=bool(mid.strict_success.any());growth=bool(high.final_G_um.gt(TARGET_G_UM).any())
        cls="universal_density_success" if g.final_rho.ge(TARGET_RHO).all() else "finite_window_candidate" if success and not loses and growth else "helps_density_preserves_lower_boundary" if g.Delta_rho_closed.max()>1e-4 and not loses else "helps_density_loses_lower_boundary" if loses else "still_growth_limited" if growth else "no_effect" if g.Delta_rho_closed.max()<1e-8 else "no_window"
        lower.append({"config_id":tag(*c[:3]),"law":c[0],"m":c[1],"prefactor_factor":c[2],"low_T2_density_exhaustion_limited":bool(low.classification.isin(["DENSIFICATION_EXHAUSTION_FAILURE","MIXED_FAILURE"]).any()),"intermediate_T2_attains_bounded_growth":success,"high_T2_growth_limited":growth,"lower_boundary_disappears":loses,"universal_density_success":bool(g.final_rho.ge(TARGET_RHO).all()),"preservation_class":cls})
    pd.DataFrame(lower).to_csv(OUT/"closed_law_lower_boundary_preservation.csv",index=False)
    comp=[]
    for c in configs:
        g=cs[cs.config_id.eq(tag(*c[:3]))];s=g[g.strict_success]
        comp.append({"model":tag(*c[:3]),"law":c[0],"closed_fraction_at_switch":.649,"A_closed_at_switch":.152,"closed_density_contribution":g.Delta_rho_closed.max(),"second_step_T2_success_interval":f"{s.T2_C.min()}-{s.T2_C.max()}" if len(s) else "none","lower_boundary_present":bool((g.T2_C<s.T2_C.min()).any()) if len(s) else False,"upper_boundary_present":bool((g.T2_C>s.T2_C.max()).any()) if len(s) else False,"conditional_comparator":False,"calibrated":False})
    comp.append({"model":"candidate_693168","law":"conditional comparator","closed_fraction_at_switch":.649,"A_closed_at_switch":.152,"closed_density_contribution":.244,"second_step_T2_success_interval":"conditional interval (inherited)","lower_boundary_present":True,"upper_boundary_present":True,"conditional_comparator":True,"calibrated":False});pd.DataFrame(comp).to_csv(OUT/"candidate693168_vs_physical_closed_laws.csv",index=False)
    best=gap.sort_values(["finite_window_count","strict_success_groups"],ascending=False).iloc[0]
    pd.DataFrame([{"decision":"retain_default","recommendation":"retain resolved_proxy_current until closed-pore size, pressure, and accommodation data constrain a physical prefactor"},{"decision":"best_screened_physical_configuration","recommendation":best.config_id},{"decision":"empirical_mode","recommendation":"use only to bound missing rate magnitude; Q_closed_emp is not a material property"},{"decision":"validation","recommendation":"none of these comparisons validates a law"}]).to_csv(OUT/"closed_law_recommended_next_action.csv",index=False)
    state={"branch":"codex/zro2-forward-closed-channel-physical-law-comparison","source_commit":"ebf82941db9b8789d018e5e7986b8e6750587e20","barrier_sha256":hashlib.sha256(BARRIER.read_bytes()).hexdigest(),"laws":list(ALLOWED_CLOSED_CHANNEL_LAWS),"mapped_configurations":[tag(*c[:3]) for c in plausible],"candidate_injection_diagnostic_only":True,"Q_closed_app_physical_input":False,"validation":False};(OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n")
    print(fixed.groupby("law").final_rho.agg(["min","max"]));print(gap.to_string(index=False))
if __name__=="__main__":main()
