#!/usr/bin/env python3
"""Focused open/closed handoff audit. Conditional diagnostics, not validation."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
import hashlib,json,os
import numpy as np
import pandas as pd

from zro2_forward.conditioned_950c import ConditionedTwoStep,run_path,matched
from zro2_forward.resolved_rules import ResolvedRuleModel,conservative_adjacent_PR
from zro2_forward.schedules import RampNoHold
from run_zro2_forward_resolved_rules import params,initial,classify,window_rows,TARGET_RHO,TARGET_G

OUT=Path("results/zro2_forward_open_closed_rate_handoff_audit")

def design():
    rows=[("resolved_default","resolved_default",1,1,True),
          ("diagnostic_open_recovery","diagnostic_open_recovery",1,1,True)]
    rows += [(f"closed_rate_{x}x","closed_rate_boost_only",x,1,True) for x in (3,10,30,100)]
    rows += [(f"closed_inventory_{x}x","closed_inventory_boost_only",1,x,True) for x in (3,10,30)]
    rows += [("balanced_handoff","balanced_handoff",1,1,True),
             ("candidate_state_injection_diagnostic","candidate_state_injection_diagnostic",1,1,False)]
    return pd.DataFrame(rows,columns=["mode_id","open_closed_handoff_mode","closed_rate_factor","closed_inventory_factor","forward_prediction_eligible"])

class InjectionModel(ResolvedRuleModel):
    """Externally supplied switch state; diagnostic only and never a prediction."""
    def __init__(self,*args,injection_rho=None,**kwargs):super().__init__(*args,**kwargs);self.injection_rho=injection_rho;self.injected=False
    def step(self,state,T_K,dt_s):
        if self.injection_rho is not None and state.rho>=self.injection_rho and not self.injected:
            p=state.pores.copy();total=p.total;weights=p.phi_open/max(p.phi_open.sum(),1e-300)
            p.phi_closed=.65*total*weights;p.phi_open=.35*total*weights;p.phi_iso[:]=0
            state=replace(state,rho=1-p.total,pores=p,A_closed=.152178631,PR_memory=1.)
            self.injected=True
        return super().step(state,T_K,dt_s)

def make_model(row,injection_rho=None,**extra):
    q=params(open_closed_handoff_mode=row["open_closed_handoff_mode"],closed_rate_factor=float(row["closed_rate_factor"]),closed_inventory_factor=float(row["closed_inventory_factor"]),**extra)
    return InjectionModel(parameters=q,injection_rho=injection_rho) if row["mode_id"]=="candidate_state_injection_diagnostic" else ResolvedRuleModel(parameters=q)

def integrate(g,col):return float(np.trapezoid(g[col].fillna(0),g.t_s)) if len(g)>1 else 0.
def enrich(h,mode,path):
    h=h.copy();h.insert(0,"mode_id",mode);h.insert(1,"path",path)
    h["phi_open"]=h.open_fraction;h["phi_precursor"]=h.isolated_fraction;h["phi_closed"]=h.closed_fraction
    return h
def fixed(row):
    m=make_model(row);frames=[];fast=[];flux=[];two=[];twoflux=[]
    hs={}
    for rate,dt in ((5,120),(50,60)):
        f,h=run_path(m,RampNoHold(rate,1500,start_C=950),initial(),dt,dt,f"{row.mode_id}_{rate}C");h=enrich(h,row.mode_id,f"{rate}C");hs[rate]=h;frames.append(h)
        collapse=h[h.open_eligibility_eff<.1];onset=h[h.rho_dot_closed_sinv>1e-10]
        fast.append({"mode_id":row.mode_id,"rate_C_min":rate,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"final_D50_m":h.iloc[-1].pore_D50_m,"final_D90_m":h.iloc[-1].pore_D90_m,"final_fine_pore_fraction":h.iloc[-1].fine_pore_fraction,"open_eligibility_collapse_T_C":collapse.iloc[0].T_C if len(collapse) else np.nan,"closed_shrinkage_onset_T_C":onset.iloc[0].T_C if len(onset) else np.nan,"forward_prediction_eligible":row.forward_prediction_eligible})
        flux.append({"mode_id":row.mode_id,"rate_C_min":rate,"Delta_rho_open":integrate(h,"rho_dot_open_sinv"),"Delta_rho_closed":integrate(h,"rho_dot_closed_sinv"),"Delta_rho_total":integrate(h,"rho_dot_total_sinv"),"state_density_gain":h.iloc[-1].rho-h.iloc[0].rho})
    mm=matched(hs[5],hs[50]);fast[-1].update({"matched_smaller_grain_fraction":float((mm.G_5_over_G_50>1).mean()),"matched_smaller_D90_fraction":float((mm.pore_D90_m_50<mm.pore_D90_m_5).mean())})
    for tag,t2 in (("low_T2",900),("near_best",1200),("high_T2",1300)):
        m=make_model(row,injection_rho=.8)
        f,h=run_path(m,ConditionedTwoStep(1400,t2,.8,40),initial(),600,600,f"{row.mode_id}_{tag}");h=enrich(h,row.mode_id,tag);frames.append(h);last=h.iloc[-1]
        two.append({"mode_id":row.mode_id,"path":tag,"T2_C":t2,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"density_ok":f.rho>=TARGET_RHO,"grain_ok":f.G_m*1e6<=TARGET_G,"classification":classify(f,h[h.T_C.eq(t2)].iloc[0] if h.T_C.eq(t2).any() else last,h.T_C.eq(t2).any()),"final_phi_open":last.phi_open,"final_phi_precursor":last.phi_precursor,"final_phi_closed":last.phi_closed,"final_A_closed":last.A_closed,"final_PR_memory":last.PR_memory})
        twoflux.append({"mode_id":row.mode_id,"path":tag,"T2_C":t2,"Delta_rho_open":integrate(h,"rho_dot_open_sinv"),"Delta_rho_closed":integrate(h,"rho_dot_closed_sinv"),"Delta_rho_total":integrate(h,"rho_dot_total_sinv"),"cumulative_closure":integrate(h,"closure_rate"),"cumulative_PR":integrate(h,"PR_coarsening_flux")})
    return pd.DataFrame(fast),pd.DataFrame(flux),pd.DataFrame(two),pd.DataFrame(twoflux),pd.concat(frames)

def chen_task(task):
    row,T1,switch,T2,hold=task;m=make_model(row,injection_rho=switch);f,h=run_path(m,ConditionedTwoStep(T1,T2,switch,hold),initial(),1800,1e9,row["mode_id"])
    second=h[h.T_C.eq(T2)];att=len(second)>0;first=second.iloc[0] if att else h.iloc[-1];c=classify(f,first,att)
    return {"mode_id":row["mode_id"],"open_closed_handoff_mode":row["open_closed_handoff_mode"],"forward_prediction_eligible":row["forward_prediction_eligible"],"T1_C":T1,"T2_C":T2,"switch_density":switch,"hold_h":hold,"classification":c,"first_step_attained":att,"rho1":first.rho,"G1_um":first.G_um,"final_rho":f.rho,"final_G_um":f.G_m*1e6,"density_ok":f.rho>=TARGET_RHO,"grain_ok":f.G_m*1e6<=TARGET_G,"strict_success":c=="SUCCESS","closed_fraction_at_switch":first.closed_fraction,"A_closed_at_switch":first.A_closed,"PR_memory_at_switch":first.PR_memory,"map_resolution_C":25}

def boundaries(points):
    out=[]
    for mode,g in points.groupby("mode_id"):
        b=window_rows(g);b.insert(0,"mode_id",mode);out.append(b)
    return pd.concat(out,ignore_index=True)

def main():
    OUT.mkdir(parents=True,exist_ok=True);d=design();d.to_csv(OUT/"handoff_mode_design.csv",index=False)
    inherited=pd.DataFrame([
      ("fast_density_loss","reduced open shrinkage","0.0587 density loss; open integral -0.0634"),("direct_transfer","too small to explain loss","2.7e-4 pore fraction"),("fast_closed_shrinkage","negligible","1.6e-8 density"),("low_T2","barrier/closed-rate limited","inherited"),("high_T2","intrinsic-growth dominated","inherited"),("mobility","cannot rescue pathway","best gap -100 C; zero all-gate cases"),("candidate_693168","conditional comparator, not calibration","closed fraction 0.649 at switch")],columns=["evidence","interpretation","value"]);inherited["non_validation"]=True;inherited.to_csv(OUT/"inherited_rate_balance_summary.csv",index=False)
    F=[];X=[];T=[];TX=[];H=[]
    for _,r in d.iterrows():
        a,b,c,e,h=fixed(r);F.append(a);X.append(b);T.append(c);TX.append(e);H.append(h);print("fixed",r.mode_id,flush=True)
    fast=pd.concat(F);flux=pd.concat(X);two=pd.concat(T);twoflux=pd.concat(TX);hist=pd.concat(H)
    fast.to_csv(OUT/"fast_rate_handoff_summary.csv",index=False);flux.to_csv(OUT/"fast_rate_flux_integrals.csv",index=False);two.to_csv(OUT/"two_step_handoff_path_summary.csv",index=False);twoflux.to_csv(OUT/"two_step_handoff_flux_integrals.csv",index=False);hist.to_csv(OUT/"handoff_diagnostic_histories.csv",index=False)
    tasks=[(r.to_dict(),T1,s,T2,h) for _,r in d.iterrows() for T1 in (1250,1300,1350,1400,1450,1500) for s in (.75,.80,.85,.88,.90) for T2 in range(900,1301,25) if T2<T1 for h in (20,40)]
    workers=min(8,os.cpu_count() or 2);rows=[]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i,z in enumerate(ex.map(chen_task,tasks,chunksize=8),1):
            rows.append(z)
            if i%500==0:print("map",i,"/",len(tasks),flush=True)
    pts=pd.DataFrame(rows)
    # Local 10 C refinement only around observed strict points.
    ref=[]
    for (mode,T1,s,h),g in pts[pts.strict_success].groupby(["mode_id","T1_C","switch_density","hold_h"]):
        rr=d[d.mode_id.eq(mode)].iloc[0].to_dict();lo=max(900,int(g.T2_C.min()-25));hi=min(T1-10,int(g.T2_C.max()+25))
        ref.extend(chen_task((rr,T1,s,t,h))|{"map_resolution_C":10} for t in range(lo,hi+1,10))
    if ref:pts=pd.concat([pts,pd.DataFrame(ref)],ignore_index=True).drop_duplicates(["mode_id","T1_C","T2_C","switch_density","hold_h"],keep="last")
    b=boundaries(pts);pts.to_csv(OUT/"chen_handoff_classification_points.csv",index=False);b.to_csv(OUT/"chen_handoff_window_boundaries.csv",index=False)
    gap=[];gates=[]
    for mode,g in pts.groupby("mode_id"):
        bb=b[b.mode_id.eq(mode)];eligible=bool(d.loc[d.mode_id.eq(mode),"forward_prediction_eligible"].iloc[0]);wins=int(bb.finite_window.sum());succ=int(g.strict_success.sum())
        gap.append({"mode_id":mode,"strict_success_count":succ,"finite_window_count":wins,"min_gap_C":bb.boundary_gap_C.min(),"median_gap_C":bb.boundary_gap_C.median(),"max_gap_C":bb.boundary_gap_C.max(),"forward_prediction_eligible":eligible})
        f50=fast[(fast.mode_id.eq(mode))&(fast.rate_C_min.eq(50))].iloc[0];fx=flux[(flux.mode_id.eq(mode))&(flux.rate_C_min.eq(50))].iloc[0]
        gates.append({"mode_id":mode,"fast_density_substantially_recovered":f50.final_rho>=.94,"smaller_grain_sign_preserved":f50.get("matched_smaller_grain_fraction",0)>.5,"named_flux_identity":abs(fx.Delta_rho_open+fx.Delta_rho_closed-fx.Delta_rho_total)<1e-10,"closed_shrinkage_non_negligible":fx.Delta_rho_closed>1e-4,"finite_bracketed_window":wins>0,"forward_prediction_eligible":eligible,"all_gates":eligible and f50.final_rho>=.94 and f50.get("matched_smaller_grain_fraction",0)>.5 and fx.Delta_rho_closed>1e-4 and wins>0})
    gaps=pd.DataFrame(gap);gaps.to_csv(OUT/"chen_handoff_boundary_gap_summary.csv",index=False);pd.DataFrame(gates).to_csv(OUT/"pathway_gate_summary.csv",index=False)
    inj=pd.concat([fast[fast.mode_id.str.contains("candidate")],two[two.mode_id.str.contains("candidate")]],ignore_index=True,sort=False);inj["injected_closed_fraction"]=.65;inj["injected_PR_memory"]=1.;inj["injected_A_closed"]=.152178631;inj["diagnostic_only"]=True;inj["forward_prediction_eligible"]=False;inj.to_csv(OUT/"candidate_state_injection_diagnostic.csv",index=False)
    ab=[]
    for _,r in d.iterrows():
        win=bool(gaps.loc[gaps.mode_id.eq(r.mode_id),"finite_window_count"].iloc[0]>0 and r.forward_prediction_eligible)
        for name in ("no_PR_memory","no_closed_transition","no_closed_shrinkage","infinite_accommodation","no_handoff_coupling","no_migration_suppression"):
            ab.append({"mode_id":r.mode_id,"ablation":name,"parent_window_present":win,"ablation_interpretable":False,"ablation_result":"requires_followup_exact_ablation" if win else "not_interpretable_parent_has_no_window","reason_not_interpretable":"followup_required" if win else "not_interpretable_parent_has_no_window"})
    pd.DataFrame(ab).to_csv(OUT/"handoff_ablation_summary.csv",index=False)
    pd.DataFrame([("barrier extrapolation","unchanged"),("handoff factors","bounded diagnostic, not calibrated"),("candidate injection","external state sufficiency diagnostic only"),("mini-map","bounded, not broad calibration"),("validation","none")],columns=["limitation","status"]).to_csv(OUT/"unresolved_limitations.csv",index=False)
    state={"status":"complete_not_validated","strict_target_rho":TARGET_RHO,"strict_target_G_um":TARGET_G,"mode_count":len(d),"strict_success_count_forward_modes":int(gaps[gaps.forward_prediction_eligible].strict_success_count.sum()),"finite_window_count_forward_modes":int(gaps[gaps.forward_prediction_eligible].finite_window_count.sum()),"candidate_injection_successes":int(gaps.loc[gaps.mode_id.eq("candidate_state_injection_diagnostic"),"strict_success_count"].iloc[0]),"all_gate_mode_count":int(pd.DataFrame(gates).all_gates.sum()),"barrier_sha256":hashlib.sha256(Path("data/zro2/bicrystal_creep_barrier_export.json").read_bytes()).hexdigest(),"validation_claim":False}
    (OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n")
    print(gaps.to_string(index=False));print(pd.DataFrame(gates).to_string(index=False))
if __name__=="__main__":main()
