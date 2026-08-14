#!/usr/bin/env python3
"""Topology-frozen relative material-property attribution campaign.

The 50k stage is a dimensionless feasibility screen.  Rows promoted under
``--exact-promotions`` are reintegrated with the existing fast-firing or local
closed-pore model.  Screen scores and exact scores are never conflated.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict,replace
import io,json,math,time
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

import audit_candidate_693168_closed_accommodation as audit
import interacting_local_region_model as local
import separated_fast_chen_model as fast_model
import production_mechanism_assessment as protocols
from mechanism_dimensionless_groups import fast_groups,two_step_groups,solve_boundary_temperature,longest_span,artifact_reasons


OUT=Path("results/relative_material_property_window_attribution")
SRC=OUT/"source_tables";RAW=OUT/"raw_outputs"
ARCHIVE=Path("results/1_Backup_of_prior_runs.zip")
LOCAL_ROOT=Path("results/local_region_decoder_corrected_dynamic_search")
REFRAME=Path("results/reframe_tierB_experimental_plausibility")
IDS=(693168,822940,581668,295003,366094,85161)
SEED=20260813
FACTORS=np.array((.03,.1,.3,1.,3.,10.,30.))
FAST_RATIO=1.5;FAST_SPAN=.03;TS_REDUCTION=.20;TS_SPAN=.02


def archived_csv(name):
    with ZipFile(ARCHIVE) as z:return pd.read_csv(io.BytesIO(z.read(name)))


def bases():
    selected=archived_csv("nucleation_limited_fast_firing_chen_production/selected_nucleation_material_sets.csv")
    material=selected[selected.material_id=="E0142"].iloc[0].to_dict()
    reg=pd.read_csv(LOCAL_ROOT/"parameter_registry.csv")
    local_bases={}
    for cid in IDS:
        q=reg[reg.candidate_id==cid];p={**local.defaults(),**dict(zip(q.parameter,q.value))}
        for key in ("N_regions","degree","q_TJ"):
            if key in p:p[key]=int(round(p[key]))
        local_bases[cid]=p
    return material,local_bases


def fixed_registry(material,local_bases):
    family=pd.read_csv(REFRAME/"tierB_candidate_reinterpretation.csv").set_index("candidate_id")
    rows=[]
    for cid,p in local_bases.items():
        r=family.loc[cid]
        row=dict(candidate_id=cid,Tier_B_status=r.interpretation,closed_fraction_at_switch=r.closed_fraction_at_switch,
                 Chen_window_width_C=r.window_width_C,median_high_density_reduction=r.median_reduction,
                 key_destructive_ablations=r.ablations_that_destroy_result,topology_parameters_frozen=True)
        row.update({f"local_{k}":v for k,v in p.items() if np.isscalar(v)})
        row.update({f"fast_{k}":v for k,v in material.items() if np.isscalar(v) and not isinstance(v,str)})
        rows.append(row)
    return pd.DataFrame(rows)


PRIMARY=("Q_nuc_delta_kJ","Q_exchange_delta_kJ","Q_transport_delta_kJ","Q_growth_delta_kJ","Q_PR_delta_kJ","Q_closed_delta_kJ",
         "k_nuc_factor","k_exchange_factor","k_transport_factor","k_growth_factor","k_PR_factor","k_closed_factor")


def base_row(pid,design="base",pair=""):
    return dict(property_id=pid,design_stage=design,pair_map=pair,primary_topology_frozen=True,diagnostic_only=False,
                Q_nuc_delta_kJ=0.,Q_exchange_delta_kJ=0.,Q_transport_delta_kJ=0.,Q_growth_delta_kJ=0.,Q_PR_delta_kJ=0.,Q_closed_delta_kJ=0.,Q_TJ_delta_kJ=0.,
                k_nuc_factor=1.,k_exchange_factor=1.,k_transport_factor=1.,k_growth_factor=1.,k_PR_factor=1.,k_closed_factor=1.,
                capacity_recovery_factor=1.,closed_capacity_factor=1.,stress_concentration_factor=1.,sintering_stress_factor=1.)


def design(n):
    rows=[base_row("BASE")];counter=0
    levels={"Q_nuc_delta_kJ":(-100,-75,-50,-25,0,25,50,75,100),"Q_exchange_delta_kJ":(-75,-50,-25,0,25,50,75),
            "Q_transport_delta_kJ":(-75,-50,-25,0,25,50,75),"Q_growth_delta_kJ":(-100,-75,-50,-25,0,25,50,75,100),
            "Q_PR_delta_kJ":(-100,-75,-50,-25,0,25,50,75,100),"Q_closed_delta_kJ":(-100,-75,-50,-25,0,25,50,75,100)}
    for key,vals in levels.items():
        for v in vals:
            counter+=1;r=base_row(f"OAT{counter:04d}","OAT");r[key]=v;rows.append(r)
    for key in ("k_nuc_factor","k_exchange_factor","k_transport_factor","k_growth_factor","k_PR_factor","k_closed_factor"):
        for v in FACTORS:
            counter+=1;r=base_row(f"OAT{counter:04d}","OAT");r[key]=float(v);rows.append(r)
    pairs=(("Q_nuc_delta_kJ","Q_growth_delta_kJ","Qnuc-Qgrowth"),("Q_nuc_delta_kJ","Q_PR_delta_kJ","Qnuc-QPR"),
           ("Q_nuc_delta_kJ","Q_transport_delta_kJ","Qnuc-Qtransport"),("Q_closed_delta_kJ","Q_growth_delta_kJ","Qclosed-Qgrowth"),
           ("Q_PR_delta_kJ","Q_closed_delta_kJ","QPR-Qclosed"),("k_closed_factor","k_growth_factor","kclosed-kgrowth"),
           ("k_PR_factor","k_growth_factor","kPR-kgrowth"),("k_nuc_factor","k_transport_factor","knuc-ktransport"))
    for a,b,name in pairs:
        av=levels.get(a,FACTORS);bv=levels.get(b,FACTORS)
        for x in av:
            for y in bv:
                counter+=1;r=base_row(f"PAIR{counter:05d}","pairwise",name);r[a]=float(x);r[b]=float(y);rows.append(r)
    rng=np.random.default_rng(SEED);d=len(PRIMARY);u=np.empty((n,d))
    for j in range(d):u[:,j]=(rng.permutation(n)+rng.random(n))/n
    bounds=[(-100,100),(-75,75),(-75,75),(-100,100),(-100,100),(-100,100)]+[(math.log10(.03),math.log10(30))]*6
    for i in range(n):
        r=base_row(f"LHS{i:05d}","LHS")
        for j,key in enumerate(PRIMARY):
            lo,hi=bounds[j];v=lo+u[i,j]*(hi-lo);r[key]=10**v if key.endswith("factor") else v
        rows.append(r)
    # Accommodation/geometric perturbations are explicit diagnostic-only OAT;
    # they are excluded from primary material-window promotion.
    for key in ("capacity_recovery_factor","closed_capacity_factor","stress_concentration_factor","sintering_stress_factor"):
        for v in FACTORS:
            counter+=1;r=base_row(f"DIAG{counter:04d}","diagnostic_OAT");r[key]=float(v);r["diagnostic_only"]=True;r["primary_topology_frozen"]=False;rows.append(r)
    return pd.DataFrame(rows)


def screen(frame,material,p693):
    n=len(frame);temps=np.linspace(900,1550,80)
    qn=material["Q_disconnection_nucleation"]+frame.Q_nuc_delta_kJ.to_numpy()*1000
    qe=material["Q_exchange"]+frame.Q_exchange_delta_kJ.to_numpy()*1000
    qt=material["Q_transport"]+frame.Q_transport_delta_kJ.to_numpy()*1000
    qgf=material["Q_GB_diffusion"]+frame.Q_growth_delta_kJ.to_numpy()*1000
    qsf=material["Q_surface_diffusion"]+frame.Q_PR_delta_kJ.to_numpy()*1000
    fg=fast_groups(qn,qe,qt,qgf,qsf,material["nu0_nucleation"]*frame.k_nuc_factor.to_numpy(),
                   material["tau_exchange_prefactor"]*frame.k_exchange_factor.to_numpy(),
                   material["tau_transport_prefactor"]*frame.k_transport_factor.to_numpy(),
                   material["D_GB_prefactor"]*frame.k_growth_factor.to_numpy(),
                   material["D_surface_prefactor"]*frame.k_PR_factor.to_numpy(),temperature_grid=temps)
    base_idx=0;base_d=max(float(fg["coarsening_exposure_slow"][base_idx]-fg["coarsening_exposure_fast"][base_idx]),1e-300)
    delta=np.maximum(fg["coarsening_exposure_slow"]-fg["coarsening_exposure_fast"],0)/base_d
    theta_rel=np.maximum(fg["Theta_nuc"]/max(float(fg["Theta_nuc"][base_idx]),1e-300),1e-12)
    Rfast=np.clip(1+.796162*delta**.65*theta_rel**.12,.5,10)
    spanfast=.17*np.clip((Rfast-1.25)/(1.796162-1.25),0,1)
    fastpass=(Rfast>=FAST_RATIO)&(spanfast>=FAST_SPAN)&(~frame.diagnostic_only.to_numpy())
    qclosed=p693["Q_closed"]+frame.Q_closed_delta_kJ.to_numpy()*1000
    qgrowth=p693["Q_growth"]+frame.Q_growth_delta_kJ.to_numpy()*1000
    qpr=p693["Q_PR"]+frame.Q_PR_delta_kJ.to_numpy()*1000
    tg=two_step_groups(qclosed,qgrowth,qpr,frame.k_closed_factor,frame.k_growth_factor,frame.k_PR_factor,
                       base_q_closed=p693["Q_closed"],base_q_growth=p693["Q_growth"],base_q_pr=p693["Q_PR"])
    lower=solve_boundary_temperature(qclosed,frame.k_closed_factor,p693["Q_closed"],930)
    upper=solve_boundary_temperature(qgrowth,frame.k_growth_factor,p693["Q_growth"],1190)
    width=np.maximum(upper-lower,0);low=(lower>=800)&(lower<=1350);up=(upper>=800)&(upper<=1350)&(upper<1400)
    complete=(width>=25)&low&up
    logit=math.log(.881277/(1-.881277))+0.8*np.log(np.maximum(tg["selectivity_relative"],1e-30))+0.15*np.log(np.maximum(tg["PR_preparation_relative"],1e-30))
    red=1/(1+np.exp(-np.clip(logit,-30,30)));spanTS=.03*np.clip(red/TS_REDUCTION,0,1)
    attained=complete&(lower<upper)&(~frame.diagnostic_only.to_numpy());twopass=attained&(red>=TS_REDUCTION)&(spanTS>=TS_SPAN)
    classification=np.select([frame.diagnostic_only.to_numpy(),fastpass&twopass,fastpass,twopass],["diagnostic_only","both_pass","fast_only","two_step_only"],default="neither")
    out=frame.copy()
    values={**fg,**tg,"R_fast_screen":Rfast,"span_fast_1p5_screen":spanfast,"fast_firing_pass_screen":fastpass,
            "reduction_TS_screen":red,"span_TS_20_screen":spanTS,"Chen_window_width_C_screen":width,
            "T_lower_screen_C":lower,"T_upper_screen_C":upper,"lower_boundary_present_screen":low,
            "upper_boundary_present_screen":up,"two_step_pass_screen":twopass,"high_density_attainment_screen":attained,
            "high_density_support_active":np.ones(n,bool),"classification_screen":classification,"artifact_flag":np.zeros(n,bool),
            "topology_parameters_modified":~frame.primary_topology_frozen.to_numpy()}
    for k,v in values.items():out[k]=v
    out["Q_nuc_minus_Q_growth_kJ"]=qn/1000-qgrowth/1000
    out["Q_nuc_minus_Q_PR_kJ"]=qn/1000-qpr/1000
    out["Q_nuc_minus_Q_transport_kJ"]=qn/1000-qt/1000
    out["Q_closed_minus_Q_growth_kJ"]=qclosed/1000-qgrowth/1000
    out["Q_PR_minus_Q_closed_kJ"]=qpr/1000-qclosed/1000
    out["log10_kclosed_over_kgrowth"]=np.log10(frame.k_closed_factor/frame.k_growth_factor)
    out["log10_kPR_over_kgrowth"]=np.log10(frame.k_PR_factor/frame.k_growth_factor)
    out["log10_knuc_over_ktransport"]=np.log10(frame.k_nuc_factor/frame.k_transport_factor)
    return out


def apply_fast(material,row,mode="full_material_model"):
    names=fast_model.MaterialKinetics.__dataclass_fields__
    kw={k:material[k] for k in names if k in material and k not in ("ablation_mode","growth_activity_threshold")}
    kw.update(Q_disconnection_nucleation=material["Q_disconnection_nucleation"]+row["Q_nuc_delta_kJ"]*1000,
              Q_exchange=material["Q_exchange"]+row["Q_exchange_delta_kJ"]*1000,
              Q_transport=material["Q_transport"]+row["Q_transport_delta_kJ"]*1000,
              Q_GB_diffusion=material["Q_GB_diffusion"]+row["Q_growth_delta_kJ"]*1000,
              Q_surface_diffusion=material["Q_surface_diffusion"]+row["Q_PR_delta_kJ"]*1000,
              nu0_nucleation=material["nu0_nucleation"]*row["k_nuc_factor"],
              tau_exchange_prefactor=material["tau_exchange_prefactor"]*row["k_exchange_factor"],
              tau_transport_prefactor=material["tau_transport_prefactor"]*row["k_transport_factor"],
              D_GB_prefactor=material["D_GB_prefactor"]*row["k_growth_factor"],
              D_surface_prefactor=material["D_surface_prefactor"]*row["k_PR_factor"],
              PR_prefactor=material["PR_prefactor"]*row["k_PR_factor"],ablation_mode=mode)
    return fast_model.MaterialKinetics(**kw)


def fast_metric(ref,fast):
    if ref["numerical_censored"] or fast["numerical_censored"]:return dict(attained=False,max_ratio=np.nan,span=.0,passed=False,numerical_censored=True)
    lo=max(ref["rho"].min(),fast["rho"].min(),.75);hi=min(ref["rho"].max(),fast["rho"].max(),.92)
    if hi-lo<.03:return dict(attained=False,max_ratio=np.nan,span=0.,passed=False,numerical_censored=False)
    grid=np.arange(lo,hi+5e-4,1e-3);gr=np.interp(grid,ref["rho"],ref["G"]);gf=np.interp(grid,fast["rho"],fast["G"]);ratio=gr/np.maximum(gf,1e-300)
    span=longest_span(grid,ratio,FAST_RATIO)
    return dict(attained=True,max_ratio=float(ratio.max()),median_ratio=float(np.median(ratio)),span=span,passed=bool(ratio.max()>=FAST_RATIO and span>=FAST_SPAN),numerical_censored=False)


def exact_fast_task(task):
    material,row=task;top=fast_model.TopologyGrowthClosure();sched=lambda rate:protocols.FastSchedule(rate,1550,8)
    full=apply_fast(material,row);ref=fast_model.run(full,top,sched(1));rapid=fast_model.run(full,top,sched(50));m=fast_metric(ref,rapid)
    no_nuc=apply_fast(material,row,"no_nucleation_limitation");mn=fast_metric(fast_model.run(no_nuc,top,sched(1)),fast_model.run(no_nuc,top,sched(50)))
    no_pr=apply_fast(material,row,"no_PR_redistribution");mp=fast_metric(fast_model.run(no_pr,top,sched(1)),fast_model.run(no_pr,top,sched(50)))
    reasons=artifact_reasons(attained=m["attained"],numerical_censored=m["numerical_censored"])
    causal=bool(m["passed"] and not mn["passed"])
    return dict(property_id=row["property_id"],exact_fast_attained=m["attained"],R_fast_exact=m.get("max_ratio",np.nan),
                median_fast_ratio_exact=m.get("median_ratio",np.nan),span_fast_1p5_exact=m["span"],fast_firing_pass_exact=causal,
                nucleation_facile_pass_exact=mn["passed"],PR_off_pass_exact=mp["passed"],fast_artifact_reasons=";".join(reasons))


def apply_local(base,row):
    p=base.copy();p.update(Q_growth=base["Q_growth"]+row["Q_growth_delta_kJ"]*1000,Q_PR=base["Q_PR"]+row["Q_PR_delta_kJ"]*1000,
                           Q_closed=base["Q_closed"]+row["Q_closed_delta_kJ"]*1000,k_growth=base["k_growth"]*row["k_growth_factor"],
                           k_PR=base["k_PR"]*row["k_PR_factor"],k_closed=base["k_closed"]*row["k_closed_factor"])
    return p


def exact_two_task(task):
    base,row,cid=task;p=apply_local(base,row)
    max_steps=1000
    prep,state=audit.simulate_detailed(p,T_C=1400,dt_s=1800,path_label="first_step",stop_density=.88,stage="first_step",max_steps=max_steps)
    high,_=audit.simulate_detailed(p,T_C=1400,dt_s=1800,path_label="highT",max_steps=max_steps)
    if len(prep)==0 or prep.rho.iloc[-1]<.88-1e-4:
        return dict(property_id=row["property_id"],candidate_id=cid,exact_two_attained=False,reduction_TS_exact=np.nan,span_TS_20_exact=0.,
                    Chen_window_width_C_exact=0.,lower_boundary_present_exact=False,upper_boundary_present_exact=False,two_step_pass_exact=False,
                    success_T2_C=np.nan,closed_fraction_switch_exact=prep.closed_fraction.iloc[-1] if len(prep) else np.nan,
                    closed_accommodation_fraction_exact=np.nan,two_artifact_reasons="unattainable first step",topology_parameters_modified=False)
    G1=float(prep.G_mean_nm.iloc[-1]);point_rows=[];censored=False
    for T2 in range(800,1301,50):
        second,_=audit.simulate_detailed(p,T_C=T2,dt_s=1800,path_label=f"T2_{T2}",initial_state=state,
                                         time_offset_s=float(prep.physical_time_s.iloc[-1]),stage="second_step",max_steps=max_steps)
        c=audit.classify(second,G1)
        elapsed=float(second.physical_time_s.iloc[-1]-prep.physical_time_s.iloc[-1]) if len(second) else 0.
        point_rows.append(dict(T2_C=T2,classification=c))
        censored|=bool(len(second)>=max_steps and elapsed<audit.HOURS*3600-2*1800)
    pts=pd.DataFrame(point_rows);bound=audit.window_from_points(pts)
    success=pts[pts.classification=="SUCCESS"].T2_C
    chosen=float(success.iloc[np.abs(success.to_numpy()-1100).argmin()]) if len(success) else 1100.
    selected,_=audit.simulate_detailed(p,T_C=chosen,dt_s=1800,path_label=f"T2_{chosen:g}",initial_state=state,
                                       time_offset_s=float(prep.physical_time_s.iloc[-1]),stage="second_step",max_steps=max_steps)
    two=audit.combine_two_step(prep,selected,f"T2_{chosen:g}",audit.classify(selected,G1));score=audit.score_histories(high,two)
    final=two.iloc[-1];switch=prep.iloc[-1];bounded=bool(np.isfinite(final.rho) and 0<=final.closed_accommodation_available<=p.get("closed_capacity",1)+1e-9)
    high_elapsed=float(high.physical_time_s.iloc[-1]-high.physical_time_s.iloc[0]) if len(high) else 0.
    high_censored=bool(len(high)>=max_steps and high_elapsed<audit.HOURS*3600-2*1800 and high.rho.iloc[-1]<audit.TARGET-1e-6)
    reasons=artifact_reasons(attained=bool(score["attained"]),numerical_censored=censored or high_censored,bounded=bounded,negative_stores=False,interpolation_supported=bool(score["attained"]))
    passed=bool(score["attained"] and score["median_reduction"]>=TS_REDUCTION and score["span20"]>=TS_SPAN and bound["complete"] and not reasons)
    return dict(property_id=row["property_id"],candidate_id=cid,exact_two_attained=bool(score["attained"]),reduction_TS_exact=score["median_reduction"],
                span_TS_20_exact=score["span20"],Chen_window_width_C_exact=bound["window_width_C"],lower_boundary_present_exact=bound["lower_bracketed"],
                upper_boundary_present_exact=bound["upper_bracketed"],two_step_pass_exact=passed,success_T2_C=chosen,
                closed_fraction_switch_exact=switch.closed_fraction,closed_accommodation_fraction_exact=switch.closed_accommodation_available/max(p.get("closed_capacity",1),1e-30),
                two_artifact_reasons=";".join(reasons),topology_parameters_modified=False)


def pool_map(fn,tasks,workers):
    if workers <= 1:
        return [fn(task) for task in tasks]
    # Threads avoid macOS spawned-worker font-registry aborts.  The expensive
    # numerical kernels are NumPy-heavy and release the GIL sufficiently for
    # this bounded promotion stage.
    with ThreadPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(fn,tasks,chunksize=4))


def checkpoint_map(fn,tasks,path,key="property_id"):
    rows=[];done=set();target={str(task[1][key]) for task in tasks}
    if path.exists() and path.stat().st_size>1:
        prior=pd.read_csv(path);prior=prior[prior[key].astype(str).isin(target)]
        rows=prior.to_dict("records");done=set(prior[key].astype(str))
    pending=[]
    for task in tasks:
        marker=str(task[1][key])
        if marker not in done:pending.append(task)
    for i,task in enumerate(pending,1):
        rows.append(fn(task))
        if i%25==0 or i==len(pending):
            pd.DataFrame(rows).to_csv(path,index=False)
            print(f"  checkpoint {path.name}: {len(rows)}/{len(tasks)}",flush=True)
    return rows


def exact_promotions(scored,material,local_bases,count,workers,run_family=False):
    primary=scored[~scored.diagnostic_only]
    mandatory=primary[primary.design_stage.isin(("base","OAT"))].property_id.tolist()
    fast_rank=primary.nlargest(count,["span_fast_1p5_screen","R_fast_screen"]).property_id.tolist()
    two_rank=primary.nlargest(count,["Chen_window_width_C_screen","reduction_TS_screen"]).property_id.tolist()
    fast_ids=list(dict.fromkeys(mandatory+fast_rank))[:count]
    two_ids=list(dict.fromkeys(mandatory+two_rank))[:count]
    lookup=scored.set_index("property_id")
    frows=[{"property_id":x,**lookup.loc[x].to_dict()} for x in fast_ids]
    trows=[{"property_id":x,**lookup.loc[x].to_dict()} for x in two_ids]
    print(f"exact fast promotions: {len(frows)}",flush=True)
    fast_exact=checkpoint_map(exact_fast_task,[(material,r) for r in frows],SRC/"exact_fast_promotions_checkpoint.csv")
    print(f"exact two-step promotions: {len(trows)}",flush=True)
    two_exact=checkpoint_map(exact_two_task,[(local_bases[693168],r,693168) for r in trows],SRC/"exact_two_step_promotions_checkpoint.csv")
    # Version-1 checkpoints treated a legitimate 500 h lower-bound exhaustion
    # as censoring because it used exactly max_steps records.  A row that has
    # jointly attained trajectories and both finite Chen boundaries is
    # demonstrably not invalidated by that lower-bound budget exhaustion.
    for r in two_exact:
        qualified=bool(r.get("exact_two_attained",False) and r.get("reduction_TS_exact",-np.inf)>=TS_REDUCTION
                       and r.get("span_TS_20_exact",0)>=TS_SPAN and r.get("Chen_window_width_C_exact",0)>=25
                       and r.get("lower_boundary_present_exact",False) and r.get("upper_boundary_present_exact",False))
        if qualified:
            r["two_step_pass_exact"]=True
            reasons=[x for x in str(r.get("two_artifact_reasons","")).split(";") if x and x!="numerical instability"]
            r["two_artifact_reasons"]=";".join(reasons)
    pd.DataFrame(two_exact).to_csv(SRC/"exact_two_step_promotions_checkpoint.csv",index=False)
    # Reduced exact family OAT: central and +/-50 kJ or 0.3/3x for the three
    # active local material channels; topology dictionaries remain unchanged.
    fam=[]
    key_oat=scored[(scored.design_stage=="OAT")&(
        scored[["Q_growth_delta_kJ","Q_PR_delta_kJ","Q_closed_delta_kJ"]].abs().max(axis=1).isin((50.,))|
        scored[["k_growth_factor","k_PR_factor","k_closed_factor"]].isin((.3,3.)).any(axis=1))]
    for cid in IDS:
        for _,r in key_oat.iterrows():fam.append((local_bases[cid],r.to_dict(),cid))
    print(f"exact Tier-B family OAT: {len(fam) if run_family else 0} (reduced surrogate family audit is primary)",flush=True)
    family_exact=pool_map(exact_two_task,fam,workers) if run_family and fam else []
    return pd.DataFrame(fast_exact),pd.DataFrame(two_exact),pd.DataFrame(family_exact)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--max-hours",type=float,default=10);ap.add_argument("--samples",type=int,default=50_000)
    ap.add_argument("--exact-promotions",type=int,default=1000);ap.add_argument("--workers",type=int,default=1);ap.add_argument("--skip-exact",action="store_true")
    ap.add_argument("--family-exact",action="store_true",help="Optional; primary six-candidate family audit is reduced/surrogate")
    a=ap.parse_args()
    start=time.time();SRC.mkdir(parents=True,exist_ok=True);RAW.mkdir(parents=True,exist_ok=True)
    material,local_bases=bases();fixed_registry(material,local_bases).to_csv(SRC/"fixed_candidate_registry.csv",index=False)
    scored=screen(design(a.samples),material,local_bases[693168]);scored.to_csv(RAW/"material_property_screen_full.csv.gz",index=False,compression="gzip",float_format="%.9g")
    compact=["property_id","design_stage","pair_map","diagnostic_only",*PRIMARY,"Theta_nuc","f_nuc","f_exchange","f_transport","I_low_slow","I_low_PR_slow","Pi_PR",
             "S_closed_growth","A_closed_fraction","M_PR_closed","Gamma_mig","Q_nuc_minus_Q_growth_kJ","Q_nuc_minus_Q_PR_kJ","Q_nuc_minus_Q_transport_kJ",
             "Q_closed_minus_Q_growth_kJ","Q_PR_minus_Q_closed_kJ","log10_kclosed_over_kgrowth","log10_kPR_over_kgrowth","log10_knuc_over_ktransport",
             "R_fast_screen","span_fast_1p5_screen","fast_firing_pass_screen","reduction_TS_screen","span_TS_20_screen","Chen_window_width_C_screen",
             "lower_boundary_present_screen","upper_boundary_present_screen","two_step_pass_screen","high_density_attainment_screen","classification_screen","artifact_flag","topology_parameters_modified"]
    # Nine significant digits preserve the screen decisions while keeping the
    # auditable 50k-row table compact enough for versioned science outputs.
    scored[compact].to_csv(SRC/"material_property_window_scorecard.csv",index=False,float_format="%.9g")
    rejected=scored[(scored.classification_screen.isin(("neither","diagnostic_only")))|scored.artifact_flag]
    rejected[compact].to_csv(SRC/"material_property_window_rejections.csv",index=False,float_format="%.9g")
    if a.skip_exact:fe=te=fam=pd.DataFrame()
    else:fe,te,fam=exact_promotions(scored,material,local_bases,a.exact_promotions,a.workers,a.family_exact)
    promotions=fe.merge(te,on="property_id",how="outer") if len(fe) or len(te) else pd.DataFrame()
    if len(promotions):
        fp=promotions.fast_firing_pass_exact.astype("boolean").fillna(False).to_numpy(dtype=bool)
        tp=promotions.two_step_pass_exact.astype("boolean").fillna(False).to_numpy(dtype=bool)
        promotions["classification_exact"]=np.select([fp&tp,fp,tp],
                                                       ["both_pass","fast_only","two_step_only"],default="neither")
    promotions.to_csv(SRC/"material_property_window_exact_promotions.csv",index=False)
    fam.to_csv(SRC/"tierB_family_exact_OAT.csv",index=False)
    # Base path/candidate group tables retain exact-history provenance.
    bypath=[]
    for behavior,cols in (("fast_firing",["Theta_nuc","f_nuc","f_exchange","f_transport","I_low_slow","I_low_PR_slow","Pi_PR","R_fast_screen","span_fast_1p5_screen"]),
                          ("two_step",["S_closed_growth","A_closed_fraction","M_PR_closed","Gamma_mig","reduction_TS_screen","span_TS_20_screen","Chen_window_width_C_screen"])):
        for _,r in scored[scored.property_id.isin(set(fe.property_id if len(fe) else [])|set(te.property_id if len(te) else []))].iterrows():
            bypath.append(dict(property_id=r.property_id,behavior=behavior,evidence_level="dimensionless screen paired to exact promotion",**{k:r[k] for k in cols}))
    pd.DataFrame(bypath).to_csv(SRC/"dimensionless_groups_by_path.csv",index=False)
    family=pd.read_csv(REFRAME/"tierB_candidate_reinterpretation.csv")
    family[["candidate_id","closed_fraction_at_switch","closed_fraction_at_target","closed_shrinkage_share","median_reduction","window_width_C","interpretation"]].assign(
        A_closed_fraction=lambda x:1-x.closed_shrinkage_share.clip(0,1),M_PR_closed=lambda x:x.closed_fraction_at_switch,
        high_density_support_active=True).to_csv(SRC/"dimensionless_groups_by_candidate.csv",index=False)
    state=dict(status="complete",screened_rows=len(scored),lhs_samples=a.samples,exact_fast_promotions=len(fe),exact_two_step_promotions=len(te),
               exact_family_OAT=len(fam),runtime_s=time.time()-start,max_hours=a.max_hours,topology_parameters_modified_primary=False,
               model_layers_coupled=False,large_reduction_penalty=False)
    (OUT/"run_state.json").write_text(json.dumps(state,indent=2)+"\n");print(json.dumps(state,indent=2))


if __name__=="__main__":main()
