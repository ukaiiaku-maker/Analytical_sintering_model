#!/usr/bin/env python3
"""Fixed-candidate dense histories and filled presentation maps for 693168.

This is deterministic map filling, not parameter search.  It imports the
frozen audit runner and never modifies model physics or candidate parameters.
"""
from pathlib import Path
from dataclasses import replace
from io import TextIOWrapper
import json, time, zipfile
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import pandas as pd

import audit_candidate_693168_closed_accommodation as audit
import interacting_local_region_model as local_model
import nucleation_fast_chen_production as nuc
import separated_fast_chen_model as fast_model
import production_mechanism_assessment as protocols

OUT=Path("results/publication_style_sintering_figures_693168")
SRC=OUT/"source_tables";DENSE=OUT/"dense_histories"
ARCHIVE=Path("results/1_Backup_of_prior_runs.zip")
TARGET=.98;TOL=.20;DT=1800

def write(path,frame):
    path.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(frame).to_csv(path,index=False,lineterminator="\n")

def prep(p,T1,switch,G0=None):
    q=dict(p)
    if G0 is not None:q["G0_nm"]=float(G0)
    h,s=audit.simulate_detailed(q,T_C=float(T1),dt_s=DT,path_label="first_step",
        stop_density=float(switch),stage="first_step")
    return q,h,s

def second(q,h,s,T2,switch,T1,G0):
    if h.rho.iloc[-1]<switch-1e-4:
        return dict(classification="UNATTAINABLE_FIRST_STEP",rho_final=h.rho.iloc[-1],G_final_nm=h.G_mean_nm.iloc[-1],growth_fraction=np.nan,target_attained=False,numerical_censor=False,closed_fraction_final=h.closed_fraction.iloc[-1],closed_shrinkage_contribution=0.,open_shrinkage_contribution=0.)
    if h.rho.iloc[-1]>=TARGET:
        return dict(classification="TARGET_REACHED_DURING_FIRST_STEP",rho_final=h.rho.iloc[-1],G_final_nm=h.G_mean_nm.iloc[-1],growth_fraction=0.,target_attained=True,numerical_censor=False,closed_fraction_final=h.closed_fraction.iloc[-1],closed_shrinkage_contribution=0.,open_shrinkage_contribution=0.)
    G1=float(h.G_mean_nm.iloc[-1])
    z,_=audit.simulate_detailed(q,T_C=float(T2),dt_s=DT,path_label=f"T2_{T2}",initial_state=s,time_offset_s=float(h.physical_time_s.iloc[-1]),stop_density=TARGET,stage="second_step")
    c=audit.classify(z,G1)
    return dict(classification=c,rho_final=float(z.rho.iloc[-1]),G_final_nm=float(z.G_mean_nm.iloc[-1]),growth_fraction=float(z.G_mean_nm.iloc[-1]/G1-1),target_attained=bool(z.rho.iloc[-1]>=TARGET-1e-6),numerical_censor=bool(not np.isfinite(z.rho.iloc[-1])),closed_fraction_final=float(z.closed_fraction.iloc[-1]),closed_shrinkage_contribution=float(z.cumulative_closed_pore_removed.iloc[-1]),open_shrinkage_contribution=float(z.cumulative_open_pore_removed.iloc[-1]))

def common(q,h,T1,T2,switch,G0):
    return dict(candidate_id=693168,T1_C=T1,T2_C=T2,rho_switch=switch,G0_nm=G0,
        rho1=float(h.rho.iloc[-1]),G1_nm=float(h.G_mean_nm.iloc[-1]),
        first_step_growth_fraction=float(h.G_mean_nm.iloc[-1]/h.G_mean_nm.iloc[0]-1),
        closed_fraction_at_switch=float(h.closed_fraction.iloc[-1]),
        accommodation_at_switch=float(h.closed_accommodation_factor.iloc[-1]),
        candidate_parameters_frozen=True,target_density=TARGET,time_budget_h=audit.HOURS,
        growth_tolerance=TOL)

def map_A(base):
    path=SRC/"chen_map_T1_T2_classification.csv"
    if path.exists():
        old=pd.read_csv(path)
        if len(old)>=779:return old
    rows=[]
    for T1 in range(1200,1551,25):
        q,h,s=prep(base,T1,.88,100)
        for T2 in range(800,min(T1-10,1350)+1,10):
            rows.append({**common(q,h,T1,T2,.88,100),**second(q,h,s,T2,.88,T1,100)})
        write(path,rows)
    return pd.DataFrame(rows)

def _route_B(args):
    base,G0,T1,switch=args;q,h,s=prep(base,T1,switch,G0);meta=common(q,h,T1,np.nan,switch,G0);rows=[]
    for T2 in range(800,1351,10):
        if T2<T1:rows.append({**meta,"T2_C":T2,**second(q,h,s,T2,switch,T1,G0)})
    return rows

def map_B(base):
    path=SRC/"chen_map_G1_T2_classification.csv";rows=[];done=set()
    if path.exists():
        old=pd.read_csv(path);rows=old.to_dict('records');done=set(zip(old.G0_nm,old.T1_C,old.rho_switch))
    routes=[(base,G0,T1,switch) for G0 in (50,75,100,125,150,200,250,300) for T1 in (1300,1350,1400,1450,1500) for switch in (.80,.84,.88,.90,.92) if (G0,T1,switch) not in done]
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures=[pool.submit(_route_B,x) for x in routes]
        for count,f in enumerate(as_completed(futures),1):
            rows.extend(f.result())
            if count%8==0:write(path,rows)
    write(path,rows);return pd.DataFrame(rows)

def _route_C(args):
    base,switch=args;q,h,s=prep(base,1400,float(switch),100);return [{**common(q,h,1400,T2,float(switch),100),**second(q,h,s,T2,float(switch),1400,100)} for T2 in range(800,1351,10)]
def map_C(base):
    path=SRC/"chen_map_switch_density_T2_classification.csv";rows=[]
    with ProcessPoolExecutor(max_workers=8) as pool:
        for result in pool.map(_route_C,[(base,x) for x in np.round(np.arange(.75,.941,.01),2)]):rows.extend(result);write(path,rows)
    return pd.DataFrame(rows)

def resample(g,n=2000):
    g=g.sort_values("physical_time_s").drop_duplicates("physical_time_s")
    if len(g)<2:return g
    grid=np.unique(np.r_[np.linspace(g.physical_time_s.min(),g.physical_time_s.max(),n),g.physical_time_s,
                         g[g.record_kind.eq('density_crossing')].physical_time_s])
    out={}
    for c in g.columns:
        if pd.api.types.is_numeric_dtype(g[c]) and not pd.api.types.is_bool_dtype(g[c]):out[c]=np.interp(grid,g.physical_time_s,g[c])
        else:out[c]=[g.iloc[np.clip(np.searchsorted(g.physical_time_s,x),0,len(g)-1)][c] for x in grid]
    out=pd.DataFrame(out);out["physical_time_s"]=grid;out["physical_time_h"]=grid/3600
    return out

def local_dense(base):
    p=dict(base);prep_h,prepared=audit.prepare_state(p,900)
    high,_=audit.simulate_detailed(p,T_C=1400,dt_s=900,path_label="highT_reference",stop_density=TARGET,stage="single_step")
    paths={"highT_reference":high}
    for T2,label in ((900,"lower_failure"),(1100,"success"),(1220,"upper_failure")):
        h,_=audit.run_two_step(p,900,T2,prepared,prep_h);h["path_label"]=label;paths[label]=h
    low,_=audit.simulate_detailed(p,T_C=1100,dt_s=900,path_label="lowT_isothermal",stop_density=TARGET,stage="single_step");paths["lowT_isothermal"]=low
    dense=pd.concat([resample(audit.augment_density_landmarks(h)) for h in paths.values()],ignore_index=True)
    ratio=audit.matched_density(high,paths["success"])
    fine=pd.read_csv("results/audit_candidate_693168_closed_accommodation/candidate_693168_T2_classification_points_fine.csv")
    scan=pd.read_csv("results/audit_candidate_693168_closed_accommodation/dense_candidate_693168_T2_scan_histories.csv")
    write(SRC/"dense_time_histories.csv",dense);write(SRC/"dense_matched_density_curves.csv",ratio);write(SRC/"dense_T2_scan_histories.csv",scan)
    write(DENSE/"local_region_presentation_histories.csv",dense)
    diag=fine.rename(columns={"rho2":"final_density","G2_nm":"final_grain_size_nm","closed_shrinkage_contribution":"closed_shrinkage_integral","open_shrinkage_contribution":"open_shrinkage_integral"})
    write(SRC/"T2_diagnostic_curves.csv",diag)
    return dense,ratio,diag

def archive_csv(name):
    with zipfile.ZipFile(ARCHIVE) as z:
        return pd.read_csv(z.open("nucleation_limited_fast_firing_chen_production/"+name))

def fast_dense():
    selected=archive_csv("selected_nucleation_material_sets.csv")
    row=selected[selected.material_id=="E0021"].iloc[0].to_dict();mat=nuc.material(row)
    rows=[]
    for rate in (1,20,50,100):
        h=fast_model.run(mat,fast_model.TopologyGrowthClosure(),protocols.FastSchedule(rate,row["peak_T_C"],row["hold_h"]))
        t=np.asarray(h["t"]);grid=np.unique(np.r_[t,np.linspace(t.min(),t.max(),2000)])
        for x in grid:
            z=dict(material_id="E0021",heating_rate_C_min=rate,path_label=f"{rate} C/min",physical_time_s=x,physical_time_h=x/3600)
            for key in ("T_C","rho","G","tau_nuc","tau_exchange","tau_transport","activity","rho_dot","G_dot"):
                z[key]=float(np.interp(x,t,h[key]))
            z["G_nm"]=z.pop("G")*1e9;z["G_dot_nm_s"]=z.pop("G_dot")*1e9;rows.append(z)
    frame=pd.DataFrame(rows);write(SRC/"dense_fast_firing_histories.csv",frame);write(DENSE/"fast_firing_presentation_histories.csv",frame)
    curves=archive_csv("fast_firing_ablation_ratio_curves.csv");curves=curves[(curves.material_id=="E0021")&(curves.ablation_mode=="full_material_model")];write(SRC/"fast_firing_ratio_curves.csv",curves)
    # Deterministic display map: peak-temperature/hold variants, fixed material.
    m=[]
    for peak in (1450,1500,1550,1600):
      for hold in (4,8,12,20):
       ref=fast_model.run(mat,fast_model.TopologyGrowthClosure(),protocols.FastSchedule(1,peak,hold))
       for rate in (20,50,100):
        fast=fast_model.run(mat,fast_model.TopologyGrowthClosure(),protocols.FastSchedule(rate,peak,hold));
        lo=max(ref["rho"].min(),fast["rho"].min(),.75);hi=min(ref["rho"].max(),fast["rho"].max(),.92)
        attained=hi>=lo+.03
        if attained:
            rho=np.linspace(lo,hi,200);ratio=np.interp(rho,ref["rho"],ref["G"])/np.interp(rho,fast["rho"],fast["G"]);value=float(np.median(ratio))
        else:value=np.nan
        m.append(dict(heating_rate_C_min=rate,peak_T_C=peak,hold_time_h=hold,matched_density_ratio=value,comparison_attained=attained))
    write(SRC/"fast_firing_heating_rate_map.csv",m);return frame,curves,pd.DataFrame(m)

def main():
    start=time.time();SRC.mkdir(parents=True,exist_ok=True);DENSE.mkdir(parents=True,exist_ok=True)
    base,_=audit.candidate_parameters();local_dense(base);fast_dense()
    A=map_A(base);B=map_B(base);C=map_C(base)
    diagnostic=A[A.T1_C==1400].copy().rename(columns={
        "rho_final":"final_density","G_final_nm":"final_grain_size_nm",
        "closed_shrinkage_contribution":"closed_shrinkage_integral",
        "open_shrinkage_contribution":"open_shrinkage_integral",
    })
    write(SRC/"T2_diagnostic_curves.csv",diagnostic)
    dictionary="""# Chen map plot data dictionary\n\n- `UNATTAINABLE_FIRST_STEP`: requested switch state was not reached within the common 500 h budget.\n- `TARGET_REACHED_DURING_FIRST_STEP`: density 0.98 was reached before a valid second step.\n- `DENSIFICATION_EXHAUSTION_FAILURE`: target missed and second-step growth <=20%.\n- `SUCCESS`: target reached and second-step growth <=20%.\n- `GRAIN_GROWTH_FAILURE`: target reached but growth >20%.\n- `MIXED_FAILURE`: target missed and growth >20%.\n- `NUMERICAL_CENSOR`: nonfinite integration output.\n\nAll tiles retain failures; no schedule-specific parameters or filtering are used.\n"""
    (SRC/"chen_map_plot_data_dictionary.md").write_text(dictionary)
    state=dict(status="complete",candidate_id=693168,physics_changed=False,runtime_s=time.time()-start,T1_T2_points=len(A),G1_T2_points=len(B),switch_T2_points=len(C),total_map_points=len(A)+len(B)+len(C))
    (OUT/"map_run_state.json").write_text(json.dumps(state,indent=2)+"\n");print(json.dumps(state,indent=2))

if __name__=="__main__":main()
