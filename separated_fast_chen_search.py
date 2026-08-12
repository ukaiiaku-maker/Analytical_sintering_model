#!/usr/bin/env python3
"""Staged fast-material, Chen-topology, and coexistence discovery."""
from dataclasses import asdict,replace
from pathlib import Path
import argparse,csv,inspect
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import separated_fast_chen_model as model
import observable_trajectory_effect_audit as effect
import production_mechanism_assessment as protocols
import plot_style as ps

OUT=Path("results/separate_fast_firing_and_chen_mechanisms")

def write(path,rows,fields=None):
    rows=list(rows);fields=fields or list(dict.fromkeys(k for r in rows for k in r)) or ["status"]
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def materials(n=256):
    qgb=(350,425,500,575,650);qs=(250,350,450,550,650);qn=(300,400,500,600,700);vs=(3e-29,8e-29,2e-28);sc=(1,3,10,30);pr=(2e-5,2e-4,2e-3);parts=("smoothing","GB_to_TJ","isolation","balanced");ze=(.25,.5,1,2);out=[]
    for i in range(n):out.append((f"M{i:03d}",model.MaterialKinetics(Q_GB_diffusion=qgb[i%5]*1e3,Q_surface_diffusion=qs[(i//5)%5]*1e3,Q_disconnection_nucleation=qn[(i//25)%5]*1e3,v_star=vs[(i//7)%3],stress_concentration=sc[(i//11)%4],PR_prefactor=pr[(i//13)%3],PR_partition=parts[(i//17)%4],zeta_eta_ratio=ze[(i//19)%4],rho0=(.65,.70,.75)[(i//23)%3],G0=(50,75,100,150,225,300)[(i//29)%6]*1e-9)))
    return out

def topology_sets():
    out=[]
    for i in range(24):out.append((f"T{i:02d}",model.TopologyGrowthClosure(mode="persistent_tj_multihit_q0" if i%2==0 else "persistent_tj_multihit_q1",TJ_drag_strength=(2,8,24)[i%3],pore_drag_strength=(0,3,10)[(i//3)%3],XJ_capacity=(.25,.5,.75)[(i//5)%3],lambda_ref=(1,2,4)[(i//7)%3],K_ref=(1,2,4)[(i//11)%3],q_TJ=i%2,pore_relax_fraction=(0,.5,1)[(i//2)%3],pore_drag_fraction=(0,.5,1)[(i//4)%3])))
    return out

def curve(a,b):
    if a["numerical_censored"] or b["numerical_censored"] or not len(a["rho"]) or not len(b["rho"]):
        return pd.DataFrame(columns=("rho","reference_G_nm","comparison_G_nm","ratio"))
    return effect.matched_curve(pd.DataFrame({"rho":a["rho"],"G_nm":a["G"]*1e9}),pd.DataFrame({"rho":b["rho"],"G_nm":b["G"]*1e9}),.75,.92,.001)

def fast_screen(n):
    summary=[];curves=[];diag=[]
    schedules=((1350,2),(1450,8),(1550,20))
    for j,(mid,p) in enumerate(materials(n),1):
        best=None
        for peak,hold in schedules:
            ref=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(1,peak,hold))
            for rate in (20,50,100):
                fast=model.run(p,model.TopologyGrowthClosure(),protocols.FastSchedule(rate,peak,hold));c=curve(ref,fast);span=effect.longest_span(c.rho.to_numpy(),c.ratio.to_numpy(),1.5) if len(c) else 0;attained=len(c)>1 and c.rho.min()<=.751 and c.rho.max()>=.919;meaningful=attained and span>=.03
                censored=bool(ref["numerical_censored"] or fast["numerical_censored"])
                row=dict(material_id=mid,peak_T_C=peak,hold_h=hold,fast_rate_C_min=rate,both_paths_attained=attained,numerical_censored=censored,max_ratio=float(c.ratio.max()) if len(c) else np.nan,span_ge_1p2=effect.longest_span(c.rho.to_numpy(),c.ratio.to_numpy(),1.2) if len(c) else 0,span_ge_1p5=span,meaningful=meaningful,rejection_reason="" if meaningful else ("numerical_censor" if censored else ("unattained_interval" if not attained else "ratio_or_span_below_threshold")),**asdict(p));summary.append(row)
                for x in c.to_dict("records"):curves.append({"material_id":mid,"peak_T_C":peak,"hold_h":hold,"fast_rate_C_min":rate,**x})
                for path,h in (("reference",ref),("fast",fast)):
                    if best is None or row["max_ratio"]>best[0]:best=(row["max_ratio"],path,h,row)
        if best:
            _,path,h,row=best;stride=max(1,len(h["rho"])//100)
            for i in range(0,len(h["rho"]),stride):diag.append(dict(material_id=mid,path=path,rho=h["rho"][i],T_C=h["T_C"][i],activity=h["activity"][i],PR_exposure=h["PR_exposure"][i],pore_D90_nm=h["pore_D90"][i]*1e9,connected_fine=h["connected_fine"][i]))
        if j%25==0:print("materials",j,"/",n,flush=True);write(OUT/"fast_firing_material_parameter_screen.csv",summary)
    return summary,curves,diag

def chen_screen(selected):
    rows=[];bounds=[];diag=[]
    for mid,p in selected:
        # Fixed material kinetics; topology alters only Gdot. A compact adaptive
        # T2 grid exposes lower/upper boundaries without retuning material data.
        for tid,tp in topology_sets():
            for G0 in (75,100,150,225,300,450,600):
                for T1 in (1250,1300,1350,1400,1450,1500):
                    for sw in (.75,.80,.84,.88,.90):
                        pm=replace(p,G0=G0*1e-9,rho0=min(p.rho0,sw-.01))
                        base=model.initial_state(pm);remove=max(sw-pm.rho0,0);weights=base["phi"]*(pm.pore_radius0/base["radii"])**2;weights/=weights.sum();base["phi"]=np.maximum(base["phi"]-remove*weights,0);base["phi"]*=max(1-sw,0)/max(base["phi"].sum(),1e-300)
                        prep=model.material_rates(sw,pm.G0,base["phi"],base["radii"],T1,pm);prep_time=min(4*3600,max(sw-pm.rho0,0)/max(prep["rho_dot"],1e-300));G1=pm.G0+prep["growth_base"]*prep_time
                        points=[]
                        for T2 in np.arange(950,min(T1,1450),25):
                            d=model.material_rates(sw,G1,base["phi"],base["radii"],T2,pm);state={"G":G1,"X_J":min(tp.XJ_capacity,prep["PR_propensity"]*prep_time),"connected_coverage":d["connected_fine"]};gf,td=model.topology_growth_factor(state,T2,tp);hold=20*3600;att=sw+d["rho_dot"]*hold>=.90;growth=d["growth_base"]*gf*hold/G1;kind="success" if att and growth<=.10 else ("growth" if growth>.10 else "density")
                            points.append((T2,kind,att,growth));diag.append(dict(material_id=mid,topology_id=tid,G0_nm=G0,G1_nm=G1*1e9,T1_C=T1,rho_switch=sw,T2_C=T2,classification=kind,X_J=state["X_J"],Lambda_over_K=td["Lambda_over_K"],pore_drag=td["pore_drag"],rho_dot=d["rho_dot"],growth_fraction=growth))
                        success=[x[0] for x in points if x[1]=="success"];lower=bool(success and any(x[1]=="density" and x[0]<min(success) for x in points));upper=bool(success and any(x[1]=="growth" and x[0]>max(success) for x in points));complete=bool(success and lower and upper and max(success)<T1)
                        row=dict(material_id=mid,topology_id=tid,G0_nm=G0,G1_nm=G1*1e9,T1_C=T1,rho_switch=sw,complete_practical_window=complete,lower_bracketed=lower,upper_bracketed=upper,T_lower=min(success) if success else np.nan,T_upper=max(success) if success else np.nan,q_TJ=tp.q_TJ,material_frozen=True,rejection_reason="" if complete else "incomplete_boundaries_or_no_success");rows.append(row)
                        if success:bounds.append(row)
    return rows,bounds,diag

def plots(fast,curves,chen,overlap,diag):
    ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);f=pd.DataFrame(fast);c=pd.DataFrame(curves);ch=pd.DataFrame(chen);o=pd.DataFrame(overlap);d=pd.DataFrame(diag)
    for x,y,name in (("Q_surface_diffusion","Q_disconnection_nucleation","Qsurface_Qnuc_map"),("Q_GB_diffusion","Q_surface_diffusion","QGB_Qsurface_map")):
        fig,ax=plt.subplots(figsize=(5,4));sc=ax.scatter(f[x]/1e3,f[y]/1e3,c=f.max_ratio,cmap="viridis",s=9);ax.set(xlabel=x+" [kJ/mol]",ylabel=y+" [kJ/mol]");plt.colorbar(sc,ax=ax,label="Max ratio");ps.finish(fig,fd/name)
    top=f.sort_values("max_ratio").iloc[-1];q=c[(c.material_id==top.material_id)&(c.fast_rate_C_min==top.fast_rate_C_min)&(c.peak_T_C==top.peak_T_C)&(c.hold_h==top.hold_h)];fig,ax=plt.subplots(figsize=(5,3.5));ax.plot(q.rho,q.ratio);[ax.axhline(v,color="#777",ls="--") for v in (1.2,1.5,2)];ax.set(xlabel="Density",ylabel="Grain ratio");ps.clean(ax);ps.finish(fig,fd/"best_fast_ratio")
    fig,ax=plt.subplots(figsize=(5,4));ok=ch[ch.complete_practical_window==True];ax.scatter(ok.G0_nm,ok.T_lower,s=8);ax.set(xlabel="G0 [nm]",ylabel="T2 lower [C]");ps.clean(ax);ps.finish(fig,fd/"Chen_map")
    fig,ax=plt.subplots(figsize=(5,3.5));z=ch[ch.complete_practical_window==True].groupby("G0_nm").apply(lambda x:(x.T_upper-x.T_lower).median(),include_groups=False);ax.plot(z.index,z.values,"o-");ax.set(xlabel="G0 [nm]",ylabel="Median window width [C]");ps.clean(ax);ps.finish(fig,fd/"Chen_window_width")
    fig,ax=plt.subplots(figsize=(5,3.5));ax.imshow([[sum(o.response_class==k)] for k in ("fast_only","Chen_only","both","neither")],aspect="auto");ax.set_yticks(range(4),("fast only","Chen only","both","neither"));ax.set_xticks([]);ps.finish(fig,fd/"overlap_map")

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--materials",type=int,default=256);args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    fast,curves,pdiag=fast_screen(args.materials);success_ids={r["material_id"] for r in fast if r["meaningful"]};rank=pd.DataFrame(fast).groupby("material_id").max_ratio.max().nlargest(4).index.tolist();selected_ids=list(success_ids) or rank;mdict=dict(materials(args.materials));selected=[(x,mdict[x]) for x in selected_ids];chen,bounds,tdiag=chen_screen(selected);chen_ids={r["material_id"] for r in chen if r["complete_practical_window"]};over=[]
    for mid in selected_ids:
        f=mid in success_ids;c=mid in chen_ids;klass="both" if f and c else ("fast_only" if f else ("Chen_only" if c else "neither"));over.append(dict(material_id=mid,material_parameters=repr(mdict[mid]),fast_success=f,Chen_success=c,response_class=klass,material_unchanged=True))
    write(OUT/"fast_firing_material_parameter_screen.csv",fast);write(OUT/"fast_firing_successful_material_sets.csv",[r for r in fast if r["meaningful"]],list(fast[0]));write(OUT/"fast_firing_rejected_material_sets.csv",[r for r in fast if not r["meaningful"]]);write(OUT/"fast_firing_ratio_curves.csv",curves);write(OUT/"fast_firing_pore_memory_diagnostics.csv",pdiag);write(OUT/"chen_topology_parameter_screen.csv",chen);write(OUT/"chen_successful_topology_sets.csv",[r for r in chen if r["complete_practical_window"]],list(chen[0]));write(OUT/"chen_rejected_topology_sets.csv",[r for r in chen if not r["complete_practical_window"]]);write(OUT/"chen_window_boundaries.csv",bounds);write(OUT/"topology_diagnostics.csv",tdiag);write(OUT/"overlap_scorecard.csv",over);write(OUT/"common_parameter_sets.csv",[r for r in over if r["response_class"]=="both"],list(over[0]));
    for klass,file in (("fast_only","fast_only_cases.csv"),("Chen_only","Chen_only_cases.csv"),("both","both_response_cases.csv")):
        write(OUT/file,[r for r in over if r["response_class"]==klass],list(over[0]))
    write(OUT/"conflict_cases.csv",[r for r in over if "conflict" in r["response_class"]],list(over[0]))
    write(OUT/"material_vs_topology_parameter_map.csv",over)
    plots(fast,curves,chen,over,tdiag)
    print("DONE fast",len(success_ids),"Chen materials",len(chen_ids),"both",sum(r["response_class"]=="both" for r in over))
if __name__=="__main__":main()
