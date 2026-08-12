#!/usr/bin/env python3
"""Generate the consolidated, evidence-backed reduced-model scorecard."""
from pathlib import Path
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plot_style as ps

OUT=Path("results/reduced_model_falsification_summary")

ROWS=[
dict(mechanism="Aggregate fixed model / PR #2 control",branch="codex/topology-constrained-mechanisms",commit="PR #2",physics="Topology-constrained mean-field fluxes",criterion="HR and Chen signs",result="negative fast-firing control",max_ratio=1.00,span_ge_1p2=0.,span_ge_1p5=0.,density_range="0.75-0.92",rho_ref_max=.92,rho_fast_max=.92,chen_windows=False,trajectory_meaningful=False,retention="baseline negative control",q_visible=False,tj_separated=False),
dict(mechanism="Pore-placement topology",branch="codex/pore-placement-topology-search",commit="historical",physics="GB/TJ/isolated pore placement",criterion="finite nanoscale Chen window",result="negative",max_ratio=1.00,span_ge_1p2=0.,span_ge_1p5=0.,density_range="to 0.90",rho_ref_max=.90,rho_fast_max=.90,chen_windows=False,trajectory_meaningful=False,retention="negative control",q_visible=False,tj_separated=True),
dict(mechanism="Pore-location action layer",branch="codex/pore-location-action-layer",commit="historical",physics="Local competing location actions",criterion="Chen map",result="negative",max_ratio=1.00,span_ge_1p2=0.,span_ge_1p5=0.,density_range="to 0.90",rho_ref_max=.90,rho_fast_max=.90,chen_windows=False,trajectory_meaningful=False,retention="abandoned mechanism",q_visible=False,tj_separated=True),
dict(mechanism="Persistent junction + TJ multihit",branch="codex/agentic-mechanism-discovery",commit="source-grounded family",physics="Persistent junction drag and multihit completion",criterion="bounded practical Chen windows",result="Chen success",max_ratio=1.00,span_ge_1p2=0.,span_ge_1p5=0.,density_range="to 0.90",rho_ref_max=.90,rho_fast_max=.90,chen_windows=True,trajectory_meaningful=False,retention="retained Chen baseline",q_visible=True,tj_separated=True),
dict(mechanism="Practical preparation-window search",branch="codex/agentic-mechanism-discovery",commit="6da2500",physics="Censor-aware preparation plus adaptive T2",criterion="T2<T1 with bracketed boundaries",result="150-300 nm windows",max_ratio=1.00,span_ge_1p2=0.,span_ge_1p5=0.,density_range="to 0.90",rho_ref_max=.90,rho_fast_max=.90,chen_windows=True,trajectory_meaningful=False,retention="retained processing map",q_visible=True,tj_separated=True),
dict(mechanism="Production PR/de-sintering memory",branch="codex/pr-desintering-fast-firing-memory",commit="3ead9c9/62394cf",physics="Conservative PR pore memory",criterion="former HR_pct>1 internal",result="internal memory; weak observable effect",max_ratio=1.301,span_ge_1p2=.00,span_ge_1p5=0.,density_range="0.85-0.92",rho_ref_max=.92,rho_fast_max=.92,chen_windows=True,trajectory_meaningful=False,retention="mechanistic baseline, not paper success",q_visible=True,tj_separated=True),
dict(mechanism="Observable trajectory audit",branch="codex/pr-desintering-fast-firing-memory",commit="4631112",physics="Matched-density effect-size analysis",criterion="ratio>=1.5 over Drho>=0.03",result="falsified observable claim",max_ratio=1.301,span_ge_1p2=.00,span_ge_1p5=0.,density_range="0.85-0.92",rho_ref_max=.92,rho_fast_max=.92,chen_windows=True,trajectory_meaningful=False,retention="controlling negative result",q_visible=True,tj_separated=True),
dict(mechanism="Heterogeneous initial state + stress",branch="codex/heterogeneous-initial-state-residual-stress",commit="40b7bb9",physics="Weighted cohorts and residual stress",criterion="observable finite-span ratio",result="transient 1.662 over Drho .008",max_ratio=1.662,span_ge_1p2=.008,span_ge_1p5=.008,density_range="0.85-0.92",rho_ref_max=.859,rho_fast_max=.880,chen_windows=False,trajectory_meaningful=False,retention="negative control",q_visible=True,tj_separated=True),
dict(mechanism="Persistent defect/topology-stress",branch="codex/persistent-defect-topology-stress-memory",commit="4f80b19",physics="Slowly relaxing defect and stored-work state",criterion="observable finite-span ratio",result="ratio 3.73; interval unattained",max_ratio=3.730,span_ge_1p2=0.,span_ge_1p5=0.,density_range="0.85-0.92",rho_ref_max=.85,rho_fast_max=.88,chen_windows=False,trajectory_meaningful=False,retention="negative control",q_visible=True,tj_separated=True),
dict(mechanism="Local connected-sink defect memory",branch="codex/local-connected-sink-defect-memory",commit="69d61e1",physics="Independent matrix/defect regions",criterion="observable finite-span ratio",result="ratio 1.758 only censored",max_ratio=1.758,span_ge_1p2=0.,span_ge_1p5=0.,density_range="0.85-0.92",rho_ref_max=.92,rho_fast_max=.92,chen_windows=False,trajectory_meaningful=False,retention="negative control",q_visible=True,tj_separated=True),
dict(mechanism="Late-stage closed-pore audit",branch="codex/late-stage-closed-pore-trajectory",commit="8ce78ef",physics="Closed store, vacancy transport, gas and detachment",criterion="observable ratio through 0.99",result="five fast paths to .95; no .98",max_ratio=2.024,span_ge_1p2=0.,span_ge_1p5=0.,density_range="0.85-0.99 targets",rho_ref_max=.972,rho_fast_max=.969,chen_windows=True,trajectory_meaningful=False,retention="late-stage negative control",q_visible=True,tj_separated=True),
]

def write(path,rows):
    fields=list(rows[0]);
    with path.open("w",newline="") as f:w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader();w.writerows(rows)

def failure_rows():
    out=[]
    for r in ROWS:
        name=r["mechanism"];out.append(dict(mechanism=name,no_Chen_window=not r["chen_windows"],Chen_but_no_fast=bool(r["chen_windows"] and not r["trajectory_meaningful"]),internal_pore_memory_weak_Grho=name in ("Production PR/de-sintering memory","Observable trajectory audit","Heterogeneous initial state + stress","Late-stage closed-pore audit"),high_ratio_short_span=name=="Heterogeneous initial state + stress",high_ratio_unattained_censored=name in ("Persistent defect/topology-stress","Local connected-sink defect memory","Late-stage closed-pore audit"),density_attainment_destroyed=name in ("Persistent defect/topology-stress","Local connected-sink defect memory"),high_density_unsupported=name!="Late-stage closed-pore audit",Chen_boundary_lost=not r["chen_windows"],parameter_nonidentifiability=name in ("Aggregate fixed model / PR #2 control","Production PR/de-sintering memory"),mean_field_cohort_averaging_suspected=name in ("Heterogeneous initial state + stress","Persistent defect/topology-stress","Local connected-sink defect memory")))
    return out

def figures():
    ps.apply_style();fd=OUT/"figures";fd.mkdir(exist_ok=True);names=[r["mechanism"] for r in ROWS];x=np.arange(len(ROWS))
    fig,ax=plt.subplots(figsize=(7.2,4));ax.plot(x,np.arange(len(x)),"o-");ax.set_xticks(x,[str(i+1) for i in x]);ax.set_yticks(np.arange(len(x)),names,fontsize=6);ax.set(xlabel="Mechanism ladder step",ylabel="");ps.clean(ax);ps.finish(fig,fd/"mechanism_ladder")
    vals=np.array([[r["chen_windows"],r["trajectory_meaningful"]] for r in ROWS]);fig,ax=plt.subplots(figsize=(6.5,4.2));ax.imshow(vals,aspect="auto",cmap=plt.matplotlib.colors.ListedColormap(["#D95F5F","#2A9D8F"]),vmin=0,vmax=1);ax.set_xticks([0,1],["Chen window","Meaningful G(rho)"]);ax.set_yticks(range(len(names)),names,fontsize=6);ps.finish(fig,fd/"Chen_vs_trajectory_matrix")
    fig,ax=plt.subplots(figsize=(6,4));ax.scatter([r["max_ratio"] for r in ROWS],[r["span_ge_1p5"] for r in ROWS]);[ax.text(r["max_ratio"],r["span_ge_1p5"],str(i+1),fontsize=6) for i,r in enumerate(ROWS)];ax.axvline(1.5,color="#555",ls="--");ax.axhline(.03,color="#555",ls="--");ax.set(xlabel="Maximum attempted mean-grain ratio",ylabel="Attained span with ratio >=1.5");ps.clean(ax);ps.finish(fig,fd/"ratio_vs_density_span")
    h=pd.read_csv("results/observable_trajectory_effect_audit/fast_firing_ratio_curves.csv");top=h.groupby("comparison_id").ratio.max().idxmax();q=h[h.comparison_id==top];fig,ax=plt.subplots(figsize=(5,3.5));ax.plot(q.rho,q.G_reference_nm,label="reference");ax.plot(q.rho,q.G_comparison_nm,label="fast");ax.set(xlabel="Density",ylabel="Mean grain size [nm]");ax.legend();ps.clean(ax);ps.finish(fig,fd/"representative_G_rho")
    mh=pd.read_csv("results/rejected_case_failure_decomposition/path_histories.csv");fig,axs=plt.subplots(1,2,figsize=(7.2,3),constrained_layout=True)
    for path,z in mh.groupby("path"):axs[0].plot(z.rho,z.cumulative_PR_desintering_work,label=path);axs[1].plot(z.rho,z.G_mean*1e9,label=path)
    axs[0].set(xlabel="Density",ylabel="PR work");axs[1].set(xlabel="Density",ylabel="Mean grain size [nm]");axs[0].legend();ps.panel_labels(axs);ps.finish(fig,fd/"internal_memory_vs_observable")
    fig,ax=plt.subplots(figsize=(6,3.5));targets=[.95,.98,.99];support=[[r["rho_ref_max"]>=t and r["rho_fast_max"]>=t for t in targets] for r in ROWS];ax.imshow(support,aspect="auto",cmap=plt.matplotlib.colors.ListedColormap(["#DDDDDD","#2A9D8F"]));ax.set_xticks(range(3),targets);ax.set_yticks(range(len(names)),names,fontsize=6);ax.set(xlabel="Jointly supported density");ps.finish(fig,fd/"high_density_support")
    fig,ax=plt.subplots(figsize=(7.2,3));ax.axis("off");nodes=[("Experimental\neffect confirmed?",.1),("No: retain\nfalsification",.35),("Yes: minimal spatial\nnetwork/ensemble",.62),("Calibrate against\n3D topology + G(rho)",.88)];
    for text,xp in nodes:ax.text(xp,.55,text,ha="center",va="center",bbox=dict(boxstyle="round",fc="white"));
    for i in range(3):ax.annotate("",xy=(nodes[i+1][1]-.09,.55),xytext=(nodes[i][1]+.09,.55),arrowprops=dict(arrowstyle="->"));ps.finish(fig,fd/"next_model_decision_tree")

def main():
    OUT.mkdir(parents=True,exist_ok=True);write(OUT/"mechanism_scorecard.csv",ROWS);write(OUT/"failure_mode_matrix.csv",failure_rows());figures();print("wrote",len(ROWS),"mechanism outcomes")
if __name__=="__main__":main()
