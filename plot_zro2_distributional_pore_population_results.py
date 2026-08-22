#!/usr/bin/env python3
"""Publication-style figures for the bounded distributional model-form test."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch

import build_zro2_distributional_pore_population_model as model

OUT=model.OUT;FIG=model.FIG;SRC=OUT/"source_tables"
COL={"open_connected":"#277da1","precursor":"#f9c74f","isolated":"#f9844a","closed":"#6a4c93"}

def save(fig,name,source):
    fig.tight_layout();fig.savefig(FIG/f"{name}.pdf",bbox_inches="tight");fig.savefig(FIG/f"{name}.png",dpi=220,bbox_inches="tight");plt.close(fig)
    return {"figure":name,"pdf":str(FIG/f"{name}.pdf"),"png":str(FIG/f"{name}.png"),"source_table":str(source),"placeholder":False}

def main():
    FIG.mkdir(parents=True,exist_ok=True);SRC.mkdir(exist_ok=True);inventory=[]
    # 1 architecture schematic
    rows=[];fig,ax=plt.subplots(figsize=(11,5));ax.axis("off")
    xs=[.08,.34,.60,.84]
    for x,t in zip(xs,model.TOPO):
        box=FancyBboxPatch((x-.09,.58),.18,.18,boxstyle="round,pad=.02",fc=COL[t],alpha=.28,ec=COL[t],lw=2);ax.add_patch(box);ax.text(x,.67,t.replace('_','\n'),ha='center',va='center',weight='bold');rows.append({"node":t,"changes_density":False})
    for a,b,label in [(0,1,"PR/pinch"),(1,2,"isolation"),(2,3,"closure")]:
        ax.annotate("",xy=(xs[b]-.1,.67),xytext=(xs[a]+.1,.67),arrowprops=dict(arrowstyle="->",lw=1.6));ax.text((xs[a]+xs[b])/2,.75,label,ha='center',fontsize=8)
    ax.annotate("precursor → closed",xy=(xs[3]-.08,.58),xytext=(xs[1]+.08,.52),ha='center',fontsize=8,arrowprops=dict(arrowstyle="->",lw=1.3,connectionstyle="arc3,rad=.22"))
    ax.text(.22,.28,"open shrinkage",ha='center',weight='bold');ax.text(.72,.28,"closed shrinkage",ha='center',weight='bold')
    ax.annotate("density",xy=(.48,.12),xytext=(.22,.25),arrowprops=dict(arrowstyle="->",lw=2,color='#d62728'));ax.annotate("density",xy=(.52,.12),xytext=(.72,.25),arrowprops=dict(arrowstyle="->",lw=2,color='#d62728'))
    ax.text(.5,.05,"Only named shrinkage changes density",ha='center',color='#d62728',weight='bold');ax.text(.5,.92,"Conservative distribution evolution → topology state → shrinkage + distributional Zener",ha='center',weight='bold')
    src=SRC/"distribution_model_schematic_source.csv";pd.DataFrame(rows).to_csv(src,index=False);inventory.append(save(fig,"distribution_model_schematic",src))

    # 2 initial families
    d=pd.read_csv(OUT/"initial_distribution_families.csv");q=d[(d.D50_input_nm==24.5)&(d.sigma_ln_input==.45)&(((d.family=='bimodal')&(d.tail_weight==.10))|((d.family!='bimodal')&(d.tail_weight==0)))]
    src=SRC/"initial_distribution_families_source.csv";q.to_csv(src,index=False);fig,ax=plt.subplots(figsize=(8,5))
    for fam,g in q.groupby('family'):ax.plot(g.radius_nm,g.phi_open,'o-',label=fam)
    ax.set(xscale='log',xlabel='Pore radius (nm)',ylabel='Open pore-volume fraction per bin');ax.legend();inventory.append(save(fig,"initial_distribution_families",src))

    # 3 evolution examples
    h=pd.read_csv(OUT/"distribution_ablation_histories.csv");q=h[h.ablation.eq('baseline')]
    cols=['time_s','D50_nm','D90_nm','D90_over_D50','sigma_ln_r','large_pore_tail_fraction','connected_fine_pore_fraction'];src=SRC/"distribution_evolution_examples_source.csv";q[cols].to_csv(src,index=False)
    fig,axs=plt.subplots(2,2,figsize=(10,7));x=q.time_s/3600
    axs[0,0].plot(x,q.D50_nm,label='D50');axs[0,0].plot(x,q.D90_nm,label='D90');axs[0,0].legend();axs[0,0].set_ylabel('Diameter (nm)')
    axs[0,1].plot(x,q.D90_over_D50);axs[0,1].set_ylabel('D90/D50');axs[1,0].plot(x,q.sigma_ln_r);axs[1,0].set_ylabel('sigma ln r')
    axs[1,1].plot(x,q.large_pore_tail_fraction,label='large tail');axs[1,1].plot(x,q.connected_fine_pore_fraction,label='connected fine');axs[1,1].legend()
    for ax in axs.flat:ax.set_xlabel('Time (h)')
    inventory.append(save(fig,"distribution_evolution_examples",src))

    # 4 phase diagram from fixed state scan
    d=pd.read_csv(OUT/"distribution_evolution_fixed_state_scan.csv");q=d[(d.D50_nm==24.5)&(d.tail_weight==.10)&(d.lambda_seg_over_r==10)].groupby(['T_C','sigma_ln','activity'],as_index=False)[['regularization_rate','damage_rate']].mean();q['dominant']='regularization';q.loc[q.damage_rate>q.regularization_rate,'dominant']='damage'
    src=SRC/"PR_regularization_damage_phase_diagram_source.csv";q.to_csv(src,index=False);fig,axs=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
    for ax,act in zip(axs,[.01,.5]):
        z=q.iloc[(q.activity-act).abs().argsort()].groupby(['T_C','sigma_ln'],as_index=False).first();piv=z.pivot(index='sigma_ln',columns='T_C',values='damage_rate')-z.pivot(index='sigma_ln',columns='T_C',values='regularization_rate');im=ax.contourf(piv.columns,piv.index,np.sign(piv.values),levels=[-1.5,0,1.5],colors=['#43aa8b','#f94144'],alpha=.75);ax.set(title=f'activity={act}',xlabel='T (°C)')
    axs[0].set_ylabel('sigma ln r');inventory.append(save(fig,"PR_regularization_damage_phase_diagram",src))

    # 5 mean versus distributional Zener
    rows=[]
    for width in (.3,.45,.65,.85,1.1):
      for tail in (0,.05,.10,.25):
       s=model.initial_state(family='bimodal',D50_nm=24.5,sigma_ln=width,tail_weight=tail);phi=s.phi[0];mean=np.sum(phi*s.radii_m)/phi.sum();rows.append({'sigma_ln':width,'tail_weight':tail,'mean_radius_Zener':phi.sum()/mean,'distributional_Zener':np.sum(phi/s.radii_m)})
    q=pd.DataFrame(rows);src=SRC/"distributional_Zener_metric_source.csv";q.to_csv(src,index=False);fig,ax=plt.subplots(figsize=(7,5));sc=ax.scatter(q.mean_radius_Zener,q.distributional_Zener,c=q.sigma_ln,s=50+250*q.tail_weight,cmap='viridis');ax.plot([q.mean_radius_Zener.min(),q.mean_radius_Zener.max()],[q.mean_radius_Zener.min(),q.mean_radius_Zener.max()],'k--');ax.set(xlabel='Mean-radius Zener metric (1/m)',ylabel='Distributional Zener metric (1/m)');fig.colorbar(sc,ax=ax,label='sigma ln r');inventory.append(save(fig,"distributional_Zener_metric",src))

    # 6 synthetic classification summary
    d=pd.read_csv(OUT/"synthetic_distribution_boundary_test.csv");q=d.groupby(['state_id','T2_C','classification'],as_index=False).size();src=SRC/"synthetic_boundary_test_source.csv";q.to_csv(src,index=False)
    order=list(d.state_id.unique());colors={'DENSIFICATION_EXHAUSTION_FAILURE':'#457b9d','SUCCESS':'#2a9d8f','GRAIN_GROWTH_FAILURE':'#e76f51','MIXED_FAILURE':'#6d597a'};fig,ax=plt.subplots(figsize=(11,6))
    for yi,state in enumerate(order):
      z=d[(d.state_id==state)&(d.representation=='discrete_bin')&(d.closed_kernel=='GB_diffusion')&(d.m==3)]
      ax.scatter(z.T2_C,np.full(len(z),yi),c=[colors[x] for x in z.classification],s=42)
    ax.set(yticks=range(len(order)),yticklabels=[x.replace('_',' ') for x in order],xlabel='T2 (°C)',title='Synthetic topology classification (discrete bins, GB diffusion, m=3)')
    ax.legend(handles=[Patch(facecolor=c,label=k.replace('_',' ').title()) for k,c in colors.items()],loc='upper center',bbox_to_anchor=(.5,-.12),ncol=2,frameon=False)
    inventory.append(save(fig,"synthetic_boundary_test",src))

    # Failure-mode map because no finite window was found.
    q=d[d.state_id.eq('narrow_fine_connected')];src=SRC/"distribution_Chen_failure_modes_source.csv";q.to_csv(src,index=False);fig,axs=plt.subplots(1,3,figsize=(12,4),sharey=True)
    for ax,(rep,g) in zip(axs,q.groupby('representation')):
      z=g[(g.closed_kernel=='GB_diffusion')&(g.m==3)];ax.scatter(z.T2_C,z.final_rho,c=[colors[x] for x in z.classification],s=45);ax.axhline(.9,color='k',ls='--');ax.set(title=rep,xlabel='T2 (°C)')
    axs[0].set_ylabel('Final density');inventory.append(save(fig,"distribution_Chen_G1_T2_failure_modes",src))

    # 8 heating-rate response
    d=pd.read_csv(OUT/"distribution_heating_rate_histories.csv");q=d[(d.family=='lognormal')&(d.peak_C==1350)&(d.hold_h==2)&d.rate_C_min.isin([.2,1,5,20,100])];src=SRC/"distribution_heating_rate_G_vs_rho_source.csv";q.to_csv(src,index=False);fig,ax=plt.subplots(figsize=(7,5))
    for rate,g in q.groupby('rate_C_min'):ax.plot(g.rho,g.G_nm,label=f'{rate:g} °C/min')
    ax.set(xlabel='Relative density',ylabel='Grain size (nm)',xlim=(.66,1));ax.legend();inventory.append(save(fig,"distribution_heating_rate_G_vs_rho",src))

    # 9 ablation matrix
    d=pd.read_csv(OUT/"distribution_ablation_matrix.csv");src=SRC/"ablation_matrix_source.csv";d.to_csv(src,index=False);fig,axs=plt.subplots(1,2,figsize=(12,6));y=np.arange(len(d));axs[0].barh(y,d.final_rho);axs[0].axvline(.9,color='k',ls='--');axs[0].set(yticks=y,yticklabels=d.ablation,xlabel='Final density');axs[1].barh(y,d.growth_fraction);axs[1].axvline(.1,color='k',ls='--');axs[1].set(yticks=y,yticklabels=[],xlabel='Growth fraction');inventory.append(save(fig,"ablation_matrix",src))
    pd.DataFrame(inventory).assign(validation=False).to_csv(OUT/"figure_inventory.csv",index=False)
    print({'figures':len(inventory),'directory':str(FIG)})

if __name__=='__main__':main()
