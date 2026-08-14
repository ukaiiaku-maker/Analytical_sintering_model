#!/usr/bin/env python3
"""Experimental-style main and supplemental figures for fixed candidate 693168."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch,Circle,FancyArrowPatch
import publication_plot_style as sty

OUT=Path('results/publication_style_sintering_figures_693168');SRC=OUT/'source_tables'
MAIN=OUT/'main_figures';SUP=OUT/'supplement_figures';CHEN=OUT/'chen_maps';FAST=OUT/'fast_firing';TWO=OUT/'two_step'
INV=[]
def save(fig,path,title,sources,category):
    sty.save(fig,path);INV.append(dict(figure_id=f'F{len(INV)+1:02d}',filename_pdf=str(path.with_suffix('.pdf').relative_to(OUT)),filename_png=str(path.with_suffix('.png').relative_to(OUT)),title=title,category=category,source_tables=sources,candidate_status='conditional Tier B; not calibrated; not validation'))
def finish_axes(ax):
    for a in np.asarray(ax).flat:sty.clean(a)
    sty.letters(ax)
def plot_paths(ax,d,x,y,labels=None):
    for lab,g in d.groupby('path_label',sort=False):
        if labels is not None and lab not in labels:continue
        ax.plot(g[x],g[y],label=lab.replace('_',' '),color=sty.COLORS.get(lab))
    ax.legend(fontsize=7)
def class_legend(ax):ax.legend(handles=[Patch(color=v,label=k.replace('_',' ').title()) for k,v in sty.CLASSES.items()],fontsize=6,ncol=2)
def categorical_grid(ax,d,x,y,title,annotate=None):
    cats=list(sty.CLASSES);code={c:i for i,c in enumerate(cats)}
    p=d.pivot_table(index=y,columns=x,values='classification',aggfunc=lambda s:s.mode().iat[0])
    z=np.vectorize(lambda q:code.get(q,np.nan),otypes=[float])(p.values)
    cmap=ListedColormap([sty.CLASSES[c] for c in cats]);norm=BoundaryNorm(np.arange(-.5,len(cats)+.5),len(cats))
    ax.pcolormesh(p.columns,p.index,z,shading='nearest',cmap=cmap,norm=norm,rasterized=True)
    ax.contour(p.columns,p.index,(z==code['SUCCESS']).astype(float),levels=[.5],colors='black',linewidths=1.2)
    ax.set(title=title);class_legend(ax)
    if annotate:ax.scatter(*annotate,marker='*',s=100,color='black',zorder=5,label='Candidate 693168')
    return p
def binned_G_map(ax,d):
    q=d.copy();edges=np.geomspace(max(q.G1_nm.quantile(.001),1),q.G1_nm.quantile(.999),36);q['G1_bin']=pd.cut(q.G1_nm,edges,labels=np.sqrt(edges[:-1]*edges[1:]),include_lowest=True)
    q=q.dropna(subset=['G1_bin']);q['G1_bin']=q.G1_bin.astype(float)
    categorical_grid(ax,q,'G1_bin','T2_C',r'Prepared grain size–$T_2$ map',annotate=(117.1,1100));ax.set(xscale='log',xlabel=r'Prepared grain size, $G_1$ [nm]',ylabel=r'Second-step temperature, $T_2$ [°C]')

def fast_figures(fast,ratio,hmap):
    q=fast
    fig,ax=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for rate,g in q.groupby('heating_rate_C_min'):
        lab=f'{rate:g} °C/min';c=sty.RATES[rate];ax[0,0].plot(g.physical_time_h,g.T_C,label=lab,color=c);ax[0,1].plot(g.physical_time_h,g.rho,label=lab,color=c);ax[1,0].plot(g.rho,g.G_nm,label=lab,color=c)
    ref=q[q.heating_rate_C_min==1]
    for rate in (20,50,100):
        g=q[q.heating_rate_C_min==rate];lo=max(ref.rho.min(),g.rho.min());hi=min(ref.rho.max(),g.rho.max());rho=np.linspace(lo,hi,400);ax[1,1].plot(rho,np.interp(rho,ref.rho,ref.G_nm)/np.interp(rho,g.rho,g.G_nm),label=f'1/{rate} °C min$^{{-1}}$',color=sty.RATES[rate])
    for y in (1.2,1.5,2):ax[1,1].axhline(y,color='.5',ls='--',lw=.8)
    labs=[('Time, $t$ [h]','Temperature, $T$ [°C]'),('Time, $t$ [h]',r'Relative density, $\rho$'),(r'Relative density, $\rho$','Grain size, $G$ [nm]'),(r'Relative density, $\rho$',r'$G_{1 °C/min}/G_{fast}$')]
    for a,(x,y) in zip(ax.flat,labs):a.set(xlabel=x,ylabel=y);a.legend()
    finish_axes(ax);fig.suptitle('E0021 nucleation-limited fast-firing envelope (not re-optimized)')
    save(fig,MAIN/'Figure_1_fast_firing_heating_rate_effect','Fast-firing heating-rate effect','dense_fast_firing_histories.csv; fast_firing_ratio_curves.csv','main')
    # Required clean derivatives.
    fig,ax=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for rate,g in q.groupby('heating_rate_C_min'):
        for a,y in zip(ax,('T_C','rho','G_nm')):a.plot(g.physical_time_h,g[y],label=f'{rate:g} °C/min',color=sty.RATES[rate])
    for a,l in zip(ax,('Temperature, $T$ [°C]',r'Relative density, $\rho$','Grain size, $G$ [nm]')):a.set(xlabel='Time, $t$ [h]',ylabel=l);a.legend()
    finish_axes(ax);save(fig,FAST/'fast_firing_T_rho_G_time','Fast-firing time histories','dense_fast_firing_histories.csv','fast_firing')
    fig,ax=plt.subplots(figsize=(6,4));
    for rate,g in q.groupby('heating_rate_C_min'):ax.plot(g.rho,g.G_nm,label=f'{rate:g} °C/min',color=sty.RATES[rate])
    ax.set(xlabel=r'Relative density, $\rho$',ylabel='Grain size, $G$ [nm]');ax.legend();sty.clean(ax);save(fig,FAST/'fast_firing_G_vs_rho','Fast-firing grain-density trajectories','dense_fast_firing_histories.csv','fast_firing')
    fig,ax=plt.subplots(figsize=(6,4));
    for rate in (20,50,100):
        g=q[q.heating_rate_C_min==rate];lo=max(ref.rho.min(),g.rho.min());hi=min(ref.rho.max(),g.rho.max());rho=np.linspace(lo,hi,400);ax.plot(rho,np.interp(rho,ref.rho,ref.G_nm)/np.interp(rho,g.rho,g.G_nm),label=f'1/{rate} °C min$^{{-1}}$',color=sty.RATES[rate])
    for y in (1.2,1.5,2):ax.axhline(y,color='.5',ls='--',lw=.8)
    ax.set(xlabel=r'Relative density, $\rho$',ylabel=r'$G_{1 °C/min}/G_{fast}$');ax.legend();sty.clean(ax);save(fig,FAST/'fast_firing_ratio_vs_rho','Fast-firing matched-density ratio','dense_fast_firing_histories.csv','fast_firing')
    fig,axs=plt.subplots(2,2,figsize=(10,8),constrained_layout=True)
    for a,(hold,g) in zip(axs.flat,hmap.groupby('hold_time_h')):
        p=g.pivot(index='peak_T_C',columns='heating_rate_C_min',values='matched_density_ratio');im=a.imshow(p.values,origin='lower',aspect='auto',cmap='viridis',vmin=1,vmax=np.nanmax(hmap.matched_density_ratio));a.set(xticks=range(len(p.columns)),xticklabels=p.columns,yticks=range(len(p.index)),yticklabels=p.index,xlabel='Heating rate [°C/min]',ylabel='Peak temperature [°C]',title=f'{hold:g} h hold');fig.colorbar(im,ax=a,label='Median matched-density grain ratio')
    sty.letters(axs);save(fig,FAST/'fast_firing_heating_rate_map','Fast-firing heating-rate map','fast_firing_heating_rate_map.csv','fast_firing')

def two_step_figures(d,ratio,diag):
    labels=['highT_reference','lower_failure','success','upper_failure']
    fig,ax=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    for a,y in zip(ax.flat[:4],('T_C','rho','G_mean_nm','G_mean_nm')):plot_paths(a,d,'physical_time_h' if a is not ax.flat[3] else 'rho',y,labels)
    ax[0,0].set(xlabel='Time, $t$ [h]',ylabel='Temperature, $T$ [°C]');ax[0,1].set(xlabel='Time, $t$ [h]',ylabel=r'Relative density, $\rho$');ax[0,2].set(xlabel='Time, $t$ [h]',ylabel='Grain size, $G$ [nm]');ax[1,0].set(xlabel=r'Relative density, $\rho$',ylabel='Grain size, $G$ [nm]')
    ax[1,1].plot(ratio.rho,ratio.reduction_TS,color=sty.COLORS['success']);ax[1,1].set(xlabel=r'Relative density, $\rho$',ylabel=r'Grain-size reduction, $1-G_{two-step}/G_{highT}$')
    ax[1,2].plot(diag.T2_C,diag.final_density,label='Final density');a2=ax[1,2].twinx();a2.plot(diag.T2_C,diag.growth_fraction,color='#D55E00',label='Growth fraction');ax[1,2].set(xlabel=r'Second-step temperature, $T_2$ [°C]',ylabel=r'Final density, $\rho$');a2.set_ylabel('Growth fraction')
    finish_axes(ax);fig.suptitle('Candidate 693168: conditional Tier B prototype; not calibrated')
    save(fig,MAIN/'Figure_2_candidate_693168_two_step_trajectory','Two-step trajectory','dense_time_histories.csv; dense_matched_density_curves.csv; T2_diagnostic_curves.csv','main')
    # Required separate plots.
    for stem,ys in [('two_step_T_rho_G_time',('T_C','rho','G_mean_nm'))]:
        fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
        for a,y,l in zip(axs,ys,('Temperature, $T$ [°C]',r'Relative density, $\rho$','Grain size, $G$ [nm]')):plot_paths(a,d,'physical_time_h',y,labels);a.set(xlabel='Time, $t$ [h]',ylabel=l)
        finish_axes(axs);save(fig,TWO/stem,'Two-step physical-time histories','dense_time_histories.csv','two_step')
    fig,ax=plt.subplots(figsize=(6,4));plot_paths(ax,d,'rho','G_mean_nm',labels);ax.set(xlabel=r'Relative density, $\rho$',ylabel='Grain size, $G$ [nm]');save(fig,TWO/'two_step_G_vs_rho','Two-step G-rho','dense_time_histories.csv','two_step')
    fig,ax=plt.subplots(figsize=(6,4));ax.plot(ratio.rho,ratio.reduction_TS,color='#009E73');ax.set(xlabel=r'Relative density, $\rho$',ylabel=r'$1-G_{two-step}/G_{highT}$');sty.clean(ax);save(fig,TWO/'two_step_reduction_vs_rho','Two-step reduction','dense_matched_density_curves.csv','two_step')
    fig,ax=plt.subplots(1,2,figsize=(10,4),constrained_layout=True);ax[0].plot(diag.T2_C,diag.final_density);ax[1].plot(diag.T2_C,diag.growth_fraction,color='#D55E00');ax[0].set(xlabel=r'$T_2$ [°C]',ylabel=r'Final density, $\rho$');ax[1].set(xlabel=r'$T_2$ [°C]',ylabel='Second-step growth fraction');finish_axes(ax);save(fig,TWO/'two_step_final_density_growth_vs_T2','Final density and growth versus T2','T2_diagnostic_curves.csv','two_step')
    fig,ax=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for a,lab in zip(ax,('lower_failure','success','upper_failure')):g=d[d.path_label==lab];a.plot(g.rho,g.G_mean_nm,color=sty.COLORS[lab]);a.set(xlabel=r'Relative density, $\rho$',ylabel='Grain size, $G$ [nm]',title=lab.replace('_',' ').title())
    finish_axes(ax);save(fig,TWO/'two_step_success_failure_triplet','Two-step outcome triplet','dense_time_histories.csv','two_step')

def chen_figures(A,B,C,diag):
    fig,ax=plt.subplots(figsize=(7,6));categorical_grid(ax,A,'T1_C','T2_C',r'$T_1$–$T_2$ fixed-candidate map',annotate=(1400,1100));ax.plot([1200,1550],[1200,1550],color='black',ls='--');ax.set(xlabel=r'First-step temperature, $T_1$ [°C]',ylabel=r'Second-step temperature, $T_2$ [°C]');save(fig,CHEN/'chen_map_T1_T2_experimental_style','T1-T2 Chen map','chen_map_T1_T2_classification.csv','chen_map')
    fig,ax=plt.subplots(figsize=(7,6));binned_G_map(ax,B);save(fig,CHEN/'chen_map_G1_T2_experimental_style','G1-T2 Chen map','chen_map_G1_T2_classification.csv','chen_map')
    fig,ax=plt.subplots(figsize=(7,6));categorical_grid(ax,C,'rho_switch','T2_C',r'Switch-density–$T_2$ map',annotate=(.88,1100));ax.set(xlabel=r'Switch density, $\rho_1$',ylabel=r'Second-step temperature, $T_2$ [°C]');save(fig,CHEN/'chen_map_switch_density_T2','Switch-density Chen map','chen_map_switch_density_T2_classification.csv','chen_map')
    fig,ax=plt.subplots(2,2,figsize=(10,8),constrained_layout=True);spec=[('final_density',r'Final density, $\rho$'),('final_grain_size_nm','Final grain size, $G$ [nm]'),('growth_fraction','Second-step growth fraction'),('closed_fraction_final',r'Closed-pore fraction, $f_{closed}$')]
    for a,(k,l) in zip(ax.flat,spec):a.plot(diag.T2_C,diag[k]);a.set(xlabel=r'$T_2$ [°C]',ylabel=l)
    finish_axes(ax);save(fig,CHEN/'T2_diagnostic_curves_experimental_style','T2 diagnostic curves','T2_diagnostic_curves.csv','chen_map')
    fig,ax=plt.subplots(1,3,figsize=(16,5),constrained_layout=True);categorical_grid(ax[0],A,'T1_C','T2_C',r'$T_1$–$T_2$',annotate=(1400,1100));ax[0].plot([1200,1550],[1200,1550],'k--');ax[0].set(xlabel=r'$T_1$ [°C]',ylabel=r'$T_2$ [°C]');binned_G_map(ax[1],B);ax[2].plot(diag.T2_C,diag.final_density,label='Final density');a2=ax[2].twinx();a2.plot(diag.T2_C,diag.growth_fraction,color='#D55E00',label='Growth');ax[2].set(xlabel=r'$T_2$ [°C]',ylabel=r'Final $\rho$');a2.set_ylabel('Growth fraction');sty.letters(ax)
    save(fig,MAIN/'Figure_3_candidate_693168_complete_Chen_map','Complete filled Chen map','three classification tables; T2 diagnostics','main')

def mechanism_figures(d,abl,six,rob):
    labels=['highT_reference','success','lower_failure','upper_failure'];fig,ax=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    for k,lab in zip(('phi_open','phi_connected','phi_iso','phi_closed'),('Open','Connected','Isolated','Closed')):ax[0,0].plot(d[d.path_label=='success'].rho,d[d.path_label=='success'][k],label=lab)
    ax[0,0].legend();ax[0,0].set(xlabel=r'$\rho$',ylabel='Pore-volume fraction')
    plot_paths(ax[0,1],d,'rho','closed_accommodation_factor',labels);ax[0,1].set(xlabel=r'$\rho$',ylabel=r'Closed accommodation, $A_{closed}$')
    for k,lab in (('cumulative_open_pore_removed','Open removal'),('cumulative_closed_pore_removed','Closed removal')):g=d[d.path_label=='success'];ax[0,2].plot(g.rho,g[k],label=lab)
    ax[0,2].legend();ax[0,2].set(xlabel=r'$\rho$',ylabel='Cumulative pore-volume removal')
    plot_paths(ax[1,0],d,'rho','PR_damage_memory',labels);ax[1,0].set(xlabel=r'$\rho$',ylabel='PR-prepared memory')
    a=abl.copy();a.median_reduction=a.median_reduction.fillna(0);ax[1,1].barh(a.ablation.str.replace('_',' '),a.median_reduction,color=np.where(a.complete,'#009E73','#D55E00'));ax[1,1].tick_params(axis='y',labelsize=6);ax[1,1].set(xlabel='Median reduction')
    ax[1,2].axis('off');xs=[.1,.32,.55,.77];names=['PR preparation','Closed transition','Finite accommodation','T₂ outcome'];
    for x,n in zip(xs,names):ax[1,2].add_patch(Circle((x,.5),.075,fill=False,lw=2));ax[1,2].text(x,.72,n,ha='center',fontsize=8)
    for a,b in zip(xs[:-1],xs[1:]):ax[1,2].add_patch(FancyArrowPatch((a+.08,.5),(b-.08,.5),arrowstyle='->',mutation_scale=15));ax[1,2].set_xlim(0,1);ax[1,2].set_ylim(0,1)
    finish_axes(ax[:,:2]);sty.letters(ax);fig.suptitle('Candidate 693168 mechanism interpretation—closed topology is the falsification target')
    save(fig,MAIN/'Figure_4_candidate_693168_mechanism_interpretation','Mechanism interpretation','dense histories; ablation audit','main')
    # supplements S1-S7
    fig,ax=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    for a,y,l,x in zip(ax.flat,('T_C','rho','G_mean_nm','G_mean_nm','rho_dot','G_dot'),('T [°C]',r'$\rho$','G [nm]','G [nm]',r'$d\rho/dt$',r'$dG/dt$'),('physical_time_h','physical_time_h','physical_time_h','rho','physical_time_h','physical_time_h')):plot_paths(a,d,x,y,labels);a.set(xlabel='Time [h]' if x!='rho' else r'$\rho$',ylabel=l)
    finish_axes(ax);save(fig,SUP/'S1_dense_time_histories_all_paths','Dense histories','dense_time_histories.csv','supplement')
    fig,ax=plt.subplots(2,3,figsize=(14,8),constrained_layout=True);g=d[d.path_label=='success'];spec=(('phi_open','Open pore fraction'),('phi_GBseg','GB-segment fraction'),('phi_TJ','TJ fraction'),('phi_iso','Isolated fraction'),('phi_closed','Closed fraction'),('N_closed','Closed-pore number proxy'))
    for a,(k,l) in zip(ax.flat,spec):a.plot(g.rho,g[k]);a.set(xlabel=r'$\rho$',ylabel=l)
    finish_axes(ax);save(fig,SUP/'S2_pore_store_evolution','Pore-store evolution','dense_time_histories.csv','supplement')
    fig,ax=plt.subplots(2,3,figsize=(14,8),constrained_layout=True);spec=(('closed_accommodation_capacity','Capacity'),('closed_accommodation_used','Used capacity'),('P_comp_closed','Completion proxy'),('closed_accommodation_factor','Accommodation factor'),('closed_shrinkage_flux','Closed shrinkage flux'),('closed_pore_contribution_to_rho_dot',r'Closed $d\rho/dt$'))
    for a,(k,l) in zip(ax.flat,spec):plot_paths(a,d,'rho',k,labels);a.set(xlabel=r'$\rho$',ylabel=l)
    finish_axes(ax);save(fig,SUP/'S3_closed_accommodation_diagnostics','Closed accommodation diagnostics','dense_time_histories.csv','supplement')
    fig,ax=plt.subplots(2,2,figsize=(10,8),constrained_layout=True);spec=(('PR_redistribution_rate','PR redistribution rate'),('rho_dot','Densification rate'),('cumulative_PR_redistributed_volume','Cumulative PR redistribution'),('cumulative_closed_pore_removed','Cumulative densifying removal'))
    for a,(k,l) in zip(ax.flat,spec):plot_paths(a,d,'rho',k,labels);a.set(xlabel=r'$\rho$',ylabel=l)
    finish_axes(ax);save(fig,SUP/'S4_PR_damage_energy_partition','PR/densification volume-equivalent partition','dense_time_histories.csv','supplement')
    fig,ax=plt.subplots(1,3,figsize=(14,5),constrained_layout=True)
    for a,k,l in zip(ax,('median_reduction','window_width_C','high_density_attainment'),('Median reduction','Window width [°C]','High-density attainment')):a.barh(abl.ablation.str.replace('_',' '),abl[k].fillna(0));a.tick_params(axis='y',labelsize=6);a.set(xlabel=l)
    finish_axes(ax);save(fig,SUP/'S5_ablation_matrix','Ablation matrix','candidate_693168_ablation_audit.csv','supplement')
    fig,ax=plt.subplots(1,4,figsize=(15,4),constrained_layout=True);x=np.arange(len(six))
    for a,k,l in zip(ax,('median_reduction','first_step_growth','closed_fraction_at_switch','window_width_C'),('Median reduction','First-step growth','Closed fraction','Window [°C]')):a.bar(x,six[k]);a.set(xticks=x,xticklabels=six.candidate_id.astype(str),ylabel=l);a.tick_params(axis='x',rotation=40)
    finish_axes(ax);save(fig,SUP/'S6_six_TierB_candidate_comparison','Six Tier-B candidates','six_tierB_candidate_comparison.csv','supplement')
    fig,ax=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for a,k,l in zip(ax,('median_reduction','window_width_C','complete'),('Median reduction','Window [°C]','Complete window')):
        p=rob.pivot(index='rho0',columns='G0_nm',values=k);im=a.imshow(p.values,origin='lower',aspect='auto');a.set(xticks=range(len(p.columns)),xticklabels=p.columns,yticks=range(len(p.index)),yticklabels=p.index,xlabel=r'$G_0$ [nm]',ylabel=r'$\rho_0$');fig.colorbar(im,ax=a,label=l)
    sty.letters(ax);save(fig,SUP/'S7_robustness_rho0_G0_heatmap','Initial-condition robustness','candidate_693168_extended_robustness.csv','supplement')

def dashboard(d,ratio,A,abl):
    fig,ax=plt.subplots(2,3,figsize=(14,8),constrained_layout=True);ax[0,0].axis('off');ax[0,0].text(.05,.95,'Candidate 693168\nconditional Tier B prototype\nnot calibrated or validated\nG₀≈103 nm; G₁≈117 nm\nclosed fraction at switch≈65%\nmain falsification target: closed topology',va='top',fontsize=11,bbox=dict(facecolor='#F3F3F3'))
    plot_paths(ax[0,1],d,'rho','G_mean_nm',['highT_reference','success']);ax[0,1].set(xlabel=r'$\rho$',ylabel='G [nm]');ax[0,2].plot(ratio.rho,ratio.reduction_TS);ax[0,2].set(xlabel=r'$\rho$',ylabel='Reduction');categorical_grid(ax[1,0],A,'T1_C','T2_C','Filled Chen map',annotate=(1400,1100));g=d[d.path_label=='success'];ax[1,1].plot(g.rho,g.phi_open,label='Open');ax[1,1].plot(g.rho,g.phi_closed,label='Closed');ax[1,1].set(xlabel=r'$\rho$',ylabel='Pore fraction');ax[1,1].legend();a=abl.copy();a.median_reduction=a.median_reduction.fillna(0);ax[1,2].barh(a.ablation.str.replace('_',' '),a.median_reduction);ax[1,2].tick_params(axis='y',labelsize=5);ax[1,2].set(xlabel='Median reduction');finish_axes(ax);save(fig,MAIN/'Figure_5_candidate_693168_dashboard','Candidate dashboard','all presentation source tables','main')

def main():
    sty.apply();
    fast=pd.read_csv(SRC/'dense_fast_firing_histories.csv');ratio_fast=pd.read_csv(SRC/'fast_firing_ratio_curves.csv');hmap=pd.read_csv(SRC/'fast_firing_heating_rate_map.csv');d=pd.read_csv(SRC/'dense_time_histories.csv');ratio=pd.read_csv(SRC/'dense_matched_density_curves.csv');diag=pd.read_csv(SRC/'T2_diagnostic_curves.csv');A=pd.read_csv(SRC/'chen_map_T1_T2_classification.csv');B=pd.read_csv(SRC/'chen_map_G1_T2_classification.csv');C=pd.read_csv(SRC/'chen_map_switch_density_T2_classification.csv');abl=pd.read_csv('results/audit_candidate_693168_closed_accommodation/candidate_693168_ablation_audit.csv');six=pd.read_csv('results/audit_candidate_693168_closed_accommodation/six_tierB_candidate_comparison.csv');rob=pd.read_csv('results/audit_candidate_693168_closed_accommodation/candidate_693168_extended_robustness.csv')
    fast_figures(fast,ratio_fast,hmap);two_step_figures(d,ratio,diag);chen_figures(A,B,C,diag);mechanism_figures(d,abl,six,rob);dashboard(d,ratio,A,abl)
    pd.DataFrame(INV).to_csv(SRC/'publication_style_figure_inventory.csv',index=False);print(f'figure_stems={len(INV)} files={2*len(INV)}')
if __name__=='__main__':main()
