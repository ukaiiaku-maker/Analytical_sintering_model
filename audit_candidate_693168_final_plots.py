#!/usr/bin/env python3
"""Publication-style visual audit for conditional Tier-B candidate 693168."""
from pathlib import Path
import csv
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Rectangle

import plot_style as ps

ROOT=Path("results/audit_candidate_693168_closed_accommodation")
FIG=ROOT/"final_figures"; TABLES=ROOT/"final_tables"
COL={"highT_reference":"#CC79A7","success":"#009E73","lower_failure":"#0072B2","upper_failure":"#D55E00"}
CLASS={"DENSIFICATION_EXHAUSTION_FAILURE":"#0072B2","SUCCESS":"#009E73","GRAIN_GROWTH_FAILURE":"#D55E00","MIXED_FAILURE":"#CC79A7","UNATTAINABLE_FIRST_STEP":"#777777"}
INVENTORY=[]

def save(fig, rel, title, source, purpose):
    path=FIG/rel;ps.finish(fig,path)
    INVENTORY.append(dict(figure_id=len(INVENTORY)+1,filename_pdf=str(path.with_suffix('.pdf').relative_to(ROOT)),filename_png=str(path.with_suffix('.png').relative_to(ROOT)),short_title=title,source_table_or_script=source,purpose=purpose,manuscript_location="Candidate 693168 audit"))

def lines(ax,d,x,y,paths=None,ylabel=None,legend=True):
    use=d if paths is None else d[d.path_label.isin(paths)]
    for label,g in use.groupby("path_label",sort=False):
        ax.plot(g[x],g[y],label=label.replace('_',' '),color=COL.get(label),alpha=.95)
    ax.set(xlabel={"rho":r"Relative density, $\rho$","physical_time_h":r"Time, $t$ [h]","T2_C":r"Second-step temperature, $T_2$ [°C]"}.get(x,x),ylabel=ylabel or y)
    if legend and use.path_label.nunique()>1:ax.legend(fontsize=7)
    ps.clean(ax)

def panels(shape,size):
    fig,axes=plt.subplots(*shape,figsize=size,constrained_layout=True);return fig,np.asarray(axes).reshape(-1)

def dashboard(d,ratio,points,abl,repro):
    fig,ax=panels((2,3),(15,8))
    ax[0].axis('off');r=repro.iloc[0]
    text=("Candidate 693168\nconditional Tier B — not validation\n"
          f"G₁ = {r.G1_nm:.1f} nm; prep growth = {100*r.first_step_growth_fraction:.1f}%\n"
          f"reduction median = {100*r.median_reduction:.1f}%\n"
          f"fine window = {r.first_success_C:.0f}–{r.last_success_C:.0f} °C\n"
          "causal core: PR + closed transition + finite accommodation")
    ax[0].text(.03,.95,text,va='top',fontsize=11,bbox=dict(boxstyle='round',facecolor='#F2F2F2'))
    lines(ax[1],d,'rho','G_mean_nm',['highT_reference','success'],r"Grain size, $G$ [nm]")
    ax[2].plot(ratio.rho,100*ratio.reduction_TS,color=COL['success']);ax[2].axhline(20,ls='--',color='.4');ax[2].set(xlabel=r"Relative density, $\rho$",ylabel="Two-step reduction [%]");ps.clean(ax[2])
    good=points[points.classification=='SUCCESS'];ax[3].axvspan(good.T2_C.min(),good.T2_C.max(),color=CLASS['SUCCESS'],alpha=.3)
    ax[3].scatter(points.T2_C,points.rho2,c=[CLASS.get(x,'#777777') for x in points.classification],s=12);ax[3].axhline(.98,ls='--',color='.4');ax[3].set(xlabel=r"$T_2$ [°C]",ylabel=r"Final $\rho$");ps.clean(ax[3])
    lines(ax[4],d,'rho','closed_accommodation_factor',['lower_failure','success','upper_failure'],r"Closed accommodation, $A_{closed}$")
    a=abl.copy();a['median_reduction']=a.median_reduction.fillna(0);ax[5].barh(a.ablation.str.replace('_',' '),100*a.median_reduction,color=np.where(a.complete,'#009E73','#D55E00'));ax[5].set(xlabel="Median reduction [%]");ax[5].tick_params(axis='y',labelsize=6);ps.clean(ax[5])
    ps.panel_labels(ax);save(fig,'candidate_693168_dashboard','Candidate overview dashboard','dense histories; fine T2; ablations','Joint visual summary')

def full_time(d,ratio):
    fig,ax=panels((4,2),(13,14));vars=[('T_C',r"Temperature, $T$ [°C]"),('rho',r"Relative density, $\rho$"),('G_mean_nm',r"Grain size, $G$ [nm]"),('rho_dot',r"Densification rate, $\dot{\rho}$ [s$^{-1}$]"),('G_dot',r"Growth rate, $\dot{G}$ [nm s$^{-1}$]"),('activity','Renewal activity')]
    for a,(v,l) in zip(ax[:6],vars):lines(a,d,'physical_time_h',v,None,l)
    lines(ax[6],d,'rho','G_mean_nm',None,r"Grain size, $G$ [nm]")
    ax[7].plot(ratio.rho,100*ratio.reduction_TS,color=COL['success']);ax[7].set(xlabel=r"Relative density, $\rho$",ylabel="Reduction [%]");ps.clean(ax[7])
    for a in ax[:6]:
        for _,g in d[d.stage=='second_step'].groupby('path_label'):
            a.axvline(g.physical_time_h.min(),color='.5',lw=.7,ls=':')
    ps.panel_labels(ax);save(fig,'time_evolution/candidate_693168_full_time_evolution','Full physical-time evolution','dense_candidate_693168_histories.csv','Compare high-T, lower failure, success, and upper failure')

def pore_store(d):
    fig,ax=panels((4,2),(13,14)); specs=[('phi_connected','Connected/open pore fraction'),('phi_GBseg','GB-segment pore fraction'),('phi_TJ','TJ pore fraction'),('phi_iso','Isolated pore fraction'),('phi_closed',r"Closed-pore fraction, $f_{closed}$"),('N_open','Open-pore number proxy'),('N_closed','Closed-pore number proxy')]
    for a,(v,l) in zip(ax[:7],specs):lines(a,d,'rho',v,None,l)
    for label,g in d.groupby('path_label'):ax[7].plot(g.rho,g.pore_radius_proxy_D50,label=f'D50 {label}',color=COL.get(label));ax[7].plot(g.rho,g.pore_radius_proxy_D90,ls='--',color=COL.get(label),label=f'D90 {label}')
    ax[7].set(xlabel=r"Relative density, $\rho$",ylabel="Pore-radius proxy");ax[7].legend(fontsize=6,ncol=2);ps.clean(ax[7]);ps.panel_labels(ax)
    save(fig,'time_evolution/candidate_693168_pore_store_evolution','Pore-store and number evolution','dense histories','Audit closed-pore population and topology transition')

def accommodation(d):
    fig,ax=panels((4,2),(13,14));specs=[('closed_fraction',r"Closed-pore fraction, $f_{closed}$"),('closed_shrinkage_flux','Closed shrinkage flux'),('open_shrinkage_flux','Open shrinkage flux'),('closed_accommodation_capacity','Accommodation capacity'),('closed_accommodation_used','Accommodation used'),('closed_accommodation_factor',r"$\Lambda_{closed}/K_{closed}$ proxy"),('P_comp_closed',r"$P_{comp,closed}$"),('closed_pore_contribution_to_rho_dot',r"Closed contribution to $\dot{\rho}$")]
    for a,(v,l) in zip(ax,specs):lines(a,d,'rho',v,None,l)
    ps.panel_labels(ax);save(fig,'time_evolution/candidate_693168_closed_accommodation','Closed-pore accommodation','dense histories','Expose finite accommodation and lower boundary')

def pr_panel(d):
    fig,ax=panels((4,2),(13,14));specs=[('PR_redistribution_rate','PR redistribution rate'),('rho_dot','Densification rate'),('PR_damage_memory','PR damage memory'),('cumulative_PR_redistributed_volume','Cumulative PR redistribution'),('cumulative_open_pore_removed','Cumulative open-pore removal'),('cumulative_closed_pore_removed','Cumulative closed-pore removal'),('closed_fraction','Closed fraction'),('removable_pore_fraction','Removable pore fraction')]
    for a,(v,l) in zip(ax,specs):lines(a,d,'rho',v,None,l)
    ps.panel_labels(ax);save(fig,'time_evolution/candidate_693168_PR_energy_topology_memory','PR and topology memory','dense histories','Distinguish named PR redistribution from densification')

def migration_panel(d):
    fig,ax=panels((4,2),(13,14));specs=[('migration_factor',r"Migration factor, $\Gamma_{mig}$"),('X_J',r"$X_J$"),('C_TJ',r"$C_{TJ}$"),('C_GBseg',r"$C_{GBseg}$"),('f_clean_GB','Clean-GB fraction'),('Lambda_over_K_TJ',r"$\Lambda_{TJ}/K_{TJ}$"),('P_comp_TJ',r"$P_{comp,TJ}$"),('pore_drag','Pore drag')]
    for a,(v,l) in zip(ax,specs):lines(a,d,'rho',v,None,l)
    ps.panel_labels(ax);save(fig,'time_evolution/candidate_693168_migration_topology_channels','Migration/topology channels','dense histories','Audit secondary migration mechanisms')

def chen(points,bound):
    present=list(points.classification.unique())
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
    for i,c in enumerate(present):
        g=points[points.classification==c];ax.scatter(g.T2_C,np.full(len(g),i),label=c.replace('_',' '),color=CLASS.get(c,'#777777'),s=24)
    ax.axvline(1400,color='.4',ls=':',label=r'$T_2=T_1$');ax.set(xlabel=r"Second-step temperature, $T_2$ [°C]",ylabel="Classification",yticks=[]);ax.legend(fontsize=7,ncol=2);ps.clean(ax)
    save(fig,'chen_maps/candidate_693168_complete_Chen_classification_map','Complete Chen classification map','fine T2 classification','Show all represented boundary classes')
    fig,ax=plt.subplots(figsize=(8,4),constrained_layout=True);good=points[points.classification=='SUCCESS'];ax.axvspan(good.T2_C.min(),good.T2_C.max(),color=CLASS['SUCCESS'],alpha=.32,label='Tier B success band');ax.scatter(points.T2_C,points.rho2,c=[CLASS.get(x,'#777777') for x in points.classification],s=18);ax.axhline(.98,color='.4',ls='--');ax.annotate(f"window = {float(bound.window_width_C):.0f} °C",xy=(good.T2_C.mean(),.982),ha='center');ax.set(xlabel=r"Second-step temperature, $T_2$ [°C]",ylabel=r"Final density, $\rho$");ax.legend();ps.clean(ax)
    save(fig,'chen_maps/candidate_693168_complete_Chen_filled_window','Complete Chen filled window','fine T2 classification','Show finite Tier-B success band')
    fig,ax=panels((4,2),(13,13));specs=[('rho2',r"Final $\rho$"),('G2_nm',r"Final $G$ [nm]"),('growth_fraction','Second-step growth fraction'),('closed_shrinkage_contribution','Integrated closed shrinkage'),('open_shrinkage_contribution','Integrated open shrinkage'),('closed_accommodation_factor','Final accommodation factor'),('PR_damage_state','Final PR memory')]
    for a,(v,l) in zip(ax[:7],specs):a.plot(points.T2_C,points[v],color='#333333');a.set(xlabel=r"$T_2$ [°C]",ylabel=l);ps.clean(a)
    codes={c:i for i,c in enumerate(present)};ax[7].scatter(points.T2_C,[codes[c] for c in points.classification],c=[CLASS.get(c,'#777') for c in points.classification],s=18);ax[7].set(xlabel=r"$T_2$ [°C]",ylabel='Classification code',yticks=list(codes.values()),yticklabels=[x.replace('_',' ') for x in codes]);ax[7].tick_params(axis='y',labelsize=6);ps.clean(ax[7]);ps.panel_labels(ax)
    save(fig,'chen_maps/candidate_693168_T2_diagnostics','Fine T2 diagnostics','fine T2 classification','Locate lower and upper boundaries')

def comparisons(abl,six,robust,fast):
    a=abl.copy();a['median_reduction']=a.median_reduction.fillna(0);fig,ax=panels((3,2),(14,13));metrics=[('median_reduction','Median reduction'),('span20',r'Span $\geq20\%$'),('window_width_C','Window width [°C]')]
    for z,(v,l) in zip(ax[:3],metrics):z.barh(a.ablation.str.replace('_',' '),a[v],color=np.where(a.complete,'#009E73','#D55E00'));z.set(xlabel=l);z.tick_params(axis='y',labelsize=6);ps.clean(z)
    for z,v,l in ((ax[3],'lower_bracketed','Lower boundary'),(ax[4],'upper_bracketed','Upper boundary'),(ax[5],'high_density_attainment','Both paths attain')):z.scatter(a[v].astype(int),np.arange(len(a)),c=np.where(a[v],'#009E73','#D55E00'));z.set(xlabel=l,yticks=np.arange(len(a)),yticklabels=a.ablation.str.replace('_',' '));z.tick_params(axis='y',labelsize=6);ps.clean(z)
    ps.panel_labels(ax);save(fig,'ablations/candidate_693168_ablation_waterfall','Ablation waterfall','candidate_693168_ablation_audit.csv','Separate attainment, boundary, and reduction causality')
    fig,ax=panels((4,2),(13,14));specs=[('median_reduction','Median reduction'),('min_reduction','Minimum reduction'),('window_width_C','Window [°C]'),('first_step_growth','First-step growth'),('closed_fraction_at_switch','Closed fraction at switch'),('span20',r'Span $\geq20\%$'),('robustness_cases_passed','Robustness cases')]
    x=np.arange(len(six));labs=six.candidate_id.astype(str)
    for z,(v,l) in zip(ax[:7],specs):
        z.bar(x,six[v],color=np.where(six.candidate_id==693168,'#009E73','#56B4E9'))
        z.set(xticks=x,xticklabels=labs,ylabel=l);z.tick_params(axis='x',rotation=35);ps.clean(z)
    ax[7].axis('off');ax[7].text(.05,.9,'Conditional Tier-B family\nCandidate 693168 is an extreme member\nAll require calibration',va='top');ps.panel_labels(ax)
    save(fig,'comparison/six_TierB_candidate_comparison','Six Tier-B candidates','six_tierB_candidate_comparison.csv','Determine whether 693168 is an outlier')
    med=robust.pivot(index='rho0',columns='G0_nm',values='median_reduction');win=robust.pivot(index='rho0',columns='G0_nm',values='window_width_C');fig,ax=panels((1,2),(12,5))
    for z,table,title,cmap in ((ax[0],med,'Median reduction','viridis'),(ax[1],win,'Window width [°C]','magma')):
        im=z.imshow(table.values,origin='lower',aspect='auto',cmap=cmap);z.set(xticks=np.arange(len(table.columns)),xticklabels=[f'{x:g}' for x in table.columns],yticks=np.arange(len(table.index)),yticklabels=[f'{x:.2f}' for x in table.index],xlabel=r'$G_0$ [nm]',ylabel=r'$\rho_0$',title=title);fig.colorbar(im,ax=z)
    ps.panel_labels(ax);save(fig,'comparison/candidate_693168_robustness_heatmap','Candidate robustness heatmap','extended robustness','Map bounded initial-condition domain')
    fig,ax=panels((2,2),(11,8));
    for z,(mid,g) in zip(ax[:2],fast.groupby('material_id')):
        z.bar(g['mode'],g.G_ref_over_G_fast,color=['#009E73' if x else '#D55E00' for x in g.fast_firing_retained]);z.axhline(1.5,ls='--',color='.4');z.set(ylabel=r'$G_{ref}/G_{fast}$',title=mid);z.tick_params(axis='x',rotation=20);ps.clean(z)
    for z,(mid,g) in zip(ax[2:],fast.groupby('material_id')):
        z.bar(g['mode'],g.density_span_ge_1p5,color='#56B4E9');z.set(ylabel=r'Density span $\geq1.5$',title=f'{mid}: nucleation-limited envelope');z.tick_params(axis='x',rotation=20);ps.clean(z)
    ps.panel_labels(ax);save(fig,'comparison/fast_firing_preservation_candidate_693168','Fast-firing preservation','frozen E0021/E0142 audit','Confirm nucleation-limited, PR-independent envelope')

def schematic():
    fig,ax=plt.subplots(figsize=(13,4),constrained_layout=True);ax.set_xlim(0,13);ax.set_ylim(0,4);ax.axis('off')
    titles=['First-step PR preparation','Closed-store transition','Low $T_2$: exhaustion','Intermediate $T_2$: success','High $T_2$: growth']
    colors=['#56B4E9','#CC79A7','#0072B2','#009E73','#D55E00']
    for i,(title,color) in enumerate(zip(titles,colors)):
        x=.2+2.55*i;ax.add_patch(Rectangle((x,.5),2.15,2.8,facecolor=color,alpha=.14,edgecolor=color));ax.text(x+1.075,3.05,title,ha='center',va='top',fontsize=9,fontweight='bold')
        for j in range(4):ax.add_patch(Circle((x+.4+.42*j,1.7+.3*(j%2)),.13+.04*(i==1),facecolor='none',edgecolor=color,lw=2))
        if i<4:ax.add_patch(FancyArrowPatch((x+2.17,1.9),(x+2.52,1.9),arrowstyle='->',mutation_scale=15,color='.3'))
    ax.text(6.5,.1,'Conceptual mechanism schematic — not a calibrated micrograph',ha='center',fontsize=9,style='italic')
    save(fig,'candidate_693168_mechanism_schematic','Mechanism schematic','conceptual synthesis','Explain PR-prepared closed accommodation and three T2 outcomes')

def final_panel(d,ratio,points,abl,six,robust):
    fig,ax=panels((2,4),(18,8));lines(ax[0],d,'rho','G_mean_nm',['highT_reference','success'],r'$G$ [nm]');ax[1].plot(ratio.rho,100*ratio.reduction_TS,color=COL['success']);ax[1].set(xlabel=r'$\rho$',ylabel='Reduction [%]');ps.clean(ax[1]);good=points[points.classification=='SUCCESS'];ax[2].axvspan(good.T2_C.min(),good.T2_C.max(),color='#009E73',alpha=.3);ax[2].scatter(points.T2_C,points.rho2,c=[CLASS.get(x,'#777') for x in points.classification],s=8);ax[2].set(xlabel=r'$T_2$ [°C]',ylabel=r'Final $\rho$');ps.clean(ax[2]);lines(ax[3],d,'rho','closed_accommodation_factor',['lower_failure','success','upper_failure'],r'$A_{closed}$');lines(ax[4],d,'rho','closed_fraction',['highT_reference','success'],r'$f_{closed}$');a=abl.copy();a['median_reduction']=a.median_reduction.fillna(0);ax[5].barh(a.ablation.str.replace('_',' '),a.median_reduction,color=np.where(a.complete,'#009E73','#D55E00'));ax[5].tick_params(axis='y',labelsize=5);ax[5].set(xlabel='Median reduction');ps.clean(ax[5]);ax[6].bar(six.candidate_id.astype(str),six.median_reduction,color=np.where(six.candidate_id==693168,'#009E73','#56B4E9'));ax[6].tick_params(axis='x',rotation=45);ax[6].set(ylabel='Median reduction');ps.clean(ax[6]);med=robust.pivot(index='rho0',columns='G0_nm',values='median_reduction');im=ax[7].imshow(med.values,origin='lower',aspect='auto');ax[7].set(xlabel=r'$G_0$',ylabel=r'$\rho_0$');fig.colorbar(im,ax=ax[7],shrink=.7);ps.panel_labels(ax);fig.suptitle('Candidate 693168: conditional Tier B, prototype-scale and not calibrated',fontweight='bold')
    ax[7].set(xticks=np.arange(len(med.columns)),xticklabels=[f'{x:g}' for x in med.columns],
              yticks=np.arange(len(med.index)),yticklabels=[f'{x:.2f}' for x in med.index],
              xlabel=r'$G_0$ [nm]',ylabel=r'$\rho_0$')
    save(fig,'candidate_693168_final_panel','Final manuscript-style panel','all final audit tables','Compact visual inspection package')

def main():
    ps.apply_style();FIG.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(ROOT/'dense_candidate_693168_histories.csv');ratio=pd.read_csv(ROOT/'dense_candidate_693168_matched_density_curves.csv');points=pd.read_csv(ROOT/'candidate_693168_T2_classification_points_fine.csv');bound=pd.read_csv(ROOT/'candidate_693168_T2_window_boundaries_fine.csv').iloc[0];abl=pd.read_csv(ROOT/'candidate_693168_ablation_audit.csv');six=pd.read_csv(ROOT/'six_tierB_candidate_comparison.csv');robust=pd.read_csv(ROOT/'candidate_693168_extended_robustness.csv');fast=pd.read_csv(ROOT/'fast_firing_preservation_audit.csv');repro=pd.read_csv(ROOT/'candidate_693168_reproduction_summary.csv')
    dashboard(d,ratio,points,abl,repro);full_time(d,ratio);pore_store(d);accommodation(d);pr_panel(d);migration_panel(d);chen(points,bound);comparisons(abl,six,robust,fast);schematic();final_panel(d,ratio,points,abl,six,robust)
    inv=pd.DataFrame(INVENTORY);inv.to_csv(ROOT/'final_figure_inventory.csv',index=False);inv.to_csv(TABLES/'final_figure_inventory.csv',index=False)
    aliases={
      'candidate_693168_G_rho_absolute':'candidate_693168_final_panel',
      'candidate_693168_reduction_vs_density':'candidate_693168_dashboard',
      'candidate_693168_physical_time_histories':'time_evolution/candidate_693168_full_time_evolution',
      'candidate_693168_closed_pore_fraction_vs_density':'time_evolution/candidate_693168_pore_store_evolution',
      'candidate_693168_closed_accommodation_vs_density':'time_evolution/candidate_693168_closed_accommodation',
      'candidate_693168_open_vs_closed_shrinkage_flux':'time_evolution/candidate_693168_closed_accommodation',
      'candidate_693168_PR_damage_closed_transition_history':'time_evolution/candidate_693168_PR_energy_topology_memory',
      'candidate_693168_fine_Chen_classification_map':'chen_maps/candidate_693168_complete_Chen_classification_map',
      'candidate_693168_fine_Chen_filled_window':'chen_maps/candidate_693168_complete_Chen_filled_window',
      'candidate_693168_T2_diagnostics':'chen_maps/candidate_693168_T2_diagnostics',
      'candidate_693168_ablation_waterfall':'ablations/candidate_693168_ablation_waterfall',
      'six_candidates_ablation_summary':'comparison/six_TierB_candidate_comparison',
      'candidate_693168_robustness_heatmap':'comparison/candidate_693168_robustness_heatmap',
      'fast_firing_preservation_candidate_693168':'comparison/fast_firing_preservation_candidate_693168',
    }
    for alias,source in aliases.items():
        for suffix in ('.pdf','.png'):
            shutil.copyfile(FIG/(source+suffix),ROOT/(alias+suffix))
    print(f'figure_stems={len(inv)} files={2*len(inv)}')

if __name__=='__main__':main()
