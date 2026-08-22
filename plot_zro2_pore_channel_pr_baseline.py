"""Render diagnostic pore-channel/PR baseline figures."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
import numpy as np, pandas as pd

ROOT=Path(__file__).resolve().parent; OUT=ROOT/"results/zro2_forward_pore_channel_pr_baseline_test"; FIG=OUT/"figures"
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150})
COL={"DENSIFICATION_EXHAUSTION_FAILURE":"#4c78a8","SUCCESS":"#59a14f","GRAIN_GROWTH_FAILURE":"#e15759","MIXED_FAILURE":"#b07aa1"}
def save(fig,n):
 FIG.mkdir(parents=True,exist_ok=True)
 for e in ("pdf","png"): fig.savefig(FIG/f"{n}.{e}",dpi=220,bbox_inches="tight")
 plt.close(fig)

def schematic():
 fig,ax=plt.subplots(figsize=(10,5)); ax.axis('off')
 boxes=[(.03,.65,'Open channels'),(.28,.65,'Precursor'),(.52,.65,'Isolated'),(.76,.65,'Closed'),(.08,.18,'Surface coarsening'),(.34,.18,'PR pinch-off'),(.58,.18,'Regularization ↔ damage'),(.82,.18,'Zener / migration')]
 for x,y,s in boxes: ax.add_patch(FancyBboxPatch((x,y),.15,.13,boxstyle='round,pad=.02',fc='#eef4f7',ec='#31576b')); ax.text(x+.075,y+.065,s,ha='center',va='center')
 for a,b in [((.18,.715),(.28,.715)),((.43,.715),(.52,.715)),((.67,.715),(.76,.715)),((.155,.65),(.155,.31)),((.415,.31),(.355,.65)),((.655,.31),(.595,.65)),((.735,.25),(.82,.25))]: ax.add_patch(FancyArrowPatch(a,b,arrowstyle='->',mutation_scale=12,color='#555'))
 ax.text(.5,.93,'Conservative pore-channel topology audit',ha='center',fontsize=14,weight='bold'); ax.text(.5,.05,'Topology transfers conserve pore volume; only named open/closed shrinkage changes density.',ha='center',style='italic'); save(fig,'pore_channel_model_schematic')

def surface_scan(s):
 q=s[(s.lambda_seg_over_r==10)&(s.gamma_GB_over_gamma_s==.5)&(s.activity==.1)&(s.W_p==.8)&(s.connected_fraction==.6)]
 fig,axs=plt.subplots(1,2,figsize=(9,4))
 for r,z in q.groupby('radius_nm'): axs[0].semilogy(z.T_C,z.tau_s_s,label=f'{r:g} nm'); axs[1].semilogy(z.T_C,np.maximum(z.J_coars_sinv,1e-30),label=f'{r:g} nm')
 axs[0].set(xlabel='Temperature (°C)',ylabel=r'$\tau_s$ (s)'); axs[1].set(xlabel='Temperature (°C)',ylabel=r'$J_{coars}$ (s$^{-1}$)'); axs[1].legend(fontsize=7,ncol=2); fig.suptitle('Surface-diffusion r⁴ scaling'); fig.tight_layout(); save(fig,'surface_diffusion_r4_rate_scan')

def instability(s):
 q=s[(s.T_C==1100)&(s.radius_nm==25)&(s.activity==.1)&(s.W_p==.8)&(s.connected_fraction==.6)]
 fig,axs=plt.subplots(1,2,figsize=(9,4))
 for val,col in [('I_PR','PR instability index'),('P_pinch','Pinch probability')]:
  piv=q.pivot(index='gamma_GB_over_gamma_s',columns='lambda_seg_over_r',values=val); im=axs[0 if val=='I_PR' else 1].imshow(piv,origin='lower',aspect='auto',extent=[piv.columns.min(),piv.columns.max(),piv.index.min(),piv.index.max()]); axs[0 if val=='I_PR' else 1].set(title=col,xlabel=r'$\lambda_{seg}/r$',ylabel=r'$\gamma_{GB}/\gamma_s$'); fig.colorbar(im,ax=axs[0 if val=='I_PR' else 1])
 fig.tight_layout(); save(fig,'PR_instability_map')

def phase(s):
 q=s[(s.T_C==1100)&(s.radius_nm==25)&(s.lambda_seg_over_r==10)&(s.gamma_GB_over_gamma_s==.5)&(s.connected_fraction==.6)]
 piv=q.pivot(index='W_p',columns='activity',values='dominant_branch'); a=(piv=='damage').astype(int)
 fig,ax=plt.subplots(figsize=(7,4)); im=ax.imshow(a,origin='lower',aspect='auto',cmap=plt.matplotlib.colors.ListedColormap(['#4c78a8','#e15759']),extent=[min(piv.columns),max(piv.columns),min(piv.index),max(piv.index)]); ax.set(xlabel='Renewal activity',ylabel=r'$W_p=\sigma_{\ln r}$',title='Regularization (blue) versus damage (red)'); fig.colorbar(im,ax=ax,ticks=[0,1],label='dominant branch'); save(fig,'regularization_vs_damage_phase_diagram')

def heating(h):
 q=h[(h['mode']=='PR_regularization_damage_v1')&(h.peak_C==1500)&(h.hold_h==2)&(h.rate_C_min.isin([.2,10,100]))]
 fig,ax=plt.subplots(figsize=(7,4.5));
 for r,z in q.groupby('rate_C_min'): ax.plot(z.rho,z.G_nm,label=f'{r:g} °C/min')
 ax.set(xlabel='Relative density',ylabel='Grain size (nm)',title='Heating-rate trajectory separation'); ax.legend(); save(fig,'heating_rate_G_vs_rho')
 fig,axs=plt.subplots(2,2,figsize=(9,6),sharex=True)
 for r,z in q.groupby('rate_C_min'):
  axs[0,0].plot(z.rho,2*z.r_nm,label=f'{r:g}'); axs[0,1].plot(z.rho,2*z.r_nm*np.exp(1.2816*z.W)); axs[1,0].plot(z.rho,z.conn); axs[1,1].plot(z.rho,z.damage_memory-z.reg_memory)
 axs[0,0].set(ylabel='D50 (nm)'); axs[0,1].set(ylabel='D90 (nm)'); axs[1,0].set(ylabel='Connected fine fraction',xlabel='Density'); axs[1,1].set(ylabel='Damage − regularization memory',xlabel='Density'); axs[0,0].legend(title='°C/min',fontsize=7); fig.suptitle('Heating-rate pore memory'); fig.tight_layout(); save(fig,'heating_rate_pore_memory')

def boundary(b):
 q=b[b.state_id=='natural_selected']; fig,axs=plt.subplots(1,3,figsize=(12,3.5),sharex=True)
 for ax,(mode,z) in zip(axs,q.groupby('mode')):
  ax.scatter(z.T2_C,z.G_nm,c=[COL[x] for x in z.classification],marker='s',s=35); ax.set(title=mode.replace('_',' '),xlabel='T2 (°C)',ylabel='G1 (nm)')
 fig.suptitle('Naturally prepared failure modes — no bounded process map run'); fig.tight_layout(); save(fig,'twostep_Chen_G1_T2_failure_modes')

def trajectories(t):
 fig,axs=plt.subplots(1,3,figsize=(11,3.5))
 labels={850:'low-T hold',1100:'two-step diagnostic',1300:'high-T hold'}
 for T,z in t.groupby('T2_C'):
  axs[0].plot(z.time_h,z.rho,label=labels.get(T,str(T))); axs[1].plot(z.time_h,z.G_nm); axs[2].plot(z.rho,z.G_nm)
 axs[0].set(xlabel='Time (h)',ylabel='Density'); axs[1].set(xlabel='Time (h)',ylabel='Grain size (nm)'); axs[2].set(xlabel='Density',ylabel='Grain size (nm)'); axs[0].legend(fontsize=7); fig.suptitle('Fixed-state hold trajectory comparison'); fig.tight_layout(); save(fig,'twostep_trajectory_comparison')

def ablations(a):
 z=a.sort_values('density_gain'); fig,axs=plt.subplots(1,2,figsize=(11,6)); axs[0].barh(z.ablation,z.density_gain,color='#4c78a8'); axs[1].barh(z.ablation,z.growth_fraction,color='#e15759'); axs[0].set(xlabel='Density gain'); axs[1].set(xlabel='Growth fraction'); fig.suptitle('Pore-channel mechanism ablations'); fig.tight_layout(); save(fig,'ablation_matrix')

def inventory():
 items={'pore_channel_model_schematic':'pore_channel_parameter_registry.csv','surface_diffusion_r4_rate_scan':'pore_channel_fixed_state_scan.csv','PR_instability_map':'pore_channel_fixed_state_scan.csv','regularization_vs_damage_phase_diagram':'pore_channel_fixed_state_scan.csv','heating_rate_G_vs_rho':'pore_channel_heating_rate_histories.csv','heating_rate_pore_memory':'pore_channel_heating_rate_histories.csv','twostep_Chen_G1_T2_failure_modes':'pore_channel_boundary_preservation_test.csv','twostep_trajectory_comparison':'pore_channel_twostep_histories.csv','ablation_matrix':'pore_channel_ablation_matrix.csv'}
 rows=[]
 for n,src in items.items(): rows.append({'figure_id':n,'pdf':f'figures/{n}.pdf','png':f'figures/{n}.png','source_table':src,'pdf_nonempty':(FIG/f'{n}.pdf').stat().st_size>1000,'png_nonempty':(FIG/f'{n}.png').stat().st_size>1000,'placeholder':False,'success_colored_map':False,'validation_claim':False})
 pd.DataFrame(rows).to_csv(OUT/'figure_inventory.csv',index=False)

def main():
 s=pd.read_csv(OUT/'pore_channel_fixed_state_scan.csv'); h=pd.read_csv(OUT/'pore_channel_heating_rate_histories.csv'); b=pd.read_csv(OUT/'pore_channel_boundary_preservation_test.csv'); t=pd.read_csv(OUT/'pore_channel_twostep_histories.csv'); a=pd.read_csv(OUT/'pore_channel_ablation_matrix.csv')
 schematic(); surface_scan(s); instability(s); phase(s); heating(h); boundary(b); trajectories(t); ablations(a); inventory(); print(FIG)
if __name__=='__main__': main()
