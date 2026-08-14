#!/usr/bin/env python3
"""Dynamic latent topology states coupled to migration, with frozen rho_dot."""
from dataclasses import dataclass
import math,numpy as np
import separated_fast_chen_model as base

R=8.31446261815324
@dataclass
class LatentState:
    C_rem_GBseg:float=.7;C_rem_TJ:float=.3;f_large_tail:float=.05;f_isolated:float=.05;pore_percolation_state:float=.8;removable_pore_memory:float=.8
    X_J:float=0.;junction_density:float=1.;segment_length:float=1.;C_TJ_constraint:float=.2;C_TJ_pore:float=.2;C_TJ_pinned:float=.1;C_TJ_relaxed:float=.1;structural_TJ_memory:float=0.;pore_filled_TJ_memory:float=0.
    sigma_res_GBseg:float=0.;sigma_res_TJ:float=0.;sigma_res_cleanGB:float=0.;stored_PR_work:float=0.;stored_shear_work:float=0.;stress_memory:float=0.;stress_relaxation_state:float=0.

def derivatives(x,rates,T,p):
    """Instantaneous latent evolution law from local state."""
    pr=rates['PR_propensity'];a=rates['activity'];rd=rates['rho_dot'];gb=rates['growth_base'];thermal=math.exp(-p['Q_relax']/R/(T+273.15))
    drem=p['recovery_sweep']*a*(1-x.removable_pore_memory)-p['generation_PR']*pr*x.removable_pore_memory-p['loss_dens']*rd*x.removable_pore_memory-p['isolation_rate']*pr*x.removable_pore_memory
    dxj=p['XJ_prod_TJ']*rd*(p['XJ_capacity']-x.X_J)+p['XJ_prod_sweep']*gb*(p['XJ_capacity']-x.X_J)-x.X_J/max(p['tau_J'],1)
    dstress=p['PR_stress_generation']*pr+p['shear_generation']*gb-x.stress_memory*(p['stress_relax_dens']*rd+thermal/max(p['tau_stress'],1))
    return drem,dxj,dstress

def advance(x,rates,T,p,dt):
    drem,dxj,ds=derivatives(x,rates,T,p);x.removable_pore_memory=float(np.clip(x.removable_pore_memory+drem*dt,0,1));x.C_rem_GBseg=.7*x.removable_pore_memory;x.C_rem_TJ=.3*x.removable_pore_memory;x.f_isolated=float(np.clip(x.f_isolated+p['isolation_rate']*rates['PR_propensity']*dt,0,1));x.f_large_tail=float(np.clip(x.f_large_tail+.5*p['generation_PR']*rates['PR_propensity']*dt,0,1));x.pore_percolation_state=x.removable_pore_memory;x.X_J=float(np.clip(x.X_J+dxj*dt,0,p['XJ_capacity']));x.structural_TJ_memory=x.X_J;x.junction_density=1+x.X_J;x.segment_length=1/max(x.junction_density,1e-9);x.C_TJ_constraint=x.X_J*(1-p['pore_relax_fraction']*x.C_TJ_pore);x.C_TJ_pinned=x.C_TJ_pore*p['pore_drag_fraction'];x.C_TJ_relaxed=x.C_TJ_pore*p['pore_relax_fraction'];x.pore_filled_TJ_memory=x.C_TJ_pinned;x.stress_memory=float(np.clip(x.stress_memory+ds*dt,0,p['stress_cap']));x.sigma_res_GBseg=.5*x.stress_memory;x.sigma_res_TJ=.35*x.stress_memory;x.sigma_res_cleanGB=.15*x.stress_memory;x.stored_PR_work+=rates['PR_propensity']*dt;x.stored_shear_work+=rates['growth_base']*dt;x.stress_relaxation_state=ds;return x

def migration_factor(x,T,G,p):
    K=p['K0']*(G/150e-9)**p['q_TJ'];thermal=math.exp(float(np.clip(-p['Q_event']/R*(1/(T+273.15)-1/1473.15),-30,30)));comp=1-math.exp(-p['lambda_ref']*thermal/max(K,1e-9));drag=p['A_J']*x.X_J+p['pore_drag_strength']*x.C_TJ_pinned+p['stress_coupling']*x.stress_memory
    return float(np.clip(comp/(1+drag),0,1)),comp,K

def shared_rates(material,state,T):return base.material_rates(state['rho'],state['G'],state['phi'],state['radii'],T,material)
