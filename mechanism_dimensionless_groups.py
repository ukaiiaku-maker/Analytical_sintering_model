#!/usr/bin/env python3
"""Dimensionless groups for topology-frozen material-property attribution.

The fast-firing and local closed-pore models are deliberately kept separate.
These functions compare their relative timescales without creating a hidden
cross-model density channel.
"""
from __future__ import annotations

import numpy as np

R = 8.31446261815324
TREF_K = 1473.15


def thermal_relative(T_C, Q):
    T=np.asarray(T_C,float)+273.15
    return np.exp(np.clip(-np.asarray(Q,float)/R*(1/T-1/TREF_K),-60,60))


def fast_groups(q_nuc,q_exchange,q_transport,q_growth,q_surface,k_nuc,k_exchange,k_transport,k_growth,k_surface,
                *,T_C=1300.,G_m=75e-9,stress_term=0.,temperature_grid=None):
    """Vectorized serial-time and low-activity groups for material rows."""
    q_nuc=np.asarray(q_nuc);T=T_C+273.15
    tau_nuc=np.exp(np.clip(q_nuc/(R*T)-stress_term,-50,50))/np.maximum(k_nuc,1e-300)
    tau_exchange=k_exchange*np.exp(np.clip(np.asarray(q_exchange)/(R*T),-50,50))
    tau_transport=k_transport*G_m**2*np.exp(np.clip(np.asarray(q_transport)/(R*T),-50,50))
    total=tau_nuc+tau_exchange+tau_transport
    activity=(tau_exchange+tau_transport)/np.maximum(total,1e-300)
    theta=tau_nuc/np.maximum(tau_exchange+tau_transport,1e-300)
    out=dict(Theta_nuc=theta,f_nuc=tau_nuc/total,f_exchange=tau_exchange/total,
             f_transport=tau_transport/total,activity=activity,tau_nuc_s=tau_nuc,
             tau_exchange_s=tau_exchange,tau_transport_s=tau_transport)
    if temperature_grid is None:return out
    temps=np.asarray(temperature_grid,float);TK=temps[:,None]+273.15
    tn=np.exp(np.clip(q_nuc[None,:]/(R*TK)-stress_term,-50,50))/np.maximum(np.asarray(k_nuc)[None,:],1e-300)
    te=np.asarray(k_exchange)[None,:]*np.exp(np.clip(np.asarray(q_exchange)[None,:]/(R*TK),-50,50))
    tt=np.asarray(k_transport)[None,:]*G_m**2*np.exp(np.clip(np.asarray(q_transport)[None,:]/(R*TK),-50,50))
    act=(te+tt)/np.maximum(tn+te+tt,1e-300)
    surface=np.asarray(k_surface)[None,:]*np.exp(np.clip(-np.asarray(q_surface)[None,:]/(R*TK),-50,50))
    growth=np.asarray(k_growth)[None,:]*np.exp(np.clip(-np.asarray(q_growth)[None,:]/(R*TK),-50,50))
    dT=np.gradient(temps)[:,None]
    low=(1-act)
    out.update(I_low_slow=np.sum(low*dT*60,axis=0),I_low_fast=np.sum(low*dT*60/50,axis=0),
               I_low_PR_slow=np.sum(low*surface*dT*60,axis=0),I_low_PR_fast=np.sum(low*surface*dT*60/50,axis=0),
               coarsening_exposure_slow=np.sum(low*growth*dT*60,axis=0),
               coarsening_exposure_fast=np.sum(low*growth*dT*60/50,axis=0),
               Pi_PR=np.mean(low*surface/np.maximum(act,1e-30),axis=0))
    return out


def solve_boundary_temperature(Q,k_factor,base_Q,base_T_C):
    """Temperature giving the base Arrhenius rate after Q/prefactor change."""
    Q=np.asarray(Q,float);k=np.maximum(np.asarray(k_factor,float),1e-300)
    base_log=-base_Q/R*(1/(base_T_C+273.15)-1/TREF_K)
    invT=1/TREF_K-(R/Q)*(base_log-np.log(k))
    return 1/np.maximum(invT,1e-12)-273.15


def two_step_groups(q_closed,q_growth,q_pr,k_closed,k_growth,k_pr,*,base_q_closed,base_q_growth,
                    base_q_pr,T2_C=1100.,T1_C=1400.,closed_fraction=.649415,
                    accommodation_fraction=.1674,lambda_over_k=.0,p_comp_tj=.0):
    q_closed=np.asarray(q_closed);q_growth=np.asarray(q_growth);q_pr=np.asarray(q_pr)
    kc=np.asarray(k_closed)*thermal_relative(T2_C,q_closed)
    kg=np.asarray(k_growth)*thermal_relative(T2_C,q_growth)
    kp=np.asarray(k_pr)*thermal_relative(T1_C,q_pr)
    selectivity=kc/np.maximum(kg,1e-300)
    base_sel=thermal_relative(T2_C,base_q_closed)/thermal_relative(T2_C,base_q_growth)
    prep_rel=kp/np.maximum(thermal_relative(T1_C,base_q_pr),1e-300)
    return dict(S_closed_growth=selectivity,selectivity_relative=selectivity/base_sel,
                M_PR_closed=closed_fraction*prep_rel/(1+prep_rel),PR_preparation_relative=prep_rel,
                A_closed_fraction=np.full_like(selectivity,accommodation_fraction,dtype=float),
                closed_accommodation_used_fraction=np.full_like(selectivity,1-accommodation_fraction,dtype=float),
                closed_accommodation_recovery_number=np.asarray(k_closed,dtype=float),
                P_comp_closed=np.full_like(selectivity,accommodation_fraction,dtype=float),
                Gamma_mig=1/np.maximum(1+kg,1e-300),Lambda_TJ_over_K_TJ=np.full_like(selectivity,lambda_over_k,dtype=float),
                P_comp_TJ=np.full_like(selectivity,p_comp_tj,dtype=float),growth_activation_number=kg)


def longest_span(x,y,threshold):
    x=np.asarray(x,float);mask=np.asarray(y,float)>=threshold;best=0.;start=None
    for i,ok in enumerate(mask):
        if ok and start is None:start=i
        if start is not None and (not ok or i==len(mask)-1):
            end=i if ok and i==len(mask)-1 else i-1
            best=max(best,float(x[end]-x[start]));start=None
    return best


def artifact_reasons(*,attained,numerical_censored=False,hidden_density_channel=False,
                     bounded=True,negative_stores=False,timestep_collapse=False,
                     diagnostic_only=False,interpolation_supported=True):
    checks=((not attained,"unattained interval"),(numerical_censored,"numerical instability"),
            (hidden_density_channel,"hidden density channel"),(not bounded,"unbounded state variable"),
            (negative_stores,"negative stores"),(timestep_collapse,"timestep collapse"),
            (diagnostic_only,"diagnostic-only mechanism"),(not interpolation_supported,"unsupported interpolation"))
    return [reason for flag,reason in checks if flag]
