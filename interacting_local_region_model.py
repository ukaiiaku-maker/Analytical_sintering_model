#!/usr/bin/env python3
"""Minimal interacting local-region pore network with named local fluxes."""
from dataclasses import dataclass
import numpy as np
from massive_latent_topology_models import R,sigmoid

@dataclass
class NetworkState:
 rho:np.ndarray;G:np.ndarray;phi_GBseg:np.ndarray;phi_TJ:np.ndarray;phi_iso:np.ndarray;phi_closed:np.ndarray;N_GBseg:np.ndarray;N_TJ:np.ndarray;N_iso:np.ndarray;N_closed:np.ndarray;connected_removable_fraction:np.ndarray;damaged_connected_fraction:np.ndarray;sweep_coalescence_seed:np.ndarray;large_attached_fraction:np.ndarray;large_TJ_fraction:np.ndarray;isolated_fraction:np.ndarray;closed_fraction:np.ndarray;X_J:np.ndarray;C_TJ:np.ndarray;C_GBseg:np.ndarray;f_clean_GB:np.ndarray;residual_stress:np.ndarray;PR_damage_memory:np.ndarray;sweep_memory:np.ndarray;closed_accommodation:np.ndarray;migration_factor:np.ndarray;densification_eligibility:np.ndarray;weights:np.ndarray
def initial(n=8,rho0=.70,G0=100.,seed=1):
 rng=np.random.default_rng(seed);w=rng.lognormal(0,.25,n);w/=w.sum();p=(1-rho0)*np.clip(rng.normal(1,.12,n),.7,1.3);p*=((1-rho0)/(w@p));gb=.6*p;tj=.2*p;iso=.2*p;z=lambda:np.zeros(n);one=lambda:np.ones(n);return NetworkState(1-p,np.full(n,G0),gb,tj,iso,z(),gb.copy(),tj.copy(),iso.copy(),z(),.8*one(),.05*one(),z(),z(),z(),iso/p,z(),.1*one(),.2*one(),.6*one(),.2*one(),z(),z(),z(),one(),one(),one(),w)
def global_observables(s):
 pore=s.phi_GBseg+s.phi_TJ+s.phi_iso+s.phi_closed;return dict(rho_global=1-float(s.weights@pore),G_mean=float(s.weights@s.G),topology_variance=float(np.average((s.connected_removable_fraction-np.average(s.connected_removable_fraction,weights=s.weights))**2,weights=s.weights)),closed_fraction=float(s.weights@s.phi_closed/max(s.weights@pore,1e-30)))
def local_fluxes(s,T_C,p):
 T=T_C+273.15;th=lambda Q:np.exp(-Q/R*(1/T-1/1473.15));activity=sigmoid((T_C-p['activity_mid'])/p['activity_width']);open_flux=p['k_open']*th(p['Q_density'])*s.phi_GBseg*s.connected_removable_fraction*s.densification_eligibility;closed_flux=p['k_closed']*th(p['Q_closed'])*s.phi_closed*s.closed_accommodation*(1-p['gas_ratio']);pr=p['k_PR']*th(p['Q_PR'])*(1-activity)**2*s.phi_GBseg;drag=p['attached_drag']*(s.large_attached_fraction+s.large_TJ_fraction)+p['junction_drag']*s.X_J+p['stress_drag']*s.residual_stress;mf=1/(1+drag);growth=p['k_growth']*th(p['Q_growth'])*mf/np.maximum(s.G,1);sweep=p['k_sweep']*(growth/np.maximum(s.G,1))**p['sweep_exp']*(s.damaged_connected_fraction+s.sweep_coalescence_seed);return dict(rho_dot_open=open_flux,rho_dot_closed=closed_flux,PR_damage=pr,sweep=sweep,G_dot=growth,migration_factor=mf,activity=activity)
def advance(s,T,p,dt,adj=None):
 f=local_fluxes(s,T,p);pr=np.minimum(f['PR_damage']*dt,s.phi_GBseg);sw=np.minimum(f['sweep']*dt,s.phi_GBseg-pr);op=np.minimum(f['rho_dot_open']*dt,s.phi_GBseg-pr-sw);cl=np.minimum(f['rho_dot_closed']*dt,s.phi_closed);s.phi_GBseg-=pr+sw+op;s.phi_TJ+=.2*pr+.2*sw;s.phi_iso+=.3*pr+.3*sw;s.phi_closed+=.1*pr+.1*sw-cl;s.phi_GBseg+=.4*pr+.4*sw;s.N_GBseg=np.maximum(s.N_GBseg-pr-p['number_loss']*sw,0);s.N_TJ+=.05*sw;s.N_iso+=.05*sw;s.N_closed=np.maximum(s.N_closed+.1*sw-cl,0);s.damaged_connected_fraction+=pr-sw;s.sweep_coalescence_seed+=.3*pr-.3*sw;s.large_attached_fraction+=.4*sw;s.large_TJ_fraction+=.2*sw;s.PR_damage_memory+=pr;s.sweep_memory+=sw;s.closed_accommodation=np.maximum(s.closed_accommodation-.1*cl,0);s.G+=f['G_dot']*dt;s.rho=1-(s.phi_GBseg+s.phi_TJ+s.phi_iso+s.phi_closed);s.migration_factor=f['migration_factor'];s.connected_removable_fraction=s.phi_GBseg/np.maximum(1-s.rho,1e-30)
 if adj is not None:
  mean=adj@s.connected_removable_fraction/np.maximum(adj.sum(axis=1),1);s.connected_removable_fraction+=p['exchange_rate']*(mean-s.connected_removable_fraction)*dt
 return f
def defaults():return dict(k_open=1.5e-5,Q_density=475e3,k_closed=2e-6,Q_closed=475e3,k_PR=1e-5,Q_PR=250e3,k_sweep=20.,sweep_exp=1.,k_growth=9e3,Q_growth=500e3,activity_mid=1180.,activity_width=70.,gas_ratio=.25,attached_drag=30.,junction_drag=10.,stress_drag=1.,number_loss=2.,exchange_rate=1e-7)
LOCAL_FUNCTIONS=(local_fluxes,)
