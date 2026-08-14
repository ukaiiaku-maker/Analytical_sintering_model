#!/usr/bin/env python3
"""Reduced, state-resolved high-density topology model for staged searches.

The implementation deliberately separates conservative pore-store transfers from
the two density-changing fluxes: open-pore and closed-pore shrinkage.
"""
from dataclasses import dataclass
import numpy as np

R = 8.31446261815324

@dataclass
class State:
    rho: float; G_nm: float; phi_open: float; phi_connected: float
    phi_isolated: float; phi_closed: float; XJ: float; stress: float

def sigmoid(x): return 1.0/(1.0+np.exp(-np.clip(x,-60,60)))

def rates(s, T_C, p):
    """Instantaneous local evolution laws."""
    T=T_C+273.15
    close_gate=sigmoid((s.rho-p['rho_close_mid'])/p['rho_close_width'])
    mobility=np.exp(-p['Q_growth']/R*(1/T-1/1623.15))
    diff=np.exp(-p['Q_density']/R*(1/T-1/1623.15))
    percol=sigmoid((s.phi_connected/(max(1-s.rho,1e-12))-p['percolation_threshold'])*p['percolation_exponent'])
    open_shrink=p['k_open']*diff*s.phi_open*percol
    closure=p['closure_rate']*close_gate*s.phi_isolated
    detach=p['detachment_rate']*close_gate*s.phi_connected
    gas_drive=max(1-p['gas_pressure_ratio'],0)
    closed_shrink=p['k_closed']*diff*s.phi_closed*gas_drive
    loss=p['connected_loss']*mobility*s.phi_connected
    recover=p['connected_recovery']*diff*s.phi_open*(1-percol)
    xprod=p['XJ_prod']*(loss+detach); xrel=s.XJ/max(p['tau_J'],1)
    stress_prod=p['stress_prod']*mobility*(s.XJ+loss); stress_rel=s.stress/max(p['tau_stress'],1)
    constraint=s.XJ*(1-p['pore_relax_fraction']*s.phi_connected/max(1-s.rho,1e-12))
    drag=p['A_J']*constraint+p['pore_drag']*s.phi_connected/max(1-s.rho,1e-12)+p['stress_coupling']*s.stress
    completion=p['lambda_TJ']/max(p['K_TJ']*(s.G_nm/150)**p['q_TJ'],1e-12)
    migration_factor=(1-np.exp(-completion))/(1+max(drag,0))
    growth=p['k_growth']*mobility*max(migration_factor,0)/max(s.G_nm,1)
    return dict(open_shrink=open_shrink,closed_shrink=closed_shrink,closure=closure,
                detach=detach,loss=loss,recover=recover,XJ_dot=xprod-xrel,
                stress_dot=stress_prod-stress_rel,G_dot=growth,
                migration_factor=migration_factor,percolation=percol)

def step(s,T_C,p,dt):
    r=rates(s,T_C,p)
    # Conservative transfers are capped by their source stores.
    loss=min(r['loss']*dt,s.phi_connected); recover=min(r['recover']*dt,s.phi_open)
    detach=min(r['detach']*dt,s.phi_connected-loss)
    closure=min(r['closure']*dt,s.phi_isolated)
    do=min(r['open_shrink']*dt,s.phi_open+recover)
    dc=min(r['closed_shrink']*dt,s.phi_closed+closure+detach)
    s.phi_open=max(s.phi_open-recover-do,0)
    s.phi_connected=max(s.phi_connected-loss+recover-detach,0)
    s.phi_isolated=max(s.phi_isolated+loss-closure,0)
    s.phi_closed=max(s.phi_closed+closure+detach-dc,0)
    s.rho=1-(s.phi_open+s.phi_connected+s.phi_isolated+s.phi_closed)
    s.G_nm=max(s.G_nm+r['G_dot']*dt,1e-9)
    s.XJ=float(np.clip(s.XJ+r['XJ_dot']*dt,0,p['XJ_capacity']))
    s.stress=float(np.clip(s.stress+r['stress_dot']*dt,0,p['stress_cap']))
    return r

def initial_state(rho0,G0,connected=.6,isolated=.1):
    pore=1-rho0
    return State(rho0,G0,pore*(1-connected-isolated),pore*connected,pore*isolated,0.,0.,0.)

def default_parameters():
    return dict(rho_close_mid=.92,rho_close_width=.02,Q_growth=500e3,Q_density=475e3,
        k_open=1.5e-5,k_closed=2e-6,k_growth=9e3,closure_rate=2e-5,
        detachment_rate=2e-6,gas_pressure_ratio=.25,percolation_threshold=.2,
        percolation_exponent=5.,connected_loss=2e-6,connected_recovery=2e-7,
        XJ_prod=8.,tau_J=2e6,tau_stress=2e6,stress_prod=.2,A_J=15.,
        pore_drag=4.,pore_relax_fraction=.5,stress_coupling=2.,lambda_TJ=2.,K_TJ=3.,q_TJ=1,
        XJ_capacity=2.,stress_cap=5.)

LOCAL_FUNCTIONS=(rates,)

def simulate(p,rho0=0.70,G0=100.,T1=1400.,T2=None,rho_switch=.88,budget_h=500.,dt=1800.):
    """Integrate a high-T or two-step path with exact state transfer."""
    s=initial_state(rho0,G0,p.get('initial_connected',.6),p.get('initial_isolated',.1))
    hist={k:[] for k in ('t','T_C','rho','G_nm','phi_open','phi_connected','phi_isolated','phi_closed','XJ','stress','rho_dot_open','rho_dot_closed','G_dot','migration_factor')}
    switched=False;t=0.;limit=budget_h*3600
    while t<=limit and s.rho<.995:
        T=T2 if switched and T2 is not None else T1
        r=rates(s,T,p)
        for k in ('t','T_C','rho','G_nm','phi_open','phi_connected','phi_isolated','phi_closed','XJ','stress'):hist[k].append(t if k=='t' else T if k=='T_C' else getattr(s,k))
        hist['rho_dot_open'].append(r['open_shrink']);hist['rho_dot_closed'].append(r['closed_shrink']);hist['G_dot'].append(r['G_dot']);hist['migration_factor'].append(r['migration_factor'])
        if T2 is not None and not switched and s.rho>=rho_switch:switched=True
        step(s,T,p,min(dt,limit-t));t+=dt
        if t>=limit:break
    out={k:np.asarray(v,float) for k,v in hist.items()};out['switched']=switched;out['final_state']=s
    return out
