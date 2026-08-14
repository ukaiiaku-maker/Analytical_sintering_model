#!/usr/bin/env python3
"""Compact analysis helpers for the staged high-density search."""
import pandas as pd
import numpy as np
import massive_latent_topology_models as model
import massive_latent_topology_objectives as objective

def best_candidate(table):
    d=pd.read_csv(table)
    if d.empty:return None
    return d.sort_values(['tier','median_reduction'],ascending=[True,False]).iloc[0].to_dict()

def count_tiers(table):
    d=pd.read_csv(table);return d.tier.value_counts(dropna=False).to_dict() if 'tier' in d else {}

def component_ablations(candidate_id,params):
    modes={'full':{},'no_connected_memory':{'connected_loss':0.,'connected_recovery':0.},
      'no_persistent_junction':{'XJ_prod':0.,'A_J':0.},'no_stress_memory':{'stress_prod':0.,'stress_coupling':0.},
      'no_multihit':{'lambda_TJ':1e6,'K_TJ':1.},'no_closed_pore':{'k_closed':0.},
      'no_detachment_closure':{'closure_rate':0.,'detachment_rate':0.},
      'topology_disabled':{'connected_loss':0.,'connected_recovery':0.,'XJ_prod':0.,'A_J':0.,'stress_prod':0.,'stress_coupling':0.,'pore_drag':0.}}
    rows=[]
    for name,change in modes.items():
      p={**params,**change};h=model.simulate(p,T1=1400,T2=None,dt=1800);t=model.simulate(p,T1=1400,T2=1200,dt=1800);rows.append({'candidate_id':candidate_id,'ablation':name,**objective.trajectory_score(h,t),'rho_high_final':h['rho'][-1],'rho_two_final':t['rho'][-1]})
    return rows
