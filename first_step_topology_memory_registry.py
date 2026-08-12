#!/usr/bin/env python3
"""Auditable registry of local first-step topology-memory hypotheses."""
from dataclasses import dataclass,asdict

@dataclass(frozen=True)
class MemoryMechanism:
    mechanism_id:str;family:str;physical_rationale:str;state_variables:str;evolution_law:str
    changes_density:bool;changes_migration_only:bool;conservative_pore_flux:bool
    observable_signatures:str;parameters:str;rejection_criteria:str

def registry():
    common="no matched-density divergence; no second-step persistence; changes shared-state rho_dot"
    return [
      MemoryMechanism("A_PR","PR_topology_generator","low-renewal exposure redistributes pores","PR_work,D50,D90,f_fine","dW=PR_propensity dt; conservative adjacent-bin transfer",False,False,True,"D90 increase; fine-pore loss","rate,Q_PR,activity_gate,partition,tau_memory",common),
      MemoryMechanism("B_REMOVE","connected_removable_inventory","track percolating removable pore area","C_remove_GB,C_remove_TJ,f_fine","loss by smoothing/isolation; recovery by local renewal",False,False,True,"connected-fine and coverage differences","loss,recovery,percolation threshold",common),
      MemoryMechanism("C_JUNCTION","persistent_junction_segment","first-step boundary sweep stores junction constraints","X_J,junction_density,segment_length,C_constraint","production by PR/TJ/sweep; thermal/activity relaxation",False,True,False,"X_J and migration-factor difference","production,relaxation,capacity,A_J",common),
      MemoryMechanism("D_STRESS","residual_stress_memory","constrained migration stores shear/back stress","sigma_res,stored_shear_work","generation by growth/PR work; Arrhenius relaxation",False,True,False,"restart migration and residual-stress difference","generation,Q_relax,coupling",common),
      MemoryMechanism("E_TJ_PARTITION","pore_TJ_partition","separate pore relaxation from structural constraint","C_TJ_pore,C_constraint,C_pinned,C_relaxed","instantaneous partition plus persistent structural state",False,True,False,"pore drag versus multihit separation","relax_fraction,drag_fraction,A_drag",common),
      MemoryMechanism("F_VAC","vacancy_multihit_memory","store accommodation packet demand","Lambda_vac,K_vac,mismatch","packet production/relaxation; Poisson completion",False,True,False,"Lambda/K and completion difference","lambda_ref,K0,q_vac,Q_event",common),
    ]

def rows():return [asdict(x) for x in registry()]
