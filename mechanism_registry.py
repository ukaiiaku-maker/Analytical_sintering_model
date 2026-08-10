#!/usr/bin/env python3
"""Source-grounded registry for sintering and grain-growth mechanisms."""
from dataclasses import dataclass
from typing import Mapping, Tuple


@dataclass(frozen=True)
class MechanismSpec:
    name: str
    mechanism_class: str
    affected_process: str
    state_variables: Tuple[str, ...]
    parameter_bounds: Mapping[str, Tuple[float, float]]
    conservation: str
    changes_density_directly: bool
    migration_only: bool
    power_channels: Tuple[str, ...]
    stress_diagnostics: Tuple[str, ...]
    source_rationale: str
    rejection_criteria: Tuple[str, ...]


COMMON_REJECTIONS=("conservation_violation","negative_store","universal_success",
                   "missing_lower_boundary","missing_upper_boundary","target_specific_tuning")


def _spec(name,cls,process,state,bounds,conservation,density,migration,power=(),stress=(),rationale=""):
    return MechanismSpec(name,cls,process,tuple(state),bounds,conservation,density,migration,
                         tuple(power),tuple(stress),rationale,COMMON_REJECTIONS)


REGISTRY={s.name:s for s in (
 _spec("GBseg_direct_densification","mixed","densification",("phi_GBseg","C_GBseg","renewal_activity"),{},"removes GBseg pore volume consistently with rho",True,False,("P_GBseg_dens",),rationale="Pore-connected boundaries alone are geometrically eligible for pore removal."),
 _spec("TJ_assisted_densification","mixed","densification",("phi_TJ","C_TJ","renewal_activity"),{"eta_TJ_dens":(0,1)},"removes TJ pore volume consistently with rho",True,False,("P_TJ_dens",),rationale="TJ-connected sinks may assist interface-normal strain but have distinct eligibility."),
 _spec("clean_GB_migration","mixed","GB migration",("f_clean_GB","G","T"),{},"no pore-volume change",False,True,("P_clean_GB",),rationale="Pore-disconnected boundaries may migrate but ordinarily do not densify."),
 _spec("GBseg_pore_drag","A_drag","GB migration",("C_GBseg","G","T"),{},"no pore-volume change",False,True,("P_GBseg_drag",),rationale="Continuous pore drag is a Class-A series resistance."),
 _spec("TJ_drag_series","A_drag","TJ migration",("C_TJ","G","T"),{},"no pore-volume change",False,True,("P_TJ_drag",),rationale="TJ drag acts in series with intrinsic GB mobility."),
 _spec("persistent_junction_drag","A_drag","GB migration",("X_J","C_TJ","f_clean_GB","G","T"),{"A_J":(1,80),"tau_J_ref_s":(1e3,2e6),"Q_relax_J_mol":(250e3,550e3)},"X_J bounded; no pore-volume change",False,True,("P_persistent_junction_drag",),rationale="Junction obstacles can persist after connected pore coverage is lost."),
 _spec("TJ_multihit_reaction","B_multihit","TJ migration",("Lambda_TJ","K_TJ","G","T","C_TJ"),{"lambda_ref":(.05,50),"K0":(1,12),"q_TJ":(0,1)},"no pore-volume change",False,True,("P_TJ_multihit",),rationale="TJ reactions are Class-B packet-completion events; fixed-packet and accommodation-demand limits are both unresolved."),
 _spec("vacancy_accommodation_multihit","B_multihit","GB migration",("Lambda_vac","K_vac","G","T"),{"q_vac":(1,2)},"no pore-volume change",False,True,("P_vacancy_accommodation",),rationale="Vacancy accommodation can require a multihit conserved-volume packet."),
 _spec("GBseg_to_TJ_relocation","mixed","pore relocation",("phi_GBseg","C_TJ","T"),{},"conservative",False,False,("P_GB_to_TJ_relocation",),rationale="Observable conservative relocation among connected stores."),
 _spec("TJ_to_GBseg_capture","mixed","pore relocation",("phi_TJ","f_clean_GB","T"),{},"conservative",False,False,("P_TJ_to_GBseg_capture",),rationale="Migrating boundaries may recapture TJ pores."),
 _spec("TJ_to_isolated_conversion","mixed","pore relocation",("phi_TJ","rho","C_TJ"),{},"conservative",False,False,("P_TJ_iso",),rationale="Connected TJ pores may become isolated late in open-pore sintering."),
 _spec("boundary_sweep_repopulation","mixed","pore relocation",("f_clean_GB","C_TJ","G_dot"),{},"conservative",False,False,("P_boundary_sweep",),rationale="Boundary sweep can repopulate pore-bearing segments."),
 _spec("stress_accumulation_release","B_multihit","stress activation",("sigma_accum_GBseg","sigma_accum_TJ","G_dot","rho_dot"),{},"no pore-volume change",False,False,("P_stress_storage","P_stress_release"),("sigma_accum_GBseg","sigma_accum_TJ"),"Coarsening may store local activation stress while densification is arrested."),
 _spec("closed_pore_placeholder","C_exchange","late-stage pore removal",("phi_iso","rho","tau_exchange"),{},"not implemented",True,False,rationale="Explicit closed-pore removal is reserved for densities beyond the open-pore model."),
)}


def candidates_for(process): return tuple(s for s in REGISTRY.values() if process in s.affected_process)

