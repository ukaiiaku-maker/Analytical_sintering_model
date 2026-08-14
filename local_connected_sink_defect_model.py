#!/usr/bin/env python3
"""Local connected-matrix/defect-rich mixture with independent evolution."""
from __future__ import annotations
from dataclasses import dataclass,replace
import numpy as np

import heterogeneous_initial_state_model as hetero
import residual_stress_memory_model as residual
import persistent_defect_topology_stress_model as persistent

MODES=("disabled_local_mixture","matrix_only","static_defect_mixture",
       "evolving_defect_memory","evolving_defect_memory_with_matrix_densification",
       "stress_retention_high")


@dataclass
class LocalMixtureParams:
    base:object;mode:str="disabled_local_mixture";G0_nm:float=150.;rho0:float=.65
    defect_weight:float=.1;defect_large_pore_factor:float=8.
    matrix_connected_fraction:float=.7;stress_scale:float=1.
    tau_defect_over_matrix:float=10.;pore_ln_sigma:float=.95


def normalized_specs(p):
    if p.mode not in MODES:raise ValueError("invalid local mixture mode")
    if p.mode=="matrix_only":return [("matrix_connected",hetero.CohortSpec(1.,G_factor=1.,pore_factor=.7,pore_ln_sigma=.45,location=(p.matrix_connected_fraction,1-p.matrix_connected_fraction-.05,.05)))]
    wd=float(np.clip(p.defect_weight,0,1));wm=1-wd
    return [("matrix_connected",hetero.CohortSpec(wm,G_factor=.9,pore_factor=.7,pore_ln_sigma=.45,location=(p.matrix_connected_fraction,max(1-p.matrix_connected_fraction-.05,0),.05))),
            ("defect_rich_large_pore",hetero.CohortSpec(wd,G_factor=1.5,pore_factor=p.defect_large_pore_factor,pore_ln_sigma=p.pore_ln_sigma,location=(.25,.25,.50),stress_sign=1,defect=True))]


def cohort_param(p,spec):
    loc=p.base.base.action.location;base=replace(loc.base,rho0=p.rho0,G0=p.G0_nm*1e-9*spec.G_factor,pore_radius0=loc.base.pore_radius0*spec.pore_factor,pore_ln_sigma=spec.pore_ln_sigma)
    lp=replace(loc,base=base,f_GBseg_init=spec.location[0],f_TJ_init=spec.location[1],f_iso_init=spec.location[2]);return replace(p.base,base=replace(p.base.base,action=replace(p.base.base.action,location=lp)))


def local_modes(p,name):
    defect=name=="defect_rich_large_pore"
    # Persistent memory is never applied to the matrix in this local model.
    # The disabled ablation recovers the heterogeneous (non-persistent)
    # baseline, which is one of the explicitly permitted parent controls.
    pmode="persistent_defect_memory" if defect and p.mode in ("evolving_defect_memory","evolving_defect_memory_with_matrix_densification","stress_retention_high") else "disabled"
    if p.mode in ("disabled_local_mixture","matrix_only"):rmode="disabled"
    elif p.mode=="static_defect_mixture":rmode="initial_only"
    elif p.mode=="evolving_defect_memory_with_matrix_densification" and not defect:rmode="disabled"
    else:rmode="mixed_evolving"
    rp=residual.ResidualStressParams(mode=rmode,sigma_res_scale=p.stress_scale if defect else .25*p.stress_scale,stress_sign="tensile" if defect else "compressive",tau_res_ref_s=2e4*(p.tau_defect_over_matrix if defect else 1))
    tau=8e5*(p.tau_defect_over_matrix if defect else 1);pp=persistent.PersistentParams(mode=pmode,defect_decay_time_s=tau)
    return rp,pp


def run(p,protocol):
    items=[];locals_={}
    for name,spec in normalized_specs(p):
        cp=cohort_param(p,spec);rp,pp=local_modes(p,name);h=persistent.run(cp,protocol,rp,pp,spec.stress_sign);items.append((spec,h));locals_[name]=h
    aggregate=hetero.aggregate_histories(items);aggregate["numerical_censored"]=any(bool(h.get("numerical_censored",False)) for _,h in items);return aggregate,locals_


def global_density_identity(p,locals_,index=0):
    specs=normalized_specs(p);w=np.array([s.weight for _,s in specs]);w/=w.sum();return float(sum(w[i]*locals_[name]["rho"][min(index,len(locals_[name]["rho"])-1)] for i,(name,_) in enumerate(specs)))


LOCAL_FUNCTIONS=(normalized_specs,local_modes)
