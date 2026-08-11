#!/usr/bin/env python3
"""Observable weighted-cohort initial microstructure layer.

Every cohort sees the same thermal protocol.  Heterogeneity is encoded only in
initial physical state and never in a schedule label.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import math
import numpy as np

import pr_desintering_memory_model as memory

INITIAL_MODES=("baseline_narrow","broad_pore_distribution","large_pore_tail",
 "bimodal_pore_distribution","broad_grain_distribution",
 "correlated_pore_grain_distribution","defect_rich_large_pore")


@dataclass(frozen=True)
class CohortSpec:
    weight: float; G_factor: float=1.; pore_factor: float=1.; pore_ln_sigma: float=.65
    location: tuple[float,float,float]=(.60,.30,.10); stress_sign: float=0.; defect: bool=False


@dataclass
class HeterogeneousParams:
    base: memory.PRMemoryParams
    initial_microstructure_mode: str="baseline_narrow"
    G0_mean_nm: float=75.; G0_ln_sigma: float=0.
    coarse_grain_tail_fraction: float=0.; coarse_grain_radius_factor: float=4.
    pore_ln_sigma: float=.65; large_pore_tail_fraction: float=.05
    large_pore_radius_factor: float=4.; bimodal_fine_fraction: float=.7
    pore_location_bias: str="baseline"; grain_pore_correlation: str="none"


def _locations(name):
    return {"baseline":(.60,.30,.10),"GBseg_fine_rich":(.80,.15,.05),
            "TJ_large_rich":(.35,.55,.10),"isolated_large_rich":(.35,.15,.50),
            "mixed":(.50,.45,.05)}[name]


def cohort_specs(p:HeterogeneousParams):
    if p.initial_microstructure_mode not in INITIAL_MODES:raise ValueError("invalid initial_microstructure_mode")
    loc=_locations(p.pore_location_bias);m=p.initial_microstructure_mode
    if m=="baseline_narrow":return [CohortSpec(1.,pore_ln_sigma=.65,location=loc)]
    if m=="broad_pore_distribution":return [CohortSpec(1.,pore_ln_sigma=p.pore_ln_sigma,location=loc)]
    if m=="large_pore_tail":
        f=p.large_pore_tail_fraction;return [CohortSpec(1-f,pore_ln_sigma=p.pore_ln_sigma,location=loc),CohortSpec(f,pore_factor=p.large_pore_radius_factor,pore_ln_sigma=.35,location=_locations("isolated_large_rich"),stress_sign=1)]
    if m=="bimodal_pore_distribution":
        f=p.bimodal_fine_fraction;return [CohortSpec(f,pore_factor=.65,pore_ln_sigma=.35,location=_locations("GBseg_fine_rich")),CohortSpec(1-f,pore_factor=p.large_pore_radius_factor,pore_ln_sigma=.45,location=loc)]
    if m=="broad_grain_distribution":
        s=max(p.G0_ln_sigma,.25);f=max(p.coarse_grain_tail_fraction,.10);return [CohortSpec((1-f)/2,G_factor=math.exp(-s),location=loc),CohortSpec((1-f)/2,G_factor=1.,location=loc),CohortSpec(f,G_factor=p.coarse_grain_radius_factor,location=loc)]
    if m=="correlated_pore_grain_distribution":
        corr=p.grain_pore_correlation
        if corr in ("large_pores_on_small_grains","TJ_pores_on_small_grains"):glarge,gsmall=1.5,.65
        else:glarge,gsmall=.75,1.35
        lloc=_locations("TJ_large_rich" if corr=="TJ_pores_on_small_grains" else ("isolated_large_rich" if corr=="isolated_pores_on_large_grains" else "mixed"))
        return [CohortSpec(.65,G_factor=gsmall,pore_factor=.75,pore_ln_sigma=.45,location=_locations("GBseg_fine_rich")),CohortSpec(.35,G_factor=glarge,pore_factor=p.large_pore_radius_factor,pore_ln_sigma=.65,location=lloc,stress_sign=1)]
    f=max(p.large_pore_tail_fraction,.05);return [CohortSpec(1-f,pore_ln_sigma=p.pore_ln_sigma,location=loc),CohortSpec(f,G_factor=1.5,pore_factor=p.large_pore_radius_factor,pore_ln_sigma=.95,location=_locations("isolated_large_rich"),stress_sign=1,defect=True)]


def cohort_params(p:HeterogeneousParams):
    out=[]
    for spec in cohort_specs(p):
        loc=p.base.base.action.location;base=loc.base
        b=replace(base,G0=p.G0_mean_nm*1e-9*spec.G_factor,pore_radius0=base.pore_radius0*spec.pore_factor,pore_ln_sigma=spec.pore_ln_sigma)
        lp=replace(loc,base=b,f_GBseg_init=spec.location[0],f_TJ_init=spec.location[1],f_iso_init=spec.location[2])
        bp=replace(p.base.base,action=replace(p.base.base.action,location=lp));out.append((spec,replace(p.base,base=bp)))
    return out


def weighted_quantile(values,weights,q):
    values=np.asarray(values,float);weights=np.asarray(weights,float);order=np.argsort(values);values=values[order];weights=weights[order]
    return float(values[np.searchsorted(np.cumsum(weights)/max(weights.sum(),1e-300),q,side="left")])


def aggregate_histories(items):
    """Interpolate cohort histories on common time and report full distributions."""
    end=max(h["t"][-1] for _,h in items);dt=max(min(np.median(np.diff(h["t"])) for _,h in items if len(h["t"])>1),1.)
    t=items[0][1]["t"].copy() if len(items)==1 else np.arange(0,end+dt*.25,dt)
    weights=np.array([s.weight for s,_ in items]);weights/=weights.sum();out={"t":t}
    scalar=["T_C","rho","G","connected_fine_pore_fraction","pore_mean_radius","large_pore_fraction","cumulative_PR_desintering_work"]
    scalar += [k for k in ("sigma_res_GBseg","sigma_res_TJ","sigma_res_large_pore","sigma_res_crack_like","residual_defect_flux") if all(k in h for _,h in items)]
    vals={k:np.vstack([np.interp(t,h["t"],h[k]) for _,h in items]) for k in scalar}
    for k in scalar:out[k]=np.average(vals[k],axis=0,weights=weights)
    out["G_mean"]=out["G"]
    out["G50"]=np.array([weighted_quantile(vals["G"][:,i],weights,.5) for i in range(len(t))])
    out["G90"]=np.array([weighted_quantile(vals["G"][:,i],weights,.9) for i in range(len(t))])
    # Cohort pore radii and pore-volume weights retain distribution information.
    d50=[];d90=[];number=[];iso=[];fg=[];ft=[]
    for i,ti in enumerate(t):
        radii=[];vol=[];nums=[];philoc=[]
        for w,(spec,h) in zip(weights,items):
            j=min(np.searchsorted(h["t"],ti),len(h["t"])-1);r=h["pore_radii"]
            stores=[h[k][j] for k in ("phi_GBseg","phi_TJ","phi_iso")];total=sum(stores)
            radii.extend(np.tile(r,3));vv=np.concatenate(stores)*w;vol.extend(vv);nums.extend(vv/np.maximum(4*np.pi*np.tile(r,3)**3/3,1e-300));philoc.append([np.sum(x)*w for x in stores])
        d50.append(weighted_quantile(radii,vol,.5));d90.append(weighted_quantile(radii,vol,.9));number.append(sum(nums));z=np.sum(philoc);fg.append(np.sum(philoc,axis=0)[0]/z);ft.append(np.sum(philoc,axis=0)[1]/z);iso.append(np.sum(philoc,axis=0)[2]/z)
    out.update(pore_D50=np.array(d50),pore_D90=np.array(d90),pore_number_proxy=np.array(number),f_GBseg=np.array(fg),f_TJ=np.array(ft),isolated_pore_fraction=np.array(iso),cohort_weights=weights)
    return out


def run(p:HeterogeneousParams,protocol,runner=memory.run):
    items=[(spec,runner(cp,protocol)) for spec,cp in cohort_params(p)]
    return aggregate_histories(items)


LOCAL_FUNCTIONS=(cohort_specs,cohort_params)
