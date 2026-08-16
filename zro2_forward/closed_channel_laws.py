"""Physically interpretable closed-pore law candidates.

The functions here are local, path-independent state laws.  They return a removal
rate for each closed-pore bin and transparent factors used to construct it.
None of the candidate constants is calibrated or validated for ZrO2.
"""
from __future__ import annotations

import numpy as np

from .material_zro2 import R

ALLOWED_CLOSED_CHANNEL_LAWS = (
    "resolved_proxy_current",
    "GB_diffusion_closed_shrinkage",
    "renewal_limited_closed_shrinkage",
    "surface_diffusion_accommodation_only",
    "gas_accommodation_limited",
    "empirical_closed_rate_scale",
)


def closed_channel_rates(state, T_K, barrier, material, q, renewal):
    """Return ``(bin removal rates, diagnostics)`` for the selected law."""
    law=q.closed_channel_law
    if law not in ALLOWED_CLOSED_CHANNEL_LAWS:
        raise ValueError(f"unknown closed_channel_law: {law}")
    p=state.pores; r=np.maximum(p.radii_m,1e-30); phi=p.phi_closed
    A=float(np.clip(state.A_closed,0,q.accommodation_max))
    activity=float(np.clip(renewal,0,1))
    capillary=2*material.gamma_s_J_m2/r
    gas=np.clip(q.closed_gas_pressure_fraction,0,1)*capillary
    gas_factor=np.clip(1-gas/np.maximum(capillary,1e-300),0,1)
    pressure=np.maximum(capillary-gas,0)
    Dgb=material.D_GB(T_K); Ds=material.D_s(T_K)
    sigma=q.C_sigma_closed*pressure
    Gstar=np.array([barrier.Gstar(float(s),T_K) for s in sigma])
    r_nuc=material.nu0_sinv*np.exp(-Gstar/(material.kB*T_K))
    tau_nuc=1/np.maximum(r_nuc,1e-300)
    tau_exchange=np.full_like(r,q.closed_exchange_tau_s)
    L=q.closed_transport_length_factor*r
    tau_transport=q.C_transport_closed*(material.kB*T_K/np.maximum(sigma*material.Omega_m3,1e-300))*L**2/max(Dgb,1e-300)
    tau_cycle=tau_nuc+tau_exchange+tau_transport
    mexp=int(q.closed_radius_exponent)
    Cgb=q.C_closed_GB*q.closed_prefactor_factor
    gb_kernel=Cgb*pressure*Dgb*material.Omega_m3/(material.kB*T_K)*r**(-mexp)
    if law == "resolved_proxy_current":
        closed_activity=activity*np.exp(np.clip(-.35*q.Q_PR_J_mol/R*(1/T_K-1/q.T_PR_ref_K),-30,30))
        tau=q.closed_tau0_s*(r/25e-9)**4/max(closed_activity*A,1e-12)
        shrink=phi/tau
        transport=1/tau
        activity_used=closed_activity
    elif law == "GB_diffusion_closed_shrinkage":
        activity_used=activity if q.closed_GB_use_renewal_activity else 1.
        shrink=phi*A*activity_used*gb_kernel
        transport=gb_kernel
    elif law == "renewal_limited_closed_shrinkage":
        activity_used=(tau_exchange+tau_transport)/np.maximum(tau_cycle,1e-300)
        shrink=phi*A*q.closed_event_strain/np.maximum(tau_cycle,1e-300)
        transport=1/np.maximum(tau_transport,1e-300)
    elif law == "surface_diffusion_accommodation_only":
        activity_used=activity; shrink=np.zeros_like(phi); transport=np.zeros_like(phi)
    elif law == "gas_accommodation_limited":
        activity_used=activity
        shrink=phi*A*gas_factor*gb_kernel
        transport=gb_kernel
    else:
        activity_used=1.
        size=(25e-9/r)**q.closed_empirical_size_exponent
        k=q.k0_closed_emp_sinv*q.closed_prefactor_factor*np.exp(-q.Q_closed_emp_J_mol/(R*T_K))
        shrink=phi*A*size*k
        transport=np.full_like(phi,k)
    shrink=np.minimum(np.maximum(shrink,0),phi/max(q.closed_rate_cap_time_s,1e-30))
    accommodation_rate=(q.C_surface_accommodation*Ds/r**4*(q.accommodation_max-A)) if law=="surface_diffusion_accommodation_only" else np.zeros_like(r)
    diag={
        "closed_channel_law":law,"closed_inventory":float(phi.sum()),"shrinkable_closed_inventory":float(phi.sum()),
        "closed_pore_radius_mean_m":float(np.sum(phi*r)/max(phi.sum(),1e-300)),"accommodation_factor":A,
        "closed_activity_factor":float(activity_used) if np.isscalar(activity_used) else float(np.mean(activity_used)),
        "closed_pressure_factor_mean_Pa":float(np.mean(pressure)),"closed_gas_factor_mean":float(np.mean(gas_factor)),
        "closed_transport_factor_mean_sinv":float(np.mean(transport)),"closed_dimensional_prefactor":Cgb if law in ("GB_diffusion_closed_shrinkage","gas_accommodation_limited") else q.k0_closed_emp_sinv if law=="empirical_closed_rate_scale" else q.closed_event_strain if law=="renewal_limited_closed_shrinkage" else q.closed_tau0_s,
        "closed_dimensional_prefactor_units":f"m^{mexp-3}" if law in ("GB_diffusion_closed_shrinkage","gas_accommodation_limited") else "s^-1" if law=="empirical_closed_rate_scale" else "dimensionless event strain" if law=="renewal_limited_closed_shrinkage" else "s",
        "D_GB_closed_m2_s":Dgb,"D_surface_closed_m2_s":Ds,"F_pressure_closed_mean_Pa":float(np.mean(pressure)),
        "F_activity_closed":float(activity_used) if np.isscalar(activity_used) else float(np.mean(activity_used)),
        "closed_radius_exponent":mexp,"C_closed_GB":Cgb,"rho_dot_closed_law_sinv":float(shrink.sum()),
        "sigma_closed_mean_Pa":float(np.mean(sigma)),"Gstar_closed_mean_J":float(np.mean(Gstar)),
        "r_nuc_closed_mean_sinv":float(np.mean(r_nuc)),"tau_nuc_closed_mean_s":float(np.mean(tau_nuc)),
        "tau_exchange_closed_mean_s":float(np.mean(tau_exchange)),"tau_transport_closed_mean_s":float(np.mean(tau_transport)),
        "tau_cycle_closed_mean_s":float(np.mean(tau_cycle)),"Lambda_closed_mean":float(np.mean(r_nuc*(tau_exchange+tau_transport))),
        "A_dot_closed_sinv":float(np.sum(accommodation_rate*phi)/max(phi.sum(),1e-300)) if phi.sum()>0 else 0.,
        "empirical_diagnostic":law=="empirical_closed_rate_scale","Q_closed_emp_not_material_property":law=="empirical_closed_rate_scale",
    }
    return shrink,diag
