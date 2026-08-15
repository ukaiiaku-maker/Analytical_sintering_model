from __future__ import annotations
import math
import numpy as np


def kinetic_state(sigma_Pa: float, T_K: float, G_m: float, connectivity: float,
                  barrier, material, C1: float = 1.0, site_density_multiplier: float = 1.0,
                  activity_multiplier: float = 1.0) -> dict[str, float]:
    geo = material.triple_line_geometry(G_m)
    Gstar = barrier.Gstar(sigma_Pa, T_K)
    rate = material.nu0_sinv * math.exp(-Gstar/(material.kB*T_K))
    tau = C1*(material.kB*T_K/max(sigma_Pa*material.Omega_m3, 1e-300))*G_m**2/max(material.D_GB(T_K), 1e-300)
    Lambda = rate*tau
    activity = np.clip(activity_multiplier*Lambda/(1+Lambda),0,1)
    edot = site_density_multiplier*np.clip(connectivity, 0, 1)*geo["eps_event"]*activity/max(tau, 1e-300)
    return {"Gstar_J": Gstar, "r_nuc_sinv": rate, "tau_sink_s": tau,
            "Lambda": Lambda, "activity": activity, "edot_sinv": edot,
            "effective_eps_event": site_density_multiplier*geo["eps_event"], **geo}


def density_rate(rho: float, edot: float) -> float:
    return max(0., (1-rho)*edot)


def connectivity(phi_open: float, phi_all: float) -> float:
    return float(np.clip(phi_open/max(phi_all, 1e-300), 0, 1))
