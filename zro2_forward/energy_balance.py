from __future__ import annotations
from dataclasses import dataclass


@dataclass
class PowerBalance:
    sigma_eff_Pa: float
    P_surf_W_m3: float
    P_dens_W_m3: float
    P_excess_W_m3: float
    efficiency: float
    stress_bound_hit: bool


def solve_effective_stress(rho: float, surface_area_rate_m2_m3_s: float, gamma_s: float,
                           rate_at_stress, zeta: float = 1., eta_geom: float = 1.,
                           sigma_min: float = 1e3, sigma_max: float = 2.5e8) -> PowerBalance:
    ps = zeta*eta_geom*gamma_s*abs(surface_area_rate_m2_m3_s)
    def pd(s): return (1-rho)*s*max(rate_at_stress(s), 0.)
    lo, hi = sigma_min, sigma_max
    if pd(lo) >= ps: sigma, hit = lo, True
    elif pd(hi) < ps: sigma, hit = hi, True
    else:
        for _ in range(32):
            mid = .5*(lo+hi)
            if pd(mid) < ps: lo = mid
            else: hi = mid
        sigma, hit = .5*(lo+hi), False
    power = pd(sigma); excess = max(0., ps-power)
    return PowerBalance(sigma, ps, power, excess, power/(ps+1e-300), hit)
