from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np


@dataclass
class PorePopulation:
    radii_m: np.ndarray
    phi_open: np.ndarray
    phi_iso: np.ndarray
    phi_closed: np.ndarray

    def copy(self):
        return PorePopulation(*(x.copy() for x in (self.radii_m, self.phi_open, self.phi_iso, self.phi_closed)))

    @property
    def total(self): return float(np.sum(self.phi_open + self.phi_iso + self.phi_closed))
    @property
    def number_open(self): return self.phi_open / ((4/3)*math.pi*self.radii_m**3)


def initial_population(rho0: float = 2.88/5.95, center_m: float = 24.5e-9,
                       bins: int = 12, ln_sigma: float = 0.38) -> PorePopulation:
    radii = np.geomspace(center_m/3, center_m*4, bins)
    weights = np.exp(-0.5*(np.log(radii/center_m)/ln_sigma)**2)
    phi = (1-rho0)*weights/weights.sum()
    z = np.zeros_like(phi)
    return PorePopulation(radii, phi, z.copy(), z.copy())


def transfer_fluxes(pop: PorePopulation, rho: float, D_s: float, activity: float,
                    excess_fraction: float, C_PR: float = 1e-23, C_iso: float = 2e-25,
                    C_close: float = 1e-25) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    r = pop.radii_m
    open_dot = np.zeros_like(r); iso_dot = np.zeros_like(r); closed_dot = np.zeros_like(r)
    rates = C_PR*D_s/r**4*max(excess_fraction, 0.)*(1-np.clip(activity, 0, 1))
    crossings = rates[:-1]*pop.phi_open[:-1]
    open_dot[:-1] -= crossings; open_dot[1:] += crossings
    large = (r/r[-1])**2
    gi = 1/(1+np.exp(-(rho-.82)/.04))*large
    gc = 1/(1+np.exp(-(rho-.90)/.03))
    ji = C_iso*D_s/r**4*gi*pop.phi_open
    jc_o = C_close*D_s/r**4*gc*pop.phi_open
    jc_i = C_close*D_s/r**4*gc*pop.phi_iso
    open_dot -= ji+jc_o; iso_dot += ji-jc_i; closed_dot += jc_o+jc_i
    return open_dot, iso_dot, closed_dot, {"bin_crossing_rate": float(crossings.sum()), "isolation_rate": float(ji.sum()), "closure_rate": float((jc_o+jc_i).sum())}


def removal_weights(pop: PorePopulation) -> np.ndarray:
    x = pop.phi_open / np.maximum(pop.radii_m, 1e-30)**4
    return x/max(float(x.sum()), 1e-300)


def diagnostics(pop: PorePopulation, fine_radius_m: float = 25e-9) -> dict[str, float]:
    po = pop.phi_open; total = max(float(po.sum()), 1e-300)
    order = np.argsort(pop.radii_m); c = np.cumsum(po[order])/total
    def q(frac): return float(pop.radii_m[order][min(np.searchsorted(c, frac), len(c)-1)])
    return {"pore_D50_m": 2*q(.5), "pore_D90_m": 2*q(.9),
            "fine_pore_fraction": float(po[pop.radii_m <= fine_radius_m].sum()/total),
            "open_fraction": float(po.sum()), "isolated_fraction": float(pop.phi_iso.sum()),
            "closed_fraction": float(pop.phi_closed.sum())}
