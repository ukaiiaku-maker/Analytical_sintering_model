"""Bin-resolved closed-pore inventory, gas, and finite accommodation state."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class ClosedPoreEvolution:
    N_closed: np.ndarray
    r_reference_m: np.ndarray
    chi_shrink: np.ndarray
    A_max: np.ndarray
    A_used: np.ndarray
    A_recovered: np.ndarray
    C_geom: np.ndarray
    P_gas_initial_Pa: np.ndarray
    V_initial_per_pore_m3: np.ndarray
    T_initial_K: float
    gas_enabled: bool = False

    @classmethod
    def initialize(cls, pores, T_K, accommodation_max=1., C_geom=1., chi_shrink=1.,
                   gas_enabled=False, gas_initial_fraction=.25, gamma_s=1., available_initial=None):
        r=np.asarray(pores.radii_m,float);v=4*np.pi*r**3/3
        N=np.divide(pores.phi_closed,v,out=np.zeros_like(r),where=v>0)
        cap=2*gamma_s/np.maximum(r,1e-30)
        amax=np.full_like(r,max(accommodation_max,0));avail=np.clip(accommodation_max if available_initial is None else available_initial,0,accommodation_max)
        return cls(N,r.copy(),np.full_like(r,np.clip(chi_shrink,0,1)),amax,
                   amax-avail,np.zeros_like(r),np.full_like(r,np.clip(C_geom,0,1)),
                   cap*np.clip(gas_initial_fraction,0,1) if gas_enabled else np.zeros_like(r),v.copy(),float(T_K),gas_enabled)

    def copy(self):
        return ClosedPoreEvolution(*(x.copy() if isinstance(x,np.ndarray) else x for x in self.__dict__.values()))

    def radii_m(self, phi_closed):
        volume=np.divide(phi_closed,self.N_closed,out=np.zeros_like(phi_closed),where=self.N_closed>0)
        dynamic=np.cbrt(np.maximum(3*volume/(4*np.pi),0))
        return np.where(self.N_closed>0,np.maximum(dynamic,1e-30),self.r_reference_m)

    def A_available(self):
        return np.clip(self.A_max-self.A_used+self.A_recovered,0,self.A_max)

    def gas_pressure_Pa(self, phi_closed, T_K):
        if not self.gas_enabled:return np.zeros_like(phi_closed)
        volume=np.divide(phi_closed,self.N_closed,out=self.V_initial_per_pore_m3.copy(),where=self.N_closed>0)
        return self.P_gas_initial_Pa*self.V_initial_per_pore_m3/np.maximum(volume,1e-300)*(T_K/self.T_initial_K)

    def advance(self, old_phi, new_phi, shrink_rate, dt_s, T_K, D_s, beta_A, recovery_rate, shape_coefficient):
        out=self.copy();removed=np.minimum(np.maximum(shrink_rate*dt_s,0),old_phi)
        closure=np.maximum(new_phi-(old_phi-removed),0)
        vref=4*np.pi*out.r_reference_m**3/3;out.N_closed+=np.divide(closure,vref,out=np.zeros_like(closure),where=vref>0)
        out.A_used=np.minimum(out.A_max,out.A_used+max(beta_A,0)*removed)
        available=out.A_available();shape=shape_coefficient*D_s/np.maximum(out.radii_m(new_phi),1e-30)**4*(out.A_max-available)
        recovery=np.maximum(recovery_rate,0)*(out.A_max-available)+np.maximum(shape,0)
        out.A_recovered=np.minimum(out.A_used,out.A_recovered+dt_s*recovery)
        return out
