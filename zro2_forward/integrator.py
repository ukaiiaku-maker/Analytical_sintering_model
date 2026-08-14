from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .material_zro2 import MaterialParameters
from .pore_population import PorePopulation, initial_population, removal_weights, transfer_fluxes, diagnostics
from .densification import kinetic_state, density_rate, connectivity
from .energy_balance import solve_effective_stress
from .grain_growth import growth_state


@dataclass
class ModelState:
    t_s: float
    T_K: float
    rho: float
    G_m: float
    pores: PorePopulation
    A_closed: float = 1.


class ForwardModel:
    def __init__(self, barrier, material=None): self.barrier=barrier; self.material=material or MaterialParameters()

    def initial_state(self):
        p=initial_population(); return ModelState(0., 298.15, 1-p.total, 10.20e-9, p)

    def rates(self, state: ModelState, T_K: float):
        p, m = state.pores, self.material
        conn = connectivity(float(p.phi_open.sum()), p.total)
        # Surface-area loss scale from diffusion; this is an instantaneous state law.
        area = float(np.sum(p.number_open*4*np.pi*p.radii_m**2))
        area_rate = -m.D_s(T_K)*area/max(state.G_m**2, 1e-300)
        def erate(s): return kinetic_state(s,T_K,state.G_m,conn,self.barrier,m)["edot_sinv"]
        power=solve_effective_stress(state.rho,area_rate,m.gamma_s_J_m2,erate)
        kin=kinetic_state(power.sigma_eff_Pa,T_K,state.G_m,conn,self.barrier,m)
        rd=density_rate(state.rho,kin["edot_sinv"])
        od=-removal_weights(p)*rd
        excess=power.P_excess_W_m3/(power.P_surf_W_m3+1e-300)
        po,pi,pc,flux=transfer_fluxes(p,state.rho,m.D_s(T_K),kin["activity"],excess)
        # Bounded accommodation proxy; only the named closed reservoir shrinks.
        tau0=1e5*(p.radii_m/25e-9)**4/(max(kin["activity"]*state.A_closed,1e-12))
        shrink=p.phi_closed/tau0
        pc-=shrink
        growth=growth_state(state.G_m,p.radii_m,p.phi_open,T_K,m)
        return od+po,pi,pc,float(shrink.sum()),growth,{**kin,**growth,**flux,**power.__dict__,**diagnostics(p)}

    def step(self,state,T_K,dt_s):
        od,ii,cc,closed_rate,growth,diag=self.rates(state,T_K)
        p=state.pores.copy(); p.phi_open=np.maximum(0,p.phi_open+dt_s*od); p.phi_iso=np.maximum(0,p.phi_iso+dt_s*ii); p.phi_closed=np.maximum(0,p.phi_closed+dt_s*cc)
        A=np.clip(state.A_closed+dt_s*(-0.2*closed_rate+(1-state.A_closed)/1e6),0,1)
        nxt=ModelState(state.t_s+dt_s,T_K,1-p.total,state.G_m+dt_s*growth["G_dot_m_s"],p,float(A))
        return nxt,diag

    def run(self, thermal_path, dt_s=1.):
        state=self.initial_state(); rows=[]
        while state.t_s < thermal_path.t_end_s:
            T=thermal_path.temperature_K(state.t_s,state.rho)
            state,d=self.step(state,T,min(dt_s,thermal_path.t_end_s-state.t_s))
            rows.append({"t_s":state.t_s,"T_K":T,"rho":state.rho,"G_m":state.G_m,**d})
        return state,rows
