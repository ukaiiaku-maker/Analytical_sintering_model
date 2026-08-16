from __future__ import annotations
from dataclasses import dataclass
import json
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
    PR_memory: float = 0.
    cumulative_PR_work: float = 0.


@dataclass
class ModelParameters:
    site_density_multiplier: float = 39.5
    sink_time_factor: float = 1.0
    surface_power_length2_m2: float = 1e-19
    C_PR_m2: float = 1e-23
    C_iso_m2: float = 2e-25
    C_close_m2: float = 1e-25
    closed_tau0_s: float = 1e5
    rho_close_mid: float = .90
    rho_close_width: float = .03
    zener_strength: float = 1.
    mobile_drag_scale: float = 1.
    stress_min_Pa: float = 1e3
    stress_max_Pa: float = 2.5e8


class ForwardModel:
    def __init__(self, barrier, material=None, parameters=None):
        self.barrier=barrier; self.material=material or MaterialParameters(); self.parameters=parameters or ModelParameters()

    def initial_state(self):
        p=initial_population(); return ModelState(0., 298.15, 1-p.total, 10.20e-9, p)

    def rates(self, state: ModelState, T_K: float):
        p, m, q = state.pores, self.material, self.parameters
        conn = connectivity(float(p.phi_open.sum()), p.total)
        # Surface-area loss scale from diffusion; this is an instantaneous state law.
        area = float(np.sum(p.number_open*4*np.pi*p.radii_m**2))
        r_eff=float(np.sum(p.phi_open*p.radii_m)/max(float(p.phi_open.sum()),1e-300))
        area_rate = -q.surface_power_length2_m2*m.D_s(T_K)*area/max(r_eff**4,1e-300)
        def erate(s): return kinetic_state(s,T_K,state.G_m,conn,self.barrier,m,q.sink_time_factor,q.site_density_multiplier)["edot_sinv"]
        power=solve_effective_stress(state.rho,area_rate,m.gamma_s_J_m2,erate,
                                     sigma_min=q.stress_min_Pa,sigma_max=q.stress_max_Pa)
        kin=kinetic_state(power.sigma_eff_Pa,T_K,state.G_m,conn,self.barrier,m,q.sink_time_factor,q.site_density_multiplier)
        rd=density_rate(state.rho,kin["edot_sinv"])
        od=-removal_weights(p)*rd
        excess=power.P_excess_W_m3/(power.P_surf_W_m3+1e-300)
        po,pi,pc,flux=transfer_fluxes(p,state.rho,m.D_s(T_K),kin["activity"],excess,q.C_PR_m2,q.C_iso_m2,q.C_close_m2,q.rho_close_mid,q.rho_close_width)
        # Bounded accommodation proxy; only the named closed reservoir shrinks.
        tau0=q.closed_tau0_s*(p.radii_m/25e-9)**4/(max(kin["activity"]*state.A_closed,1e-12))
        shrink=p.phi_closed/tau0
        pc-=shrink
        growth=growth_state(state.G_m,p.radii_m,p.phi_open,T_K,m,zener_strength=q.zener_strength,mobile_drag_scale=q.mobile_drag_scale)
        open_total=od+po
        tau=np.full_like(p.radii_m,np.inf)
        mask=(p.phi_open>0)&(open_total<0); tau[mask]=p.phi_open[mask]/(-open_total[mask])
        arrays={"pore_radii_m_json":json.dumps(p.radii_m.tolist()),
                "phi_open_json":json.dumps(p.phi_open.tolist()),
                "phi_iso_json":json.dumps(p.phi_iso.tolist()),
                "phi_closed_json":json.dumps(p.phi_closed.tolist()),
                "phi_open_dot_json":json.dumps(open_total.tolist()),
                "phi_iso_dot_json":json.dumps(pi.tolist()),
                "phi_closed_dot_json":json.dumps(pc.tolist()),
                "tau_remove_s_json":json.dumps(tau.tolist())}
        return open_total,pi,pc,float(shrink.sum()),growth,{**kin,**growth,**flux,
               "rho_dot_open_sinv":rd,"rho_dot_closed_sinv":float(shrink.sum()),
               **power.__dict__,**diagnostics(p),**arrays}

    def step(self,state,T_K,dt_s):
        od,ii,cc,closed_rate,growth,diag=self.rates(state,T_K)
        p=state.pores.copy(); p.phi_open=np.maximum(0,p.phi_open+dt_s*od); p.phi_iso=np.maximum(0,p.phi_iso+dt_s*ii); p.phi_closed=np.maximum(0,p.phi_closed+dt_s*cc)
        A=np.clip(state.A_closed+dt_s*(-0.2*closed_rate+(1-state.A_closed)/1e6),0,1)
        nxt=ModelState(state.t_s+dt_s,T_K,1-p.total,state.G_m+dt_s*growth["G_dot_m_s"],p,float(A))
        return nxt,diag

    def run(self, thermal_path, dt_s=10., record_every_s=60., initial_state=None):
        state=initial_state if initial_state is not None else self.initial_state(); rows=[]
        next_record=0.
        while state.t_s < thermal_path.t_end_s:
            T=thermal_path.temperature_K(state.t_s,state.rho)
            trial=min(dt_s,thermal_path.t_end_s-state.t_s)
            od,ii,cc,_,growth,_=self.rates(state,T)
            maximum=max(float(np.max(np.abs(od))),float(np.max(np.abs(ii))),float(np.max(np.abs(cc))),1e-300)
            trial=min(trial,2e-3/maximum,.01*state.G_m/max(growth["G_dot_m_s"],1e-300))
            trial=max(min(trial,thermal_path.t_end_s-state.t_s),1e-6)
            state,d=self.step(state,T,trial)
            if state.t_s+1e-9 >= next_record or state.t_s >= thermal_path.t_end_s:
                rows.append({"t_s":state.t_s,"T_K":T,"rho":state.rho,"G_m":state.G_m,
                             "A_closed":state.A_closed,**d})
                next_record=state.t_s+record_every_s
        return state,rows
