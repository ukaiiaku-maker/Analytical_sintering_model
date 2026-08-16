from __future__ import annotations
from dataclasses import dataclass, replace
import json
import numpy as np

from .conditioned_950c import BARRIER
from .barrier_json import BarrierModel
from .densification import kinetic_state, density_rate, connectivity
from .energy_balance import solve_effective_stress
from .grain_growth import growth_state
from .integrator import ForwardModel, ModelParameters, ModelState
from .material_zro2 import MaterialParameters, R
from .pore_population import diagnostics, removal_weights
from .closed_channel_laws import closed_channel_rates


@dataclass
class ResolvedRuleParameters(ModelParameters):
    mechanism_mode: str = "resolved_rules"
    open_closed_handoff_mode: str = "resolved_default"
    closed_rate_factor: float = 1.
    closed_inventory_factor: float = 1.
    handoff_beta: float = .8
    handoff_closed_beta: float = 9.
    closed_tau0_s: float = 3.0e3
    closed_channel_law: str = "resolved_proxy_current"
    closed_radius_exponent: int = 3
    C_closed_GB: float = 1.0
    closed_prefactor_factor: float = 1.0
    closed_GB_use_renewal_activity: bool = True
    C_sigma_closed: float = 1.0
    C_transport_closed: float = 1.0
    closed_transport_length_factor: float = 1.0
    closed_exchange_tau_s: float = 1.0
    closed_event_strain: float = 1.0e-3
    closed_gas_pressure_fraction: float = .25
    C_surface_accommodation: float = 1.0e-24
    k0_closed_emp_sinv: float = 1.0e-3
    Q_closed_emp_J_mol: float = 130000.
    closed_empirical_size_exponent: float = 4.
    closed_rate_cap_time_s: float = 10.
    activity_mid: float = .15
    activity_width: float = .06
    activity_power: float = 1.5
    Q_PR_J_mol: float = 180000.
    T_PR_ref_K: float = 1573.15
    PR_to_iso_fraction: float = .12
    PR_memory_phi_scale: float = 2.0e-3
    precursor_prepare_tau_s: float = 5.0e4
    closed_transition_tau_s: float = 2.0e4
    rho_iso_mid: float = .78
    rho_closed_mid: float = .84
    rho_transition_width: float = .035
    accommodation_max: float = 1.
    accommodation_prepare_tau_s: float = 5.0e3
    accommodation_capacity_phi: float = .20
    PR_memory_decay_s: float = 2.0e6
    closed_migration_drag: float = 5.
    PR_migration_drag: float = 1.5
    gb_mobility_mode: str = "CS_endpoint_calibrated"
    M0_factor: float = 1.
    Q_M_J_mol_override: float | None = None
    no_PR_damage: bool = False
    no_closed_transition: bool = False
    no_closed_shrinkage: bool = False
    infinite_closed_accommodation: bool = False
    no_TJ_multihit: bool = True
    no_residual_stress: bool = True
    no_pore_drag: bool = False
    no_persistent_junction_state: bool = True
    no_sweep_coalescence: bool = True
    no_network_heterogeneity: bool = True


def resolved_material(base: MaterialParameters | None = None, M0_factor=1., Q_M_J_mol=None):
    m=base or MaterialParameters()
    return replace(m,M0_m4_J_s=m.M0_m4_J_s*M0_factor,Q_M_J_mol=m.Q_M_J_mol if Q_M_J_mol is None else Q_M_J_mol,
                   mobility_prefactor_status="global intrinsic growth uncertainty; never enters densification")

def conservative_adjacent_PR(phi_open,rate):
    crossing=rate[:-1]*phi_open[:-1];flux=np.zeros_like(phi_open);flux[:-1]-=crossing;flux[1:]+=crossing
    return flux,crossing


class ResolvedRuleModel(ForwardModel):
    def __init__(self, barrier=None, material=None, parameters=None):
        q=parameters or ResolvedRuleParameters()
        if q.mechanism_mode != "resolved_rules": raise ValueError("resolved model requires mechanism_mode='resolved_rules'")
        m=resolved_material(material,q.M0_factor,q.Q_M_J_mol_override)
        super().__init__(barrier or BarrierModel.load(BARRIER),m,q)

    def rates(self,state:ModelState,T_K:float):
        p,m,q=state.pores,self.material,self.parameters
        allowed={"resolved_default","diagnostic_open_recovery","closed_rate_boost_only","closed_inventory_boost_only","balanced_handoff","candidate_state_injection_diagnostic"}
        if q.open_closed_handoff_mode not in allowed: raise ValueError(f"unknown open/closed handoff mode: {q.open_closed_handoff_mode}")
        open_phi=float(p.phi_open.sum());conn=connectivity(open_phi,p.total)
        area=float(np.sum(p.number_open*4*np.pi*p.radii_m**2));r_eff=float(np.sum(p.phi_open*p.radii_m)/max(open_phi,1e-300))
        area_rate=-q.surface_power_length2_m2*m.D_s(T_K)*area/max(r_eff**4,1e-300)
        def erate(s):return kinetic_state(s,T_K,state.G_m,conn,self.barrier,m,q.sink_time_factor,q.site_density_multiplier)["edot_sinv"]
        power=solve_effective_stress(state.rho,area_rate,m.gamma_s_J_m2,erate,sigma_min=q.stress_min_Pa,sigma_max=q.stress_max_Pa)
        kin=kinetic_state(power.sigma_eff_Pa,T_K,state.G_m,conn,self.barrier,m,q.sink_time_factor,q.site_density_multiplier)
        tau_nuc=1/max(kin["r_nuc_sinv"],1e-300);tau_exchange=.35*kin["tau_sink_s"];tau_transport=.65*kin["tau_sink_s"];tau_cycle=tau_nuc+tau_exchange+tau_transport
        renewal=(tau_exchange+tau_transport)/tau_cycle
        fine=float(p.phi_open[p.radii_m<=25e-9].sum()/max(open_phi,1e-300));close_gate=1/(1+np.exp(-(state.rho-q.rho_closed_mid)/q.rho_transition_width));open_eligibility_base=1-.98*close_gate
        open_path_eligibility=open_eligibility_base;removable=conn*(.1+.9*fine)*open_path_eligibility
        geo=m.triple_line_geometry(state.G_m);edot=q.site_density_multiplier*geo["eps_event"]*renewal/max(tau_exchange+tau_transport,1e-300)*removable
        if q.open_closed_handoff_mode=="diagnostic_open_recovery":
            open_path_eligibility=1.;removable=1.;edot=kin["edot_sinv"]
        rho_open=density_rate(state.rho,edot);open_shrink=-removal_weights(p)*rho_open
        low=1/(1+np.exp(-(q.activity_mid-renewal)/q.activity_width))*max(1-renewal,0.)**q.activity_power
        theta=np.exp(np.clip(-q.Q_PR_J_mol/R*(1/T_K-1/q.T_PR_ref_K),-50,50));topology=.2+.8*fine
        rate=(0. if q.no_PR_damage else q.C_PR_m2)*m.D_s(T_K)/np.maximum(p.radii_m,1e-30)**4*low*theta*topology
        pr,crossing=conservative_adjacent_PR(p.phi_open,rate)
        iso_gate=1/(1+np.exp(-(state.rho-q.rho_iso_mid)/q.rho_transition_width));large=(p.radii_m/p.radii_m[-1])**2
        to_iso=q.PR_to_iso_fraction*rate*iso_gate*large*p.phi_open+p.phi_open*state.PR_memory/q.precursor_prepare_tau_s*iso_gate*large
        transition_on=0. if q.no_closed_transition else 1.;inventory_factor=q.closed_inventory_factor if q.open_closed_handoff_mode=="closed_inventory_boost_only" else 1.;to_closed_iso=inventory_factor*transition_on*close_gate*(.02+.98*state.PR_memory)*p.phi_iso/q.closed_transition_tau_s;to_closed_open=inventory_factor*transition_on*close_gate*state.PR_memory*p.phi_open/q.closed_transition_tau_s;to_closed=to_closed_iso+to_closed_open
        open_dot=open_shrink+pr-to_iso-to_closed_open;iso_dot=to_iso-to_closed_iso;closed_dot=to_closed.copy()
        A=1. if q.infinite_closed_accommodation else np.clip(state.A_closed,0,q.accommodation_max)
        closed_activation=renewal*np.exp(np.clip(-.35*q.Q_PR_J_mol/R*(1/T_K-1/q.T_PR_ref_K),-30,30))
        shrink_law,closed_law_diag=closed_channel_rates(state,T_K,self.barrier,m,q,renewal)
        closed_availability=float(shrink_law.sum());open_availability=max(rho_open,0.);handoff_readiness=closed_availability/(closed_availability+open_availability+1e-300)
        closed_factor=q.closed_rate_factor if q.open_closed_handoff_mode=="closed_rate_boost_only" else 1.
        if q.open_closed_handoff_mode=="balanced_handoff":
            open_path_eligibility=open_eligibility_base*(1-q.handoff_beta*handoff_readiness)
            removable=conn*(.1+.9*fine)*open_path_eligibility
            edot=q.site_density_multiplier*geo["eps_event"]*renewal/max(tau_exchange+tau_transport,1e-300)*removable
            rho_open=density_rate(state.rho,edot);open_shrink=-removal_weights(p)*rho_open
            open_dot=open_shrink+pr-to_iso-to_closed_open
            closed_factor=1+q.handoff_closed_beta*state.PR_memory*A
        shrink=np.zeros_like(p.phi_closed) if q.no_closed_shrinkage else closed_factor*shrink_law;closed_dot-=shrink;rho_closed=float(shrink.sum())
        base=growth_state(state.G_m,p.radii_m,p.phi_open,T_K,m,zener_strength=q.zener_strength,mobile_drag_scale=q.mobile_drag_scale)
        intrinsic=m.M_GB(T_K)*m.gamma_GB_J_m2/max(state.G_m,1e-30);pore_gamma=1. if q.no_pore_drag else base["Gamma_growth"]
        closed_drag=1/(1+q.closed_migration_drag*float(p.phi_closed.sum())/max(1-state.rho,1e-12));pr_drag=1/(1+q.PR_migration_drag*state.PR_memory);event_gamma=.1+.9*np.sqrt(np.clip(renewal,0,1));Gamma=float(np.clip(pore_gamma*closed_drag*pr_drag*event_gamma,0,1));actual=intrinsic*Gamma
        growth={**base,"M_GB_intrinsic":m.M_GB(T_K),"Gamma_migration":Gamma,"G_dot_intrinsic_m_s":intrinsic,"G_dot_actual_m_s":actual,"G_dot_m_s":actual,"pore_Zener_drag_contribution":pore_gamma,"closed_accommodation_migration_contribution":closed_drag,"persistent_TJ_contribution":1.}
        total_open=open_shrink+pr-to_iso-to_closed_open;ta=np.full_like(p.radii_m,np.inf);mask=(p.phi_open>0)&(total_open<0);ta[mask]=p.phi_open[mask]/(-total_open[mask])
        arrays={"pore_radii_m_json":json.dumps(p.radii_m.tolist()),"phi_open_json":json.dumps(p.phi_open.tolist()),"phi_iso_json":json.dumps(p.phi_iso.tolist()),"phi_closed_json":json.dumps(p.phi_closed.tolist()),"phi_open_dot_json":json.dumps(open_dot.tolist()),"phi_iso_dot_json":json.dumps(iso_dot.tolist()),"phi_closed_dot_json":json.dumps(closed_dot.tolist()),"tau_remove_s_json":json.dumps(ta.tolist())}
        diag={**kin,**growth,**power.__dict__,**diagnostics(p),**arrays,**closed_law_diag,"tau_nuc_s":tau_nuc,"tau_exchange_s":tau_exchange,"tau_transport_s":tau_transport,"tau_cycle_s":tau_cycle,"activity":renewal,"activity_open":renewal,"activity_closed":closed_activation,"connected_removable_factor":removable,"open_path_eligibility":open_path_eligibility,"open_eligibility_base":open_eligibility_base,"open_eligibility_eff":open_path_eligibility,"closed_availability":closed_availability,"handoff_readiness":handoff_readiness,"tau_open_s":open_phi/max(rho_open,1e-300),"tau_closed_s":float(p.phi_closed.sum())/max(rho_closed,1e-300),"local_activation_stress_Pa":power.sigma_eff_Pa,"rho_dot_open_sinv":rho_open,"rho_dot_closed_sinv":rho_closed,"rho_dot_total_sinv":rho_open+rho_closed,"PR_coarsening_flux":float(crossing.sum()),"PR_relocation_flux":float(crossing.sum()),"PR_to_isolated_flux":float(to_iso.sum()),"PR_to_closed_precursor_flux":float(to_closed.sum()),"bin_crossing_rate":float(crossing.sum()),"isolation_rate":float(to_iso.sum()),"closure_rate":float(to_closed.sum()),"closed_shrinkage_flux":rho_closed,"PR_low_activity_gate":low,"PR_thermal_factor":theta,"PR_memory":state.PR_memory,"cumulative_PR_work":state.cumulative_PR_work,"A_closed":A,"mechanism_mode":"resolved_rules","open_closed_handoff_mode":q.open_closed_handoff_mode,"gb_mobility_mode":q.gb_mobility_mode,"M0_factor":q.M0_factor,"Q_M_kJ_mol":m.Q_M_J_mol/1000}
        return open_dot,iso_dot,closed_dot,rho_closed,growth,diag

    def step(self,state,T_K,dt_s):
        od,ii,cc,closed_rate,growth,diag=self.rates(state,T_K);p=state.pores.copy();p.phi_open=np.maximum(0,p.phi_open+dt_s*od);p.phi_iso=np.maximum(0,p.phi_iso+dt_s*ii);p.phi_closed=np.maximum(0,p.phi_closed+dt_s*cc)
        pr_rate=diag["PR_coarsening_flux"]/self.parameters.PR_memory_phi_scale;memory=np.clip(state.PR_memory+dt_s*(pr_rate*(1-state.PR_memory)-state.PR_memory/self.parameters.PR_memory_decay_s),0,1)
        if self.parameters.infinite_closed_accommodation:A=1.
        else:
            preparation=memory*(self.parameters.accommodation_max-state.A_closed)/self.parameters.accommodation_prepare_tau_s
            if self.parameters.closed_channel_law=="surface_diffusion_accommodation_only": preparation+=diag["A_dot_closed_sinv"]
            A=np.clip(state.A_closed+dt_s*(preparation-closed_rate/self.parameters.accommodation_capacity_phi),0,self.parameters.accommodation_max)
        work=state.cumulative_PR_work+dt_s*diag["PR_coarsening_flux"]
        diag.update(diagnostics(p));diag.update({"pore_radii_m_json":json.dumps(p.radii_m.tolist()),"phi_open_json":json.dumps(p.phi_open.tolist()),"phi_iso_json":json.dumps(p.phi_iso.tolist()),"phi_closed_json":json.dumps(p.phi_closed.tolist()),"A_closed":float(A),"PR_memory":float(memory),"cumulative_PR_work":work})
        return ModelState(state.t_s+dt_s,T_K,1-p.total,state.G_m+dt_s*growth["G_dot_m_s"],p,float(A),float(memory),work),diag


def resolved_initial_state(state):
    return replace(state,A_closed=0.,PR_memory=0.,cumulative_PR_work=0.)


ABLATIONS=("no_PR_damage","no_closed_transition","no_closed_shrinkage","infinite_closed_accommodation","no_TJ_multihit","no_residual_stress","no_pore_drag","no_persistent_junction_state","no_sweep_coalescence","no_network_heterogeneity")
