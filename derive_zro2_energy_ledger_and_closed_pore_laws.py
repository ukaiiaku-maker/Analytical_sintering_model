"""Diagnostic ZrO2 energy ledger and physical closed-pore law candidates.

This module does not modify the accepted forward integrator.  It evaluates
schedule-independent constitutive candidates and reconstructs an energy ledger
from already-generated histories.  Results are diagnostic, not validation.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd

from zro2_forward.barrier_json import BarrierModel
from zro2_forward.material_zro2 import MaterialParameters

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results/zro2_forward_energy_ledger_closed_pore_derivation"
FIG = OUT / "figures"
SOURCE = ROOT / "results/zro2_forward_final_summary_figures"
BARRIER_PATH = ROOT / "data/zro2/bicrystal_creep_barrier_export.json"
EPS = 1.0e-30
MAT = MaterialParameters()
BARRIER = BarrierModel.load(BARRIER_PATH)
R_REF_M = 25.0e-9


def barrier_sha256() -> str:
    return hashlib.sha256(BARRIER_PATH.read_bytes()).hexdigest()


def capillary_stress(radius_m: float, gas_pressure_pa: float = 0.0,
                     geometry_factor: float = 1.0) -> float:
    if radius_m <= 0:
        return 0.0
    return max(geometry_factor * 2.0 * MAT.gamma_s_J_m2 / radius_m - gas_pressure_pa, 0.0)


def pore_area_density(phi: np.ndarray, radius_m: np.ndarray) -> float:
    """Spherical pore area per bulk volume, sum(3 phi_i/r_i), m^-1."""
    p = np.asarray(phi, float)
    r = np.asarray(radius_m, float)
    return float(np.sum(3.0 * p[r > 0] / r[r > 0]))


def gb_area_density(grain_m: float) -> float:
    return MAT.C_GB / max(grain_m, EPS)


def external_area_density(rho: float, specimen_length_m: float = 1.0e-3) -> float:
    """Compact-scale external area proxy; geometry is explicitly unresolved."""
    return 6.0 * max(rho, 0.0) ** (2.0 / 3.0) / specimen_length_m


def closed_state_rate(mode: str, T_K: float, radius_m: float, phi_closed: float,
                      shrinkability: float, accommodation: float,
                      gas_fraction: float = 0.0, exponent: int = 3,
                      geometry_factor: float = 1.0) -> dict[str, float | bool | str]:
    """Evaluate one bin without any processing-path identifiers.

    The GB coefficient r_ref**(m-2) supplies the units required by the proposed
    expression.  It is an unresolved dimensional coefficient, not a fitted energy.
    """
    pcap = geometry_factor * 2.0 * MAT.gamma_s_J_m2 / max(radius_m, EPS)
    pgas = np.clip(gas_fraction, 0.0, 1.5) * pcap
    sigma = capillary_stress(radius_m, pgas, geometry_factor)
    active_inventory = max(phi_closed, 0.0) * max(shrinkability, 0.0) * max(accommodation, 0.0)
    Dgb, Ds = MAT.D_GB(T_K), MAT.D_s(T_K)
    gstar = float(BARRIER.Gstar(sigma, T_K))
    rnuc = MAT.nu0_sinv * math.exp(-gstar / (MAT.kB * T_K))
    tau_sink = math.inf if sigma <= 0 else (
        MAT.kB * T_K / (sigma * MAT.Omega_m3)
        * radius_m**2 / max(Dgb, EPS)
    )
    tau_exchange = radius_m**2 / max(Ds, EPS)
    lam = 0.0 if not np.isfinite(tau_sink) else rnuc * tau_sink
    activity = lam / (1.0 + lam)
    tau_cycle = 1.0 / max(rnuc, EPS) + tau_exchange + tau_sink
    scale = (R_REF_M / max(radius_m, EPS)) ** exponent

    shape_rate = Ds / max(radius_m, EPS) ** 2 * max(1.0 - accommodation, 0.0)
    if mode == "renewal_limited_closed_shrinkage":
        rate = active_inventory * scale / max(tau_cycle, EPS)
        status = "physical_inputs_plus_semi_phenomenological_event_scale"
    elif mode == "GB_diffusion_closed_shrinkage":
        dimensional_prefactor = R_REF_M ** (exponent - 2)
        rate = (dimensional_prefactor * active_inventory * Dgb * MAT.Omega_m3
                * sigma / (MAT.kB * T_K) * max(radius_m, EPS) ** (-exponent))
        status = "physical_transport_with_unresolved_dimensional_geometry"
    elif mode == "gas_limited_closed_shrinkage":
        dimensional_prefactor = R_REF_M ** (exponent - 2)
        rate = (dimensional_prefactor * active_inventory * Dgb * MAT.Omega_m3
                * sigma / (MAT.kB * T_K) * max(radius_m, EPS) ** (-exponent))
        status = "physical_counterpressure_with_unresolved_geometry"
    elif mode == "surface_diffusion_accommodation_only":
        rate = 0.0
        status = "physical_shape_change_non_densifying"
    elif mode == "empirical_reduced_closure":
        q_emp = 300_000.0
        rate = active_inventory * 1.0e5 * math.exp(-q_emp / (8.31446261815324 * T_K))
        status = "empirical_diagnostic_only"
    else:
        raise ValueError(mode)
    return {
        "mode": mode, "T_K": T_K, "radius_m": radius_m, "phi_closed": phi_closed,
        "shrinkability": shrinkability, "accommodation": accommodation,
        "gas_fraction": gas_fraction, "pcap_Pa": pcap, "Pgas_Pa": pgas,
        "sigma_Pa": sigma, "D_GB_m2_s": Dgb, "D_s_m2_s": Ds,
        "Gstar_J": gstar, "temperature_extrapolated": not BARRIER.temperature_in_fit_range(T_K),
        "r_nuc_sinv": rnuc, "tau_sink_s": tau_sink, "tau_exchange_s": tau_exchange,
        "Lambda": lam, "activity": activity, "shape_rate_sinv": shape_rate,
        "rho_dot_closed_sinv": max(float(rate), 0.0), "exponent": exponent,
        "physical_status": status,
    }


def conservative_transfer(phi: np.ndarray, flux: float, dt: float,
                          source: int, destination: int) -> np.ndarray:
    out = np.asarray(phi, float).copy()
    amount = min(max(flux, 0.0) * max(dt, 0.0), out[source])
    out[source] -= amount
    out[destination] += amount
    return out


def density_identity(phi_open: np.ndarray, phi_iso: np.ndarray,
                     phi_closed: np.ndarray) -> float:
    return 1.0 - float(np.sum(phi_open) + np.sum(phi_iso) + np.sum(phi_closed))


def _write_csv(name: str, rows) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / name, index=False)
    return frame


def write_source_registries() -> None:
    rules = [
        ("R1", "docs/ZRO2_FORWARD_ACTIVE_MODEL_EQUATIONS_AND_PARAMETERS.md", "physical interpretation", "Fast firing is primarily nucleation-limited onset.", "Retain the stress-resolved fitted barrier.", "already_implemented"),
        ("R2", "docs/FINAL_FAST_FIRING_AND_TWO_STEP_MECHANISM_SYNTHESIS.md", "mechanism separation", "Two-step response uses PR-prepared closed/accommodation memory.", "Audit preparation and shrinkage separately.", "implemented_as_proxy"),
        ("R3", "docs/ZRO2_FORWARD_CLOSED_CHANNEL_PHYSICAL_MAPPING.md", "density identity", "Only named pore-volume shrinkage changes density.", "Keep redistribution conservative.", "already_implemented"),
        ("R4", "docs/ZRO2_FORWARD_CLOSED_CHANNEL_PROPERTY_GAP.md", "property gap", "The closed rate lacks a fully derived dimensional prefactor.", "Expose the prefactor as a calibration target.", "missing_physical_mapping"),
        ("R5", "docs/ZRO2_FORWARD_QCLOSED_APPARENT_PROPERTY_CORRECTION.md", "apparent slope", "Q_closed_app is post-run and not a material property.", "Never introduce physical Q_closed.", "already_implemented"),
        ("R6", "docs/ZRO2_FORWARD_FINAL_SUMMARY_FIGURES.md", "energy closure", "The strict GB-area-loss equality generated the displayed diagnostic trajectories.", "Retain only as an ablation until a full ledger is audited.", "conflicting_current_model"),
        ("R7", "docs/ZRO2_FORWARD_NEXT_DECISION_FROM_PROPERTY_BOUNDS.md", "next decision", "Physical property mapping precedes new mechanism search.", "Run fixed-state analytical tests first.", "already_implemented"),
        ("R8", "1_2-step_discussion5(2).docx", "unavailable", "Requested manuscript/SI discussion was not supplied or found.", "Do not infer its contents; use repository-derived windows only.", "literature_needed"),
        ("R9", "results/final_mechanism_synthesis_and_property_windows", "candidate family", "Candidate 693168 is conditional mechanism evidence only.", "Do not copy its values into ZrO2 physical parameters.", "already_implemented"),
    ]
    _write_csv("source_rule_manifest.csv", [dict(zip(
        ["rule_id", "source_file", "source_section", "statement", "implementation_implication", "status"], r)) for r in rules])

    lit = [
        ("L01", "Coble, J. Appl. Phys. 32 (1961), doi:10.1063/1.1736107", "final-stage GB/volume diffusion", "closed-pore shrinkage from curvature-driven diffusion", "D_GB,gamma,r,T,geometry", "r,phi,T", "geometry dependent", "GB diffusivity", "geometry", True, False, False, False, "full text not locally available", "primary bibliographic/abstract record; equation normalization requires full text"),
        ("L02", "Nichols and Mullins, J. Appl. Phys. 36 (1965), doi:10.1063/1.1714360", "surface diffusion", "surface transport relaxes curvature without center approach", "D_s,gamma,r,T", "shape,r,T", "r^-4 characteristic shape time", "surface diffusivity", "geometry", False, True, False, False, "geometry coefficient", "primary bibliographic record; used only for non-densifying accommodation"),
        ("L03", "Entrapped-gas sintering study (1971), doi:10.1016/0025-5416(71)90060-7", "closed gas opposition", "sigma=max(2 gamma/r-Pgas,0)", "gamma,r,Pgas", "r,V,T,Pgas", "r^-1 stress", "none", "ideal-gas inventory", False, False, False, True, "initial gas amount", "primary record; abstract-level access"),
        ("L04", "Densification of ceramics containing entrapped gases (1989), doi:10.1016/0955-2219(89)90020-4", "gas-limited densification", "gas pressure retards or arrests pore shrinkage", "transport,gas,geometry", "r,V,T", "model dependent", "transport law", "geometry", True, False, False, False, "full coefficients", "primary record; abstract-level access"),
        ("L05", "Watanabe and Masuda, J. Jpn. Soc. Powder Metall. 29 (1982), doi:10.2497/jjspm.29.151", "pore pinning", "pinning strength scales with pore fraction divided by radius", "r_p,f_v,gamma_GB", "pore bins,G", "R_Z proportional to r_p/f_v", "none", "geometry", False, False, True, True, "shape factor", "primary paper record"),
        ("L06", "Chen and Wang, Nature 404 (2000), doi:10.1038/35004548", "two-step kinetic window", "densification retained while boundary migration is suppressed", "GB diffusion,GB mobility,T", "rho,G,T", "not applicable", "separate diffusion/migration barriers", "experiment", True, False, True, False, "material-specific kinetics", "primary abstract/full summary accessible"),
        ("L07", "Wang, Chen and Chen, JACS 89 (2006), doi:10.1111/j.1551-2916.2005.00763.x", "two-step boundary", "lower densification and upper migration boundaries", "diffusion,mobility,junction state", "rho,G,T,time", "not applicable", "mechanism dependent", "experiment", True, False, True, False, "ZrO2-specific mapping", "primary abstract record; Y2O3 system"),
        ("L08", "Mazaheri et al., Mater. Sci. Eng. A 492 (2008), repository PDF", "3Y-TZP densification and growth", "closed gas retards densification; pore closure changes pinning", "Q_s,Q_g,pore state", "rho,G,T", "not extracted", "Q_s=485+/-12 kJ/mol", "experiment/MSC", True, False, True, False, "closed-bin prefactor", "local primary PDF; reported growth slope 546+/-23 kJ/mol"),
        ("L09", "Pouchly, Maca and Shen (2013), repository PDF", "two-step MSC", "kinetic changes accompany closed-porosity stage", "T,rho,G", "rho,G,T,time", "not extracted", "apparent MSC slopes", "experiment", True, False, True, False, "state-resolved pore data", "local primary PDF"),
        ("L10", "1_2-step_discussion5(2).docx", "requested discussion", "not accessible", "unknown", "unknown", "unknown", "unknown", "unknown", False, False, False, False, "entire source", "file absent; no claims extracted"),
    ]
    cols = ["literature_id","source","mechanism","governing_equation","physical_inputs","state_variables","radius_exponent","activation_energy_source","prefactor_source","changes_density_directly","changes_shape_only","changes_migration_only","directly_implementable","missing_inputs","notes"]
    _write_csv("literature_seed_registry.csv", [dict(zip(cols, r)) for r in lit])


def write_law_registries() -> None:
    laws = [
        ("renewal_limited_closed_shrinkage", "phi chi A (r_ref/r)^m/tau_cycle", "G*(sigma,T),D_GB,D_s", "event strain/geometry", True, False, "semi-phenomenological"),
        ("GB_diffusion_closed_shrinkage", "C phi chi A D_GB Omega sigma/(kBT) r^-m", "D_GB,Omega,sigma,T,r", "C has units m^(m-2)", True, False, "dimensional calibration target"),
        ("surface_diffusion_accommodation_only", "A_dot proportional D_s r^-4(Amax-A)", "D_s,r,shape", "shape coefficient", False, False, "physical scaling"),
        ("gas_limited_closed_shrinkage", "GB law with sigma=max(2gamma/r-Pgas,0)", "D_GB,gamma,r,Pgas", "gas inventory and geometry C", True, False, "physical plus geometry target"),
        ("empirical_reduced_closure", "phi chi A k0 exp(-Qemp/RT)", "state and empirical pair", "k0,Qemp", True, True, "empirical diagnostic only"),
    ]
    cols = ["law_id","equation","physical_inputs","unresolved_terms","changes_density_directly","uses_empirical_Q","physical_status"]
    _write_csv("closed_pore_law_registry.csv", [dict(zip(cols, r)) for r in laws])
    mappings = [
        ("gamma_s","MaterialParameters.gamma_s_J_m2","physical",1.0,"J m^-2"),
        ("D_GB","0.056 exp(-380000/RT)","physical/current fixed",0.056,"m2 s^-1 prefactor"),
        ("D_s","0.10 exp(-380000/RT)","physical/current fixed",0.10,"m2 s^-1 prefactor"),
        ("Omega","MaterialParameters.Omega_m3","physical",MAT.Omega_m3,"m3"),
        ("Gstar","fitted barrier JSON","physical fit; extrapolated below fit",np.nan,"J"),
        ("phi_closed","evolving closed inventory","state-derived",np.nan,"1"),
        ("chi_shrink","closed-bin shrinkability","state-derived/proxy",np.nan,"1"),
        ("A_closed","available accommodation","state-derived/proxy",np.nan,"1"),
        ("C_GB_c", "r_ref^(m-2)", "semi-phenomenological dimensional target", np.nan, "m^(m-2)"),
        ("Q_closed_emp","empirical comparator only","empirical",300.0,"kJ mol^-1"),
    ]
    _write_csv("closed_pore_parameter_mapping.csv", [dict(zip(["parameter","mapping","category","value","units"], r)) for r in mappings])
    units = []
    for m in (3,4):
        units += [
            {"law_id":"renewal_limited_closed_shrinkage","exponent":m,"rate_units":"s^-1","prefactor_units":"1","unit_check_pass":True,"note":"dimensionless inventory and radius ratio divided by cycle time"},
            {"law_id":"GB_diffusion_closed_shrinkage","exponent":m,"rate_units":"s^-1","prefactor_units":f"m^{m-2}","unit_check_pass":True,"note":"required dimensional coefficient is not yet geometry-derived"},
            {"law_id":"gas_limited_closed_shrinkage","exponent":m,"rate_units":"s^-1","prefactor_units":f"m^{m-2}","unit_check_pass":True,"note":"same transport audit with explicit counterpressure"},
        ]
    units += [{"law_id":"surface_diffusion_accommodation_only","exponent":4,"rate_units":"shape s^-1","prefactor_units":"1","unit_check_pass":True,"note":"D_s/r^2 rate form; density contribution identically zero"},
              {"law_id":"empirical_reduced_closure","exponent":0,"rate_units":"s^-1","prefactor_units":"s^-1","unit_check_pass":True,"note":"diagnostic empirical comparator; not promotable"}]
    _write_csv("closed_pore_unit_audit.csv", units)
    groups = [
        ("capillary_to_gas","Pgas/(2 gamma/r)","gas arrest as group approaches one"),
        ("renewal_transport","r_nuc tau_sink","activity=Lambda/(1+Lambda)"),
        ("relative_radius","r/r_ref","controls m=3 or m=4 size sensitivity"),
        ("closed_availability","phi_closed chi_shrink A_closed","required multiplicative inventory"),
        ("zener_ratio","G f_v/(C_Z r_p)","migration increasingly constrained above unity"),
    ]
    _write_csv("closed_pore_dimensionless_groups.csv", [dict(zip(["group_id","definition","interpretation"],r)) for r in groups])
    growth = [
        ("intrinsic","M_GB gamma_GB/G","physical mobility","changes G only"),
        ("pore_Zener","R_Z=C_Z r_p/f_v","geometry-derived with C_Z unresolved","changes migration only"),
        ("junction","Gamma_TJ","state/proxy","changes migration only"),
        ("preparation","Gamma_PR","state/proxy","changes migration only"),
        ("accommodation","Gamma_accommodation","state/proxy","changes migration only"),
    ]
    _write_csv("growth_pinning_law_registry.csv", [dict(zip(["term","law","status","coupling"],r)) for r in growth])


def reconstruct_energy_ledger() -> pd.DataFrame:
    paths = [SOURCE / "final_heating_rate_panel_source.csv", SOURCE / "final_twostep_panel_source.csv"]
    frames = []
    for path in paths:
        if path.exists():
            q = pd.read_csv(path)
            if "run_id" in q:
                chosen = list(q["run_id"].drop_duplicates().astype(str))
                chosen = chosen[:1] + chosen[-1:] if len(chosen) > 1 else chosen
                q = q[q["run_id"].astype(str).isin(chosen)].copy()
            q["source_history"] = path.name
            frames.append(q)
    h = pd.concat(frames, ignore_index=True)
    rows = []
    for (src, run), g in h.groupby(["source_history","run_id"], sort=False):
        g = g.sort_values("physical_time_s").drop_duplicates("physical_time_s").copy()
        t = g.physical_time_s.to_numpy(float)
        rho = g.rho.to_numpy(float); grain = g.G_nm.to_numpy(float) * 1e-9
        phi = (g.phi_open_total.fillna(0)+g.phi_iso_total.fillna(0)+g.phi_closed_total.fillna(0)).to_numpy(float)
        radius = np.maximum(g.pore_D50_nm.fillna(g.pore_D90_nm).fillna(25).to_numpy(float)*0.5e-9, 1e-12)
        Eext = MAT.gamma_s_J_m2*np.array([external_area_density(x) for x in rho])
        Egb = MAT.gamma_GB_J_m2*np.array([gb_area_density(x) for x in grain])
        Epore = MAT.gamma_s_J_m2*3.0*phi/radius
        pgas = 0.25*2*MAT.gamma_s_J_m2/radius
        Egas = pgas*np.maximum(g.phi_closed_total.fillna(0).to_numpy(float),0)
        def release(x):
            return -np.gradient(x,t,edge_order=1) if len(x)>1 else np.zeros_like(x)
        Ps, Pg, Pp = release(Eext), release(Egb), release(Epore)
        Pgas = np.gradient(Egas,t,edge_order=1) if len(Egas)>1 else np.zeros_like(Egas)
        Pav = np.maximum(Ps,0)+np.maximum(Pg,0)+np.maximum(Pp,0)-Pgas
        ro = g.rho_dot_open.fillna(0).to_numpy(float); rc = g.rho_dot_closed.fillna(0).to_numpy(float)
        so = g.sigma_eff_open.fillna(g.sigma_eff).fillna(0).to_numpy(float)
        sc = g.sigma_eff.fillna(0).to_numpy(float)
        Popen=np.maximum(so*ro,0); Pclosed=np.maximum(sc*rc,0)
        Ppr=np.maximum(g.PR_work.fillna(0).to_numpy(float),0)*np.maximum(Pp,0)
        Pgrowth=np.maximum(Pg,0)
        Pdrag=np.maximum(g.G_dot_intrinsic.fillna(0).to_numpy(float)-g.G_dot_actual.fillna(0).to_numpy(float),0)*MAT.gamma_GB_J_m2*MAT.C_GB/np.maximum(grain**2,EPS)
        Psum=Popen+Pclosed+Ppr+Pgrowth+Pdrag+np.maximum(Pgas,0)
        for j,(_,x) in enumerate(g.iterrows()):
            rows.append({"source_history":src,"run_id":run,"physical_time_s":t[j],"T_C":x.T_C,"rho":rho[j],"G_nm":x.G_nm,
                         "area_reconstruction":"spherical_D50_plus_compact_scale_external_proxy",
                         "P_surf_release_W_m3":Ps[j],"P_GB_release_W_m3":Pg[j],"P_pore_release_W_m3":Pp[j],"P_gas_cost_W_m3":Pgas[j],
                         "P_available_W_m3":Pav[j],"P_open_dens_W_m3":Popen[j],"P_closed_dens_W_m3":Pclosed[j],
                         "P_surface_smooth_W_m3":0.0,"P_pore_coarsen_W_m3":Ppr[j],"P_GB_growth_W_m3":Pgrowth[j],"P_drag_W_m3":Pdrag[j],
                         "P_gas_W_m3":max(Pgas[j],0),"P_other_W_m3":0.0,"P_excess_W_m3":max(Pav[j]-Popen[j]-Pclosed[j],0),
                         "Pi_dens":(Popen[j]+Pclosed[j])/max(Pav[j],EPS),"Pi_total":Psum[j]/max(Pav[j],EPS),
                         "ledger_residual_W_m3":Pav[j]-Psum[j],"budget_violation":bool(Psum[j]>Pav[j]+max(1e-9,1e-6*abs(Pav[j]))),
                         "strict_balance_sigma_Pa":sc[j],"capillary_sigma_Pa":capillary_stress(radius[j],0.25*2*MAT.gamma_s_J_m2/radius[j])})
    return _write_csv("energy_ledger_diagnostic_histories.csv", rows)


def write_energy_tables(hist: pd.DataFrame) -> None:
    channels = [
        ("P_open_dens","sigma_open rho_dot_open",True,"densifying work"),
        ("P_closed_dens","sum sigma_closed_i rho_dot_closed_i",True,"densifying work"),
        ("P_surface_smooth","surface diffusion accommodation",False,"shape dissipation"),
        ("P_pore_coarsen","conservative PR/coarsening",False,"topology dissipation"),
        ("P_GB_growth","grain-boundary area release/migration",False,"migration channel"),
        ("P_drag","pore/Zener/junction resistance",False,"resistance channel"),
        ("P_gas","gas compression/pressure work",False,"counter-work"),
        ("P_other","unresolved nonnegative residual channel",False,"never silently rescaled"),
    ]
    _write_csv("energy_ledger_channel_registry.csv", [dict(zip(["channel","definition","changes_density","interpretation"],r)) for r in channels])
    audit = hist.groupby(["source_history","run_id"],as_index=False).agg(
        n_points=("physical_time_s","size"), violation_fraction=("budget_violation","mean"),
        maximum_Pi_dens=("Pi_dens","max"), maximum_Pi_total=("Pi_total","max"),
        minimum_residual_W_m3=("ledger_residual_W_m3","min"))
    audit["accepted_forward_trajectory_changed"] = False
    audit["interpretation"] = "diagnostic reconstruction; violations are recorded, not rescaled"
    audit.to_csv(OUT/"energy_ledger_balance_audit.csv",index=False)
    hist.to_csv(OUT/"energy_ledger_fixed_path_test.csv",index=False)
    closures = [
        ("capillary_force_closure","local pore curvature minus gas",False,True,True,True,False,False,False,True,"primary candidate; diagnostic only"),
        ("bounded_power_consistency_closure","local capillary force plus ledger deficit check",False,True,True,True,True,False,False,True,"primary candidate with explicit deficit"),
        ("strict_GB_area_loss_balance","root of Pdens=GB-area-loss power",True,False,False,False,False,True,True,True,"diagnostic_only"),
    ]
    cols=["closure_id","sigma_source","uses_GB_area_loss_only","uses_surface_area","uses_pore_area","uses_gas_pressure","uses_energy_ledger","changes_density_rate","changes_closed_rate","schedule_label_free","physical_status"]
    _write_csv("stress_closure_comparison.csv",[dict(zip(cols,r)) for r in closures])


def write_growth_audit(hist: pd.DataFrame) -> None:
    g=hist[["source_history","run_id","physical_time_s","rho","G_nm"]].copy()
    g["phi_total"]=1-g.rho
    g["pore_radius_nm"]=12.5
    g["R_Z_nm"]=4*g.pore_radius_nm/(3*np.maximum(g.phi_total,EPS))
    g["migration_changes_density_directly"]=False
    g["status"]="bin data unavailable in exported history; D50-scale reconstruction"
    g.to_csv(OUT/"grain_growth_dissipation_audit.csv",index=False)


def write_equations() -> None:
    text = r"""# Diagnostic energy-ledger equations

This ledger is an audit of existing trajectories; it does not replace the accepted
forward calculation and is not validation.

Per bulk volume, the stored-energy terms are

\[
E_{surf}=\gamma_s A_{surf},\quad E_{GB}=\gamma_{GB}C_{GB}/G,\quad
E_{pore}=\gamma_s\sum_i3\phi_i/r_i,\quad E_{gas}=\sum_iP_{gas,i}\phi_{c,i}.
\]

The external-area term uses an explicitly provisional compact-scale geometry.
Exported histories do not contain pore-bin arrays, so their pore area is reconstructed
from total pore fraction and D50. Candidate-law tests remain bin-resolved.

Signed rates are \(P_j=-dE_j/dt\), except gas cost \(+dE_{gas}/dt\), and

\[
P_{available}=[P_{surf}]_+ +[P_{GB}]_+ +[P_{pore}]_+-P_{gas,cost}.
\]

Named expenditures are open and closed densification, surface smoothing,
conservative pore redistribution, GB growth, drag, gas work, and an explicit
unresolved channel. The audit records rather than repairs violations of
\(\sum P_{spent}\le P_{available}+\epsilon\). It also reports
\(\Pi_{dens}=P_{dens}/P_{available}\), \(\Pi_{total}=P_{spent}/P_{available}\),
and \(P_{excess}=\max(P_{available}-P_{dens},0)\).

The old equality \(P_{dens}=\gamma_{GB}C_{GB}\dot G/G^2\) is retained only as
`strict_GB_area_loss_balance`, a diagnostic ablation. It is too narrow because it
omits pore/external surface release, gas work, smoothing, redistribution, and drag,
and it feeds the migration law back into activation stress.
"""
    (OUT/"energy_ledger_equations.md").write_text(text)


def write_property_comparison() -> None:
    rows=[
        ("Gstar_effective_Q","stress- and temperature-resolved; no scalar input","Delta Q_nuc survival 0 to +50 kJ/mol","diagnostic only","barrier slices start above all tested process temperatures"),
        ("Q_GB",380.0,"fixed 380 kJ/mol","inside","unchanged physical input"),
        ("Q_s",380.0,"fixed 380 kJ/mol","inside","unchanged physical input"),
        ("Q_growth",MAT.Q_M_J_mol/1000,"current fixed mobility law","reference","not refit here"),
        ("Q_closed_app","computed slope only","reduced Delta Q_closed -25 to +100 kJ/mol OAT","not comparable","inventory/activity must be removed before interpretation"),
        ("closed_availability_ratio","actual first-state phi*chi*A is very small","candidate-family threshold not copied","below","inventory/availability gap"),
        ("PR_rate_scale","D_s r^-4 topology law","PR prefactor threshold 0.3x base","unmapped","geometry coefficient missing"),
        ("growth_prefactor","current fixed mobility","0.1x base threshold","not changed","candidate 693168 not copied"),
        ("density_boundary","fixed-state scan result","lower boundary required","tested","see boundary table"),
        ("growth_boundary","fixed-state scan result","upper boundary required","tested","see boundary table"),
    ]
    _write_csv("physical_to_reduced_property_window_comparison.csv",[dict(zip(["quantity","physical_candidate","reduced_works_envelope","comparison","interpretation"],r)) for r in rows])


def main() -> None:
    OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(exist_ok=True)
    write_source_registries(); write_law_registries(); write_equations()
    hist=reconstruct_energy_ledger(); write_energy_tables(hist); write_growth_audit(hist)
    write_property_comparison()
    _write_csv("barrier_guardrail.csv",[{"path":str(BARRIER_PATH.relative_to(ROOT)),"sha256":barrier_sha256(),"unchanged_by_this_analysis":True,
                                         "minimum_fit_temperature_C":float(BARRIER.temperatures_K.min()-273.15),"all_candidate_temperatures_extrapolated":True}])
    print(f"wrote diagnostic registries and ledger to {OUT}")


if __name__ == "__main__":
    main()
