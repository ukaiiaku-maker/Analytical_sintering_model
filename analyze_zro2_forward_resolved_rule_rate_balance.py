#!/usr/bin/env python3
"""Read-only physics audit of resolved-rule rate balances; not validation."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pandas as pd

from zro2_forward.conditioned_950c import ConditionedTwoStep, run_path
from zro2_forward.resolved_rules import ResolvedRuleModel
from zro2_forward.schedules import RampNoHold
from run_zro2_forward_resolved_rules import params, initial, TARGET_RHO, TARGET_G

OUT = Path("results/zro2_forward_resolved_rule_rate_balance_audit")


def integ(x, column):
    return float(np.trapezoid(x[column].fillna(0), x.t_s)) if len(x) > 1 else 0.0


def tau_d90(row):
    values = np.asarray(json.loads(row.tau_remove_s_json), float)
    radii = np.asarray(json.loads(row.pore_radii_m_json), float)
    finite = np.isfinite(values)
    return float(np.interp(row.pore_D90_m / 2, radii[finite], values[finite])) if finite.any() else np.inf


def normalized(frame, case, source="resolved"):
    x = frame.copy()
    x["audit_case"] = case
    x["source_layer"] = source
    x["total_rho_dot_sinv"] = x.rho_dot_open_sinv.fillna(0) + x.rho_dot_closed_sinv.fillna(0)
    x["PR_to_precursor_flux"] = x.get("PR_to_isolated_flux", x.get("isolation_rate", 0))
    x["precursor_to_closed_flux"] = x.get("PR_to_closed_precursor_flux", x.get("closure_rate", 0))
    x["phi_open"] = x.open_fraction
    x["phi_precursor"] = x.isolated_fraction
    x["phi_closed"] = x.closed_fraction
    x["closed_accommodation_used"] = (1 - x.A_closed) * params().accommodation_capacity_phi
    x["closed_accommodation_capacity"] = params().accommodation_capacity_phi
    x["P_surf"] = x.P_surf_W_m3
    x["P_dens"] = x.P_dens_W_m3
    x["P_excess"] = x.P_excess_W_m3
    x["G_dot_intrinsic"] = x.get("G_dot_intrinsic_m_s", x.G_dot_clean_m_s)
    x["G_dot_actual"] = x.get("G_dot_actual_m_s", x.G_dot_m_s)
    x["tau_remove_D90_s"] = x.apply(tau_d90, axis=1)
    return x


def run_resolved(label, q, path, dt=60, record=60):
    _, h = run_path(ResolvedRuleModel(parameters=q), path, initial(), dt, record, label)
    return normalized(h, label)


def summary(x):
    rows = []
    for case, g in x.groupby("audit_case"):
        last = g.iloc[-1]
        rows.append({"audit_case": case, "source_layer": last.source_layer, "final_rho": last.rho,
                     "final_G_um": last.G_um, "Delta_rho_open": integ(g, "rho_dot_open_sinv"),
                     "Delta_rho_closed": integ(g, "rho_dot_closed_sinv"),
                     "Delta_rho_total_flux": integ(g, "total_rho_dot_sinv"),
                     "Delta_rho_state": last.rho-g.iloc[0].rho,
                     "cumulative_PR_coarsening": integ(g, "PR_coarsening_flux"),
                     "cumulative_closure": integ(g, "closure_rate"), "final_phi_open": last.phi_open,
                     "final_phi_precursor": last.phi_precursor, "final_phi_closed": last.phi_closed,
                     "final_A_closed": last.A_closed, "final_Gamma_migration": last.Gamma_migration,
                     "non_validation": True})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Fixed diagnostic paths only: no Chen or mobility search.
    cases = []
    cases.append(run_resolved("resolved_default", params(), ConditionedTwoStep(1400, 1100, .85, 20), 120, 120))
    cases.append(run_resolved("resolved_M0_factor_0p3", params(gb_mobility_mode="bounded_uncertainty_factor", M0_factor=.3), ConditionedTwoStep(1400, 1100, .85, 20), 120, 120))
    for qv in (500, 550):
        base = params()
        factor = np.exp((qv*1000-base.Q_M_J_mol_override if base.Q_M_J_mol_override else qv*1000-405238.395)/(8.314462618*1773.15))
        cases.append(run_resolved(f"resolved_QM_{qv}", params(gb_mobility_mode="activation_energy_envelope", M0_factor=factor, Q_M_J_mol_override=qv*1000), ConditionedTwoStep(1400, 1100, .85, 20), 120, 120))
    baseline = pd.read_csv("results/zro2_forward_pdf_conditioned_950C_comparison/pdf_conditioned_dense_histories.csv")
    baseline = baseline[baseline.case.eq("conditioned_50C_min")].copy()
    for c, default in {"tau_nuc_s":np.nan,"tau_exchange_s":np.nan,"tau_transport_s":np.nan,"tau_cycle_s":np.nan,
                       "PR_coarsening_flux":baseline.get("bin_crossing_rate",0),"PR_to_isolated_flux":baseline.get("isolation_rate",0),
                       "PR_to_closed_precursor_flux":baseline.get("closure_rate",0),"PR_memory":0,"Gamma_migration":baseline.Gamma_growth,
                       "G_dot_intrinsic_m_s":baseline.G_dot_clean_m_s,"G_dot_actual_m_s":baseline.G_dot_m_s,
                       "open_path_eligibility":np.nan}.items(): baseline[c] = default
    baseline = normalized(baseline, "previous_PDF_conditioned_50C", "previous_pdf_conditioned")
    resolved5 = run_resolved("resolved_5C", params(), RampNoHold(5,1500,start_C=950), 60, 60)
    resolved50 = run_resolved("resolved_50C", params(), RampNoHold(50,1500,start_C=950), 30, 30)
    trajectories = pd.concat([baseline, resolved5, resolved50] + cases, ignore_index=True)
    keep = ["audit_case","source_layer","t_s","T_C","rho","G_um","rho_dot_open_sinv","rho_dot_closed_sinv","total_rho_dot_sinv",
            "PR_coarsening_flux","PR_to_precursor_flux","precursor_to_closed_flux","closure_rate","closed_shrinkage_flux","phi_open","phi_precursor","phi_closed",
            "A_closed","closed_accommodation_used","closed_accommodation_capacity","activity","tau_nuc_s","tau_exchange_s","tau_transport_s","tau_cycle_s",
            "sigma_eff_Pa","P_surf","P_dens","P_excess","Gamma_migration","G_dot_intrinsic","G_dot_actual","pore_D90_m","fine_pore_fraction",
            "connected_removable_factor","open_path_eligibility","tau_remove_D90_s"]
    for c in keep:
        if c not in trajectories: trajectories[c] = np.nan
    trajectories[keep].to_csv(OUT/"density_flux_comparison.csv", index=False)
    trajectories[[c for c in keep if c in {"audit_case","source_layer","t_s","T_C","rho","G_um","activity","tau_nuc_s","tau_exchange_s","tau_transport_s","tau_cycle_s","sigma_eff_Pa","P_surf","P_dens","P_excess","Gamma_migration","G_dot_intrinsic","G_dot_actual"}]].to_csv(OUT/"trajectory_comparison_summary.csv", index=False)
    trajectories[["audit_case","t_s","T_C","rho","phi_open","phi_precursor","phi_closed","PR_coarsening_flux","PR_to_precursor_flux","precursor_to_closed_flux","closure_rate"]].to_csv(OUT/"pore_store_comparison.csv", index=False)
    trajectories[["audit_case","t_s","T_C","rho","A_closed","closed_accommodation_used","closed_accommodation_capacity","closed_shrinkage_flux","Gamma_migration"]].to_csv(OUT/"accommodation_comparison.csv", index=False)

    # Fast-rate loss decomposition.
    fast = pd.concat([baseline, resolved50], ignore_index=True)
    fs = summary(fast)
    b, r = fs.iloc[0], fs.iloc[1]
    collapse = resolved50[resolved50.open_path_eligibility < .10]
    collapse_row = collapse.iloc[0] if len(collapse) else resolved50.iloc[-1]
    loss = pd.DataFrame([{
        "previous_final_rho": b.final_rho, "resolved_final_rho": r.final_rho, "density_loss": b.final_rho-r.final_rho,
        "open_shrinkage_change": r.Delta_rho_open-b.Delta_rho_open, "density_shifted_to_closed_store": r.final_phi_closed-b.final_phi_closed,
        "closed_shrinkage_attained": r.Delta_rho_closed, "PR_coarsening_integral": r.cumulative_PR_coarsening,
        "closure_integral": r.cumulative_closure, "previous_final_D90_m": baseline.iloc[-1].pore_D90_m,
        "resolved_final_D90_m": resolved50.iloc[-1].pore_D90_m, "previous_final_fine_fraction": baseline.iloc[-1].fine_pore_fraction,
        "resolved_final_fine_fraction": resolved50.iloc[-1].fine_pore_fraction, "previous_final_tau_remove_D90_s": baseline.iloc[-1].tau_remove_D90_s,
        "resolved_final_tau_remove_D90_s": resolved50.iloc[-1].tau_remove_D90_s, "resolved_final_connected_removable_factor": resolved50.iloc[-1].connected_removable_factor,
        "open_path_collapse_T_C": collapse_row.T_C, "open_path_collapse_t_s": collapse_row.t_s,
        "diagnosis":"open_shrinkage_rate_reduced; direct_transfer_too_small_to_explain_loss; closed_shrinkage_too_slow", "non_validation":True}])
    loss.to_csv(OUT/"fast_rate_density_loss_decomposition.csv", index=False)

    # Fixed representative boundary paths.
    boundary=[]; growth=[]
    for tag,t2 in (("low_T2",900),("near_best",1200),("high_T2",1300)):
        h=run_resolved(tag,params(),ConditionedTwoStep(1400,t2,.8,40),120,120)
        last=h.iloc[-1]; dr_open=integ(h,"rho_dot_open_sinv"); dr_closed=integ(h,"rho_dot_closed_sinv")
        reason = "mixed"
        if dr_closed < .02 and last.phi_closed > .005: reason="closed_shrinkage_too_slow"
        if last.open_path_eligibility < .15 and last.phi_open > .01: reason="closure_too_early"
        if last.activity < .02: reason="barrier_activity_too_low"
        boundary.append({"path":tag,"T2_C":t2,"Delta_rho_open":dr_open,"Delta_rho_closed":dr_closed,"Delta_rho_total":dr_open+dr_closed,
                         "missed_density_target":max(0,TARGET_RHO-last.rho),"cumulative_PR_coarsening":integ(h,"PR_coarsening_flux"),
                         "cumulative_closure":integ(h,"closure_rate"),"final_A_closed":last.A_closed,"final_phi_closed":last.phi_closed,
                         "final_open_pore_volume":last.phi_open,"final_tau_remove_D90_s":last.tau_remove_D90_s,"failure_classification":reason,"non_validation":True})
        intrinsic=integ(h,"G_dot_intrinsic"); actual=integ(h,"G_dot_actual")
        growth.append({"path":tag,"T2_C":t2,"G_dot_intrinsic_integrated_m":intrinsic,"G_dot_actual_integrated_m":actual,
                       "migration_activity_factor_integrated_ratio":actual/max(intrinsic,1e-30),"Zener_pore_drag_mean":h.pore_Zener_drag_contribution.mean(),
                       "closed_accommodation_migration_factor_mean":h.closed_accommodation_migration_contribution.mean(),
                       "PR_memory_migration_factor_mean":(1/(1+params().PR_migration_drag*h.PR_memory)).mean(),
                       "event_completion_factor_mean":(h.Gamma_migration/(h.pore_Zener_drag_contribution*h.closed_accommodation_migration_contribution/(1+params().PR_migration_drag*h.PR_memory))).replace([np.inf,-np.inf],np.nan).mean(),
                       "final_G_um":last.G_um,"growth_fraction":last.G_um/h.iloc[0].G_um-1,"T_upper_growth_C":t2 if last.G_um>TARGET_G else np.nan,
                       "failure_control":"intrinsic_growth_dominates_at_density_attaining_temperature" if last.G_um>TARGET_G else "density_not_attained_before_growth_limit","non_validation":True})
    pd.DataFrame(boundary).to_csv(OUT/"lower_boundary_rate_balance.csv",index=False)
    pd.DataFrame(growth).to_csv(OUT/"upper_boundary_growth_balance.csv",index=False)
    pd.DataFrame(boundary).merge(pd.DataFrame(growth),on=["path","T2_C","non_validation"]).to_csv(OUT/"boundary_failure_comparison.csv",index=False)

    oldab=pd.read_csv("results/zro2_forward_resolved_rules/resolved_rule_ablation_summary.csv")
    oldab["parent_window_present"]=False; oldab["ablation_interpretable"]=False
    oldab["reason_not_interpretable"]="not_interpretable_parent_has_no_window"
    oldab["ablation_result"]="not_interpretable_parent_has_no_window"
    oldab.to_csv(OUT/"resolved_rule_ablation_summary_reinterpreted.csv",index=False)

    cand=pd.read_csv("results/local_region_decoder_corrected_dynamic_search/local_region_state_histories_compact.csv").query("candidate_id==693168 and path=='two_step'")
    switch=cand.loc[cand.T_C.lt(1400)].iloc[0] if cand.T_C.lt(1400).any() else cand.iloc[-1]
    rb=pd.read_csv("results/zro2_forward_resolved_rules/resolved_rule_chen_classification_points.csv")
    comp=pd.DataFrame([
        {"model":"resolved_rule","closed_fraction_at_switch":rb.closed_fraction_at_switch.max(),"A_closed_at_switch":rb.A_closed_at_switch.max(),"PR_memory_at_switch":rb.PR_memory_at_switch.max(),"first_step_growth":rb.G1_um.max()/0.05-1,"high_density_reduction":.016589,"Chen_window_width_C":0,"lower_boundary_C":np.nan,"upper_boundary_C":np.nan,"calibrated":False},
        {"model":"candidate_693168","closed_fraction_at_switch":switch.closed,"A_closed_at_switch":switch.closed_accommodation,"PR_memory_at_switch":switch.PR_memory,"first_step_growth":.137068,"high_density_reduction":.893813,"Chen_window_width_C":275,"lower_boundary_C":925,"upper_boundary_C":1200,"calibrated":False}])
    comp.to_csv(OUT/"resolved_vs_candidate693168_state_comparison.csv",index=False)
    summary(trajectories).to_csv(OUT/"trajectory_integral_summary.csv",index=False)
    print(loss.to_string(index=False)); print(pd.DataFrame(boundary).to_string(index=False))

if __name__ == "__main__": main()
