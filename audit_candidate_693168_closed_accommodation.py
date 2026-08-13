#!/usr/bin/env python3
"""Focused reproducibility and mechanism audit for conditional Tier-B candidate 693168.

This module is intentionally an analysis layer.  It does not change the frozen
local-region evolution laws, material kinetics, target, time budget, or decoded
candidate parameters.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math
import time

import numpy as np
import pandas as pd

import interacting_local_region_model as model
import interacting_local_region_decoder as decoder
import interacting_local_region_objectives as objectives
from local_region_decoder_corrected_postprocess import ablation_parameters
import massive_latent_topology_optimizers as optimizers
import interacting_local_region_search as production_search


CANDIDATE_ID = 693168
SOURCE = Path("results/local_region_decoder_corrected_dynamic_search")
OUT = Path("results/audit_candidate_693168_closed_accommodation")
TABLES = OUT / "final_tables"
TARGET = 0.98
SWITCH = 0.88
T1_C = 1400.0
SUCCESS_T2_C = 1100.0
PRODUCTION_T2_C = 1200.0
LOWER_T2_C = 900.0
UPPER_T2_C = 1220.0
HOURS = 500.0
GROWTH_TOLERANCE = 0.20


def write_csv(path: Path, rows) -> None:
    if isinstance(rows, pd.DataFrame):
        frame = rows
    else:
        frame = pd.DataFrame(list(rows))
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def candidate_parameters(candidate_id=CANDIDATE_ID):
    registry = pd.read_csv(SOURCE / "parameter_registry.csv")
    selected = registry[registry.candidate_id == candidate_id]
    if len(selected) != len(decoder.NAMES):
        raise RuntimeError(f"candidate {candidate_id}: expected {len(decoder.NAMES)} decoded parameters, found {len(selected)}")
    registry_decoded = dict(zip(selected.parameter, selected.value))
    # Candidate IDs are row indices in the deterministic 1,000,000-row LHS.
    # Reconstruct from seed/ID so CSV formatting cannot perturb the physics.
    sample = optimizers.latin_hypercube(
        1_000_000, len(decoder.NAMES), production_search.SEED
    )[candidate_id]
    decoded = decoder.decode(sample)
    for key in decoder.NAMES:
        if not np.isclose(float(decoded[key]), float(registry_decoded[key]), rtol=2e-5, atol=1e-10):
            raise RuntimeError(f"registry mismatch for {key}: {registry_decoded[key]} vs {decoded[key]}")
    return {**model.defaults(), **decoded}, decoded


def _weighted_quantile(values, weights, q):
    order = np.argsort(values)
    v, w = np.asarray(values)[order], np.asarray(weights)[order]
    c = np.cumsum(w) / max(np.sum(w), 1e-30)
    return float(np.interp(q, c, v))


def state_row(state, flux, p, t_s, T_C, stage, path_label, step_s=0.0):
    w = state.weights
    gb = float(w @ state.phi_GBseg)
    tj = float(w @ state.phi_TJ)
    iso = float(w @ state.phi_iso)
    closed = float(w @ state.phi_closed)
    pore = max(gb + tj + iso + closed, 1e-30)
    n_gb = float(w @ state.N_GBseg)
    n_tj = float(w @ state.N_TJ)
    n_iso = float(w @ state.N_iso)
    n_closed = float(w @ state.N_closed)
    local_phi = state.phi_GBseg + state.phi_TJ + state.phi_iso + state.phi_closed
    local_n = state.N_GBseg + state.N_TJ + state.N_iso + state.N_closed
    radius_proxy = np.cbrt(np.maximum(local_phi / np.maximum(local_n, 1e-30), 0.0))
    closed_radius_proxy = np.cbrt(np.maximum(state.phi_closed / np.maximum(state.N_closed, 1e-30), 0.0))
    rho_dot_open = float(w @ flux["rho_dot_open"])
    rho_dot_closed = float(w @ flux["rho_dot_closed"])
    g_dot = float(w @ flux["G_dot"])
    pr_rate = float(w @ flux["PR_damage"])
    cap = float(p.get("closed_capacity", 1.0))
    accommodation = float(w @ state.closed_accommodation)
    capacity_tau = max(float(p.get("capacity_tau", 1e30)), 1e-30)
    growth_ref = float(w @ state.G)
    return dict(
        candidate_id=CANDIDATE_ID, path_label=path_label, stage=stage,
        physical_time_s=float(t_s), physical_time_h=float(t_s / 3600.0),
        T_C=float(T_C), rho=1.0 - pore, G_mean_nm=growth_ref,
        G50_nm=_weighted_quantile(state.G, w, .5),
        G90_nm=_weighted_quantile(state.G, w, .9),
        rho_dot=rho_dot_open + rho_dot_closed, G_dot=g_dot,
        growth_fraction_from_switch=np.nan, target_attained=(1.0 - pore >= TARGET),
        classification="", solver_step_s=float(step_s), record_kind="solver_state",
        activity=float(w @ flux["activity"]), renewal_activity=float(w @ flux["activity"]),
        connected_eligibility=float(w @ state.connected_removable_fraction),
        densification_eligibility=float(w @ state.densification_eligibility),
        tau_nuc_s=np.nan, tau_exchange_s=np.nan, tau_transport_s=np.nan,
        sigma_base=np.nan, sigma_act_total=np.nan,
        stress_concentration=float(w @ state.residual_stress),
        phi_open=gb + tj + iso, phi_connected=gb + tj, phi_GBseg=gb,
        phi_TJ=tj, phi_iso=iso, phi_closed=closed,
        N_open=n_gb + n_tj + n_iso, N_connected=n_gb + n_tj,
        N_GBseg=n_gb, N_TJ=n_tj, N_iso=n_iso, N_closed=n_closed,
        pore_radius_proxy_D50=_weighted_quantile(radius_proxy, w, .5),
        pore_radius_proxy_D90=_weighted_quantile(radius_proxy, w, .9),
        closed_radius_proxy_D50=_weighted_quantile(closed_radius_proxy, w, .5),
        closed_radius_proxy_D90=_weighted_quantile(closed_radius_proxy, w, .9),
        D50_nm=np.nan, D90_nm=np.nan,
        large_pore_fraction=float(w @ state.large_attached_fraction),
        connected_fine_pore_fraction=gb / pore,
        removable_pore_fraction=gb / pore,
        closed_fraction=closed / pore, isolated_fraction=iso / pore,
        closed_accommodation_capacity=cap,
        closed_accommodation_available=accommodation,
        closed_accommodation_used=max(cap - accommodation, 0.0),
        closed_accommodation_recovery=max(cap - accommodation, 0.0) / capacity_tau,
        closed_accommodation_factor=np.clip(accommodation / max(cap, 1e-30), 0.0, 1.0),
        P_comp_closed=np.clip(accommodation / max(cap, 1e-30), 0.0, 1.0),
        Lambda_closed=np.nan, K_closed=np.nan, Lambda_over_K_closed=np.nan,
        closed_shrinkage_flux=rho_dot_closed, open_shrinkage_flux=rho_dot_open,
        rho_dot_open=rho_dot_open, rho_dot_closed=rho_dot_closed,
        closed_pore_contribution_to_rho_dot=rho_dot_closed,
        open_pore_contribution_to_rho_dot=rho_dot_open,
        PR_damage_memory=float(w @ state.PR_damage_memory),
        PR_redistribution_rate=pr_rate,
        H_PR=np.nan, H_dens=np.nan, w_PR=np.nan, w_dens=np.nan,
        cumulative_PR_surface_energy_loss=np.nan,
        cumulative_densifying_work=np.nan,
        cumulative_non_densifying_work=np.nan,
        cumulative_PR_redistributed_volume=0.0,
        cumulative_open_pore_removed=0.0,
        cumulative_closed_pore_removed=0.0,
        X_J=float(w @ state.X_J), C_TJ=float(w @ state.C_TJ),
        C_GBseg=float(w @ state.C_GBseg), f_clean_GB=float(w @ state.f_clean_GB),
        migration_factor=float(w @ flux["migration_factor"]),
        pore_drag=float(p.get("attached_drag", 0.0) * w @ (state.large_attached_fraction + state.large_TJ_fraction)),
        persistent_junction_drag=float(p.get("junction_drag", 0.0) * w @ state.X_J),
        P_TJ_multihit=float(w @ flux["P_comp_TJ"]),
        P_persistent_junction_drag=np.nan, P_clean_GB=np.nan,
        Lambda_TJ=float(w @ flux["Lambda_TJ"]), K_TJ=float(w @ flux["K_TJ"]),
        Lambda_over_K_TJ=float(w @ flux["Lambda_TJ"] / max(float(w @ flux["K_TJ"]), 1e-30)),
        P_comp_TJ=float(w @ flux["P_comp_TJ"]), residual_stress=float(w @ state.residual_stress),
    )


def simulate_detailed(p, *, T_C, dt_s, path_label, initial_state=None,
                      time_offset_s=0.0, stop_density=TARGET, stage="single_step",
                      max_steps=120_000):
    state = (model.initial(p["N_regions"], rho0=p.get("rho0", .70),
                           G0=p.get("G0_nm", 100.0), p=p)
             if initial_state is None else model.clone_state(initial_state))
    adjacency = model.network_adjacency(p["N_regions"], p)
    rows = []
    elapsed = 0.0
    stagnant = 0
    previous_rho = None
    cum_pr = cum_open = cum_closed = 0.0
    while elapsed <= HOURS * 3600 and len(rows) < max_steps:
        flux = model.local_fluxes(state, T_C, p)
        row = state_row(state, flux, p, time_offset_s + elapsed, T_C, stage, path_label)
        row.update(cumulative_PR_redistributed_volume=cum_pr,
                   cumulative_open_pore_removed=cum_open,
                   cumulative_closed_pore_removed=cum_closed)
        rows.append(row)
        rho = row["rho"]
        if rho >= stop_density:
            break
        if previous_rho is not None and rho - previous_rho < 1e-7:
            stagnant += 1
        else:
            stagnant = 0
        previous_rho = rho
        if stagnant >= 12:
            break
        remaining = HOURS * 3600 - elapsed
        if remaining <= 0:
            break
        step = min(float(dt_s), remaining)
        density_rate = max(row["rho_dot"], 0.0)
        if density_rate > 0:
            step = min(step, .002 / density_rate)
            gap = max(stop_density - rho, 0.0)
            step = min(step, max(.8 * gap / density_rate, 1e-3))
        row["solver_step_s"] = step
        # Integrate the same bounded local amounts used by model.advance;
        # raw rate*dt would over-count whenever a store is exhausted in a step.
        cum_pr += float(state.weights @ np.minimum(
            np.maximum(flux["PR_damage"], 0.0) * step,
            np.maximum(state.phi_GBseg, 0.0)))
        cum_open += float(state.weights @ np.minimum(
            np.maximum(flux["rho_dot_open"], 0.0) * step,
            np.maximum(state.phi_GBseg, 0.0)))
        cum_closed += float(state.weights @ np.minimum(
            np.maximum(flux["rho_dot_closed"], 0.0) * step,
            np.maximum(state.phi_closed, 0.0)))
        model.advance(state, T_C, p, step, adjacency)
        elapsed += step
        if not all(np.all(np.isfinite(getattr(state, key))) for key in
                   ("rho", "G", "phi_GBseg", "phi_TJ", "phi_iso", "phi_closed", "closed_accommodation")):
            break
    frame = pd.DataFrame(rows)
    return frame, model.clone_state(state)


def prepare_state(p, dt_s):
    return simulate_detailed(p, T_C=T1_C, dt_s=dt_s, path_label="first_step",
                             stop_density=SWITCH, stage="first_step")


def combine_two_step(prep, second, label, classification):
    second = second.iloc[1:].copy() if len(second) else second
    for key in ("cumulative_PR_redistributed_volume", "cumulative_open_pore_removed",
                "cumulative_closed_pore_removed"):
        if len(second):
            second[key] += float(prep[key].iloc[-1])
    frame = pd.concat([prep.copy(), second], ignore_index=True)
    frame["path_label"] = label
    frame["classification"] = classification
    switch_rows = np.flatnonzero(frame.stage.to_numpy() == "second_step")
    if len(switch_rows):
        G1 = frame.iloc[max(int(switch_rows[0]) - 1, 0)].G_mean_nm
        frame["growth_fraction_from_switch"] = frame.G_mean_nm / G1 - 1.0
    return frame


def augment_density_landmarks(frame, landmarks=(.90, .92, .95, .98)):
    """Insert linearly interpolated event rows at requested density crossings."""
    groups=[]
    for label, group in frame.groupby("path_label", sort=False):
        group=group.sort_values("physical_time_s").copy()
        additions=[]
        for target in landmarks:
            if target < group.rho.min() or target > group.rho.max() or np.isclose(group.rho, target, atol=1e-10).any():
                continue
            upper=np.flatnonzero(group.rho.to_numpy() >= target)
            if not len(upper) or upper[0] == 0:
                continue
            j=int(upper[0]);a=group.iloc[j-1];b=group.iloc[j]
            fraction=(target-a.rho)/max(b.rho-a.rho,1e-30)
            row=a.copy()
            for key in group.columns:
                if (pd.api.types.is_numeric_dtype(group[key])
                        and not pd.api.types.is_bool_dtype(group[key])
                        and key not in ("candidate_id",)):
                    row[key]=a[key]+fraction*(b[key]-a[key])
            row["rho"]=target;row["record_kind"]="density_crossing"
            row["target_attained"]=bool(target>=TARGET)
            additions.append(row)
        if additions:
            group=pd.concat([group,pd.DataFrame(additions)],ignore_index=True).sort_values("physical_time_s")
        groups.append(group)
    return pd.concat(groups,ignore_index=True)


def classify(frame, G1):
    attained = float(frame.rho.iloc[-1]) >= TARGET - 1e-6
    growth = float(frame.G_mean_nm.iloc[-1]) / max(G1, 1e-30) - 1.0
    if attained and growth <= GROWTH_TOLERANCE:
        return "SUCCESS"
    if attained:
        return "GRAIN_GROWTH_FAILURE"
    if growth <= GROWTH_TOLERANCE:
        return "DENSIFICATION_EXHAUSTION_FAILURE"
    return "MIXED_FAILURE"


def run_two_step(p, dt_s, T2_C, prepared=None, prep_frame=None, keep_history=True):
    if prepared is None or prep_frame is None:
        prep_frame, prepared = prepare_state(p, dt_s)
    if prep_frame.rho.iloc[-1] < SWITCH - 1e-4:
        return prep_frame, "UNATTAINABLE_FIRST_STEP"
    G1 = float(prep_frame.G_mean_nm.iloc[-1])
    second, _ = simulate_detailed(
        p, T_C=T2_C, dt_s=dt_s, path_label=f"T2_{T2_C:g}C",
        initial_state=prepared, time_offset_s=float(prep_frame.physical_time_s.iloc[-1]),
        stop_density=TARGET, stage="second_step",
    )
    c = classify(second, G1)
    return combine_two_step(prep_frame, second, f"T2_{T2_C:g}C", c), c


def matched_density(high, two, lo=.90, hi=.99, spacing=1e-4):
    top = min(float(high.rho.max()), float(two.rho.max()), hi)
    if top < lo:
        return pd.DataFrame()
    grid = np.arange(lo, top + spacing / 2, spacing)
    gh = np.interp(grid, high.rho, high.G_mean_nm)
    gt = np.interp(grid, two.rho, two.G_mean_nm)
    return pd.DataFrame(dict(
        rho=grid, G_highT_nm=gh, G_two_step_nm=gt,
        reduction_TS=1.0 - gt / gh, ratio_TS=gh / gt,
        highT_interpolation_supported=(grid >= high.rho.min()) & (grid <= high.rho.max()),
        two_step_interpolation_supported=(grid >= two.rho.min()) & (grid <= two.rho.max()),
        both_paths_attained=True,
    ))


def window_from_points(points):
    success = points[points.classification == "SUCCESS"].T2_C.to_numpy()
    lower = bool((points.classification == "DENSIFICATION_EXHAUSTION_FAILURE").any())
    upper = bool((points.classification == "GRAIN_GROWTH_FAILURE").any())
    width = float(success.max() - success.min()) if len(success) >= 2 else 0.0
    return dict(first_success_C=float(success.min()) if len(success) else np.nan,
                last_success_C=float(success.max()) if len(success) else np.nan,
                lower_boundary_C=float(success.min()) if len(success) and lower else np.nan,
                upper_boundary_C=float(success.max()) if len(success) and upper else np.nan,
                window_width_C=width, lower_bracketed=lower, upper_bracketed=upper,
                complete=bool(width >= 25 and lower and upper))


def scan_T2(p, dt_s, temperatures, prepared=None, prep_frame=None, keep=None):
    if prepared is None or prep_frame is None:
        prep_frame, prepared = prepare_state(p, dt_s)
    G1 = float(prep_frame.G_mean_nm.iloc[-1])
    rows, histories = [], {}
    for T2 in temperatures:
        path, c = run_two_step(p, dt_s, float(T2), prepared, prep_frame)
        second = path[path.stage == "second_step"]
        final = path.iloc[-1]
        rows.append(dict(
            candidate_id=CANDIDATE_ID, T2_C=float(T2), classification=c,
            rho2=float(final.rho), G2_nm=float(final.G_mean_nm),
            growth_fraction=float(final.G_mean_nm / G1 - 1.0),
            target_attained=bool(final.rho >= TARGET - 1e-6),
            closed_shrinkage_contribution=float(np.trapezoid(second.closed_shrinkage_flux, second.physical_time_s)) if len(second)>1 else 0.0,
            open_shrinkage_contribution=float(np.trapezoid(second.open_shrinkage_flux, second.physical_time_s)) if len(second)>1 else 0.0,
            closed_accommodation_state=float(final.closed_accommodation_available),
            closed_accommodation_factor=float(final.closed_accommodation_factor),
            PR_damage_state=float(final.PR_damage_memory),
            numerical_censored=bool(not np.isfinite(final.rho)),
        ))
        if keep is not None and float(T2) in keep:
            histories[float(T2)] = path
    points = pd.DataFrame(rows)
    boundary = window_from_points(points)
    points["lower_boundary_marker"] = np.isclose(points.T2_C, boundary["first_success_C"], equal_nan=False)
    points["upper_boundary_marker"] = np.isclose(points.T2_C, boundary["last_success_C"], equal_nan=False)
    return points, boundary, histories, prep_frame, prepared


def score_histories(high, two):
    score = objectives.trajectory_score(
        {"rho": high.rho.to_numpy(), "G_nm": high.G_mean_nm.to_numpy()},
        {"rho": two.rho.to_numpy(), "G_nm": two.G_mean_nm.to_numpy()},
    )
    return score


def timestep_audit(p):
    rows = []
    configurations = [("original_exact", 1800), ("max_30_min", 1800),
                      ("max_15_min", 900), ("max_5_min", 300),
                      ("strict_density_increment", 300)]
    for label, dt in configurations:
        prep, prepared = prepare_state(p, dt)
        high, _ = simulate_detailed(p, T_C=T1_C, dt_s=dt, path_label="highT_reference")
        two, _ = run_two_step(p, dt, PRODUCTION_T2_C, prepared, prep)
        points, boundary, _, _, _ = scan_T2(p, dt, range(800, 1301, 10), prepared, prep)
        score = score_histories(high, two)
        switch = prep.iloc[-1]
        rows.append(dict(
            run_label=label, max_timestep_min=dt/60, exact_switch_density=switch.rho,
            switch_time_h=switch.physical_time_h, T_at_switch_C=T1_C,
            G1_nm=switch.G_mean_nm,
            first_step_growth_fraction=switch.G_mean_nm/prep.G_mean_nm.iloc[0]-1,
            closed_fraction_at_switch=switch.closed_fraction,
            closed_accommodation_at_switch=switch.closed_accommodation_available,
            success_T2_interval=f"{boundary['first_success_C']:g}-{boundary['last_success_C']:g}" if boundary["complete"] else "",
            **boundary, median_reduction=score["median_reduction"],
            minimum_reduction=score["min_reduction"], maximum_reduction=score["max_reduction"],
            span_above_20pct=score["span20"], both_path_attainment=score["attained"],
            numerical_censor=bool(points.numerical_censored.any()),
        ))
    return pd.DataFrame(rows)


def ablation_audit(base):
    names = ["full_model", "no_PR_damage", "no_closed_transition", "no_closed_shrinkage",
             "infinite_closed_accommodation", "no_sweep_coalescence", "no_TJ_multihit",
             "no_residual_stress", "no_attached_pore_drag", "no_persistent_junction",
             "no_local_region_heterogeneity", "topology_disabled"]
    rows = []
    for name in names:
        p = base if name == "full_model" else ablation_parameters(base, name)
        prep, prepared = prepare_state(p, 1800)
        high, _ = simulate_detailed(p, T_C=T1_C, dt_s=1800, path_label="highT_reference")
        two, _ = run_two_step(p, 1800, SUCCESS_T2_C, prepared, prep)
        score = score_histories(high, two)
        points, boundary, _, _, _ = scan_T2(p, 1800, range(800, 1301, 10), prepared, prep)
        final_two = two.iloc[-1]
        switch = prep.iloc[-1]
        required_attainment = bool(score["attained"])
        rows.append(dict(
            candidate_id=CANDIDATE_ID, ablation=name,
            median_reduction=score["median_reduction"], min_reduction=score["min_reduction"],
            max_reduction=score["max_reduction"], span20=score["span20"],
            high_density_attainment=required_attainment, complete_Chen_window=boundary["complete"],
            **boundary, closed_fraction_at_switch=switch.closed_fraction,
            closed_accommodation_at_switch=switch.closed_accommodation_available,
            required_for_high_density_attainment=False, required_for_lower_boundary=False,
            required_for_upper_boundary=False, required_for_trajectory_reduction=False,
            fast_firing_preserved=np.nan,
            loss_reason="preserved" if boundary["complete"] and score["span20"] >= .02 else "joint trajectory/window criteria lost",
        ))
    frame = pd.DataFrame(rows)
    full = frame.iloc[0]
    mask = frame.ablation != "full_model"
    frame.loc[mask, "required_for_high_density_attainment"] = frame.loc[mask, "high_density_attainment"] < bool(full.high_density_attainment)
    frame.loc[mask, "required_for_lower_boundary"] = frame.loc[mask, "lower_bracketed"] < bool(full.lower_bracketed)
    frame.loc[mask, "required_for_upper_boundary"] = frame.loc[mask, "upper_bracketed"] < bool(full.upper_bracketed)
    frame.loc[mask, "required_for_trajectory_reduction"] = frame.loc[mask, "span20"] < float(full.span20) - 1e-8
    return frame


def extended_robustness(base):
    rows = []
    for rho0 in (.60, .65, .70, .75, .80):
        for G0 in (50., 75., 100., 150., 225., 300.):
            p = {**base, "rho0": rho0, "G0_nm": G0}
            prep, prepared = prepare_state(p, 1800)
            high, _ = simulate_detailed(p, T_C=T1_C, dt_s=1800, path_label="highT_reference")
            two, _ = run_two_step(p, 1800, SUCCESS_T2_C, prepared, prep)
            score = score_histories(high, two)
            points, boundary, _, _, _ = scan_T2(p, 1800, range(800, 1301, 20), prepared, prep)
            switch = prep.iloc[-1]
            rows.append(dict(candidate_id=CANDIDATE_ID, rho0=rho0, G0_nm=G0,
                             high_density_attainment=score["attained"],
                             median_reduction=score["median_reduction"], span20=score["span20"],
                             window_width_C=boundary["window_width_C"], complete=boundary["complete"],
                             first_step_growth=switch.G_mean_nm/G0-1,
                             closed_fraction_at_switch=switch.closed_fraction,
                             classification="Tier_B_like" if score["span20"]>=.02 and boundary["complete"] else "outside_joint_criteria"))
    return pd.DataFrame(rows)


def six_candidate_comparison():
    accepted = pd.read_csv(SOURCE / "accepted_tier_candidates.csv")
    ablation = pd.read_csv(SOURCE / "ablation_summary.csv")
    robust = pd.read_csv(SOURCE / "robustness_summary.csv")
    rows=[]
    for _, row in accepted.iterrows():
        cid=int(row.candidate_id); abl=ablation[ablation.candidate_id==cid]
        params=pd.read_csv(SOURCE/"parameter_registry.csv")
        p=dict(zip(params[params.candidate_id==cid].parameter,params[params.candidate_id==cid].value))
        losses=abl[abl.loss_from_full].ablation.tolist()
        rc=int(robust[(robust.candidate_id==cid)&(robust.complete)&(robust.attained)].shape[0])
        plausibility=("extreme_but_informative" if row.median_reduction>.75 or row.closed_fraction_at_switch>.9
                      else "plausible_Tier_B_prototype")
        rows.append(dict(candidate_id=cid, first_step_growth=row.first_step_growth_fraction,
                         closed_fraction_at_switch=row.closed_fraction_at_switch,
                         median_reduction=row.median_reduction,min_reduction=row.min_reduction,
                         max_reduction=row.max_reduction,span20=row.span20,
                         window_width_C=row.window_width_C,
                         T2_success_range="see candidate classification table",
                         closed_capacity=p.get("closed_capacity",np.nan),
                         k_closed=p.get("k_closed",np.nan),k_PR=p.get("k_PR",np.nan),
                         ablation_losses=";".join(losses),robustness_cases_passed=rc,
                         physical_plausibility_flag=plausibility))
    return pd.DataFrame(rows)


def dense_and_fine(base):
    prep, prepared = prepare_state(base, 900)
    high, _ = simulate_detailed(base, T_C=T1_C, dt_s=900, path_label="highT_reference")
    paths = {"highT_reference": high}
    for T2, label in ((LOWER_T2_C,"lower_failure"),(SUCCESS_T2_C,"success"),(UPPER_T2_C,"upper_failure")):
        path, _ = run_two_step(base, 900, T2, prepared, prep)
        path["path_label"] = label
        paths[label]=path
    dense=augment_density_landmarks(pd.concat(paths.values(),ignore_index=True))
    ratio=matched_density(high,paths["success"])
    points,boundary,scan_histories,_,_=scan_T2(base,1800,range(800,1301,5),keep={900.,1100.,1220.})
    scan_dense=[]
    for T2,h in scan_histories.items():
        x=h.copy();x["scan_T2_C"]=T2;scan_dense.append(x)
    return dense,ratio,points,pd.DataFrame([boundary]),pd.concat(scan_dense,ignore_index=True)


def fast_firing_table():
    frame=pd.read_csv(SOURCE/"fast_firing_preservation.csv")
    out=[]
    for _,r in frame.iterrows():
        for mode,meaningful in (("full",r.full_meaningful),("nucleation_facile",r.nucleation_facile_meaningful),("PR_off",r.PR_off_meaningful)):
            out.append(dict(material_id=r.material_id,mode=mode,G_ref_over_G_fast=r.full_max_ratio if mode!="nucleation_facile" else 1.0,
                            density_span_ge_1p5=r.full_span_ge_1p5 if meaningful else 0.0,
                            fast_firing_retained=bool(meaningful),source=r.source,
                            note="archived frozen material envelope; local-region audit does not modify MaterialKinetics"))
    return pd.DataFrame(out)


def parameter_table(decoded, full):
    rows=[]
    for k,v in full.items():
        rows.append(dict(candidate_id=CANDIDATE_ID,parameter=k,value=v,
                         source="decoded_frozen_registry" if k in decoded else "frozen_model_default",
                         in_dynamic_fingerprint=k in decoder.FINGERPRINT))
    return pd.DataFrame(rows)


def main():
    start=time.time(); OUT.mkdir(parents=True,exist_ok=True);TABLES.mkdir(parents=True,exist_ok=True)
    base,decoded=candidate_parameters()
    parameter_table(decoded,base).to_csv(OUT/"candidate_693168_parameter_vector.csv",index=False)
    timestep=timestep_audit(base);write_csv(OUT/"candidate_693168_reproduction_summary.csv",timestep)
    dense,ratio,points,boundaries,scan_hist=dense_and_fine(base)
    write_csv(OUT/"candidate_693168_highT_vs_twostep_dense_histories.csv",dense[dense.path_label.isin(["highT_reference","success"])])
    write_csv(OUT/"candidate_693168_ratio_curve.csv",ratio)
    write_csv(OUT/"candidate_693168_closed_accommodation_history.csv",dense)
    write_csv(OUT/"candidate_693168_T2_classification_points_fine.csv",points)
    write_csv(OUT/"candidate_693168_T2_window_boundaries_fine.csv",boundaries)
    write_csv(OUT/"dense_candidate_693168_histories.csv",dense)
    write_csv(OUT/"dense_candidate_693168_matched_density_curves.csv",ratio)
    write_csv(OUT/"dense_candidate_693168_T2_scan_histories.csv",scan_hist)
    ablation=ablation_audit(base);write_csv(OUT/"candidate_693168_ablation_audit.csv",ablation)
    robust=extended_robustness(base);write_csv(OUT/"candidate_693168_extended_robustness.csv",robust)
    six=six_candidate_comparison();write_csv(OUT/"six_tierB_candidate_comparison.csv",six)
    fast=fast_firing_table();write_csv(OUT/"fast_firing_preservation_audit.csv",fast)
    # Required final-table aliases are intentionally compact and auditable.
    index=dense.groupby(["path_label","stage"]).agg(rows=("rho","size"),time_start_s=("physical_time_s","min"),time_end_s=("physical_time_s","max"),rho_start=("rho","min"),rho_end=("rho","max")).reset_index()
    aliases={
      "candidate_693168_dense_history_index.csv":index,
      "candidate_693168_matched_density_ratio_curve.csv":ratio,
      "candidate_693168_T2_classification_fine.csv":points,
      "candidate_693168_Chen_boundaries_fine.csv":boundaries,
      "candidate_693168_closed_accommodation_summary.csv":dense.groupby("path_label").agg(closed_fraction_switch=("closed_fraction","max"),closed_fraction_final=("closed_fraction","last"),closed_removed=("cumulative_closed_pore_removed","last"),open_removed=("cumulative_open_pore_removed","last"),accommodation_min=("closed_accommodation_available","min")).reset_index(),
      "candidate_693168_ablation_summary_final.csv":ablation,
      "six_TierB_candidate_comparison_final.csv":six,
      "candidate_693168_robustness_summary_final.csv":robust,
    }
    for name,frame in aliases.items():write_csv(TABLES/name,frame)
    state=dict(candidate_id=CANDIDATE_ID,status="audit_calculations_complete",runtime_s=time.time()-start,
               frozen_model_files_changed=False,dense_history_rows=len(dense),fine_T2_points=len(points),
               exact_parameter_fingerprint=decoder.fingerprint(decoded),target=TARGET,switch=SWITCH)
    (OUT/"audit_run_state.json").write_text(json.dumps(state,indent=2)+"\n")
    print(json.dumps(state,indent=2))


if __name__ == "__main__":
    main()
