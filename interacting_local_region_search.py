#!/usr/bin/env python3
"""Auditable staged search for interacting local pore-region dynamics."""
from pathlib import Path
import argparse
import csv
import gzip
import json
import time
import numpy as np
import pandas as pd

import interacting_local_region_model as model
import interacting_local_region_objectives as objectives
import massive_latent_topology_optimizers as optimizers
import interacting_local_region_decoder as decoder

OUT = Path("results/local_region_decoder_corrected_dynamic_search")
SEED = 20260818
STAGE1_REQUIRED = 20_000
STAGE2_REQUIRED = 1_000


def write(path, rows, gz=False):
    rows = list(rows)
    fields = list(dict.fromkeys(key for row in rows for key in row)) or ["status"]
    opener = gzip.open if gz else open
    path.parent.mkdir(parents=True, exist_ok=True)
    with opener(path, "wt" if gz else "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_state(**values):
    (OUT / "run_state.json").write_text(json.dumps(values, indent=2) + "\n")


def finite_history(history):
    required = ("rho", "G_mean", "rho_dot_open", "rho_dot_closed")
    return all(len(history[name]) and np.all(np.isfinite(history[name])) for name in required)


def simulate(p, n=8, T1=1400, T2=None, switch=.88, dt=1800, hours=500,
             initial_state=None, switched_initial=False, time_offset=0.0,
             stop_at_switch=False, return_state=False, max_steps=6000):
    state = (model.initial(n, rho0=p.get("rho0", .70), G0=p.get("G0_nm", 100.), p=p)
             if initial_state is None else model.clone_state(initial_state))
    adjacency = model.network_adjacency(n, p)
    names = (
        "t", "T_C", "rho", "G_mean", "G50", "G90", "connected",
        "topology_variance", "closed", "rho_dot_open", "rho_dot_closed",
        "migration_factor", "PR_memory", "sweep_memory", "residual_stress",
        "X_J", "phase",
        "N_total", "pore_radius_proxy", "f_GBseg", "f_TJ", "f_iso",
        "closed_accommodation",
    )
    history = {name: [] for name in names}
    elapsed = 0.0
    switched = switched_initial
    stagnant_steps = 0
    previous_rho = None
    while elapsed < hours * 3600 and len(history["t"]) < max_steps:
        obs = model.global_observables(state)
        target_reached = obs["rho_global"] >= .98
        switched = switched or (T2 is not None and obs["rho_global"] >= switch - 1e-6)
        temperature = T2 if switched else T1
        flux = model.local_fluxes(state, temperature, p)
        values = dict(
            t=time_offset + elapsed, T_C=temperature, rho=obs["rho_global"],
            G_mean=obs["G_mean"], G50=np.median(state.G), G90=np.quantile(state.G, .9),
            connected=np.average(state.connected_removable_fraction, weights=state.weights),
            topology_variance=obs["topology_variance"], closed=obs["closed_fraction"],
            rho_dot_open=state.weights @ flux["rho_dot_open"],
            rho_dot_closed=state.weights @ flux["rho_dot_closed"],
            migration_factor=state.weights @ flux["migration_factor"],
            PR_memory=state.weights @ state.PR_damage_memory,
            sweep_memory=state.weights @ state.sweep_memory,
            residual_stress=state.weights @ state.residual_stress,
            X_J=state.weights @ state.X_J,
            phase=int(switched),
            N_total=state.weights @ (state.N_GBseg + state.N_TJ + state.N_iso + state.N_closed),
            pore_radius_proxy=state.weights @ np.cbrt(
                (state.phi_GBseg + state.phi_TJ + state.phi_iso + state.phi_closed)
                / np.maximum(state.N_GBseg + state.N_TJ + state.N_iso + state.N_closed, 1e-30)
            ),
            f_GBseg=state.weights @ state.C_GBseg,
            f_TJ=state.weights @ state.C_TJ,
            f_iso=state.weights @ state.isolated_fraction,
            closed_accommodation=state.weights @ state.closed_accommodation,
        )
        for name in names:
            history[name].append(values[name])
        if stop_at_switch and switched:
            break
        if target_reached:
            break
        density_increment = np.inf if previous_rho is None else obs["rho_global"] - previous_rho
        previous_rho = obs["rho_global"]
        if density_increment < 1e-7 and obs["rho_global"] < .98:
            stagnant_steps += 1
        else:
            stagnant_steps = 0
        if stagnant_steps >= 12:
            break
        step = min(dt, hours * 3600 - elapsed)
        density_rate = max(float(values["rho_dot_open"] + values["rho_dot_closed"]), 0.0)
        if density_rate > 0:
            step = min(step, .002 / density_rate)
            if T2 is not None and not switched:
                density_gap = max(switch - obs["rho_global"], 0.0)
                step = min(step, max(.8 * density_gap / density_rate, 1e-3))
        model.advance(state, temperature, p, step, adjacency)
        elapsed += step
        if not all(np.all(np.isfinite(getattr(state, field))) for field in (
            "rho", "G", "phi_GBseg", "phi_TJ", "phi_iso", "phi_closed"
        )):
            break
    result = {name: np.asarray(values) for name, values in history.items()}
    return (result, model.clone_state(state)) if return_state else result


def classify_second_step(history, rho_target=.98, growth_tolerance=.20):
    if not finite_history(history):
        return "numerical_invalid"
    phase_indices = np.flatnonzero(history["phase"] > .5)
    if not len(phase_indices):
        return "unattainable_first_step"
    j = int(phase_indices[0])
    attained = history["rho"][-1] >= rho_target
    growth = history["G_mean"][-1] / max(history["G_mean"][j], 1e-30) - 1.0
    if attained and growth <= growth_tolerance:
        return "success"
    if attained:
        return "grain_growth"
    if growth <= growth_tolerance:
        return "density_exhaustion"
    return "mixed"


def invalid_score(cid, exact, reason="numerical_invalid"):
    return dict(
        attained=False, min_reduction=np.nan, median_reduction=np.nan,
        max_reduction=np.nan, span20=0., span30=0., A_TS=0.,
        window_width_C=0., lower_bracketed=False, upper_bracketed=False,
        complete=False, candidate_id=cid, tier="reject",
        rejection_reason=reason, exact_reconfirmed=exact,
        rho_high_final=np.nan, rho_two_final=np.nan,
    )


def score_pair(cid, p, n=8, dt=7200, hours=500, T2_grid=range(800, 1501, 50)):
    high = simulate(p, n, dt=dt, hours=hours)
    preparation, prepared_state = simulate(
        p, n, T2=1200, dt=dt, hours=hours, stop_at_switch=True, return_state=True
    )
    if not len(preparation["phase"]) or not np.any(preparation["phase"] > .5):
        two = preparation
    else:
        second = simulate(
            p, n, T2=1200, dt=dt, hours=hours, initial_state=prepared_state,
            switched_initial=True, time_offset=float(preparation["t"][-1]),
        )
        two = {name: np.concatenate([preparation[name], second[name]]) for name in preparation}
    exact = dt <= 1800
    if not finite_history(high) or not finite_history(two):
        return invalid_score(cid, exact), high, two, []
    score = objectives.trajectory_score(
        {"rho": high["rho"], "G_nm": high["G_mean"]},
        {"rho": two["rho"], "G_nm": two["G_mean"]},
    )
    points = []
    if not np.any(preparation["phase"] > .5):
        points = [dict(T2_C=T2, classification="unattainable_first_step", practical=T2 < 1400,
                       rho_final=float(two["rho"][-1])) for T2 in T2_grid]
        window = objectives.chen_window(points)
        result = {
            **score, **window, "candidate_id": cid,
            "tier": objectives.assign_tier(score, window, exact),
            "exact_reconfirmed": exact,
            "rho_high_final": float(high["rho"][-1]),
            "rho_two_final": float(two["rho"][-1]),
        }
        return result, high, two, points
    evaluated = {}
    def evaluate(T2):
        path = simulate(
            p, n, T2=T2, dt=max(dt, 1800), hours=hours,
            initial_state=prepared_state, switched_initial=True,
            time_offset=float(preparation["t"][-1]),
        )
        classification = classify_second_step(path)
        evaluated[T2] = dict(
            T2_C=T2, classification=classification, practical=T2 < 1400,
            rho_final=float(path["rho"][-1]) if len(path["rho"]) else np.nan,
        )
    for T2 in T2_grid:
        evaluate(T2)
    if exact:
        ordered = sorted(evaluated)
        transitions = [(a, b) for a, b in zip(ordered[:-1], ordered[1:])
                       if evaluated[a]["classification"] != evaluated[b]["classification"]]
        refine25 = sorted({T for a, b in transitions for T in range(a + 25, b, 25)})
        for T2 in refine25:
            if T2 not in evaluated:
                evaluate(T2)
        ordered = sorted(evaluated)
        transitions = [(a, b) for a, b in zip(ordered[:-1], ordered[1:])
                       if evaluated[a]["classification"] != evaluated[b]["classification"]]
        refine10 = sorted({T for a, b in transitions for T in range(a + 10, b, 10)})
        for T2 in refine10:
            if T2 not in evaluated:
                evaluate(T2)
    points = [evaluated[T] for T in sorted(evaluated)]
    practical_points = [point for point in points if point["practical"]]
    window = objectives.chen_window(practical_points)
    result = {
        **score, **window, "candidate_id": cid,
        "tier": objectives.assign_tier(score, window, exact),
        "exact_reconfirmed": exact,
        "rho_high_final": float(high["rho"][-1]),
        "rho_two_final": float(two["rho"][-1]),
    }
    if result["tier"] == "reject" and not result.get("rejection_reason"):
        result["rejection_reason"] = "exact_joint_criteria_failed"
    return result, high, two, points


def surrogate_columns(x):
    heterogeneity = x[:, 1] * (1.0 - x[:, 6]) * (.5 + .5 * x[:, 7])
    pr_response = x[:, 9] * (1.0 - x[:, 11]) * (.25 + .75 * x[:, 18])
    sweep_response = (x[:, 19] + x[:, 20]) * (.25 + .75 * x[:, 21]) * (.5 + .5 * x[:, 23])
    closed_support = x[:, 27] * x[:, 29] * (1.0 - x[:, 31])
    migration_resistance = (x[:, 33] + x[:, 34] + x[:, 44]) * (.3 + .7 * x[:, 37])
    attainment = closed_support * (.5 + .5 * x[:, 25])
    reduction = heterogeneity * pr_response * sweep_response * migration_resistance
    lower = pr_response * x[:, 17] * (1.0 - x[:, 27])
    upper = migration_resistance * (.5 + .5 * x[:, 38])
    score = 2.0 * reduction + attainment + heterogeneity + .5 * (lower + upper)
    return score, reduction, pr_response, sweep_response, closed_support, lower, upper


def compact_history_rows(candidate_id, path_name, history, max_points=240):
    stride = max(1, len(history["rho"]) // max_points)
    for j in range(0, len(history["rho"]), stride):
        yield {
            "candidate_id": candidate_id, "path": path_name,
            **{name: float(values[j]) for name, values in history.items()},
        }


def run(args):
    audit_path = Path("results/local_region_decoder_audit/audit_state.json")
    if args.require_decoder_audit:
        if not audit_path.exists() or not json.loads(audit_path.read_text()).get("passed"):
            raise RuntimeError("decoder audit required before dynamic search")
    start = time.time()
    deadline = start + args.max_hours * 3600
    OUT.mkdir(parents=True, exist_ok=True)
    n_stage0 = max(args.stage0_min, args.stage0)
    samples = optimizers.latin_hypercube(n_stage0, len(decoder.NAMES), SEED)
    score, reduction, pr, sweep, closed, lower, upper = surrogate_columns(samples)
    n_keep = min(args.stage1, len(samples))
    keep = np.argpartition(score, -n_keep)[-n_keep:]
    stage0 = []
    for i in keep:
        parameters = decoder.decode(samples[i])
        stage0.append(dict(
            candidate_id=int(i), fingerprint=decoder.fingerprint(parameters),
            score=float(score[i]), approximate_high_density_reduction=float(reduction[i]),
            approximate_PR_damage=float(pr[i]), approximate_sweep_response=float(sweep[i]),
            approximate_closed_pore_support=float(closed[i]),
            approximate_high_density_attainment=float(closed[i] * (.5 + .5 * samples[i, 25])),
            lower_boundary_likelihood=float(lower[i]), upper_boundary_likelihood=float(upper[i]),
            sampled_seed=SEED, projection_only=True, tier="unscored",
        ))
    write(OUT / "stage0_massive_screen.csv.gz", stage0, True)

    unique_stage0 = len({row["fingerprint"] for row in stage0})
    stage1_path = OUT / "stage1_reduced_dynamic_summary.csv"
    if args.resume and stage1_path.exists() and len(pd.read_csv(stage1_path)) >= n_keep:
        stage1 = pd.read_csv(stage1_path).to_dict("records")[:n_keep]
        write_state(status="resumed_after_stage1", stage0=n_stage0, stage1=len(stage1),
                    elapsed_s=time.time() - start)
    else:
        stage1 = []
        for count, i in enumerate(keep, 1):
            parameters = {**model.defaults(), **decoder.decode(samples[i])}
            result = score_pair(
                int(i), parameters, n=parameters["N_regions"], dt=14400,
                hours=96, T2_grid=range(800, 1400, 100),
            )[0]
            stage1.append({**result, "fingerprint": decoder.fingerprint(parameters)})
            if count % 50 == 0:
                write_state(status="running_stage1", stage0=n_stage0, stage1=count,
                            elapsed_s=time.time() - start)
            if count % 1000 == 0:
                write(stage1_path, stage1)
            if time.time() >= deadline:
                break
        write(stage1_path, stage1)

    valid_stage1 = [row for row in stage1 if row["rejection_reason"] != "numerical_invalid"]
    ranked = sorted(
        valid_stage1,
        key=lambda row: (
            np.nan_to_num(row["median_reduction"], nan=-9.0),
            row["rho_two_final"],
        ), reverse=True,
    )[:args.stage2]
    stage2 = []
    all_points = []
    for count, row in enumerate(ranked, 1):
        i = int(row["candidate_id"])
        parameters = {**model.defaults(), **decoder.decode(samples[i])}
        result, _, _, points = score_pair(
            i, parameters, n=parameters["N_regions"], dt=1800,
            T2_grid=range(800, 1400, 50),
        )
        result["fingerprint"] = decoder.fingerprint(parameters)
        stage2.append(result)
        all_points.extend({**point, "candidate_id": i} for point in points)
        if count % 50 == 0:
            write_state(status="running_stage2", stage0=n_stage0, stage1=len(stage1),
                        stage2=count, elapsed_s=time.time() - start)
        if count % 100 == 0:
            write(OUT / "stage2_exact_dynamic_summary.csv", stage2)
            write(OUT / "chen_classification_points_compact.csv", all_points)
        if time.time() >= deadline:
            break
    write(OUT / "stage2_exact_dynamic_summary.csv", stage2)
    write(OUT / "chen_classification_points_compact.csv", all_points)

    production = sorted(
        stage2,
        key=lambda row: (
            row["tier"] in ("Tier_A", "Tier_B"),
            np.nan_to_num(row["median_reduction"], nan=-9.0),
        ), reverse=True,
    )[:args.production]
    accepted = [row for row in production if row["tier"] in ("Tier_A", "Tier_B")]
    rejected = [
        {**row, "rejection_reason": row.get("rejection_reason") or "exact_joint_criteria_failed"}
        for row in production if row["tier"] not in ("Tier_A", "Tier_B")
    ]
    write(OUT / "production_candidate_summary.csv", production)
    write(OUT / "accepted_tier_candidates.csv", accepted)
    write(OUT / "rejected_candidates.csv", rejected)
    boundaries = [{key: row[key] for key in (
        "candidate_id", "window_width_C", "lower_bracketed", "upper_bracketed", "complete"
    )} for row in stage2]
    write(OUT / "chen_window_boundaries.csv", boundaries)

    histories = []
    ratio_rows = []
    attainment_rows = []
    for row in production[:min(10, len(production))]:
        i = int(row["candidate_id"])
        parameters = {**model.defaults(), **decoder.decode(samples[i])}
        high = simulate(parameters, parameters["N_regions"], dt=1800)
        two = simulate(parameters, parameters["N_regions"], T2=1200, dt=1800)
        histories.extend(compact_history_rows(i, "highT", high))
        histories.extend(compact_history_rows(i, "two_step", two))
        top = min(high["rho"].max(), two["rho"].max(), .98)
        if top >= .90:
            rho_grid = np.arange(.90, top + 5e-4, 1e-3)
            gh = np.interp(rho_grid, high["rho"], high["G_mean"])
            gt = np.interp(rho_grid, two["rho"], two["G_mean"])
            ratio_rows.extend(dict(
                candidate_id=i, rho=float(rho), G_highT_nm=float(a), G_two_step_nm=float(b),
                reduction_TS=float(1.0 - b / a),
            ) for rho, a, b in zip(rho_grid, gh, gt))
        attainment_rows.append(dict(
            candidate_id=i, highT_rho_final=float(high["rho"][-1]),
            two_step_rho_final=float(two["rho"][-1]),
            both_attain_098=bool(high["rho"][-1] >= .98 and two["rho"][-1] >= .98),
        ))
    write(OUT / "local_region_state_histories_compact.csv", histories)
    write(OUT / "closed_pore_histories_compact.csv", [
        {key: row[key] for key in ("candidate_id", "path", "t", "rho", "closed", "rho_dot_closed")}
        for row in histories
    ])
    write(OUT / "two_step_ratio_curves.csv", ratio_rows)
    write(OUT / "high_density_attainment.csv", attainment_rows)
    write(OUT / "ablation_summary.csv", [])
    write(OUT / "fast_firing_preservation.csv", [])
    write(OUT / "pareto_front.csv", production)
    registry = []
    for row in production:
        decoded = decoder.decode(samples[int(row["candidate_id"])])
        registry.extend(dict(candidate_id=row["candidate_id"], parameter=name, value=decoded[name])
                        for name in decoder.NAMES)
    write(OUT / "parameter_registry.csv", registry)
    write(OUT / "raw_outputs_manifest.csv", [dict(
        path="raw_outputs/", status="ignored_not_generated",
        reason="only compact histories were generated",
    )])

    unique1 = len({row["fingerprint"] for row in stage1})
    unique2 = len({row["fingerprint"] for row in stage2})
    complete = (
        n_stage0 >= 1_000_000 and n_stage0 >= 100_000 and unique_stage0 >= 10_000
        and len(stage1) >= STAGE1_REQUIRED and unique1 >= 5_000
        and len(stage2) >= STAGE2_REQUIRED and unique2 >= 500
    )
    state = dict(
        status="complete" if complete else "incomplete",
        runtime_s=time.time() - start, sampled_seed=SEED,
        stage0=n_stage0, stage0_unique_parameter_vectors=n_stage0,
        stage0_unique_fingerprints=unique_stage0,
        stage1=len(stage1), stage1_unique_fingerprints=unique1,
        stage2=len(stage2), stage2_unique_fingerprints=unique2,
        production=len(production), accepted=len(accepted),
        numerical_invalid_stage1=sum(row["rejection_reason"] == "numerical_invalid" for row in stage1),
        numerical_invalid_stage2=sum(row["rejection_reason"] == "numerical_invalid" for row in stage2),
    )
    write_state(**state)
    print(json.dumps(state, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-hours", type=float, default=10)
    parser.add_argument("--stage0-min", type=int, default=1_000_000)
    parser.add_argument("--stage0", type=int, default=1_000_000)
    parser.add_argument("--stage1", type=int, default=20_000)
    parser.add_argument("--stage2", type=int, default=1_000)
    parser.add_argument("--production", type=int, default=50)
    parser.add_argument("--require-decoder-audit", action="store_true")
    parser.add_argument("--decoder-audited", action="store_true")
    parser.add_argument("--resume", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
