#!/usr/bin/env python3
"""Apply production gates after exact local-region confirmation."""
from pathlib import Path
import argparse
import csv
import io
import zipfile
import json
import numpy as np
import pandas as pd

import interacting_local_region_decoder as decoder
import interacting_local_region_model as model
import interacting_local_region_search as search
import massive_latent_topology_optimizers as optimizers

OUT = Path("results/local_region_decoder_corrected_dynamic_search")
ARCHIVE = Path("results/1_Backup_of_prior_runs.zip")


def write(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(key for row in rows for key in row)) or ["status"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def preparation_audit(samples, candidates):
    rows = []
    for candidate in candidates.itertuples(index=False):
        cid = int(candidate.candidate_id)
        p = {**model.defaults(), **decoder.decode(samples[cid])}
        history = search.simulate(p, p["N_regions"], T2=1200, dt=1800)
        switch = np.flatnonzero(history["phase"] > .5)
        if len(switch):
            j = int(switch[0])
            growth = history["G_mean"][j] / history["G_mean"][0] - 1.0
            rho_switch = history["rho"][j]
            closed = history["closed"][j]
        else:
            growth, rho_switch, closed = np.nan, np.nan, np.nan
        support = bool(np.max(history["rho_dot_closed"]) > 1e-12)
        rows.append(dict(
            candidate_id=cid, rho_switch_exact=rho_switch,
            switch_error=abs(rho_switch - .88) if np.isfinite(rho_switch) else np.nan,
            G0_nm=history["G_mean"][0], G1_nm=history["G_mean"][j] if len(switch) else np.nan,
            first_step_growth_fraction=growth, first_step_growth_le_20pct=bool(growth <= .20),
            closed_fraction_at_switch=closed, high_density_support_active=support,
            production_gate_passed=bool(growth <= .20 and support and abs(rho_switch - .88) <= 2e-4),
        ))
    return rows


def archived_fast_firing():
    member = "nucleation_limited_fast_firing_chen_production/fast_firing_timing_audit.csv"
    with zipfile.ZipFile(ARCHIVE) as archive:
        frame = pd.read_csv(io.BytesIO(archive.read(member)))
    rows = []
    for material_id in ("E0021", "E0142"):
        group = frame[frame.material_id == material_id]
        full = group[group.ablation_mode == "full_material_model"].iloc[0]
        facile = group[group.ablation_mode == "no_nucleation_limitation"].iloc[0]
        pr_off = group[group.ablation_mode == "no_PR_redistribution"].iloc[0]
        rows.append(dict(
            material_id=material_id, source="results/1_Backup_of_prior_runs.zip:" + member,
            full_meaningful=bool(full.meaningful), full_max_ratio=full.max_ratio,
            full_span_ge_1p5=full.span_ge_1p5,
            nucleation_facile_meaningful=bool(facile.meaningful),
            PR_off_meaningful=bool(pr_off.meaningful),
            preserved=bool(full.meaningful and not facile.meaningful),
            local_region_topology_changes_material_rho_dot=False,
            note="Frozen material baseline; local-region search does not modify MaterialKinetics",
        ))
    return rows


def ablation_parameters(base, name):
    p = dict(base)
    if name == "no_PR_damage":
        p["k_PR"] = 0.0
    elif name == "no_sweep_coalescence":
        p["k_sweep_damaged"] = p["k_sweep_connected"] = 0.0
    elif name == "no_local_region_heterogeneity":
        p.update(weight_sigma=0.0, rho_sigma=0.0, G_sigma=0.0, cluster=0.0)
    elif name == "homogenized_network":
        p.update(degree=max(int(p["N_regions"]) - 1, 1), cluster=1.0)
    elif name == "no_interregion_exchange":
        p["exchange_rate"] = 0.0
    elif name == "no_pore_detachment":
        p["detachment"] = 0.0
    elif name == "no_recapture":
        p["recapture"] = 0.0
    elif name == "no_closed_transition":
        p["closed_transition"] = 0.0
    elif name == "no_closed_shrinkage":
        p["k_closed"] = 0.0
    elif name == "infinite_closed_accommodation":
        p.update(closed_capacity=1e6, capacity_tau=1.0)
    elif name == "no_attached_pore_drag":
        p["attached_drag"] = p["pore_drag_fraction"] = 0.0
    elif name == "no_persistent_junction":
        p["junction_drag"] = p["XJ_prod"] = 0.0
    elif name == "no_TJ_multihit":
        p.update(lambda_TJ=1e6, K_TJ=1.0, q_TJ=0)
    elif name.startswith("q_TJ_"):
        p["q_TJ"] = int(name[-1])
    elif name == "no_residual_stress":
        p.update(stress_PR=0.0, stress_shear=0.0, stress_migration=0.0, stress_nucleation=0.0)
    elif name == "topology_disabled":
        p.update(attached_drag=0.0, junction_drag=0.0, pore_drag_fraction=0.0,
                 stress_migration=0.0, lambda_TJ=1e6, K_TJ=1.0, q_TJ=0)
    return p


def causal_ablations(samples, accepted):
    names = (
        "full_model", "no_PR_damage", "no_sweep_coalescence",
        "no_local_region_heterogeneity", "homogenized_network",
        "no_interregion_exchange", "no_pore_detachment", "no_recapture",
        "no_closed_transition", "no_closed_shrinkage",
        "infinite_closed_accommodation", "no_attached_pore_drag",
        "no_persistent_junction", "no_TJ_multihit", "q_TJ_0", "q_TJ_1",
        "q_TJ_2", "no_residual_stress", "topology_disabled",
    )
    rows = []
    for candidate in accepted:
        cid = int(candidate["candidate_id"])
        base = {**model.defaults(), **decoder.decode(samples[cid])}
        for name in names:
            p = ablation_parameters(base, name)
            score = search.score_pair(
                cid, p, n=p["N_regions"], dt=1800, hours=500,
                T2_grid=range(800, 1400, 50),
            )[0]
            rows.append(dict(
                candidate_id=cid, ablation=name, tier=score["tier"],
                median_reduction=score["median_reduction"], span20=score["span20"],
                window_width_C=score["window_width_C"], lower_bracketed=score["lower_bracketed"],
                upper_bracketed=score["upper_bracketed"], complete=score["complete"],
                rho_high_final=score["rho_high_final"], rho_two_final=score["rho_two_final"],
                loss_from_full=(name != "full_model" and score["tier"] not in ("Tier_A", "Tier_B")),
            ))
    return rows


def robustness(samples, accepted):
    rows = []
    for candidate in accepted[:3]:
        cid = int(candidate["candidate_id"])
        base = {**model.defaults(), **decoder.decode(samples[cid])}
        for rho0 in (.65, .70, .75):
            for G0 in (75., 100., 150.):
                p = {**base, "rho0": rho0, "G0_nm": G0}
                score = search.score_pair(cid, p, n=p["N_regions"], dt=1800,
                                          T2_grid=range(800, 1400, 50))[0]
                rows.append(dict(candidate_id=cid, rho0=rho0, G0_nm=G0,
                                 tier=score["tier"], median_reduction=score["median_reduction"],
                                 window_width_C=score["window_width_C"], complete=score["complete"],
                                 attained=score["attained"]))
    return rows


def production_histories(samples, accepted):
    histories, ratios, attainment = [], [], []
    for candidate in accepted:
        cid = int(candidate["candidate_id"])
        p = {**model.defaults(), **decoder.decode(samples[cid])}
        _, high, two, _ = search.score_pair(
            cid, p, n=p["N_regions"], dt=1800, T2_grid=range(800, 1400, 50)
        )
        for path, history in (("highT", high), ("two_step", two)):
            histories.extend(search.compact_history_rows(cid, path, history, max_points=300))
        top = min(high["rho"].max(), two["rho"].max(), .98)
        if top >= .90:
            rho = np.arange(.90, top + 5e-4, 1e-3)
            gh = np.interp(rho, high["rho"], high["G_mean"])
            gt = np.interp(rho, two["rho"], two["G_mean"])
            ratios.extend(dict(candidate_id=cid, rho=float(x), G_highT_nm=float(a),
                               G_two_step_nm=float(b), reduction_TS=float(1 - b / a))
                          for x, a, b in zip(rho, gh, gt))
        attainment.append(dict(candidate_id=cid, highT_rho_final=float(high["rho"][-1]),
                               two_step_rho_final=float(two["rho"][-1]),
                               both_attain_098=bool(high["rho"][-1] >= .98 and two["rho"][-1] >= .98)))
    write(OUT / "local_region_state_histories_compact.csv", histories)
    write(OUT / "closed_pore_histories_compact.csv", [
        {key: row[key] for key in ("candidate_id", "path", "t", "rho", "closed",
                                    "closed_accommodation", "rho_dot_closed")}
        for row in histories
    ])
    write(OUT / "two_step_ratio_curves.csv", ratios)
    write(OUT / "high_density_attainment.csv", attainment)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--histories-only", action="store_true")
    args = parser.parse_args()
    exact = pd.read_csv(OUT / "stage2_exact_dynamic_summary.csv")
    provisional = exact[exact.tier.isin(["Tier_A", "Tier_B"])].copy()
    samples = optimizers.latin_hypercube(1_000_000, len(decoder.NAMES), search.SEED)
    if args.histories_only:
        current = pd.read_csv(OUT / "accepted_tier_candidates.csv").to_dict("records")
        production_histories(samples, current)
        print(f"regenerated_histories={len(current)}")
        return
    audits = preparation_audit(samples, provisional)
    write(OUT / "production_gate_audit.csv", audits)
    audit_by_id = {row["candidate_id"]: row for row in audits}
    accepted, rejected = [], []
    for row in provisional.to_dict("records"):
        audit = audit_by_id[int(row["candidate_id"])]
        combined = {**row, **audit}
        if audit["production_gate_passed"]:
            accepted.append(combined)
        else:
            combined["rejection_reason"] = (
                "first_step_growth_above_20pct" if not audit["first_step_growth_le_20pct"]
                else "inactive_high_density_support_or_switch_mismatch"
            )
            rejected.append(combined)
    write(OUT / "accepted_tier_candidates.csv", accepted)
    other = exact[~exact.tier.isin(["Tier_A", "Tier_B"])].sort_values(
        "median_reduction", ascending=False, na_position="last"
    ).head(max(0, 50 - len(rejected))).to_dict("records")
    for row in other:
        row["rejection_reason"] = row.get("rejection_reason") or "trajectory_or_window_below_threshold"
    write(OUT / "rejected_candidates.csv", rejected + other)
    write(OUT / "fast_firing_preservation.csv", archived_fast_firing())
    if accepted:
        ablations = causal_ablations(samples, accepted)
        write(OUT / "ablation_summary.csv", ablations)
        write(OUT / "robustness_summary.csv", robustness(samples, accepted))
        production_histories(samples, accepted)
    else:
        write(OUT / "ablation_summary.csv", [dict(
            candidate_id="", ablation="not_run", status="no_candidate_passed_preparation_gate",
            reason="causal production ablations cannot promote a preparation-ineligible candidate",
        )])
    state_path = OUT / "run_state.json"
    state = json.loads(state_path.read_text())
    state.update(provisional_tier_candidates=len(provisional), production_accepted=len(accepted),
                 accepted=len(accepted), production_gate_applied=True)
    state_path.write_text(json.dumps(state, indent=2) + "\n")
    print(f"provisional={len(provisional)} production_accepted={len(accepted)}")


if __name__ == "__main__":
    main()
