#!/usr/bin/env python3
"""Audit the frozen material-property screen under the reframed Tier-B rules.

This is a recorded-data analysis, not a new parameter search.  ``--max-hours``
is an upper runtime guard retained for compatibility with longer future runs.
Topology parameters are neither loaded nor modified.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import time
from zipfile import ZipFile

import numpy as np
import pandas as pd


ARCHIVE = Path("results/1_Backup_of_prior_runs.zip")
OUT = Path("results/reframe_tierB_experimental_plausibility")
FAST_RATIO = 1.5
FAST_SPAN = 0.03
TWO_STEP_REDUCTION = 0.20
TWO_STEP_SPAN = 0.02


def archived_csv(z: ZipFile, name: str) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(z.read(name)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-hours", type=float, default=10.0)
    args = ap.parse_args()
    started = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    with ZipFile(ARCHIVE) as z:
        prefix = "separated_mechanism_production_search/"
        registry = archived_csv(z, prefix + "material_parameter_registry.csv")
        screen = archived_csv(z, prefix + "fast_firing_material_screen.csv")
        selected = archived_csv(z, "nucleation_limited_fast_firing_chen_production/selected_nucleation_material_sets.csv")
        timing = archived_csv(z, "nucleation_limited_fast_firing_chen_production/fast_firing_timing_audit.csv")
        strict = archived_csv(z, "strict_chen_window_production/strict_window_recheck.csv")

    pre = screen[screen.attained.fillna(False)].groupby("material_id").agg(
        fast_prepass_max_ratio=("max_ratio", "max"), fast_prepass_max_span=("span_ge_1p5", "max")
    ).reset_index()
    frame = registry.merge(pre, on="material_id", how="left")
    frame["fast_prepass"] = (frame.fast_prepass_max_ratio >= FAST_RATIO) & (frame.fast_prepass_max_span >= FAST_SPAN)
    frame["causal_fast_audit_available"] = frame.material_id.isin(selected.material_id)
    frame["fast_firing_pass"] = False
    frame["causal_fast_ratio"] = np.nan
    frame["causal_fast_median_ratio"] = np.nan
    frame["causal_fast_span_ge_1p5"] = np.nan
    frame["nucleation_facile_removes"] = pd.NA
    frame["PR_off_preserves_allowed"] = pd.NA
    frame["transport_or_exchange_control_rejected"] = pd.NA
    frame["two_step_evaluated"] = frame.material_id.isin(("E0021", "E0142"))
    frame["two_step_tier"] = "not evaluated"
    frame["two_step_pass"] = False
    frame["two_step_reduction_threshold"] = TWO_STEP_REDUCTION
    frame["two_step_span_threshold"] = TWO_STEP_SPAN
    for mid in ("E0021", "E0142"):
        audit = timing[timing.material_id == mid].set_index("ablation_mode")
        full = audit.loc["full_material_model"]
        no_nuc = audit.loc["no_nucleation_limitation"]
        no_pr = audit.loc["no_PR_redistribution"]
        transport = audit.loc["transport_only"]
        exchange = audit.loc["exchange_limited_variant"]
        idx = frame.material_id == mid
        frame.loc[idx,"fast_firing_pass"] = bool(full.meaningful and full.max_ratio >= FAST_RATIO and full.span_ge_1p5 >= FAST_SPAN and not bool(no_nuc.meaningful))
        frame.loc[idx,"causal_fast_ratio"] = float(full.max_ratio)
        frame.loc[idx,"causal_fast_median_ratio"] = float(full.median_ratio)
        frame.loc[idx,"causal_fast_span_ge_1p5"] = float(full.span_ge_1p5)
        frame.loc[idx,"nucleation_facile_removes"] = not bool(no_nuc.meaningful)
        frame.loc[idx,"PR_off_preserves_allowed"] = bool(no_pr.meaningful)
        frame.loc[idx,"transport_or_exchange_control_rejected"] = (not bool(exchange.meaningful)) and (not bool(transport.attained))
        tiers = set(strict[strict.material_id == mid].tier)
        best = "Tier_B" if "Tier_B" in tiers else ("Tier_C" if "Tier_C" in tiers else "reject")
        frame.loc[idx,"two_step_tier"] = best
        frame.loc[idx,"two_step_pass"] = best == "Tier_B"

    frame["joint_qualitative_behavior_pass"] = frame.fast_firing_pass & frame.two_step_pass
    frame["closed_accommodation_state_recorded"] = False
    frame["joint_full_evidence_pass"] = frame.joint_qualitative_behavior_pass & frame.closed_accommodation_state_recorded
    frame["topology_parameters_modified"] = False
    frame["Q_nuc_minus_Q_GB_kJ_mol"] = (frame.Q_disconnection_nucleation-frame.Q_GB_diffusion)/1000
    frame["Q_surface_minus_Q_GB_kJ_mol"] = (frame.Q_surface_diffusion-frame.Q_GB_diffusion)/1000
    frame["Q_transport_minus_Q_exchange_kJ_mol"] = (frame.Q_transport-frame.Q_exchange)/1000
    frame["Q_nuc_over_Q_GB"] = frame.Q_disconnection_nucleation/frame.Q_GB_diffusion
    frame["Q_surface_over_Q_GB"] = frame.Q_surface_diffusion/frame.Q_GB_diffusion
    frame["Q_transport_over_Q_exchange"] = frame.Q_transport/frame.Q_exchange
    frame["relative_ordering"] = np.where(
        (frame.Q_disconnection_nucleation>frame.Q_GB_diffusion)&(frame.Q_surface_diffusion>frame.Q_GB_diffusion),
        "Q_surface,Q_nuc > Q_GB", "other ordering")
    frame["interpretation"] = np.select(
        [frame.joint_qualitative_behavior_pass, frame.fast_firing_pass, frame.fast_prepass],
        ["joint qualitative evidence; closed-accommodation coupling not evaluated",
         "causal fast-firing evidence; no Tier-B two-step overlap",
         "fast-firing pre-screen only; causal ablations unavailable"],
        default="no causal fast-firing evidence")
    frame.to_csv(OUT/"relative_material_property_window_reframed.csv",index=False)
    sens = frame[frame.material_id.isin(("E0021","E0142"))][[
        "material_id","Q_GB_diffusion","Q_surface_diffusion","Q_disconnection_nucleation","Q_exchange","Q_transport",
        "Q_nuc_minus_Q_GB_kJ_mol","Q_surface_minus_Q_GB_kJ_mol","Q_transport_minus_Q_exchange_kJ_mol",
        "fast_firing_pass","causal_fast_ratio","causal_fast_median_ratio","causal_fast_span_ge_1p5",
        "nucleation_facile_removes","PR_off_preserves_allowed","two_step_tier","two_step_pass",
        "joint_qualitative_behavior_pass","closed_accommodation_state_recorded","joint_full_evidence_pass",
        "topology_parameters_modified","interpretation"]]
    sens.to_csv(OUT/"tierB_material_property_sensitivity.csv",index=False)
    state = dict(status="complete",max_hours=args.max_hours,runtime_s=time.monotonic()-started,
                 materials_in_frozen_registry=len(frame),causal_materials=len(sens),
                 joint_qualitative_materials=int(frame.joint_qualitative_behavior_pass.sum()),
                 joint_full_evidence_materials=int(frame.joint_full_evidence_pass.sum()),
                 topology_parameters_modified=False,new_parameter_search=False)
    (OUT/"material_property_window_run_state.json").write_text(json.dumps(state,indent=2)+"\n")
    print(json.dumps(state,indent=2))


if __name__ == "__main__":
    main()
