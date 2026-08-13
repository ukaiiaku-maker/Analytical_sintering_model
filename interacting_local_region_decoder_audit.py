#!/usr/bin/env python3
"""Decoder diversity and dynamical-influence preflight."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np

import interacting_local_region_decoder as decoder
import interacting_local_region_model as model
import massive_latent_topology_optimizers as optimizers

OUT = Path("results/local_region_decoder_audit")


def write(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(key for row in rows for key in row)) or ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def dynamic_signature(unit_vector):
    p = {**model.defaults(), **decoder.decode(unit_vector)}
    state = model.initial(p["N_regions"], p=p, seed=13)
    adjacency = model.network_adjacency(p["N_regions"], p, seed=13)
    signature = [
        p["N_regions"], np.var(state.weights), np.var(state.rho), np.var(state.G),
        adjacency.sum(), np.var(adjacency.sum(axis=1)),
        np.mean(state.connected_removable_fraction), np.mean(state.closed_fraction),
    ]
    for temperature in (1000, 1250, 1450):
        flux = model.local_fluxes(state, temperature, p)
        signature.extend(float(np.mean(flux[name])) for name in (
            "rho_dot_open", "rho_dot_closed", "PR_damage", "sweep", "detachment",
            "recapture", "closed_transition", "G_dot", "migration_factor",
            "activity", "Lambda_TJ", "K_TJ", "P_comp_TJ",
        ))
        model.advance(state, temperature, p, 600.0, adjacency)
        signature.extend(float(np.mean(getattr(state, name))) for name in (
            "phi_GBseg", "phi_TJ", "phi_iso", "phi_closed", "N_GBseg", "N_TJ",
            "N_iso", "N_closed", "PR_damage_memory", "sweep_memory", "X_J",
            "residual_stress", "closed_accommodation", "migration_factor",
        ))
    return np.nan_to_num(np.asarray(signature), nan=0., posinf=1e300, neginf=-1e300)


def influence_rows():
    baseline = np.full(len(decoder.NAMES), .37)
    base_signature = dynamic_signature(baseline)
    rows = []
    for index, name in enumerate(decoder.NAMES):
        probe = baseline.copy()
        probe[index] = .83
        delta = dynamic_signature(probe) - base_signature
        scale = np.maximum(np.abs(base_signature), 1e-15)
        effect = float(np.max(np.abs(delta) / scale))
        rows.append(dict(sampled_column=name, diagnostic_only=False,
                         dynamic_influence=effect, unused=effect <= 1e-12))
    return rows


def run(n):
    samples = optimizers.latin_hypercube(n, len(decoder.NAMES), 20260817)
    decoded = [decoder.decode(row) for row in samples]
    fingerprints = [decoder.fingerprint(p) for p in decoded]
    parameter_rows = []
    for index, name in enumerate(decoder.NAMES):
        values = np.array([p[name] for p in decoded])
        parameter_rows.append(dict(
            parameter=name, minimum=values.min(), maximum=values.max(),
            variance=values.var(), source_column=index,
            correlation=float(np.corrcoef(samples[:, index], values)[0, 1]),
        ))
    unique_vectors = len({tuple(round(float(p[name]), 12) for name in decoder.NAMES) for p in decoded})
    influence = influence_rows()
    unused = [row for row in influence if row["unused"]]
    summary = dict(
        sampled_rows=n, unique_parameter_vectors=unique_vectors,
        unique_dynamic_fingerprints=len(set(fingerprints)),
        unique_N_regions=len({p["N_regions"] for p in decoded}),
        passed=(unique_vectors >= .95 * n and len(set(fingerprints)) >= 1000 and not unused),
        unused_columns=len(unused), dynamically_influential_columns=len(influence) - len(unused),
    )
    write(OUT / "decoder_parameter_map.csv", parameter_rows)
    write(OUT / "decoder_fingerprint_summary.csv", [summary])
    write(OUT / "unused_sampled_columns.csv", influence)
    main_out = Path("results/local_region_decoder_corrected_dynamic_search")
    main_out.mkdir(parents=True, exist_ok=True)
    write(main_out / "decoder_parameter_map.csv", parameter_rows)
    write(main_out / "decoder_fingerprint_summary.csv", [summary])
    write(main_out / "unused_sampled_columns.csv", influence)
    (OUT / "audit_state.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=10_000)
    run(parser.parse_args().samples)


if __name__ == "__main__":
    main()
