# Analytical Sintering Model

Reduced-order analytical sintering models and automated search scripts for finding a common formulation that can reproduce both:

1. fast-heating-rate sintering trajectories, where rapid heating gives more densification per grain-growth increment; and
2. two-step sintering trajectories, where a short high-temperature step followed by a lower-temperature hold gives better density-vs-grain-size efficiency than a long high-temperature path.

The current code is intentionally exploratory. It is a compact 0-D model intended for mechanism search and Codex-driven refactoring, not a final calibrated materials model.

## Repository contents

- `sinter_reference_model_v3_multibin.py` — self-contained multi-bin pore-population reference model.
- `sinter_reference_model_v6_grainstress_multibin.py` — wrapper adding grain-size capillary baseline stress and grain-size activity-window controls.
- `sweep_lambda_window_priority.py` — general Lambda-window priority search harness.
- `sweep_lambda_window_priority_v4.py` — corrected v6 launcher that passes grain-stress controls into the search.
- `debug_lambda_v6_grainstress.py` — deterministic diagnostic for second-step activity windows.
- `docs/CODEX_HANDOFF.md` — detailed Codex task instructions and success criteria.

## Quick start

```bash
python3 -m pip install -r requirements.txt
python3 debug_lambda_v6_grainstress.py
python3 sweep_lambda_window_priority_v4.py \
  --model sinter_reference_model_v6_grainstress_multibin \
  --n 500 \
  --rho-target 0.92 \
  --outdir sweep_lambda_window_priority_v6_test
```

## Topology-constrained prototype

```bash
python3 -m pytest -q
python3 run_topology_diagnostics.py
python3 search_topology_initial.py --n 12
python3 stress_test_topology_memory.py
python3 stress_test_pore_bin_memory.py
python3 density_window_processing_map.py
python3 initial_condition_factorial_map.py
python3 smoothing_gate_identifiability.py
python3 topology_gate_identifiability.py
python3 two_step_window_map.py
python3 mechanism_discrimination_study.py
python3 expanded_phase_space_exploration.py --workers 4
python3 expanded_upper_size_extension.py
python3 expanded_size_onset_refinement.py
python3 expanded_phase_space_analysis.py
python3 growth_mechanism_sensitivity.py --workers 4
python3 pore_junction_pinning_sensitivity.py --workers 4
python3 pore_location_topology_sensitivity.py
```

`topology_constrained_sintering.py` separates topology, stress, serial renewal
times, event yields, mechanism fluxes, and nonnegative dissipation weights. It
does not use a scalar total-efficiency multiplier. Its default memory mode is
the conservative, observable `pore_bin_redistribution`; the former empirical
topology-damage state remains available as an explicit ablation mode.
The bounded growth-mobility audit compares the unchanged baseline with
explicit junction-limited and threshold-mobility grain-migration closures;
see `docs/NANOSCALE_GROWTH_SUPPRESSION_MECHANISM.md` before interpreting the
prototype parameters.
The follow-up pore/junction-pinning audit resolves migration resistance using
the instantaneous pore-bin and connected-topology state. Its negative result
and limitations are documented in `docs/PORE_JUNCTION_PINNING_MECHANISM.md`.
The pore-placement ladder resolves GB-segment, triple-junction, and isolated
pore stores with conservative fluxes; see
`docs/PORE_LOCATION_TOPOLOGY_COUPLING.md` for its bounded negative result.

## Development status

The current target is not parameter fitting alone. The model needs a physically credible coupling among renewal-limited densification, grain growth, pore topology, pore-size-distribution evolution, stress generation, and competing dissipation. The existing implementation provides a starting point for automated searches and staged mechanism tests.
