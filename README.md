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
python3 pore_location_agentic_sensitivity.py
python3 agentic_mechanism_search.py --workers 1
python3 adaptive_T2_boundary_search.py
python3 preparation_window_search.py
python3 production_mechanism_assessment.py
python3 production_mechanism_postprocess.py
python3 joint_pr_desintering_search.py --workers 8
python3 joint_pr_desintering_search.py --full --selected-variant PR_attrition_moderate --workers 8 --resume-chen
python3 pr_desintering_postprocess.py
python3 production_pr_desintering_assessment.py --workers 8
python3 tj_constraint_ablation.py
python3 production_pr_postprocess.py
python3 generate_paper_figures.py
python3 generate_supplement_figures.py
python3 observable_trajectory_effect_audit.py
python3 joint_heterogeneity_residual_stress_search.py --workers 4
python3 rejected_case_failure_decomposition.py
python3 persistent_defect_memory_screen.py --workers 4
python3 local_connected_sink_defect_search.py --workers 4
python3 late_stage_closed_pore_search.py --workers 4
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
The follow-up action layer allocates local GB-segment, triple-junction,
capture, isolation, and migration actions through nonnegative competing
propensities. Its selected Chen-map audit and ablations are documented in
`docs/AGENTIC_PORE_LOCATION_TOPOLOGY.md`.
The source-grounded discovery layer compares persistent junction populations
and Class-B TJ multihit completion against that negative control; see
`docs/SOURCE_MECHANISM_PRIORS.md` and
`docs/AGENTIC_MECHANISM_SEARCH_REPORT.md`.
The frozen-mechanism production campaign and its negative joint-response result
are documented in `docs/PRODUCTION_MECHANISM_ASSESSMENT.md`. Raw point tables
are intentionally ignored; compact review tables and seven figures are kept in
`results/production_mechanism_assessment/`.
The follow-up local PR/de-sintering competition test is documented in
`docs/PR_DESINTERING_FAST_FIRING_MEMORY.md`. It preserves the frozen negative
control as `early_memory_mode="disabled"` and adds only conservative,
non-densifying pore redistribution competing with renewal removal.
The production confirmation, local PR robustness audit, and pore-occupied
versus structurally constrained TJ ablation are reported in
`docs/PRODUCTION_PR_DESINTERING_ASSESSMENT.md`.
The dedicated interpretation and focused q0/q1 evidence for pore-filled versus
structurally constrained TJs are in `docs/TJ_PORE_CONSTRAINT_INTERPRETATION.md`.
The manuscript-ready figure package is reproducibly generated into
`results/paper_figures/`; its ordering, captions, and visual conventions are
documented in `docs/FIGURE_MANIFEST.md`, `docs/FIGURE_CAPTION_DRAFTS.md`, and
`docs/PLOT_STYLE_GUIDE.md`.
The subsequent observable trajectory-effect audit replaces the weak internal
`HR_pct > 1%` criterion with a matched-density ratio and finite-span test; see
`docs/OBSERVABLE_TRAJECTORY_EFFECT_AUDIT.md` before making fast-firing claims.
The bounded weighted-cohort and residual-stress follow-up is documented in
`docs/HETEROGENEITY_RESIDUAL_STRESS_TRAJECTORY_SEARCH.md`; it retains the
observable finite-density-span criterion and reports a negative result.

## Development status

The current target is not parameter fitting alone. The model needs a physically credible coupling among renewal-limited densification, grain growth, pore topology, pore-size-distribution evolution, stress generation, and competing dissipation. The existing implementation provides a starting point for automated searches and staged mechanism tests.
