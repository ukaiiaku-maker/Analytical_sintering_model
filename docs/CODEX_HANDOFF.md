# Codex handoff: automated search for a unified analytical sintering model

## Objective

Build, audit, and extend a reduced-order analytical sintering model that can reproduce, in one formulation, both:

1. **Fast-heating-rate sintering:** rapid heating gives a better density-vs-grain-size trajectory than slow heating, i.e. smaller grain size at the same target density.
2. **Two-step sintering:** a short high-temperature step followed by a lower-temperature hold gives a better density-vs-grain-size trajectory than a long high-temperature hold.

The task is not to tune a single plot. The goal is to find a physically credible coupling among renewal-limited densification, grain growth, pore topology, pore-size-distribution evolution, stress generation, and non-densifying dissipation.

## Current code

The repository seed includes:

- `sinter_reference_model_v3_multibin.py`: self-contained multi-bin pore distribution reference model.
- `sinter_reference_model_v6_grainstress_multibin.py`: wrapper adding explicit grain-size capillary baseline stress and grain-size activity-window control.
- `sweep_lambda_window_priority.py`: automated parameter search harness.
- `sweep_lambda_window_priority_v4.py`: compatibility launcher.
- `debug_lambda_v6_grainstress.py`: deterministic diagnostic script for second-step activity windows.

Run the diagnostic first:

```bash
python3 debug_lambda_v6_grainstress.py
```

Then run a small search:

```bash
python3 sweep_lambda_window_priority_v4.py \
  --model sinter_reference_model_v6_grainstress_multibin \
  --n 500 \
  --rho-target 0.92 \
  --outdir sweep_lambda_window_priority_v6_test
```

Do not trust any large sweep until the diagnostic behavior is physically sensible.

## Scientific framing

Interface-normal strain during sintering is modeled as a renewal process. Each event requires source/sink activation and point-defect transport/exchange. The active sink population is an evolving state variable, not a fixed continuum boundary condition.

The reduced model should distinguish:

- geometrically eligible boundary/pore topology;
- kinetic activity, represented by `Lambda = r_nuc * tau_sink` and `activity = Lambda/(1 + Lambda)`;
- densifying strain rate at active eligible boundaries;
- non-densifying interfacial evolution, including surface/pore smoothing, grain growth, pore drag, clean-boundary migration, and de-sintering/P-R-like pore redistribution;
- pore-size distribution evolution, because small pores close faster while large pores persist, reduce boundary coverage, and may dissipate more drag energy.

## Main research question

Find a common analytical formulation that predicts both:

```text
G_20Cmin(rho*) < G_0.2Cmin(rho*)
```

and

```text
G_two_step(rho*) < G_highT(rho*)
```

with meaningful percentage improvement at a fixed density, while still producing substantial densification on realistic time/temperature paths.

## Known problems from prior attempts

1. **Lambda alone is not enough.** High `Lambda` can occur in coarse-grain states if `tau_sink` becomes large. A useful activity window must be paired with densification efficiency, such as `d rho / d ln G`.
2. **Raw grain-size difference is a poor score.** Use percent improvements and impose maximum grain-size sanity bounds.
3. **Stress decomposition needs careful auditing.** The current v6 wrapper separates `sigma_base_grain ~ 2 gamma/R_G` from a concentration/geometric stress term. Check whether the base model already includes part of this stress to avoid double counting.
4. **Pore topology is simplified.** The multi-bin distribution helps, but topology is still described by reduced moments and projected area coverage, not a generated effective network.
5. **Power channels are partly phenomenological.** If a channel is claimed to control the trajectory, it needs diagnostic output and an ablation switch.

## Proposed Codex work plan

### Phase 0 — Repository health

- Run `python3 -m py_compile *.py`.
- Add a `tests/` directory with deterministic smoke tests for diagnostics and a small sweep.
- Ensure output directories and generated CSV/PNG files are ignored by git.

### Phase 1 — Audit current v6 behavior

- Run `debug_lambda_v6_grainstress.py`.
- Confirm whether the activity-window mechanism behaves correctly as initial `G0` is varied.
- Add output columns for `sigma_base`, `sigma_concentration`, `Lambda`, `activity`, `rho_dot`, `dGdt/G`, and `d rho / d ln G`.

### Phase 2 — Improve objective functions

Add a score based on:

```text
HR_pct = 100*(G_slow - G_fast)/G_slow
TS_pct = 100*(G_highT - G_twoStep)/G_highT
window_efficiency = median(activity * rho_dot/(dGdt/G + eps))
```

Successful cases should pass both percent-improvement tests and have positive activity/trajectory-efficiency diagnostics.

### Phase 3 — Mechanism modularization

Refactor into modules:

- renewal kinetics;
- stress generation and partitioning;
- pore topology/multi-bin distribution;
- grain growth / Zener / drag;
- objective functions and search.

All mechanisms should be switchable for ablation.

### Phase 4 — Test candidate coupling strategies

Test at least three coupling architectures:

1. **Lagrange-multiplier / constrained dissipation partition.** Allocate interfacial power among densification, clean GB growth, pore drag, PR/de-sintering, and pore redistribution subject to a thermodynamic constraint.
2. **Effective microstructure topology loop.** Construct an effective microstructure from instantaneous descriptors, evaluate topology, evolve kinetics for one adaptive step, update topology, and repeat.
3. **Brute-force / Bayesian search.** Use Latin-hypercube or Optuna-style global exploration, followed by local refinement near feasible regions.

### Phase 5 — Add triple-line drag option

Implement an optional nucleation-rate-limited triple-line drag mechanism:

```text
P_TL_drag = coefficient * triple_line_density * v_GB * drag_force
activity_TL = Lambda_TL/(1 + Lambda_TL)
```

Make it switchable and include diagnostics.

### Phase 6 — Automated report

For each candidate model, generate a Markdown report with:

- parameter table;
- fast-heating trajectory comparison;
- two-step trajectory comparison;
- activity and efficiency windows;
- pore-distribution evolution;
- ablation summary showing which mechanisms are necessary.

## Acceptance criteria

A candidate formulation is acceptable only if it satisfies all of the following:

- reaches the target density without runaway coarsening;
- fast heating improves the density-vs-grain-size trajectory;
- two-step sintering improves the density-vs-grain-size trajectory;
- success persists in held-out protocols not used in the search;
- mechanism diagnostics are internally consistent.

## Avoid

- Do not optimize only raw grain-size differences in nm.
- Do not accept high-`Lambda` solutions if `d rho/d ln G` is poor.
- Do not hide physical changes in unnamed tuning constants.
- Do not replace the renewal model with a purely empirical densification law unless the renewal diagnostics are preserved for comparison.
