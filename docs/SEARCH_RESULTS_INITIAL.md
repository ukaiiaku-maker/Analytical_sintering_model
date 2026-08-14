# Initial deterministic results

## Tests and protocols

The invariant tests cover nonnegative/conservative partition weights, exact
pore-volume closure, nonnegative pore populations, required diagnostics, and
the expected decrease of boundary coverage when pore size increases at fixed
pore volume.

Four baseline paths were compared at the unchanged target density `rho=0.90`:
0.2 and 20 C/min ramps to 1450 C, a 1350 C isothermal path, and a 1350 to
1250 C density-triggered two-step path. All four reach the target. Exact values
are recorded in `results/initial/protocol_summary.csv` and
`percent_metrics.csv`.

Baseline values are `HR_pct = -6.4503%` and `TS_pct = 11.5345%`. Grain sizes at
the target are 424.31 nm (slow), 451.68 nm (fast), 394.37 nm (high-T), and
348.88 nm (two-step).

## What works

- The pore, topology, stress, renewal-time, flux, efficiency, and power
  diagnostics are emitted on every step.
- Pore conservation is enforced without clipping density independently.
- The baseline produces a positive two-step percentage improvement.
- The second-step `rho0=0.83` grid exposes both activity and `E_G`, preventing
  coarse-grain high-Lambda states from being labeled efficient by activity
  alone.

## What fails

The baseline heating-rate percentage is negative: the rapid path reaches the
target with more grain growth. This means the current serial event-time and
topology memory are not yet sufficient to reproduce fast-heating benefit. The
small seeded search is intentionally bounded and uses `min(HR_pct, TS_pct)` plus
`E_G`; it does not optimize raw nanometers and rejects any sample where a path
misses the target.

The best target-reaching row in the 12-case seeded search has
`HR_pct = -7.3453%`, `TS_pct = 14.9760%`, and median `E_G = 7.2691`. Thus the
search strengthened only the already-correct two-step response and provides no
evidence that parameter variation inside these narrow bounds repairs the
missing slow-ramp topology-memory mechanism.

## Next experiments

1. Add an explicit accumulated intermediate-temperature surface-coarsening
   state so slow ramps degrade removable-pore topology before activation.
2. Make event growth depend on which completion step occurred, rather than a
   single per-event increment.
3. Replace normalized propensities with a convex Onsager partition constrained
   by interfacial free-energy decrease.
4. Run ablations for P-R, triple-line drag, and pore coarsening before expanding
   the search bounds.
