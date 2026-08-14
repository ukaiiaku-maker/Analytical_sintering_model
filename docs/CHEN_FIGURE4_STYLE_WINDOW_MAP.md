# Chen/Wang Figure-4-style process-window audit

## Status and scope

This is a fixed-parameter process-window prediction and model-discrimination exercise, not a fit to Chen/Wang Y2O3 data and not validation. PR #2 remains draft. The bounded subset requested for an initial calculation was executed: G0 = 50, 75, 100, 150, 225, 300, 450 nm; T1 = 1250, 1300, 1350 C; rho switch = 0.75, 0.80, 0.85, 0.88; and T2 = 1000--1300 C in 50 C increments. Density and connectivity smoothing gates use identical material parameters and 96 h first/second-step budgets.

Every second step begins from an actually simulated first-step state. No synthetic state or extended hold is used. The campaign contains 168 first-step conditions, 1,176 second-step trajectories, and 10,584 classifications (three density targets times three growth tolerances).

## Classification result

All 168 requested first-step states reached their switch density. Thus this bounded subset contains no `UNATTAINABLE_FIRST_STEP` region; the absence is reported rather than manufactured.

At rho target 0.90 and the primary 5% growth tolerance, neither gate produces a `SUCCESS` point:

| gate | densification exhaustion | grain-growth failure | mixed failure | success |
|---|---:|---:|---:|---:|
| density | 270 | 183 | 135 | 0 |
| connectivity | 278 | 182 | 128 | 0 |

At a 10% growth tolerance a very narrow window appears only for G0 = 450 nm and T2 = 1200 C: six density-gate and twelve connectivity-gate first-step states succeed. No rho target 0.95 or 0.98 case succeeds at any tolerance. These high-density failures are direct evidence of missing late-stage physics in the current model.

The map has physically distinct lower- and upper-temperature failures. T2 <= 1100 C is predominantly densification exhaustion; T2 >= 1250 C is predominantly grain-growth failure; 1150--1200 C contains mixed behavior. However, the required overlap of adequate densification and <=5% growth is absent.

## Particle-size result

The calculation does **not** support the desired nanoscale two-step domain. At 5% tolerance every size lacks a window. Small sizes fail primarily through relative grain growth/mixed behavior. As size increases, relative grain growth becomes less severe and the dominant failure shifts toward densification exhaustion; only the 450 nm class succeeds under the relaxed 10% criterion.

Because no nanoscale success-to-failure transition exists, no adaptive size point was added. Refining an interval would imply a boundary not present in the coarse result. The outcome is classified `INCONSISTENT_WITH_SWEEP` for nanoscale robustness and size-boundary behavior.

The connectivity gate modestly improves the relaxed window count (12 versus 6) but does not repair the qualitative topology: it cannot create a 5% window, nanoscale eligibility, or high-density attainment. It therefore improves some kinetics without supplying the missing scale separation.

## Boundary interpretation

For each actual first-step state, `window_boundaries.csv` reports:

- `T_lower_density_C`: lowest tested T2 reaching the density target;
- `T_upper_no_growth_C`: highest tested T2 below the growth tolerance;
- successful T2 minimum/maximum and width when the two sets overlap.

At 5% tolerance the density and no-growth intervals do not overlap. This is the central negative result: the model activates enough densification only after grain growth has already exceeded the allowed increment.

## Mechanism classification

| proposed element | assessment | reason |
|---|---|---|
| lower-T densification exhaustion | SUPPORTED | reproduced broadly at 1000--1100 C |
| upper-T grain-growth failure | SUPPORTED | reproduced broadly at 1250--1300 C |
| finite 5% two-step window | INCONSISTENT_WITH_SWEEP | no successful point |
| finite 10% two-step window | SUPPORTED, NARROW | only G0=450 nm, T2=1200 C |
| nanoscale two-step robustness | INCONSISTENT_WITH_SWEEP | relative growth is worst at small G0 |
| connectivity-controlled shift | SUPPORTED WEAKLY | changes counts but not qualitative boundary |
| rho >= 0.95 map | REQUIRED MISSING PHYSICS | no target attainment |

The simplest diagnosis is that the current grain-growth law and densification kinetics do not develop the Chen/Wang-style timescale separation at nanoscale G1. Explicit late-stage/closed-pore physics remains necessary for full-density maps, but was intentionally not added in this task. A future mechanism test must also address why fractional grain growth is excessive at small G, rather than tuning success thresholds.

## Outputs

- `two_step_window_points.csv`: all state diagnostics and classifications.
- `window_boundaries.csv`: lower/upper bounds per first-step state.
- `failure_modes.csv`: every non-successful classification.
- `size_window_summary.csv`: size-level mechanism classification.
- requested G1/rho1 classification, boundary, growth, and density maps.

Run with `python3 two_step_window_map.py`.
