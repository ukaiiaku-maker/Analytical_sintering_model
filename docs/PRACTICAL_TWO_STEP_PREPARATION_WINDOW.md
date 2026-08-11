# Practical two-step preparation-window search

## Scope

This audit uses the adaptive, censor-aware Chen-window solver without changing
densification kinetics, pore redistribution, the 0.90 target, mechanism
parameters, or the uniform 96 h budgets.  It evaluates `mech_009`, `mech_019`,
and explicit q=0 variants over all 4,320 requested combinations of T1, heating
rate, switch density, and G0.  Six switch-density states are extracted from
each shared first-step trajectory, avoiding redundant integration without
merging route-specific practical classifications.

Of 4,320 preparation states, 1,040 do not attain the switch density, 468 reach
the final target during step 1 and are ineligible, and 823 exceed 20% first-
step growth.  The remaining 1,989 routes correspond to 1,829 exactly distinct
instantaneous states and receive adaptive second-step searches.  Kinetic and
practical (`T2<T1`) boundaries remain separate.

## Main result

First-step preparation converts part of the previously kinetic-only 150–300
nm domain into complete practical windows.  Under the strict combination of
5% first-step growth and 5% second-step growth, complete practical counts are:

| mechanism | complete cases | minimum G0 | maximum width |
|---|---:|---:|---:|
| mech_009 q=1 | 106 | 150 nm | 135 C |
| mech_019 q=1 | 120 | 150 nm | 125 C |
| mech_009 q=0 | 54 | 225 nm | 25 C |
| mech_019 q=0 | 177 | 150 nm | 85 C |

Every counted case is `COMPLETE_WINDOW`: density exhaustion is bracketed
below, grain-growth failure is bracketed above, T2 is below T1, step 1 remains
below the target, and both independent growth tolerances are satisfied.
Upper-censored and kinetic-only cases are retained in separate tables.

## Preparation tradeoff

Increasing T1 is necessary for the strict q=1 practical windows: both q=1
survivors first appear at T1=1350 C and remain through 1500 C.  The largest
strict counts occur around 1400 C, then fall.  Higher T1 therefore does not
immediately destroy the window, but it increasingly destroys *admissible
preparation*: median first-step growth rises from 0.0036 at 1250 C to 0.316 at
1500 C.  The number of attained histories above 20% growth rises from zero at
1250 C to 341 at 1500 C.  The useful preparation domain is consequently a
finite compromise, not a monotonic benefit of higher T1.

The 0.2 C/min routes produce no strict 5%/5% practical successes.  Successful
routes span 1–100 C/min, but heating rate has negligible point-biserial
correlation with success once the first-step state is included (`r≈0.01`).
Heating rate acts through preparation history rather than as a label in the
mechanism.

## Which state predicts practical success?

For routes meeting 5% preparation growth and evaluated at 5% second-step
growth, simple point-biserial correlations with complete practical success
rank as follows:

| descriptor | correlation magnitude |
|---|---:|
| T1 | 0.54 |
| connected TJ coverage C_TJ | 0.51 |
| connected GB-segment coverage C_GBseg | 0.48 |
| clean-GB fraction | 0.48 (negative) |
| G1 | 0.45 |
| Lambda_TJ/K_TJ | 0.32 |
| rho1 | 0.07 |
| X_J | 0.07 |
| heating rate | 0.01 |

Thus connected topology and prepared grain size predict the practical window
better than density or the scalar persistent population alone.  T1 remains a
strong proxy because it controls both boundary bracketing and the prepared
topology.  These are univariate associations, not a causal calibration.

## q=0 alternatives

Both fixed-packet q=0 variants remain viable under `T2<T1`; they are not
silently discarded.  `mech_009_q0` is narrow and shifts its strict onset to
225 nm.  `mech_019_q0` produces the largest number of strict complete cases,
but its maximum width (85 C) is smaller than q=1.  This reverses any simple
claim that q=1 is always preferable and reinforces the need for experimental
TJ reaction-count measurements.

## Nanoscale limit

No mechanism produces a complete practical window below 150 nm when both
preparation and second-step growth are restricted to 5%.  Relaxing only the
second-step tolerance to 10% admits 100 nm cases for `mech_019` and
`mech_019_q0`, but not a strict 5% window.  At 10–20% preparation tolerance,
additional 75–100 nm cases appear.  The 75–100 nm region is therefore still
marginal and tolerance-sensitive, not a robust Chen-style domain.

## Censoring and guardrails

Practical maps contain both lower-bound-right and upper-bound-right censoring,
while kinetic maps retain the extrapolative 1550 C cap.  No censored record is
counted as a complete practical window.  Targets reached during step 1 are
explicitly rejected.  First- and second-step growth fractions are stored and
scored independently.  The 96 h budget, target, and mechanism parameters are
identical for every path.  The fixed-budget ramp also reports whether nominal
T1 was reached before the density switch, avoiding hidden schedule identity.

## Outputs

All requested tables and plots are under `results/preparation_window_search/`.
The main tables are `first_step_preparation_states.csv`,
`adaptive_second_step_boundaries.csv`, `practical_two_step_windows.csv`,
`kinetic_only_windows.csv`, `censored_preparation_cases.csv`,
`preparation_rejected_cases.csv`, and `mechanism_preparation_summary.csv`.
The complete calculation took about 18,878 s in the available single-process
environment.
