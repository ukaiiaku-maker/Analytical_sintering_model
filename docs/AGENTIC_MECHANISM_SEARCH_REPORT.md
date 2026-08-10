# Source-grounded agentic mechanism search

## Status

This is a mechanism-discovery result, not a calibrated or validated model.
The branch preserves the aggregate, pore-location, and first action-layer
negative controls.  It does not change densification kinetics, pore
redistribution, `rho_target=0.90`, or the fixed 96 h step budgets.

The bounded search finds an interesting coupled mechanism: a persistent
junction obstacle population combined with a Class-B, accommodation-demand
TJ multihit closure produces finite 5–10% Chen-style windows from 150–300 nm
and above.  Isolated successes at 75–100 nm have zero sampled width and are
not claimed as robust windows.

## Discovery decision

The source priors are recorded in `SOURCE_MECHANISM_PRIORS.md`; all candidates
and their conservation/density roles are declared in `mechanism_registry.py`.
The previous failure required migration resistance that does not disappear
instantaneously with connected pore coverage and that separates GB migration
from GB-diffusion densification.  The first iteration therefore selected:

1. `persistent_junction_drag` (Class A): a bounded state `X_J` is produced by
   instantaneous TJ-connected densification events, relocation/capture, and
   boundary sweep, and decays by thermal relaxation.  It changes migration
   only.  Unrealized migration power is `P_persistent_junction_drag`.
2. `TJ_multihit_reaction` (Class B): a packet completes with
   `P_comp=Pr[Poisson(Lambda_TJ)>=K_TJ]`.  Both fixed packet (`q=0`) and
   accommodation-demand (`q=1`, `K` increasing with grain size) are explicit
   modes.  It changes migration only and reports `Lambda`, `K`, `Lambda/K`,
   completion probability, and smooth/intermittent/stagnant regime.

Vacancy multihit, stress accumulation/release, and Class-C exchange remain
registered but were not stacked into this iteration.  This preserves
identifiability of the first positive map.

## Search and rejection accounting

The deterministic 64-set design varies only named physical parameters:
`A_J`, junction lifetime, junction production, TJ event count, required packet
count, event activation energy, and the stated `q_TJ` limit.  Every set is in
`parameter_registry.csv`.  The reduced map uses 75/150/300/600 nm,
`T1=1300 C`, switch density 0.80/0.85, and five T2 values.  Rejections are
persisted with reasons.

Eighteen sets are rejected: eleven persistent-only cases have no success
window and seven multihit/combined cases lose the upper grain-growth boundary.
No set is accepted without both a lower densification-exhaustion and an upper
grain-growth failure.  Forty-six sets pass the reduced logical screen; only
the top two are escalated, preventing a large full-map sweep.

The two escalated cases are both coupled `persistent_tj_multihit_q1`:

- `mech_009`: `A_J=12`, `tau_J_ref=1.35e6 s`, `XJ_prod_TJ=1`,
  `lambda_ref=0.5`, `K0=1`, `Q_event=340 kJ/mol`;
- `mech_019`: `A_J=48`, `tau_J_ref=4.5e5 s`, `XJ_prod_TJ=2`,
  `lambda_ref=1.5`, `K0=2`, `Q_event=380 kJ/mol`.

These are bounded hypothesis parameters, not fitted values.  The factor-four
difference in `A_J` and different event packet definitions giving similar map
topology demonstrate remaining non-identifiability.

## Full Chen-style maps

The full map uses the specified 50–600 nm initial grains, three T1 values,
three switch densities, T2 from 900–1300 C in 25 C increments, and both 5%
and 10% growth tolerances.  It contains 4,896 classifications.

`mech_009` produces 120 successes at each tolerance, 37 grain-growth failures,
and more than 1,000 densification-exhaustion failures.  Finite windows occur
for every tested first-step group from 150–600 nm; maximum widths are 50 C at
150–450 nm and 75 C at 600 nm.  There are no successes below 150 nm.

`mech_019` produces 125 successes at 5% and 126 at 10%, with 26–27
grain-growth failures and more than 1,040 densification-exhaustion failures.
Finite 50 C windows occur at 150–450 nm and a 75 C window at 600 nm.  The
75–100 nm successes have zero sampled width and occur only for a few
first-step groups; they are marginal points, not a nanoscale processing
window.

Thus the coupled closure passes the narrow scientific target in the
150–300 nm range without universal success.  It retains a lower T2 boundary
where density exhausts and an upper boundary where growth activates.  The
result is qualitatively robust to the 5% versus 10% tolerance for these two
cases—the identical or nearly identical success counts imply a sharp
activation separation rather than tolerance-driven scoring.

## Mechanism interpretation

Persistent drag alone is insufficient in most of the bounded design: eleven
of thirteen persistent-only cases have no window.  Conversely, overly strong
multihit suppression produces successes but erases the upper growth boundary
in seven cases.  The surviving `q=1` coupled cases occupy the middle regime:
stored junction obstacles persist after connected coverage declines, while
the temperature-sensitive multihit completion probability rises enough at
high T2 to reactivate migration.  This supplies the required separation:
GB-diffusion densification remains unchanged, low-T2 densification can still
fail, and high-T2 migration can still activate.

The nominal interrupted ramp demonstrates that `X_J` is a real state rather
than a label: it evolves from 0 to about 0.117 using local event, relocation,
and sweep rates.  Nominal TJ multihit histories span `Lambda/K` from nearly
zero to about 2.1–2.4 and completion probabilities from zero to 0.93–0.97,
crossing stagnant, intermittent, and near-smooth regimes.  These diagnostics
make the mechanism experimentally falsifiable.

## What remains unresolved

This is not validation.  The lower nanoscale onset is about 150 nm for finite
windows, while 50 nm has none and 75–100 nm only marginal points in the
stronger case.  `A_J`, junction lifetime, event-count prefactor, and required
packet scaling are not independently identifiable from processing maps alone.
The action model also treats `X_J` as a scalar population rather than a
measured junction-density/segment-length distribution.

The most useful calibration data are interrupted second-step measurements of
triple-junction density and GB segment length together with boundary/TJ pore
occupancy, followed by a restart measurement of grain-growth rate at the same
density and grain size.  This would constrain `X_J` production/relaxation
separately from `Lambda/K`.  If those data reject the persistent population,
the next registered mechanism should be stress accumulation/release.  If they
support persistence but not the `q=1` packet law, vacancy-accommodation
multihit and its `q_vac=1/2` limits should be tested next.  Closed-pore physics
is still required before extending claims toward full density.

## Anti-cheat and reproducibility

Local closures receive only instantaneous state, temperature, and fixed
parameters.  Tests verify exact baseline recovery, schedule-label absence,
Poisson completion limits, unchanged densification for migration-only modes,
nonnegative pore stores and numbers, exact `rho=1-sum(phi)`, conservative
relocation, exclusion of isolated pores from open-pore removal, strict
classification, ineligible already-reached targets, persisted rejection
reasons, and successful-point criteria.  The complete run took 2,171 s in one
worker after sandbox multiprocessing was unavailable.

All requested CSVs and deterministic plots are in
`results/agentic_mechanism_search/`.

## Adaptive boundary search and censored windows

The original 900–1300 C map was a classification grid, not a complete
per-state boundary search.  The correction retains that grid as stage 1, then
extends each promising state independently to 800 C when its lower edge is
left-censored and to a generic, explicitly model-extrapolative 1550 C cap when
its upper growth boundary is missing.  Only changing 25 C brackets are refined
at 10 C spacing.  Kinetic maps permit `T2>=T1` to diagnose the mathematical
window; practical maps require the literal two-step condition `T2<T1`.

### Boundary completion

The adaptive calculation evaluates 3,484 T2 points and produces 576 boundary
records.  In the kinetic map, no success window is upper-censored at 1550 C:
every reported success region has a densification-exhaustion point below it
and a grain-growth-failure point above it.  One 50 nm first-step state for
`mech_019` is lower-bound-right-censored because no density-attaining T2 is
found; it is not counted as a window.  All other non-window states are
reported as `NO_OVERLAP`, not silently discarded.

No prior 150–300 nm success window becomes censored after extension.  Instead,
the old 1300 C ceiling is shown to have truncated several upper boundaries.
For `mech_009`, complete 5% kinetic windows occur in all nine first-step groups
at 150, 225, and 300 nm, with maximum refined widths of 55, 75, and 135 C;
upper boundaries reach 1310, 1320, and 1375 C.  At 10%, the corresponding
maximum widths are 75, 90, and 150 C and upper boundaries reach 1325, 1335,
and 1395 C.  `mech_019` likewise retains all nine groups at 150–300 nm, with
5% widths up to 65/90/125 C and 10% widths up to 75/105/140 C.  Thus the
finite 150–300 nm kinetic result survives, but some of it lies above the
first-step preparation temperature and is diagnostic rather than practical.

### Practical two-step subset

The practical map is substantially narrower.  `mech_009` has complete 5%
windows only for three groups each at 150 and 225 nm; at 10% it additionally
has marginal 75–100 nm groups and three groups each at 150 and 225 nm.
`mech_019` has complete 5% practical windows in selected groups at 75, 100,
150, and 225 nm and complete 10% cases from 50–150 nm.  The 50 nm point has
zero sampled success width and remains marginal.  Across both candidates and
tolerances there are 38 complete practical boundary records, while many
otherwise successful kinetic states are `UPPER_BOUND_RIGHT_CENSORED` at the
`T2<T1` boundary or `LOWER_BOUND_RIGHT_CENSORED` because density is not reached
before that boundary.  They are not accepted as complete practical windows.

Consequently, the model produces both lower densification-exhaustion and upper
grain-growth failure over the expanded kinetic range, but only a subset is a
literal two-step processing window.  The kinetic and practical outputs are
kept in separate CSVs and plots.

### Local parameter refinement

Twenty-two one-at-a-time variants use factors 0.5/2 for persistent resistance,
production, and event count, factors 0.3/3 for junction lifetime, ±40 kJ/mol
for the event barrier, and explicit `q=0` alternatives.  All retain lower and
upper failures in the reduced map and remain `promising_reduced_only`; none is
promoted without its own adaptive full map.  The `q=1` family generally gives
more successes, but both visible `q=0` variants retain finite reduced windows
(4–7 and 7–9 points at 5–10%).  The entire local OAT neighborhood therefore
remains plausible, which reinforces parameter non-identifiability rather than
selecting a unique fit.

### Corrected answers

1. **Were upper boundaries bracketed?** Yes for every accepted kinetic
   window, after extensions reaching as high as roughly 1475 C at 600 nm.
2. **Did previous windows become censored?** No 150–300 nm kinetic window did;
   several practical subsets are censored by `T2<T1`.
3. **Are 150–300 nm windows finite?** Yes in the kinetic map, with both
   boundaries completed per state.
4. **Are they practical?** Only selected first-step groups, mainly through
   225 nm, are complete literal two-step windows.  Higher-G windows are often
   kinetic-only.
5. **Are both failure modes retained?** Yes over the adaptive kinetic domain;
   censor labels replace claims where the practical domain ends first.
6. **Which parameters remain promising?** Both survivor neighborhoods and
   both `q` limits retain reduced boundaries; `q=1` is stronger but not
   uniquely identified.

The new tables are `adaptive_T2_boundary_points.csv`,
`adaptive_window_boundaries.csv`, `censored_window_cases.csv`,
`practical_T2_less_than_T1_windows.csv`, `kinetic_window_map.csv`,
`practical_two_step_map.csv`, and `extended_parameter_refinement.csv`.
