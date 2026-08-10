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
