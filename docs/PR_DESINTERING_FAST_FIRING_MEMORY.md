# PR/De-sintering Fast-firing Memory

## Scope

This is a mechanism prototype, not validation. The production negative control
at `eaac130` was pushed before this branch was created. Its four frozen
persistent-junction/TJ-multihit bases, density targets, and 96 h budgets were
not retuned. The new closure is selected by `early_memory_mode`; `disabled`
delegates exactly to the prior solver.

The local law uses only temperature, renewal activity, connected pore-bin
volume, pore radius, and connected-topology availability. It allocates
competition between densifying removal and PR/smoothing action. PR action
transfers volume conservatively from fine GB-segment pores toward larger bins,
TJ sites, and isolated sites; it never removes pore volume. No protocol name,
heating-rate label, target, or schedule class enters the closure.

## Bounded design and escalation

The reduced screen compared disabled, PR competition, moderate connected-fine
attrition, and strong attrition for all four frozen bases. The disabled modes
retained 9–21 complete practical Chen windows and produced exactly zero
beneficial fast-firing cases, recovering the production negative result.

Twelve enabled candidates passed both gates. The moderate variant was selected
for the full fast grid because it is less aggressive than the strong variant
and was already non-universal. All four moderate candidates remained joint
positive:

| candidate | complete Chen windows | beneficial full-grid cases |
|---|---:|---:|
| mech_009 | 16 | 1,437 |
| mech_019 | 19 | 655 |
| mech_009 q0 | 8 | 1,518 |
| mech_019 q0 | 21 | 605 |

The beneficial subsets coexist with harmful, neutral, and unattainable points;
the mechanism therefore creates a finite response region rather than a
universal fast-firing bonus. Both compared paths must attain the target. The
0.2 C/min path is used when attainable, otherwise the declared 1 C/min
reference is selected independently for each target.

## Physical result

Across beneficial full-grid cases, median `HR_pct` is 3.04–3.35%, with maxima
of 10.4–23.1% depending on the frozen base. At matched density the fast path
accumulates about `1.1e6–1.6e6` less PR work and retains approximately
`0.0064–0.0070` more connected fine-pore fraction. Its connected-pore mean
radius is lower by roughly 0.7–1.0 nm. These are observable pore-distribution
memories, not inferred efficiency multipliers.

The response is strongest at 1400–1500 C, for 8–20 h holds, and for 50–100 nm
starting grains. It appears in all initial topologies, but the base matters:
`mech_009` favors GB-segment-rich states more strongly, while `mech_019`
responses are more evenly distributed among TJ-rich, baseline, and mixed
states. Benefits occur at every tested attainable density from 0.85 to 0.92.

The two processing routes remain mechanistically distinct. Fast firing wins by
crossing the low-renewal PR-active interval quickly and preserving removable
connected pores. Two-step sintering still wins through the existing prepared
`X_J`/TJ-multihit migration suppression and its lower densification-exhaustion
and upper grain-growth boundaries. PR redistribution was not inserted into the
grain-growth closure.

## Stress and energy accounting

The prototype reports instantaneous and cumulative densifying and
non-densifying work, including PR surface-energy expenditure. It deliberately
does not add a new accumulated stress state: doing so without an independent
elastic-storage closure would add non-identifiability. PR flux changes no
density directly, and tests enforce exact `rho = 1 - sum(phi)` and nonnegative
pore and number bins.

## Answers and limitations

1. **Does PR competition create fast-firing benefit?** Yes, in bounded,
   attainable regions for every frozen base; disabled modes remain negative.
2. **Are practical Chen windows preserved?** Yes, 8–21 complete reduced-map
   windows per moderate candidate, with both boundaries and `T2 < T1`.
3. **Are connected pores preserved?** Yes at matched density in successful
   cases: more fine connected volume and a smaller mean connected-pore radius.
4. **Sensitivity?** Peak temperature and hold duration are influential;
   benefits concentrate at small `G0` and vary with initial pore placement.
5. **Different pathways?** Yes: early pore-removability preservation for fast
   firing versus late migration suppression for two-step sintering.
6. **Measurements?** Interrupted ramps should measure pore-size distributions
   separated by GB segment/TJ/isolated location, connected fine-pore coverage,
   and surface-area loss at matched temperature and density. Simultaneous
   dilatometry and grain-size measurements would test whether PR exposure rises
   before renewal densification becomes active.

The main uncertainty is the absolute PR rate and its partition among smoothing,
GB-to-TJ transfer, and isolation. Those parameters are physical hypotheses,
not calibrated values. A decisive experiment should constrain the interrupted-
ramp pore redistribution before any further search or validation claim.
