# Production Mechanism Assessment

## Scope and guardrails

This is a production-scale stress test of four frozen candidates (`mech_009`,
`mech_019`, and their q0 ablations), not a calibration or validation claim.
Mechanism parameters, density targets, and the 96 h budget are common to every
schedule. Failed targets remain censored. The practical two-step definition
requires `T2 < T1`, both kinetic boundaries, at least 25 C window width, and
the stated preparation and second-step growth tolerances.

The calculation covered 1,440 first-step trajectories, 6,585 admissible
routes, 6,306 unique instantaneous states, and 34,560 fast-firing
trajectories (138,240 target records). A production integration ceiling of
900 s was used. A six-case check against 180 s gave maximum density difference
`1.64e-4`, maximum relative grain-size difference `0.223%`, and maximum state
difference about `5.66e-4`; this numerical approximation is adequate for map
classification but should be tightened before parameter inference.

## Practical Chen-style windows

| mechanism | Tier A | Tier B | minimum Tier-A G1 (nm) | maximum Tier-A width (C) |
|---|---:|---:|---:|---:|
| mech_009 | 570 | 668 | 150.90 | 145 |
| mech_019 | 641 | 853 | 151.08 | 160 |
| mech_009_q0 | 0 | 153 | -- | -- |
| mech_019_q0 | 497 | 1,132 | 182.56 | 85 |

All scored successful windows have both lower densification-exhaustion and
upper grain-growth bounds. Thus the current family can generate finite,
practical Chen-map topology. The q0/q1 comparison is not uniquely identifiable
from processing windows: removing the q term destroys Tier A for `mech_009`
but retains many Tier-A windows for `mech_019`.

## Fast-firing map and joint test

No frozen candidate produced a beneficial fast-firing comparison at any
attainable target. Harmful counts were 6,098, 6,311, 6,400, and 6,385 for
`mech_009`, `mech_019`, `mech_009_q0`, and `mech_019_q0`, respectively; the
remainder were neutral or unattainable. The 0.2 C/min reference never reached
a comparable target, so every valid comparison used the declared 1 C/min
fallback. Harmful responses occur in all four starting topologies. Therefore
`joint_positive` is false for every candidate.

This is evidence that the frozen family is sufficient for the two-step
window, but insufficient for coexistence with a fast-heating advantage. The
failure is physically interpretable: persistent junction/multihit terms can
suppress migration after preparation, yet the same accumulated resistance is
at least as strong during slower firing. There is no sufficiently strong,
observable slow-ramp degradation of subsequent densification efficiency.

## Minimum ingredients and remaining ambiguity

The present evidence supports geometric eligibility, pore-size/location
memory, separate densification and migration channels, persistent junction
state, high-temperature TJ multihit reactivation, adaptive preparation, and
censor-aware boundary classification. Redistribution, a density gate,
mean-grain junction resistance, connected pinning, or action weights alone
have each proved insufficient.

The next isolated mechanism test should couple slow intermediate-temperature
exposure to an observable loss of removable connected fine pores (or to a
measured stress accumulation/release state) while leaving densification
kinetics untouched in the ablation. It must recover this frozen negative
control when disabled and must not be an arbitrary growth multiplier.

Useful discriminating measurements are interrupted-ramp pore-size/location
distributions, connected boundary coverage, triple-junction occupancy,
grain-growth onset versus temperature, and matched-density grain size across
several heating rates. These would constrain the missing memory channel more
directly than final density alone.

## Figure guide

1. Practical Chen map: Tier A/B/C states in `T2` versus prepared `G1`.
2. Kinetic versus practical window widths.
3. First-step temperature/switch-density preparation tradeoff.
4. Fast-firing HR map by initial topology; no positive region appears.
5. Representative histories of `X_J`, multihit completion, connected stores,
   and clean-boundary fraction.
6. Densification, clean-GB, persistent-junction, multihit power channels and
   total activation stress.
7. Mechanism-ingredient chain used by the frozen candidates.

The compact CSVs are review artifacts; ignored `raw_*` files retain the full
point-level evidence locally. This result is a negative-control baseline and
does not justify validation or merging the draft PR.
