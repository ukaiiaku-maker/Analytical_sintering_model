# Local Connected-Sink Defect Memory

## Hypothesis and implementation

This branch replaces the previous cohort-wide defect eligibility penalty with
an explicit local mixture. `matrix_connected` and
`defect_rich_large_pore` subpopulations evolve independently under the same
thermal schedule. Each carries its own density, grain size, pore distribution,
location fractions, PR work, and residual/defect state. Global density is the
normalized weighted average of local densities, and global mean/median/tail
grain and pore metrics are computed from the weighted populations.

In the localized modes, persistent defect memory is never applied to the
matrix. The matrix can continue ordinary connected-sink densification while
the defect-rich cohort stores damage. Residual stress and defect state are not
pore-volume sources. Isolated stores are not removed by the open-pore
densification channels.

The disabled mode recovers the permitted non-persistent heterogeneous
baseline. A uniform 1,000-adaptive-step numerical censor applies to active
residual/persistent cohorts; ordinary 96 h paths require fewer steps. Censored
paths are explicitly rejected and cannot be scored.

## Bounded screen

The 96-case fractional design produces 192 fast/reference comparisons and
crosses all requested factors and ablations. No case changes the 96 h budget,
targets, or observable criterion.

| Mode | Comparisons | Jointly attained | Numerically censored | Maximum ratio | Meaningful |
|---|---:|---:|---:|---:|---:|
| disabled local mixture | 32 | 10 | 0 | 1.021 | 0 |
| matrix only | 32 | 14 | 0 | 1.005 | 0 |
| static defect mixture | 32 | 0 | 20 | 1.018 | 0 |
| evolving defect memory | 32 | 0 | 20 | 1.758 | 0 |
| evolving memory + matrix densification | 32 | 0 | 20 | 1.034 | 0 |
| high stress retention | 32 | 0 | 20 | 1.150 | 0 |

All 192 comparisons are rejected: 88 for reference nonattainment, 80 for
numerical censoring, and 24 for ratios below 1.5. The only local ratios above
1.5 occur in numerically censored evolving-defect cases. They are not evidence
of sustained separation. No candidate passed the trajectory gate, so no Chen
calculation was run.

## Required questions

1. **Does localization allow both paths to attain 0.85--0.92?** It improves
   attainment for matrix-only and disabled controls, but no active persistent
   mode jointly attains the full interval.
2. **Sustained mean-grain separation?** No.
3. **Mean, median, tail, or tail only?** None qualifies. Above-threshold mean
   ratios occur only on censored paths; tail metrics cannot rescue them.
4. **Does matrix densification continue?** Yes locally. However, the weighted
   global path remains limited by the defect-rich volume and fixed schedule.
5. **Finite-span or attainment artifact?** Still an attainment/numerical-censor
   artifact, not a finite jointly attained effect.
6. **Does residual stress persist?** Defect cohorts retain stress longer by
   construction, but sustained stress does not produce a valid global path.
7. **Chen preservation?** Not scored because no fast-firing candidate passed
   the prerequisite gate.
8. **Discriminating measurements?** Spatially registered nano-CT/FIB-SEM pore
   topology, local diffraction stress, grain-size distributions around pore
   clusters, and local strain/densification maps during interrupted ramps would
   distinguish local defect memory from global densification suppression.

Persistent defect memory must be either spatially stronger than this reduced
mixture, coupled to late-stage closed-pore physics, or replaced by a different
experimentally constrained mechanism.

