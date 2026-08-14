# Reduced-Model Falsification Summary

## Controlling conclusion

The source-grounded migration-suppression family reproduces finite practical
Chen-style two-step windows. No tested reduced mechanism reproduces an
experimentally meaningful fast-firing mean-grain trajectory under the
controlling criterion: `G_reference/G_fast >= 1.5` continuously over a jointly
attained `Delta rho >= 0.03`.

The earlier production `HR_pct > 1%` result is therefore mechanistically useful
but not a paper-ready fast-firing result. The matched-density observable audit
supersedes that interpretation.

## Mechanism outcomes

| Mechanism | Branch / commit | Added physics | Outcome | Maximum attempted ratio | Density support | Chen? | Meaningful fast trajectory? | Disposition |
|---|---|---|---|---:|---|:---:|:---:|---|
| Aggregate PR #2 control | `codex/topology-constrained-mechanisms` | topology-constrained fluxes | negative fast control | 1.00 | to 0.92 | no | no | baseline negative control |
| Pore-placement topology | `codex/pore-placement-topology-search` | GB/TJ/isolated stores | no nanoscale window | 1.00 | to 0.90 | no | no | negative control |
| Pore-location action | `codex/pore-location-action-layer` | competing local actions | negative | 1.00 | to 0.90 | no | no | abandon as sole mechanism |
| Persistent junction + TJ multihit | `codex/agentic-mechanism-discovery` | junction drag/multihit | bounded Chen windows | 1.00 | to 0.90 | yes | no | retain Chen baseline |
| Practical preparation map | commit `6da2500` | censor-aware first/second steps | 150–300 nm practical windows | 1.00 | to 0.90 | yes | no | retain processing map |
| Production PR memory | commits `3ead9c9`, `62394cf` | conservative pore memory | internal memory, small HR regions | 1.301 | 0.85–0.92 | yes | no | mechanism baseline only |
| Observable trajectory audit | `4631112` | matched-density effect sizes | fast-firing claim falsified | 1.301 | 0.85–0.92 | yes | no | controlling negative result |
| Heterogeneity + residual stress | `40b7bb9` | weighted cohorts/stress | 1.662 over only 0.008 density | 1.662 | reference to 0.859 | no recheck | no | negative control |
| Persistent defect memory | `4f80b19` | stored defect/work/stress | 3.73 in unattained interval | 3.730 | incomplete 0.85–0.92 | no recheck | no | negative control |
| Local connected-sink mixture | `69d61e1` | independent matrix/defect regions | 1.758 only censored | 1.758 | incomplete | no recheck | no | negative control |
| Late-stage closed pores | `8ce78ef` | closure, vacancy transport, gas | five fast paths to 0.95; none to 0.98 | 2.024 raw; 1.110 attained 0.90–0.95 | joint below 0.95 | baseline preserved | no | late-stage negative control |

Full machine-readable fields, including q0/q1 visibility and TJ
pore/constraint separation, are in `mechanism_scorecard.csv`.

## Two distinct scientific targets

### A. Chen-style two-step behavior

The successful family has a lower second-step boundary set by densification
exhaustion and an upper boundary set by grain-growth activation. It enforces
the practical `T2 < T1` condition and requires both boundaries to be bracketed.
Persistent junction drag and TJ multihit completion, followed by the
censor-aware preparation search, produce selected practical windows in the
150–300 nm range.

### B. Fast-firing `G(rho)` separation

This requires a much stronger observable: ratio at least 1.5 over density span
at least 0.03, with both paths attaining every compared density. None of the
PR-memory, heterogeneous, persistent, local-mixture, or closed-pore closures
passes. High raw ratios arise only over short spans or in unattained/censored
paths.

## Why Chen-window success is not enough

A finite two-step window asks whether a lower-temperature restart can retain
densification while suppressing migration between two kinetic boundaries. A
continuous-heating trajectory test asks whether prior thermal history creates
a large, persistent difference in grain size at the same density. These are
different counterfactuals. A model can separate diffusion from migration well
enough to create a bounded restart window while its continuous-ramp memory is
too small, too brief, or too damaging to density attainment to create a large
`G(rho)` separation.

## Failure synthesis

The sequence of negative controls is informative:

- topology and pore placement move boundaries but do not establish continuous
  trajectory memory;
- PR memory visibly changes pore D90, connected fine pores, and stored work,
  but the mean-grain response remains weak;
- heterogeneity and residual stress amplify local ratios only briefly;
- persistent/global defect penalties destroy joint density attainment;
- localizing damage lets the matrix densify but still leaves active defect
  cases censored or weak;
- explicit closed-pore transport permits limited 0.95 attainment but does not
  sustain differential grain growth, and no path supports 0.98–0.99.

This does not prove that fast-firing separation is absent in real materials.
It falsifies the tested class of compact reduced closures as a sufficient
explanation under the specified observable criterion.

## Current decision

Do not add another scalar closure. First establish experimentally whether the
large mean-grain separation exists below density 0.92, near final density, or
primarily in distribution tails. If confirmed, seek user approval for a
minimal spatial pore-network or interacting local-region model. See
`NEXT_MODEL_CLASS_DECISION.md`.

