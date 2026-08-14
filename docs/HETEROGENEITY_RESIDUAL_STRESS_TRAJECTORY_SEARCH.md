# Heterogeneous Initial State and Residual-Stress Trajectory Search

## Scope

This prototype tests two explicitly ablatable additions to the frozen
PR/de-sintering production model: weighted initial microstructure cohorts and
local residual-stress states. It retains the 96 h budget, production density
targets, densification kinetics, TJ definitions, and PR baseline
(`PR_plus_connected_fine_attrition`, `k_PR_ref_s=2e-4`). All four frozen bases
remain registered. No schedule name, ramp label, target, or response class
enters a local closure.

The cohort layer resolves narrow, broad, tailed, bimodal, broad-grain,
correlated pore/grain, and defect-rich initial states. Every cohort experiences
the same thermal history. Aggregation reports `G_mean`, `G50`, `G90`, pore
`D50/D90`, pore number, connected fine pores, and GB-segment/TJ/isolated
locations. Residual stress has separate GB-segment, TJ, large-pore, and
crack-like states. It alters activation or conservative redistribution but is
not itself a pore-volume source.

## Design and observable criterion

The bounded campaign contains 144 registered parameter sets: 28
distribution-only cases, 20 stress-only cases, and 96 combined
fractional-factorial cases. Reference paths use 1 C/min and are never silently
mixed with other references. Fast paths use 20 or 100 C/min. Targets are
0.85, 0.88, 0.90, and 0.92. A candidate qualifies only if
`G_mean,reference/G_mean,fast >= 1.5` continuously over `Delta rho >= 0.03`
with both paths attaining the interval. Tail-only qualification is recorded
but cannot establish a mean-grain result.

## Findings

No candidate qualifies.

| Stage | Cases | Fully attained | Maximum `G_mean` ratio | Meaningful |
|---|---:|---:|---:|---:|
| Initial distribution only | 28 | 5 | 1.103 | 0 |
| Residual stress only | 20 | 3 | 1.187 | 0 |
| Combined | 96 | 7 | 1.662 | 0 |

The combined maximum (`C091`) is instructive but rejected. Its ratio exceeds
1.5 only from density 0.850 to 0.858, a span of 0.008, and the slow path does
not attain the complete 0.85--0.92 interval. At density 0.85 the slow/fast pore
`D90` values are 450/218 nm, connected fine-pore fractions are 0.024/0.057,
large-pore fractions are 0.979/0.935, and cumulative PR works are
503,000/29,400 model units. Thus the combined closure creates observable pore
memory and a brief grain response, but not a sustained trajectory effect.

Maximum ratios across all sampled points are 1.663 for `G_mean`, 1.639 for
`G50`, and 1.696 for `G90`. No case qualifies only through `G90`; the failure
is the finite-span/attainment requirement, not merely the selected grain
statistic.

No case was promoted to the adaptive Chen preservation calculation because no
case passed the mandatory fast-firing gate. `Chen_preservation_summary.csv`
therefore contains no scored window rather than a fabricated or censored
success. The already committed baseline Chen maps remain unchanged.

## Required interpretation

1. **Starting distributions alone?** No. Maximum ratio 1.103.
2. **Residual stress alone?** No. Maximum ratio 1.187.
3. **Combined larger than either?** Yes locally (1.662), but only over density
   span 0.008 and therefore not meaningful.
4. **Which observables respond?** Mean, median, and tail grain metrics respond
   briefly; pore D90, large-pore fraction, connected fine pores, and PR work
   show the clearest memory.
5. **Does fast heating preserve fine connected pores?** In the strongest
   rejected case, yes, while also suppressing the large-pore tail.
6. **Does slow heating accumulate stress-assisted PR damage?** Yes in selected
   combined cases, including C091, but this does not persist over sufficient
   density span.
7. **Chen preservation?** Not re-scored because no candidate passed the
   prerequisite trajectory gate; baseline Chen results are preserved.
8. **Below density 0.92?** All scoring is at or below 0.92. Nothing is inferred
   at 0.95 or above.
9. **Calibration measurements?** Interrupted-ramp 3D pore D50/D90 and location
   fractions, grain EBSD/TEM distributions (mean/G50/G90), diffraction-based
   residual stress separated by pore-rich region, TJ strain mapping, and
   matched-density dilatometry are the most direct constraints.

The current physics produces plausible internal microstructural memory but
still does not create an experimentally meaningful grain-size--density
trajectory separation below rho=0.92.

