# Final fast-firing and two-step mechanism synthesis

## A. Executive conclusion

Fast firing and two-step sintering are compatible in the current framework, but they arise from different dominant mechanisms.

Fast firing is controlled primarily by nucleation-limited material kinetics. A slow ramp spends more time in a low-activity interval in which grain growth can proceed without equivalent density gain; a fast ramp crosses that interval more quickly. Two-step behavior in the current best candidate is controlled primarily by Plateau-Rayleigh (PR)-prepared closed-pore accommodation memory. The first step creates a high-density support state; the second-step temperature then selects between density exhaustion, useful densification with bounded growth, and thermally activated grain growth.

The overlap is a finite relative-property region, not a unique tuned point. From 1,903 unique exact-promoted cases, 485 are fast-only, 119 are two-step-only, 73 pass both, and 1,226 pass neither. The dimensionless screen predicted 19,880 both-pass cases, or 272.3 times the exact count. Surrogate screening is therefore useful only for promotion and is not final evidence.

Candidate 693168 remains a **conditional Tier B mechanism candidate, not validation or a paper-ready calibration**. Its large high-temperature/two-step grain-size separation is the intended experimental-scale signature when both paths attain the interval and states remain bounded.

No model physics, topology parameters, material parameters, or time budgets were changed for this synthesis. No new search was run.

## B. Fast-firing mechanism

The fast-firing response requires a finite nucleation-onset window.

1. Nucleation must be difficult enough to create a low-activity waiting interval during heating.
2. The slow ramp must spend enough time in that interval for non-densifying growth or redistribution to accumulate.
3. Nucleation cannot be so difficult that the fast path fails to densify.
4. Once nucleation activates, exchange and transport must remain fast enough for both paths to attain the matched-density interval.

The exact base result reaches `G_ref/G_fast = 1.796` and remains above 1.5 over a density span of 0.17. Among the 1,000 exact fast-firing promotions, 55.8% pass the full rule, only 19.3% pass after making nucleation facile, and 71.7% pass with PR disabled. Thus PR-off does not invalidate fast firing; the nucleation-facile response is the causal ablation.

The local OAT window is finite and nonmonotonic. `Delta Q_nuc = 0` through `+50 kJ/mol` survives, while `-25` and `+75 kJ/mol` fail. Lower barriers make onset too facile; sufficiently higher barriers lose the required attainable trajectory separation. This is a relative-onset result, not a claim that one absolute activation energy applies to all materials.

## C. Two-step mechanism

Candidate 693168 attributes the two-step response to a PR-prepared closed/accommodation topology:

1. PR/surface evolution during the first step moves pore inventory into a persistent closed-pore store.
2. At low T2, closed-pore shrinkage and finite accommodation are insufficient to attain the target density, producing the lower boundary.
3. At intermediate T2, closed-pore densification remains active while migration is sufficiently suppressed, producing the success band.
4. At high T2, grain growth activates strongly enough to produce the upper boundary.

The fine exact T2 map gives a complete success band from 925 to 1205 °C, bracketed below by density exhaustion and above by grain-growth failure. The base exact property audit gives a 250 °C window because it uses the coarser common promotion protocol; both sources show the same boundary topology.

The destructive ablations for 693168 are no PR damage, no closed transition, no closed shrinkage, and infinite closed accommodation. These establish the causal chain inside the current model. TJ multihit, attached pore drag, residual stress, and sweep/coalescence can matter for other Tier-B candidates but are not the primary causes of 693168's strong trajectory separation.

The local material-property constraints are:

- `Delta Q_closed = -25` through `+100 kJ/mol` survives. At `-50 kJ/mol`, the lower boundary disappears, so the case is not a complete Chen window.
- `Delta Q_growth` remains viable across the full tested `±100 kJ/mol`; its activation-energy limit was not found.
- The PR prefactor must be at least `0.3x` in the tested OAT; `0.1x` fails.
- The growth prefactor must be at least `0.1x` to preserve the upper boundary; `0.03x` loses it.

Candidate 693168 approaches a very high closed fraction and derives nearly all high-density support from closed shrinkage. That makes its accommodation state the principal calibration concern, not proof of quantitative validity.

## D. Relative-property inequalities

The useful conditions are qualitative timescale inequalities rather than unique absolute barriers.

### 1. Nucleation-onset condition

At low-to-intermediate temperature, nucleation waiting must dominate enough to create history:

`tau_nuc / (tau_exchange + tau_transport) > O(1)` over a finite heating interval,

while remaining finite enough that the fast path activates. The 73 exact joint successes occupy `Theta_nuc = 23.4` to `2.96e5`; this broad, coverage-limited range is not a universal threshold.

### 2. Attainment condition

After activation:

`tau_exchange + tau_transport < available processing time`,

so both paths reach the scored density interval. An unattained ratio is never scored.

### 3. Closed-accommodation selectivity

At a useful second-step temperature:

`tau_closed shrinkage < processing time < tau_grain growth`,

or equivalently the effective closed-shrinkage/growth selectivity must be high enough to densify before significant migration. Joint exact successes occupy `S_closed/growth = 0.116` to `126`, again as an observed promoted range rather than a universal bound.

### 4. Lower-bound condition

Closed shrinkage and accommodation must become insufficient at lower T2. If accommodation is unlimited or closed shrinkage remains successful everywhere, the density-exhaustion boundary disappears and there is no finite Chen window.

### 5. Upper-bound condition

Grain-boundary migration must reactivate strongly enough at high T2 to exceed the growth tolerance. If growth is suppressed at every T2, universal success replaces the bounded processing window.

### 6. PR-preparation condition

PR/surface redistribution must be strong enough to prepare the closed/accommodation state. The tested PR prefactor threshold is `0.3x` base. It must not destroy all densification or remove attainable high-density support.

## E. Generality to crystalline particle systems

The two effects do not require one absolute activation energy. They require relative orderings among nucleation, exchange, transport, PR/surface evolution, closed-pore shrinkage/accommodation, and grain-growth/migration.

These ingredients are plausible in crystalline systems because densification requires serial activated source/sink renewal, exchange, and transport, whereas grain growth and migration can be limited by different topology and accommodation processes. A material can therefore show fast firing through nucleation-limited onset without requiring PR as its causal fast channel, while its two-step response can still depend on PR-prepared closed topology.

The six exact Tier-B base candidates demonstrate a qualitative family, not a calibrated material class. Candidate 693168 is the strongest separation comparator; 822940 provides a lower-closed-fraction comparator; 581668 provides an intermediate topology. The family material-window comparison remains a reduced transfer from exact 693168 OAT, not six independent exact material searches.

## F. Experimental falsification

The proposed interpretation should be tested with:

- matched-density `G_mean`, `G50`, and `G90` for multiple heating rates and two-step schedules;
- densification-onset time during interrupted ramps;
- independent exchange and transport relaxation measurements;
- three-dimensional open/closed pore fraction at the first-step switch and during T2 holds;
- pore D50, D90, and large-pore-tail evolution;
- connected fine-pore fraction and percolation;
- trapped-gas pressure or another closed-pore accommodation proxy;
- high-temperature grain-growth mobility across the upper T2 boundary;
- interrupted first-step and second-step tomography to test persistence of the prepared state.

The strongest falsification would be either (a) a nucleation-facile material retaining the same fast-firing separation, or (b) candidate-like two-step behavior without a persistent schedule-dependent closed-pore/accommodation state. The primary calibration targets remain the closed-pore fraction and accommodation trajectory.

Machine-readable evidence is in `results/final_mechanism_synthesis_and_property_windows/source_tables/`. Only exact-promoted rows support final classifications; surrogate counts are labeled separately.
