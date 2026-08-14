# Manuscript mechanism narrative: fast firing and two-step sintering

## Main conclusion

The reduced analytical framework supports fast firing and two-step sintering within one mechanistic program, but it does not attribute the two outcomes to a single universal rate multiplier. The exact-promoted evidence instead resolves two dominant kinetic pathways. Fast firing is governed by nucleation-limited onset: rapid heating shortens the time spent in a low-activity interval during which grain growth can accumulate without comparable density gain. The best current two-step response is governed by Plateau-Rayleigh (PR)-prepared closed-pore accommodation memory: the first step prepares a persistent high-density pore state, and the second-step temperature selects between densification exhaustion at low temperature, useful densification with bounded growth at intermediate temperature, and grain-growth activation at high temperature.

These mechanisms are compatible over a finite relative-property region. Of 1,903 unique exact-promoted perturbations, 485 pass only the fast-firing criterion, 119 pass only the two-step criterion, 73 pass both, and 1,226 pass neither. The overlap is therefore neither universal nor confined to one parameter point. By contrast, the dimensionless feasibility screen labeled 19,880 rows as both-pass, 272.3 times the exact count. Final classifications consequently use exact-promoted trajectories only; the screen is a promotion tool, not evidence.

Candidate 693168 is the strongest current two-step comparator. It remains conditional Tier B because its first-step growth exceeds the Tier-A bound, its closed-pore/accommodation state has not been calibrated against experiment, and its broad predicted response has not been quantitatively fitted. The result is a mechanism hypothesis with explicit falsification targets, not validation or a paper-ready calibration.

## Fast-firing mechanism

The fast-firing response arises from the ordering of three serial timescales: densifying-event nucleation, exchange, and transport. Before event nucleation activates, the system can remain in a low-density-gain regime while non-densifying growth or redistribution proceeds. A slow ramp dwells in that regime. A fast ramp traverses it more quickly and reaches the active densification regime with less accumulated grain growth. The experimentally relevant output is therefore not a rate at one temperature, but the matched-density trajectory ratio,

\[
R_{\mathrm{fast}}(\rho)=\frac{G_{\mathrm{reference}}(\rho)}{G_{\mathrm{fast}}(\rho)}.
\]

The exact base trajectory reaches a maximum ratio of 1.796 and remains above 1.5 over a density span of 0.17. The causal ablations support a nucleation-onset interpretation. Among 1,000 exact fast-firing promotions, 55.8% satisfy the full trajectory rule, only 19.3% retain it when nucleation is made facile, and 71.7% retain it with PR disabled. Thus PR can coexist with fast firing but is not the required causal channel in the current material envelope. Making nucleation facile removes the waiting interval and weakens or eliminates the trajectory separation.

The onset requirement is finite on both sides. In the exact one-at-a-time scan, shifts of `Delta Q_nuc = 0` through `+50 kJ/mol` retain the effect; `-25 kJ/mol` fails because nucleation becomes too facile, while `+75 kJ/mol` fails because the useful attainable separation is lost. This nonmonotonic window is consistent with a timing mechanism: nucleation must be slow enough to create history, but not so slow that the fast path cannot activate and densify.

Exchange and transport remain part of the necessary serial framework even though they do not dominate the local sensitivity. After nucleation, their combined time must be short enough for both heating paths to reach the matched-density interval. Ratios from unattained or unsupported intervals are never scored.

## Two-step mechanism

The candidate-693168 response is best described as PR-prepared closed-pore accommodation memory. During the first step, PR/surface evolution and the connected-to-closed transition prepare a closed pore inventory. That inventory persists into the second step and supports high-density shrinkage through a finite accommodation state. The second-step temperature then selects among three regimes.

At low `T2`, the available closed-pore shrinkage and accommodation cannot complete densification within the common time budget. The path stalls below the target and defines the lower boundary. At intermediate `T2`, closed shrinkage remains active while grain-boundary migration is sufficiently limited; the path reaches the target with bounded grain growth. At high `T2`, grain growth becomes thermally active enough to exceed the fixed tolerance and defines the upper boundary. The fine exact map places the successful band between approximately 925 and 1205 °C. The coarser common property-audit protocol reports a 250 °C band; both resolutions retain the same lower-exhaustion/success/upper-growth topology.

The destructive ablations identify the causal channels within the current model. Removing PR damage, the closed transition, closed shrinkage, or finite closed accommodation destroys the joint trajectory/window result. These findings do not make TJ multihit, pore drag, residual stress, or sweep/coalescence universally irrelevant. They show only that those channels are not the primary explanation for the strong candidate-693168 response; other members of the six-candidate Tier-B family retain different ablation dependencies.

The material-property scan likewise identifies a finite local window. `Delta Q_closed = -25` through `+100 kJ/mol` retains the complete response, whereas at `-50 kJ/mol` the lower boundary disappears. `Delta Q_growth` remains viable across the full tested `±100 kJ/mol`, so its activation-energy limit was not found. The PR prefactor must reach at least `0.3x` base in the tested scan; `0.1x` fails. The growth prefactor must reach at least `0.1x` to preserve the upper boundary; `0.03x` suppresses growth so strongly that the required upper bracket is lost.

## Relative-property inequalities

The mechanism is more naturally stated as a set of relative timescale conditions than as a unique set of activation energies.

### Nucleation-onset inequality

During low-to-intermediate-temperature heating, the nucleation time must be large relative to the exchange-plus-transport time over a finite interval,

\[
\Theta_{\mathrm{nuc}}=\frac{\tau_{\mathrm{nuc}}}
{\tau_{\mathrm{exchange}}+\tau_{\mathrm{transport}}}>\mathcal{O}(1),
\]

while remaining finite enough for the fast path to activate. The 73 exact joint successes occupy an observed `Theta_nuc` range from 23.4 to `2.96e5`. This is a coverage-limited model range, not a universal threshold.

### Post-onset attainment inequality

Once nucleation activates,

\[
\tau_{\mathrm{exchange}}+\tau_{\mathrm{transport}}
< t_{\mathrm{available}},
\]

must hold strongly enough for both compared paths to attain the scoring interval. Otherwise, a large ratio is censored rather than accepted.

### Closed-shrinkage/growth selectivity

Within the useful second-step band, closed-pore shrinkage must operate before appreciable grain growth,

\[
\tau_{\mathrm{closed\ shrinkage}} < t_{\mathrm{hold}}
< \tau_{\mathrm{grain\ growth}}.
\]

The corresponding observed selectivity group, `S_closed/growth`, ranges from 0.116 to 126 among the 73 exact joint successes. Its breadth demonstrates conditional combinations rather than a single scalar cutoff.

### Boundary conditions

A complete two-step window additionally requires both inequalities to reverse at its edges. Closed shrinkage/accommodation must become insufficient at low `T2`, or no lower density-exhaustion boundary exists. Grain growth must reactivate at high `T2`, or no upper boundary exists. PR/surface redistribution must also be strong enough to prepare the closed store without eliminating attainable densification.

## Why large high-temperature/two-step grain-size separation is plausible

The purpose of two-step processing is to separate densification from grain growth. A large matched-density grain-size difference is therefore an intended physical signature, not an artifact by magnitude alone. For candidate 693168, both the high-temperature reference and two-step path attain the density interval from 0.95 to 0.98 with supported interpolation. In the timestep-controlled audit, the high-temperature path reaches approximately 858–1123 nm while the representative two-step path remains near 118 nm. The result survives 30-, 15-, and 5-minute maximum-step audits, with only modest boundary movement.

This magnitude is plausible as a qualitative experimental-scale contrast because high-temperature migration can compound rapidly after density support becomes available, whereas an intermediate second step can retain a nanoscale structure. The model does not claim that the predicted absolute grain sizes or reduction percentage are quantitatively correct for a named material. The distinction is essential: large separation is admissible when both paths attain the interval and the states are numerically bounded, but quantitative agreement requires independent calibration.

## Why candidate 693168 remains conditional Tier B

Candidate 693168 is not promoted to Tier A for four reasons.

1. Its first-step growth is approximately 13.7%, above the strict Tier-A preparation limit.
2. Its closed fraction at the switch is approximately 0.649 of the remaining model pore inventory and approaches unity at high density. These are model-store fractions, not measured stereological closed porosity.
3. Closed shrinkage supplies approximately 99% of its modeled high-density support. The mapping from this internal channel to measured pore shrinkage, gas pressure, and accommodation capacity is not calibrated.
4. Its material robustness was evaluated exactly only for the frozen 693168 topology. The six-candidate family comparison transfers reduced margins; it is not six independent exact material-property searches.

The candidate is therefore valuable as the strongest conditional mechanism comparator and as a source of discriminating experiments. It is not evidence that its parameter values or internal state magnitudes apply universally.

## Experimental falsification and calibration priorities

The fast-firing hypothesis is most directly tested by measuring densification onset and grain-size evolution during interrupted ramps. At matched density, the model predicts smaller `G_mean`, `G50`, and `G90` for the fast path after the onset interval. A nucleation-facile material should show a substantially weaker heating-rate separation. Independent exchange and transport relaxation measurements are needed to confirm that post-onset completion, rather than continued waiting, permits density attainment.

The two-step hypothesis requires direct observation of the prepared pore state. Three-dimensional open/closed pore fraction, pore D50 and D90, the large-pore tail, and connected fine-pore/percolation measures should be collected at the first-step switch and during second-step holds. A trapped-gas pressure measurement or physically justified accommodation proxy should be followed across low-failure, successful, and upper-growth `T2` conditions. Grain-growth mobility should be measured across the upper boundary. Interrupted first- and second-step tomography should test whether the prepared state persists long enough to control the matched-density grain trajectory.

The most important calibration targets are the conversion between model stores and measurable 3D pore fractions, the closed-pore shrinkage rate, accommodation capacity and recovery, and high-temperature grain-growth mobility. Calibration should be performed against observables not used to select candidate 693168, followed by a blind schedule comparison.

## Limitations and non-claims

- This work does not validate the reduced model against a specific material.
- It does not provide a paper-ready quantitative calibration.
- Candidate 693168 remains conditional Tier B and is not a Tier-A result.
- Surrogate screen rows are not mechanistic evidence; exact-promoted trajectories control final classification.
- The observed activation-energy and prefactor windows are conditional on the tested model and parameter envelope, not universal material constants.
- Fast firing and two-step behavior were evaluated with the same perturbation vector in separate existing model layers; the synthesis does not introduce a shared dynamically coupled state.
- The near-unity closed-store trajectory in candidate 693168 is a primary falsification target, not a measured fact.
- No new simulation, search, mechanism, topology change, classification change, or parameter retuning was introduced in preparing this narrative.

Quantitative claims in this document are traceable to `results/final_mechanism_synthesis_and_property_windows/source_tables/` and to the timestep-controlled candidate-693168 audit.
