# Non-claims and guardrails

## Purpose

This document defines the language and evidence limits for presenting the fast-firing and two-step mechanism synthesis. It should accompany manuscript drafting, captions, talks, and data releases so that conditional model findings are not converted into validation claims.

## Evidence hierarchy

1. **Exact-promoted trajectories control final classification.** The final exact union contains 1,903 cases: 485 fast-only, 119 two-step-only, 73 both-pass, and 1,226 neither.
2. **The surrogate screen is not evidence.** It screened 50,655 perturbations and predicted 19,880 both-pass rows, but only 73 both-pass cases appear in the exact promoted union. The 272.3-fold discrepancy must be disclosed whenever the screen is shown.
3. **OAT windows are local model observations.** `Delta Q_nuc = 0…+50 kJ/mol`, `Delta Q_closed = -25…+100 kJ/mol`, the tested `Delta Q_growth = ±100 kJ/mol`, the `0.3x` PR threshold, and the `0.1x` growth threshold apply to the frozen tested model and envelope.
4. **Candidate-family results are exact at their base states but not equally exact in material-property robustness.** The six-candidate material comparison uses a reduced transfer from exact candidate-693168 OAT.

## Required non-claims

### No validation claim

The present result is a conditional mechanism attribution. It has not been validated against an independently measured material, pore-state trajectory, and held-out processing schedule.

Do not write:

- “The model is validated.”
- “The mechanism is proven.”
- “The parameters are material constants.”

Preferred language:

- “The exact-promoted model supports a conditional mechanism hypothesis.”
- “The ablation and property-window audits are consistent with…”
- “The predicted state trajectory provides a falsification target.”

### No paper-ready calibration claim

The model has not been quantitatively calibrated for closed-pore fraction, gas accommodation, PR preparation, or absolute grain-growth mobility. Publication-style figures improve communication; they do not make the parameterization paper ready.

Preferred language:

- “Manuscript-ready mechanism narrative” refers to writing structure, not quantitative validation.
- “Publication-style figure” refers to graphical quality, not evidentiary status.

### Candidate 693168 is conditional Tier B

Candidate 693168 is not Tier A. Its first-step growth is approximately 13.7%, its closed-store trajectory is large and uncalibrated, and its material robustness has been evaluated exactly only at its frozen topology. It is the strongest current separation comparator and an informative experimental target.

Do not promote it because the high-temperature/two-step separation is large. Large separation is allowed when attained and stable, but it does not cure missing calibration.

### Large grain-size separation is not an artifact by magnitude alone

The candidate's high-temperature path reaches approximately 858–1123 nm while the representative two-step path remains near 118 nm over `rho = 0.95–0.98`. Both paths attain the interval and the result survives tighter timesteps. This is the intended qualitative two-step signature.

Do not describe the result as invalid solely because the ratio is large. Instead state that its magnitude is experimentally plausible but quantitatively uncalibrated.

### Surrogate screening is not final evidence

Surrogate/dimensionless scores are allowed for promotion, phase-space visualization, and prioritization. They must be labeled as screening evidence and cannot determine final both-pass counts, causal claims, or material-property windows.

### Closed-pore/accommodation state is the main calibration target

Candidate 693168 derives nearly all modeled high-density support from closed shrinkage and approaches a nearly fully closed remaining-pore inventory at high density. These are internal store fractions, not directly measured stereological pore fractions. Closed-pore fraction, shrinkage kinetics, trapped gas, accommodation capacity, and recovery are the primary experimental constraints.

Do not present the modeled closed fraction as a measured porosity.

### Material-property windows are not universal constants

The reported activation-energy and prefactor intervals are observed survival regions in a finite tested model envelope. They are not independent rectangular tolerances and should not be transferred directly to arbitrary materials.

Preferred language:

- “Observed exact OAT survival range.”
- “Coverage-limited joint envelope.”
- “Relative ordering consistent with the tested model.”

Avoid:

- “Universal critical activation energy.”
- “Necessary bound for all crystalline powders.”
- “Material-independent prefactor threshold.”

## Mechanism-specific guardrails

### Fast firing

- Attribute the current effect primarily to nucleation-limited onset.
- State that nucleation-facile ablation removes or weakens the effect.
- State that PR-off can preserve the effect and is not a rejection.
- Retain matched-density attainment and the continuous density-span criterion.
- Do not claim that PR can never influence fast firing in another material or closure.

### Two-step sintering

- Attribute candidate 693168 primarily to PR-prepared closed-pore accommodation memory.
- Retain the lower density-exhaustion and upper grain-growth boundaries.
- Report the fixed target, time budget, and growth tolerance with every Chen map.
- State that no PR damage, no closed transition, no closed shrinkage, and infinite accommodation destroy the candidate's joint result.
- Do not claim that TJ multihit, pore drag, residual stress, or sweep/coalescence are universally irrelevant; their roles differ across the Tier-B family.

### Joint behavior

- State that the same relative-property perturbation vector passed two separate exact model layers.
- Do not imply that the synthesis introduced a fully coupled fast-firing/two-step dynamic state.
- Do not infer exact labels for the unpromoted portion of the 50k screen.

## Numerical and data guardrails

- Do not score unattained density intervals.
- Do not lower target density to rescue a case.
- Do not extend time budgets selectively.
- Do not hide negative pore stores, unbounded states, unsupported interpolation, or timestep collapse.
- Do not replace named mechanism channels with an unexplained efficiency scalar.
- Preserve exact pore-volume conservation and nonnegative state stores.
- Retain rejection reasons and distinguish censored, unattained, and mechanistically failed cases.
- Keep archived-result deletions separate from science and writing commits.

## Figure and caption guardrails

Every figure should identify whether its evidence is exact, surrogate, reduced transfer, or schematic. Captions should state the candidate tier and unresolved calibration when candidate 693168 is shown. Screening maps should include a visible “not final evidence” label. Candidate-family figures should state that material robustness is a reduced transfer. Large separation should be described as an intended qualitative signature with uncalibrated magnitude.

## Minimum disclosure paragraph

The following statement, or its equivalent, should appear in any manuscript based on the current branch:

> The reported framework provides a conditional mechanism attribution rather than validation. Final classifications use exact-promoted trajectories; the dimensionless screen is used only for promotion. Candidate 693168 remains Tier B because its closed-pore/accommodation trajectory and absolute magnitude are not independently calibrated. Reported activation-energy and prefactor windows are observed within the tested model envelope and are not universal material constants.

## Conditions required before a validation claim

A validation claim requires all of the following:

1. independent mapping between model pore stores and 3D open/closed pore measurements;
2. calibrated nucleation onset, exchange, and transport times;
3. calibrated closed-pore shrinkage and accommodation capacity/recovery;
4. calibrated grain-growth mobility across the upper T2 boundary;
5. reproduction of matched-density grain and pore distributions, not density alone;
6. successful prediction of held-out heating rates and two-step schedules without retuning;
7. uncertainty propagation and replicate experimental support.

Until those conditions are met, use “conditional,” “consistent with,” “mechanism hypothesis,” and “falsification target,” not “validated,” “proven,” or “universal.”
