# Candidate 693168 closed-accommodation audit

## Status

Candidate 693168 reproduces as a **conditional Tier B prototype**. It is not Tier A, calibrated, validated, or paper-ready. The focused audit changes no model physics, material kinetics, topology parameters, density targets, time budgets, candidate IDs, or classification rules.

The exact parameter vector was reconstructed from the recorded Latin-hypercube seed and candidate row, not from rounded CSV text. Its fingerprint is `b46304a8d7c7055a`, matching the production search.

## Numerical reproduction

The exact 30-minute run reproduces the prepared state: switch density 0.880000, prepared mean grain size 117.067 nm, 13.71% first-step growth, and 64.94% of the remaining pore volume in the closed store. The independently timestep-controlled audit gives:

| Maximum step | Success band | Width | Median reduction, 0.95–0.98 |
|---|---:|---:|---:|
| 30 min | 930–1200 °C | 270 °C | 89.4% |
| 15 min | 930–1190 °C | 260 °C | 86.5% |
| 5 min | 930–1190 °C | 260 °C | 88.3% |

The switch density and (G_1) are invariant at reported precision. Boundaries shift by 10 °C and median reduction by at most 2.5 percentage points, within the requested tolerances. The result therefore does not collapse under tighter integration.

The original production helper enforces a 30-minute minimum for second-step calls (`max(dt, 1800)`). The audit did not change it; instead, it uses a separate runner over the same local evolution laws so 15- and 5-minute limits are genuine.

## Absolute trajectory magnitude

The magnitude is extreme. Over matched density 0.95–0.98, the 5-minute reconstruction has:

- high-temperature reference grain size 858–1123 nm;
- representative 1100 °C two-step grain size 117.7–118.4 nm;
- minimum/median/maximum reduction 85.0%/86.5%/87.3% for the 15-minute reconstruction, 87.3%/88.3%/88.8% at 5 minutes, and 88.3%/89.4%/90.4% for the production-compatible 30-minute reconstruction.

Thus the large ratio is not interpolation censoring: both paths attain the interval and all support flags are true. It is caused by very large absolute high-temperature reference growth combined with nearly arrested second-step growth. This is a prototype-scale prediction and cannot be treated as quantitatively plausible without grain-growth and pore-topology calibration.

## Closed-pore trajectory

At the switch, the remaining-pore composition is 64.94% closed and the available accommodation is 0.152 on a capacity of 0.909. The successful path reaches the target with approximately 94.5% of the remaining pore volume closed. The high-temperature reference approaches a fully closed remaining-pore population. These fractions are model store fractions, not measured stereological closed porosity.

Closed-pore shrinkage supplies most of the density gain after the switch and is essential for reaching 0.98. This makes the predicted closed population and accommodation history the primary falsification target. The trajectory is internally conservative and bounded, but its magnitude is physically aggressive.

Finite accommodation creates the lower boundary. At low (T_2), closed shrinkage becomes too slow relative to the fixed time budget and the path stalls below 0.98. The finite state remains nonnegative and no greater than capacity. Infinite accommodation does not make low temperatures universally successful in this candidate; instead, it changes the prepared topology and destroys the joint trajectory/window result. Therefore the relevant claim is that finite accommodation is required for this candidate—not the stronger claim that its removal erases only the lower boundary.

The upper boundary is thermal grain-growth activation. Above roughly 1200 °C the target is reached, but second-step growth exceeds the fixed 20% production tolerance.

## Ablations

Exact reruns confirm that removing PR damage, closed transition, closed shrinkage, or finite accommodation destroys the joint high-density trajectory/window result. No sweep/coalescence, no TJ multihit, no attached-pore drag, no persistent-junction state, no network heterogeneity, and the registered topology-disabled migration closure retain a complete window. Residual-stress removal retains a window but reduces the median trajectory advantage substantially.

This candidate is therefore best interpreted as **PR-prepared closed-store accommodation memory**, not as a universal TJ-multihit or heterogeneous-network result. The fact that `topology_disabled` preserves the effect means the decoder's closed-pore channel is causally upstream of the migration-topology ablation currently named “topology disabled”; it does not mean all topology is irrelevant.

## Robustness

The bounded (5\times6) initial-condition audit retains a complete finite window throughout ρ0=0.65–0.80 and G0=50–300 nm, without retuning. The ρ0=0.60 row is not accepted: preparation growth is extreme, the closed fraction approaches one, and a complete bracketed window is lost. Robustness is therefore broad within the original neighborhood but not universal.

## Fast firing

Frozen E0021 and E0142 material records remain meaningful with maximum ratios 1.860 and 1.796 over density span 0.17. Nucleation-facile ablation removes the effect and PR-off preserves it. The local-region audit does not alter `MaterialKinetics` or its density rate; fast firing was not reoptimized.

## Answers to audit questions

1. **Timestep reproduction:** yes; the candidate survives 30/15/5-minute limits.
2. **Absolute grain sizes:** high-T 858–1123 nm versus two-step 118 nm over 0.95–0.98.
3. **Physical plausibility:** numerically supported but physically extreme.
4. **Closed fraction at switch:** 64.94% of remaining pore volume.
5. **High-density support:** predominantly the closed-shrinkage channel.
6. **Finite accommodation:** required for the joint candidate and controls low-temperature exhaustion together with thermal shrinkage.
7. **Lower boundary:** insufficient closed-pore shrinkage/accommodation within the common time budget.
8. **Upper boundary:** thermally activated grain growth.
9. **Destructive ablations:** no PR damage, no closed transition, no closed shrinkage, infinite accommodation.
10. **Candidate-family context:** 693168 is the strongest and an extreme member, not the only Tier-B result.
11. **Best falsification experiment:** interrupted first-step and second-step measurements of open/closed pore fraction, pore-number distribution, internal gas pressure/accommodation proxy, and absolute grain size at matched density.
12. **Before validation:** calibrate the mapping from model stores to 3D connected/closed porosity, closed-pore shrinkage kinetics, gas accommodation capacity/recovery, and high-temperature grain-growth kinetics.

## Missing observables

The model does not expose physical pore diameters, explicit nucleation/exchange/transport times, absolute PR/densifying energy, or a separately derived closed-pore Poisson Λ/K event model. The final tables retain explicit `NaN` fields or clearly named proxies rather than inventing these quantities. See `final_missing_plot_inputs.md`.
