# Equation limitations and approximations

## Implemented proxies

1. The fast-model connectivity factor is a clipped linear density gate. It is not a derived pore-percolation law.
2. The local-region activity is a sigmoid temperature gate multiplied by residual-stress inhibition. It is not the serial nucleation/exchange/transport law.
3. The closed-pore radius factor uses \(\phi_{\rm closed}/N_{\rm closed}\) in code-normalized units. It is not a calibrated pore diameter.
4. Closed accommodation is a bounded capacity/recovery state consumed by closed shrinkage. It is not a derived defect-chemistry or gas-diffusion law.
5. The local-region TJ completion factor is a sigmoid proxy. Exact Poisson completion is implemented only in the earlier agentic TJ-multihit model.
6. Stress, action power, and PR work channels are reduced proxies unless a source report explicitly supplies physical units.

## Numerical regularization

Exponential arguments are clipped: generally \([-50,50]\) or \([-60,60]\) in final relative-temperature expressions, with wider \([-700,700]\) protection in earlier aggregate Arrhenius helpers. Denominators use small positive floors. Probabilities and state fractions are clipped to physical bounds. Local pore removal is capped by available volume. Grain growth in the local model uses an exact inverse-size step.

These protections prevent overflow and negative stores; they are not additional physical mechanisms.

## Evidence limitations

The 50,655-row screen is a promotion device. Its algebraic fast and two-step scores are screening-only. It predicted 19,880 both-pass rows, while the 1,903-row exact union contains 73 both-pass rows. Exact-promoted trajectories control every final classification.

The six-candidate family is exact at the frozen base states, but its material-property comparison is a reduced transfer from exact candidate-693168 OAT. It is not six independent exact property campaigns.

## Candidate-693168 limitations

Candidate 693168 is conditional Tier B:

- first-step growth is approximately 13.71%, above Tier A;
- 64.94% of remaining pore volume is in the modeled closed store at the switch;
- the successful path approaches a nearly closed remaining inventory at target;
- high-density densification is dominated by the closed-shrinkage channel;
- accommodation capacity, recovery, gas ratio, and store-to-tomography mapping are uncalibrated.

The large high-temperature/two-step grain-size separation is supported by attained trajectories and timestep checks, but its magnitude is not fitted. It is a falsifiable prediction, not validation.

## Negative controls and non-final closures

Aggregate dissipation partitioning, density-gated smoothing memory, pore placement, local action allocation, persistent-junction drag, TJ multihit, residual stress, connected-sink mixtures, and earlier late-stage closed-pore models remain in the audit because they shaped the mechanism progression. They must not be presented as the final joint mechanism.

Persistent junction plus exact Poisson TJ multihit produced earlier Chen windows. Candidate 693168 does not require those channels in destructive ablation. Conversely, the candidate's closed-store equations must not be retroactively assigned to those earlier windows.

## Classification approximations

Chen boundaries depend on a finite temperature grid, adaptive extensions, local 10 °C refinement, a fixed time budget, a specified target density, and a specified growth tolerance. A finite window means bracketing within that protocol definition, not a universal material phase boundary.

Fast firing depends on supported matched-density interpolation and the fixed ratio/span rule. Raw endpoint differences, unattained ratios, and censored trajectories are excluded.

## Non-claims

- No equation audit constitutes validation.
- No reported parameter window is a universal constant.
- No surrogate row is final evidence.
- No hidden closed-pore \(\Lambda/K\) equation is claimed.
- No model physics, parameter, time budget, or classification was changed for this documentation branch.
