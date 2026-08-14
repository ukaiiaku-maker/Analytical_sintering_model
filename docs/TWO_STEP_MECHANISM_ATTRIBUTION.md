# Two-step mechanism attribution

## Result

Candidate 693168 remains a **conditional Tier-B mechanism candidate, not a validated or paper-ready model**. Its exact base reproduction reaches the high-density interval on both paths, produces a 0.906 grain-size reduction over a 0.03 density span, and retains a 250 °C finite Chen window with a lower densification-exhaustion boundary and upper grain-growth boundary.

The inferred chain is:

1. First-step PR/topology evolution builds a closed-pore store.
2. Closed shrinkage plus finite accommodation sustains density gain at low T2.
3. Grain-boundary migration remains suppressed in the success band.
4. At still lower T2, closed shrinkage/accommodation cannot reach the target.
5. At higher T2, thermally activated grain growth creates the upper boundary.

Large high-T/two-step grain-size separation is the desired trajectory signature and is not an artifact when density support, numerical stability, and state bounds hold.

## Exact local sensitivities

The PR prefactor is the largest exact rank correlate of Chen-window width (`Spearman r = 0.954`); `M_PR_closed` follows (`r = 0.912`) and `log10(k_PR/k_growth)` is third (`r = 0.842`). These correlations are association within promoted cases and do not prove sufficiency.

OAT checks show:

- `k_PR`: 0.03× and 0.1× fail; 0.3×–30× pass.
- `Q_PR`: the full ±100 kJ/mol tested range passes, although reduction and window width change.
- `Q_closed`: -50 kJ/mol and below lose the lower boundary; -25 through +100 kJ/mol pass.
- `k_closed`: 0.03×–30× pass but move both reduction and window width.
- `Q_growth`: the full ±100 kJ/mol tested range passes.
- `k_growth`: 0.03× loses the upper boundary; 0.1×–30× pass.

Thus the most identifiable local orderings are PR production relative to growth and the presence of both closed-shrinkage and growth boundaries. Absolute barrier widths are much less constrained by this campaign.

## PR and accommodation

PR is not required for fast firing, but it is required in the current 693168 two-step interpretation as the first-step preparation route. The destructive ablations are no PR damage, no closed transition, no closed shrinkage, and infinite closed accommodation. Infinite accommodation destroys the joint result because the finite store is part of what brackets the response; it is not a harmless numerical limit.

The dominant uncertainty is whether the modeled closed fraction and accommodation trajectory are physically realizable. At target density, the base candidate approaches an almost entirely closed pore inventory and derives roughly 99% of its high-density densification flux from closed shrinkage. Independent measurements of closed-pore volume, pressure/accommodation, and shrinkage kinetics are therefore essential.

## Boundary controls

The lower boundary is controlled by effective closed shrinkage, available accommodation, and the prepared closed inventory. The upper boundary is controlled by thermally activated migration/growth. A finite success band requires both; a success point without either bracket is rejected.
