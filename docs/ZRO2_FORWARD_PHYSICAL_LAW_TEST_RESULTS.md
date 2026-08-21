# Physical closed-law analytical test results

## Guardrails

All fixed-state limiting tests passed: density rate vanishes with zero closed
inventory or nonpositive driving stress; increasing gas pressure reduces driving
stress; the tested radius laws decrease with radius; surface accommodation alone
does not densify; and conservative transfers preserve pore volume. The m=3 and m=4
GB forms pass dimensional checks only after assigning a coefficient with units
m\(^{m-2}\); its magnitude is not yet physically derived.

The scan contains 101,376 rows over 850–1350 °C, six radii, four inventories, four
accommodation levels, four gas ratios, and 20/40/96 h. All rows flag barrier
temperature extrapolation. Some bracketed states produce finite density gain inside
the reduced diagnostic band, while surface accommodation produces zero density gain
as required.

## Boundary preservation

The fixed selected first-step state produced only
`DENSIFICATION_EXHAUSTION_FAILURE` at the lowest temperatures and `MIXED_FAILURE`
once intrinsic growth exceeded tolerance. None of the renewal, GB diffusion,
gas-limited, surface-only, or empirical comparator modes produced a success region.
Thus no law preserves a lower/success/upper topology, zero laws are promoted, and
the optional bounded process map and Chen-map figure were not run.

## Energy ledger

Four representative exported paths were reconstructed. Budget-violation fractions
range from about 0.32 to 0.83. Near-zero reconstructed available power creates very
large diagnostic ratios, so the raw powers and residual must be examined alongside
`Pi_dens` and `Pi_total`. The result demonstrates that the prior GB-area-only equality
does not close a general energy budget. Because pore bins were not exported, the
historical ledger is an audit-quality reconstruction, not a precise free-energy
calculation.

No accepted forward-model physics, barrier file, diffusivity, or mobility parameter
was changed. This analytical audit is **not validation**. No validation claim is made.
