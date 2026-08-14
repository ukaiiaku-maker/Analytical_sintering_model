# Topology-constrained mechanism architecture

## Scope

`topology_constrained_sintering.py` is a parallel prototype, not a silent
rewrite of v3/v6. It makes the model's hypotheses searchable through named
parameters and switchable mechanisms.

The update loop is:

1. infer boundary coverage, connectivity, isolation, P-R eligibility, and
   triple-line activity from `rho`, `G`, and the pore bins;
2. calculate grain-scale baseline stress and topology-dependent concentration
   stress;
3. calculate separate nucleation, exchange, transport, and triple-line times;
4. evaluate mechanism flux objects with density, growth, pore-bin, and power
   contributions;
5. normalize nonnegative power propensities under topology compatibility;
6. evolve the state while projecting pore volume exactly onto
   `rho = 1 - sum(phi)`.

## Physical distinctions

Renewal activity is reported as `Lambda/(1+Lambda)`, but it does not determine
densification alone. `completion_rate`, `density_gain_per_event`,
`grain_growth_per_event`, topology yield, and `E_G` are separate. Exchange and
transport power diagnostics remain separate from their kinetic times: a large
power channel is not automatically labeled rate limiting.

The current partition is a positive-propensity normalization, the simplest
constrained prototype. Its weights are nonnegative and sum to one. It is not an
`eta_total`, because every channel is named, recorded, topology-gated, and can
be ablated. A later version should replace normalization with a convex
free-energy minimization and explicit Onsager matrix.

## Known limitations

- Topology is inferred from moments rather than a generated network.
- Bin radii are fixed; coarsening transfers volume upward rather than moving
  bin coordinates.
- Default parameters are dimensional plausibility estimates, not calibration.
- The present formulation gives a two-step advantage but the baseline
  heating-rate sign is wrong; this is retained as a scientific failure rather
  than hidden by a target change.
