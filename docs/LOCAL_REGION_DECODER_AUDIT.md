# Local-region decoder audit

The corrected decoder maps 47 Latin-hypercube coordinates into all six requested physical blocks: network structure, PR damage, sweep/coalescence, closed-pore support, migration/topology, and residual stress. A one-column-at-a-time dynamical-influence audit perturbs each coordinate and compares initialized network statistics, fluxes at 1000/1250/1450 °C, and short evolved states.

The 10,000-row preflight produced 10,000 unique decoded parameter vectors and 10,000 unique dynamic fingerprints across `N_regions = 8, 16, 32`. All 47 columns change at least one dynamical signature; none is unused or diagnostic-only. The search fingerprint contains every decoded physical parameter, not only `N_regions` or a reporting subset.

The local laws receive only temperature, the current state, graph adjacency, and material parameters. They contain no schedule, slow/fast, protocol, or target labels.
