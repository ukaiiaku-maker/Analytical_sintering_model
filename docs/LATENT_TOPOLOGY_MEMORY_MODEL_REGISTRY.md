# Latent Topology Memory Model Registry

Implemented dynamic states are:

- connected-removable inventory: C_rem_GBseg, C_rem_TJ, large-tail/isolation and percolation/removable memory;
- junction/segment memory: XJ, junction density, segment length, constraint/pore/pinned/relaxed partitions;
- residual-stress/work memory: location stresses, PR/shear work, capped stress memory and relaxation state.

All latent states are bounded. The combined migration factor uses multihit completion, junction/pore drag, and residual-stress resistance. It never changes the frozen material `rho_dot`.
