# Mechanism Options and Coupling Strategies

This document describes candidate approaches for building a unified analytical sintering model that can reproduce both fast-heating-rate sintering and two-step sintering.

The options should be treated as complementary rather than mutually exclusive. The recommended path is to combine them into a staged reduced-order framework:

1. constrained dissipation partition;
2. effective topology / microstructure state;
3. modular multi-mechanism kinetics;
4. automated search after the physical state and constraints are explicit.

## 1. Multi-mechanism renewal model

The core densification model should remain a renewal process. Densification requires:

1. nucleation or availability of climb-mediating defects;
2. local exchange / point-defect formation and annihilation;
3. point-defect transport;
4. possible mobility or drag of the defect, triple line, or boundary segment;
5. a pore-connected topology so the strain actually removes pore volume.

A general event time can be represented as:

```text
tau_event = tau_nuc + tau_exchange + tau_transport + tau_drag
```

with densification strain rate:

```text
edot = Delta_epsilon_event / tau_event
```

The existing two-step model used a simpler two-process renewal form:

```text
Lambda = r_nuc * tau_sink
activity = Lambda / (1 + Lambda)
```

That remains useful, but it should not be the only descriptor once exchange resistance, triple-line drag, topology, and pore distribution are included.

### Goal for Codex

Generalize the renewal model so the active densification rate can include:

```text
tau_nuc
tau_exchange
tau_transport
tau_triple_line_drag
```

and report diagnostics for each limiting time.

### Acceptance test

The model should identify whether a given condition is:

```text
nucleation-limited
exchange-limited
transport-limited
triple-line-drag-limited
topology-limited
```

without relying on a single scalar efficiency factor.

## 2. Nucleation-rate-limited triple-line drag

A plausible missing mechanism is triple-line or pore-edge drag that becomes important when boundaries move while pores and triple lines must reorganize.

This can be represented as an additional serial time:

```text
tau_event = tau_nuc + tau_sink + tau_TL
```

or as a competing dissipative channel:

```text
P_TL_drag = force_TL * velocity_TL * active_TL_density
```

Triple-line drag should depend on:

```text
triple-line density
pore size distribution
pore connectivity
grain boundary velocity
temperature
local stress
density / topology stage
```

### Why it may solve the problem

Fast heating can help because the compact spends less time in intermediate-temperature regimes where coarsening and topology degradation occur before densification activates.

Two-step sintering can help because the first high-temperature step activates densification, but the second lower-temperature step reduces grain-boundary migration and triple-line/pore drag.

### Implementation goal

Add a `triple_line_drag` mechanism with:

```python
MechanismFlux(
    rho_dot_contribution,
    G_dot_contribution,
    pore_bin_flux,
    power_dissipated,
    diagnostics
)
```

At first, make it phenomenological but topology-dependent.

## 3. Constrained dissipation / Lagrange-multiplier partition

A Lagrange-multiplier-style formulation is useful if it enforces physical constraints, not if it becomes another fitting parameter.

The state evolves under a finite interfacial power budget:

```text
P_available = -dF_interface/dt + P_external
```

This power is partitioned among:

```text
P_densification
P_clean_GB_migration
P_pore_connected_GB_migration
P_pore_drag
P_triple_line_drag
P_PR_desintering
P_exchange
P_transport
```

A constrained solver should enforce:

```text
P_available = sum(P_i)
pore volume conservation
solid volume conservation
nonnegative pore-bin populations
rho = 1 - total_pore_volume
0 <= f_pore, f_clean, f_PR, f_TL <= 1
```

The multiplier can be interpreted as a shadow price for available interfacial power or for topological compatibility.

### Codex goal

Implement a small constrained partition function:

```python
partition = solve_dissipation_partition(state, mechanisms, constraints)
```

The first version can be simple: allocate fractions using positive mechanism propensities and renormalize under constraints. A later version can use an optimizer.

### Important distinction

The largest dissipation term need not be the rate-limiting step. Keep kinetic times and dissipated powers separate.

## 4. Effective topology loop

This is probably the most important option.

Instead of running a spatial microstructure simulation, construct an effective topology from instantaneous reduced descriptors:

```text
rho
G
pore radii r_i
pore volume fractions phi_i
pore number densities N_i
```

Then compute topology descriptors:

```text
f_pore        = pore-connected boundary fraction
f_clean       = clean boundary fraction
f_PR          = P-R / de-sinter eligible fraction
f_TL          = triple-line-drag active fraction
connectivity  = pore network connectivity metric
isolation     = isolated pore fraction
```

Then run one analytical model step and update the topology.

### Proposed loop

```python
state = initial_state()

for step in time:
    topology = infer_effective_topology(state)
    mechanisms = evaluate_mechanisms(state, topology, T)
    partition = solve_dissipation_partition(state, mechanisms)
    state = evolve_state(state, partition, dt)
```

### Why this is better than raw parameter search

It gives the model a way to decide when boundaries remove pore volume versus when they only migrate, drag pores, or dissipate energy.

### Initial topology formulas

Start with simple monotonic functions:

```text
f_pore ~ 1 - exp(-chi * projected_pore_area / GB_area)
f_clean = 1 - f_pore
f_PR ~ low_density_factor * small_ligament_factor * pore_connectivity_factor
f_TL ~ pore_boundary_fraction * triple_line_density(G, pore bins)
```

Then refine from diagnostics.

## 5. Multi-bin pore distribution

A two-bin pore model is too weak for Plateau-Rayleigh and large-pore memory effects. Use logarithmic pore-size bins.

State variables:

```text
r_i      pore radius bins
phi_i    pore volume fraction in bin i
N_i      pore number density in bin i
```

At fixed pore volume, larger pores imply fewer pores:

```text
N_i = phi_i / ((4/3) pi r_i^3)
```

Boundary coverage should use projected area:

```text
coverage ~ sum_i N_i * pi * r_i^2
```

Densification should remove small pores faster:

```text
removal_weight_i ~ phi_i * (r_ref / r_i)^q
```

P-R/de-sintering should move volume from small bins to larger bins:

```text
phi_i -> phi_{i+1}
```

Pore drag and coalescence should also push volume and/or number toward larger bins.

### Acceptance test

At equal pore volume, a distribution with larger pores should have:

```text
lower pore-boundary coverage
lower densification efficiency
higher large-pore memory
stronger drag/coalescence penalty when mobile
```

## 6. Baseline stress plus concentration stress

The activation stress should not be a single fitted stress.

Use:

```text
sigma_local = K_microstructure * (sigma_base + sigma_concentration + sigma_external)
```

where:

```text
sigma_base ~ 2 gamma / R_characteristic
```

and `sigma_concentration` comes from geometric evolution, coarsening, triple-line constraint, or local topology.

Important: choose the length scale carefully.

Candidate baseline stresses:

```text
pore-scale capillary stress: 2 gamma_s / r_pore
grain-scale sintering stress: 2 gamma_s / R_grain
neck/triple-line local stress: K_TL * gamma / r_neck
```

Do not assume one scale is correct globally. Codex should expose them as separate diagnostics and test which is needed.

## 7. Automated search

After the architecture is explicit, search over physical parameters:

```text
Q_nuc
Q_exchange
Q_transport
Q_TL
Vact
Dgb
GB mobility
surface mobility
pore drag coefficient
triple-line drag coefficient
P-R coefficient
connectivity thresholds
pore-bin transfer exponents
stress concentration factors
```

Use percent objectives:

```text
HR_pct = 100*(G_slow - G_fast)/G_slow
TS_pct = 100*(G_highT - G_twoStep)/G_highT
```

Add trajectory efficiency:

```text
E_G = d rho / d ln G
```

and penalize:

```text
failure to reach target density
runaway grain growth
unphysical stresses
negative pore bins
hidden nonconservation
success only on training schedule
```

## Recommended next implementation order

1. Refactor state and topology.
2. Add mechanism modules.
3. Add constrained partition.
4. Add multi-bin pore tracking if not already clean.
5. Add triple-line drag mechanism.
6. Add combined percent/objective score.
7. Run small deterministic diagnostics.
8. Run limited search.
9. Run held-out validation.
10. Only then run large automated search.
