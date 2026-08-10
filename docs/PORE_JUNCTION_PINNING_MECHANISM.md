# Observable pore/junction pinning audit

## Outcome

This fixed-candidate prototype does **not** satisfy the nanoscale-window
acceptance criterion. It preserves the lower densification-exhaustion boundary,
the upper grain-growth boundary, exact baseline recovery, and exact pore
conservation, but it produces no finite 5--10% success window at or below
300 nm. This is a useful negative result, not validation and not a reason to
increase a resistance constant.

## State-resolved closure

The new `growth_mode="pore_junction_pinning"` is an ablation alongside the
unchanged `baseline` and `junction_limited` modes. It uses

```text
M_eff/M_GB = 1 / (1 + tau_pin/tau_GB)

tau_pin/tau_GB = R_pin (G_ref/G)^n
  (connected_coverage / coverage_ref)
  (pore_junction_occupancy / junction_ref)
  (1 - isolated_fraction)
  [fine_floor + (1-fine_floor) removable_fine_fraction]
  exp[(Q_pin-Q_growth)/R (1/T-1/T_ref)]

pore_junction_occupancy = 1 - exp[-L_pore G^2]
L_pore = sum_i N_i 2 pi r_i
```

Every quantity is instantaneous and observable or derived from the pore-bin
state. There are no protocol, ramp-rate, target, or schedule-class inputs.
The single fixed hypothesis uses `R_pin=24`, `G_ref=150 nm`, `n=1`,
`T_ref=1250 C`, `Q_pin=500 kJ/mol`, coverage reference 0.35, junction
occupancy reference 0.50, and fine-pore floor 0.20. These values were not swept
after observing the map.

Pinning acts only on GB migration. For this new mode the dissipation partition
retains the free-migration pore-drag propensity, so a shared state has exactly
the same instantaneous densification rate as baseline while its grain-growth
rate is reduced. The existing junction-limited implementation is preserved.

## Bounded Chen-style map

The calculation contains 4,590 second-step trajectories and 9,180
classifications at fixed `rho_target=0.90` and fixed 96 h budgets:

- `G0 = 25, 35, 50, 75, 100, 150, 225, 300, 450, 600 nm`
- `T1 = 1250, 1300, 1350 C`
- `rho_switch = 0.75, 0.80, 0.85`
- `T2 = 900...1300 C` in 25 C increments
- growth tolerances 5% and 10%

| mode | tolerance | sizes with any success | sizes with finite >=25 C window |
|---|---:|---|---|
| baseline | 5% | 600 nm | none |
| baseline | 10% | 450, 600 nm | 600 nm |
| junction limited | 5% | 450, 600 nm | none |
| junction limited | 10% | 300, 450, 600 nm | 450, 600 nm |
| pore/junction pinning | 5% | 450, 600 nm | none |
| pore/junction pinning | 10% | 450, 600 nm | 450, 600 nm |

For pinning, the minimum growth required among target-reaching paths is 15.79%
at 300 nm and 71.65% at 225 nm. Thus even isolated nanoscale successes are
absent. The mode still contains, across both tolerances, 2,048
densification-exhaustion failures, 350 grain-growth failures, 620 mixed
failures, and 42 successes. It does not turn the map into universal success.

Baseline rows reproduce the committed fixed-model control within `7.8e-16` in
density and `4.6e-13 nm` in grain size.

## Why pinning fails

The state resolution identifies a release instability hidden by the mean-grain
junction closure. Connected pore-boundary coverage decreases as densification
proceeds. The pinning time ratio therefore collapses and the migration mobility
rises late in the second step. At 300 nm, a representative `T1=1300 C`,
`rho_switch=0.85` path reaches the density boundary near `T2=1225 C`, but by
then accumulated growth is about 16.4%. Lower T2 retains pinning but exhausts
densification; higher T2 releases migration and fails the growth criterion.

The pore-junction occupancy proxy is also nearly saturated for much of the
map, so connected coverage carries most of the state dependence. This makes
the closure observable but does not supply persistent nanoscale pinning after
pores detach or become isolated.

## Implication

Migration suppression is necessary but this connected-pore form is not
sufficient. Increasing `R_pin` would obscure the diagnosed release problem.
The next mechanism should distinguish persistent solute/second-phase or
structural junction drag from removable connected-pore pinning, or resolve a
junction population whose relaxation is not identical to loss of connected
pore coverage. Independent interrupted-hold measurements of boundary mobility,
connected coverage, and grain growth would discriminate those alternatives.

## Reproduction

```bash
python3 pore_junction_pinning_sensitivity.py --workers 4
python3 -m pytest -q
```

CSV trajectories, classifications, and boundaries plus 5%/10% Chen-style
maps and the pinning-release diagnostic are in
`results/pore_junction_pinning/`.
