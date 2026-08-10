# Nanoscale grain-growth suppression mechanism audit

## Status

This is a mechanism prototype and negative/partial result, not a validation or
calibration. The fixed-model phase-space map remains the negative-control
baseline. No densification kinetics, pore redistribution law, density target,
or per-step time budget was changed in this audit.

The junction-limited closure shifts the 10% growth-window onset from 450 nm to
300 nm, but it does **not** create a 5% window below 450 nm. It therefore moves
the model in the required direction without yet reproducing a robust 50--300
nm Chen/Wang-style window.

## Closures

`growth_mode="baseline"` is an exact identity: the mobility factor is one and
the trajectories reproduce the committed expanded-phase-space control to
within `7.8e-16` in density and `4.6e-13 nm` in grain size.

The junction-limited mode uses a serial migration time:

```text
M_eff/M_GB = 1 / (1 + tau_junction/tau_GB)

tau_junction/tau_GB = R_ref (G_ref/G)^n
    exp[(Q_junction - Q_growth)/R (1/T - 1/T_ref)]
```

with one global parameter set for every size and schedule:
`R_ref=50`, `G_ref=150 nm`, `n=2`, `T_ref=1250 C`, and
`Q_junction=520 kJ/mol`. The `G^-2` dependence represents increasing junction
density at nanoscale grain size. The higher junction activation energy keeps
GB diffusion/densification active while migration is more strongly suppressed
at low T. These values are hypotheses requiring independent mobility data;
they were not fitted by initial condition or target.

The threshold alternative compares capillary drive `4 gamma_s/G` with a
thermally activated, size-dependent junction resistance and applies a smooth
activation sigmoid. At the tested default it is effectively a negative
control: its window onset is unchanged from baseline.

The migration power offered to the nonnegative dissipation partition is based
on the unsuppressed migration propensity. Junction resistance dissipates the
unrealized part as a named `junction_drag_power` channel. This is important:
using the already-suppressed growth rate as the power propensity inadvertently
changed the densification partition, violating the requested separation of
growth from densification.

All closures use only instantaneous temperature, grain size, topology, and
material parameters. They contain no protocol, ramp-rate, schedule, target, or
time-budget labels.

## Bounded map

The deterministic map contains 4,590 second-step trajectories and 9,180
classifications:

- `G0 = 25, 35, 50, 75, 100, 150, 225, 300, 450, 600 nm`
- `T1 = 1250, 1300, 1350 C`
- `rho_switch = 0.75, 0.80, 0.85`
- `T2 = 900...1300 C` in 25 C increments
- `rho_target = 0.90`
- growth tolerances 5% and 10%
- fixed 96 h first- and second-step budgets

The full records are in
`results/growth_mechanism_sensitivity/growth_mode_trajectories.csv`, with
classifications and boundaries in the adjacent CSV files.

| mode | tolerance | grain sizes with success | widest sampled window |
|---|---:|---|---:|
| baseline | 5% | 600 nm | 0 C |
| baseline | 10% | 450, 600 nm | 25 C |
| junction limited | 5% | 450, 600 nm | 0 C |
| junction limited | 10% | 300, 450, 600 nm | 50 C |
| threshold mobility | 5% | 600 nm | 0 C |
| threshold mobility | 10% | 450, 600 nm | 25 C |

At `G0=300 nm`, junction limitation gives one sampled successful temperature,
`T2=1225 C`, for each `T1` at `rho_switch=0.85`; the required growth is about
8.4%. This is a zero-width point on the 25 C grid, not yet a robust finite
window. No success appears at 225 nm or below.

The junction mode is not universal success. Across both tolerances it retains
2,133 densification-exhaustion failures, 225 grain-growth failures, 407 mixed
failures, and 238 unattainable first steps. Thus both the lower kinetic boundary
and upper migration boundary survive. The 119 unattainable first-step states
per tolerance are recorded and never scored as second-step successes.

## Interpretation

The serial junction time reduces the erroneous nanoscale fractional growth and
separates diffusion from migration in the intended direction. Its limited
success shows that migration mobility is part of the missing physics. However,
the strict window still begins too coarsely, and the apparent 300 nm result is
only a single T2 point at the 10% criterion. Increasing the resistance further
without independent evidence would be parameter forcing, so no broader search
was performed.

The remaining discrepancy suggests that a single mean-grain junction time is
insufficient. The next discriminating experiment/model step should constrain
junction mobility independently (grain-growth-only holds versus temperature
and initial G), then test a state-resolved junction population or pore/junction
pinning law. Such a law should make connected pore or triple-junction density
retard migration more strongly than boundary diffusion. Late-stage closed-pore
physics remains out of scope here.

## Reproduction

```bash
python3 growth_mechanism_sensitivity.py --workers 4
python3 -m pytest -q
```

Figures include Chen-style classifications for all modes, boundary curves,
window width versus grain size, and thermal mobility activation. Exact pore
conservation, nonnegative bins, baseline recovery, uniform budgets, and success
classification requirements are covered by automated tests.
