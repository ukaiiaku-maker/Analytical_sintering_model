# Density-window processing map

## Scope and status

This audit maps where the current conservative pore-bin memory model predicts
benefit, neutrality, harm, or unattainable densification. It does not require a
universal protocol advantage and does not constitute validation.

One kinetic/mechanism parameter set is used throughout. Three initial classes
vary only initial descriptors:

| Class | rho0 | G0 (nm) | Pore scale (nm) | Log-width | Initial mean r (nm) | Initial f_pore |
|---|---:|---:|---:|---:|---:|---:|
| loose_fine | 0.70 | 100 | 16 | 0.85 | 33.59 | 0.555 |
| baseline_intermediate | 0.75 | 150 | 22 | 0.65 | 36.76 | 0.560 |
| predensified_partially_isolated | 0.82 | 250 | 35 | 0.45 | 47.48 | 0.503 |

The descriptors are varied together as representative classes, so this map
does not identify independent causal sensitivities for `rho0`, `G0`, pore
scale, and width.

## Fixed protocol budgets

Every trajectory is run once to a fixed protocol budget and then sampled at all
targets. There is no target-specific extension or retuning.

- Heating ramps use the same 8 h hold at 1450 C. Total ramp-plus-hold budgets
  are 126.75, 55.5, 31.75, 12.75, 10.375, 9.1875, and 8.475 h for 0.2, 0.5,
  1, 5, 10, 20, and 50 C/min, respectively.
- High-T isothermal and two-step paths use 96 h.
- Targets are 0.80, 0.85, 0.88, 0.90, 0.92, 0.95, and 0.98.

The ensemble contains 210 unique trajectories and 1,470 protocol-target rows.
Of the eligible rows, 718 are explicitly recorded as unattainable.

## Canonical density-resolved comparison

The canonical heating comparison is 0.2 versus 20 C/min. The canonical
two-step comparison is 1350 C isothermal versus 1350→1250 C with switch density
0.85.

| Initial class | Target | HR_pct | TS_pct | Classification |
|---|---:|---:|---:|---|
| loose_fine | 0.80 | -21.67 | 0.00 | harmful |
| loose_fine | 0.85 | -15.01 | 0.00 | harmful |
| loose_fine | 0.88 | -6.28 | 9.93 | harmful |
| loose_fine | 0.90 | 1.84 | 11.42 | beneficial |
| baseline_intermediate | 0.80 | -7.04 | 0.00 | harmful |
| baseline_intermediate | 0.85 | -11.81 | 0.00 | harmful |
| baseline_intermediate | 0.88 | -5.31 | 9.10 | harmful |
| baseline_intermediate | 0.90 | 4.47 | 11.36 | beneficial |
| predensified_partially_isolated | 0.85 | -2.58 | 0.00 | harmful |
| predensified_partially_isolated | 0.88 | -11.60 | 4.62 | harmful |
| predensified_partially_isolated | 0.90 | -5.91 | 10.55 | harmful |

At 0.92 the canonical slow, fast, and two-step paths fail for every class while
the 1350 C high-T path reaches target. At 0.95 and 0.98 all four canonical paths
fail. These points are unscored, not assigned zero benefit.

Targets at or below `rho0` are marked ineligible for processing scores. For
example, the predensified class is already above 0.80 at time zero.

## 1. Where does fast heating help?

For the loose/fine and baseline classes, fast heating becomes beneficial only
at the late part of the attainable trajectory, `rho = 0.90`. At that target:

- loose/fine is positive for rates 1--50 C/min relative to 0.2 C/min, with
  `HR_pct` from 0.52 to 1.88%;
- baseline is positive for rates 1--50 C/min, with `HR_pct` from 2.25 to 4.54%.

No tested rate improves the trajectory at targets 0.80, 0.85, or 0.88. The
predensified class has no positive heating-rate point through its attainable
range. Thus the predicted fast-heating window requires both sufficient density
evolution and an initial state that spends time below the smoothing density
gate.

## 2. Where does two-step sintering help?

The canonical two-step path is neutral before its 0.85 switch, positive by
0.88, and more positive at 0.90 for all three classes.

The full grid retains a finite rather than universal window. Counting only
genuine two-step points where the high step is executed and the sampled target
lies above the switch:

| Class | Target | Scored | Positive | Mean TS_pct |
|---|---:|---:|---:|---:|
| loose/fine | 0.80 | 12 | 11 | 3.05 |
| loose/fine | 0.85 | 24 | 22 | 8.61 |
| loose/fine | 0.88 | 28 | 25 | 10.60 |
| loose/fine | 0.90 | 24 | 20 | 8.59 |
| baseline | 0.80 | 12 | 11 | 1.37 |
| baseline | 0.85 | 24 | 22 | 5.69 |
| baseline | 0.88 | 28 | 25 | 9.84 |
| baseline | 0.90 | 24 | 20 | 8.54 |
| predensified | 0.88 | 12 | 11 | 5.19 |
| predensified | 0.90 | 12 | 10 | 7.85 |

Switches at or below `rho0` are retained in the raw table but flagged because
they skip the high-temperature step and are low-T isothermal paths rather than
genuine two-step schedules.

## 3. Where do both effects coexist?

Under the canonical protocols, both effects coexist only at `rho = 0.90` for
the loose/fine and baseline classes. The baseline class has the larger HR
benefit (4.47%) while retaining essentially the same TS benefit (11.36%). No
coexistence appears for the predensified class before the attainable range ends.

## Density-window summary

All eligible canonical points are attainable in the early and intermediate
windows through 0.90. In the intermediate window:

- loose/fine: 33% of scored targets have positive HR and 67% positive TS;
  mean `HR_pct = -6.48%`, mean `TS_pct = 7.11%`;
- baseline: 33% positive HR and 67% positive TS; mean `HR_pct = -4.22%`,
  mean `TS_pct = 6.82%`;
- predensified: 0% positive HR and 67% positive TS; mean `HR_pct = -6.70%`,
  mean `TS_pct = 5.06%`.

Canonical fraction attainable is zero in the late and near-final windows
because a comparison is scored only when all required paths reach the target.

## 4. Where does the model fail to densify?

- Through 0.85, every eligible heating, high-T, and two-step path reaches target.
- At 0.88, all heating and high-T paths reach; 52/60 two-step paths reach for
  loose/baseline and 54/60 for predensified.
- At 0.90, all heating and high-T paths reach, but only 36/60 two-step paths
  reach for each class.
- At 0.92, no heating path and no two-step path reaches. Two of three high-T
  paths reach for each class.
- At 0.95 and 0.98, none of the 70 protocols per class reaches target.

The dedicated `unattainable_cases.csv` retains every failure with `NaN` target
time, grain size, topology, and efficiency metrics.

## 5. Which failures indicate missing late-stage physics?

The sharp loss of attainability above 0.90 is common to all three initial
classes and all heating rates. It is therefore not explained by one initial
pore distribution or schedule. The current topology loop increases isolation
while connectivity and removable coverage fall, but it has no distinct
late-stage isolated-pore shrinkage, gas-pressure balance, vacancy transport to
remote sinks, or closed-pore thermodynamics. Targets 0.95--0.98 are classified
`requires late-stage physics`, not merely unfavorable processing.

## 6. Which initial descriptors shift the optimum?

The baseline class gives the strongest HR crossover. At 0.90 its slow/fast mean
pore-radius contrast is 56.54/47.50 nm and cumulative redistribution contrast is
0.170/0.0022. The loose/fine class begins with a broader distribution and more
large-pore volume, so further redistribution produces a smaller relative
contrast and only 1.84% HR benefit. The predensified class starts above the
pre-densification smoothing gate; slow/fast cumulative redistribution is only
0.0224/0.00024, insufficient to reverse the original growth ordering.

This suggests that initial density relative to the smoothing gate and the
available removable fine-pore population shift the optimum most strongly.
Because `rho0`, `G0`, pore scale, and width were varied together, a factorial
sensitivity audit is required before assigning independent causality.

## Classification rule

- `beneficial`: both canonical percentage metrics exceed +0.5%;
- `neutral`: the target is already present initially, or both metrics are
  within ±0.5%;
- `harmful`: an attainable canonical point does not meet the beneficial or
  neutral definition;
- `unattainable`: a required path fails below 0.95;
- `requires late-stage physics`: a required path fails at or above 0.95.

Exact tables and plots are in `results/density_window_processing_map/`.
