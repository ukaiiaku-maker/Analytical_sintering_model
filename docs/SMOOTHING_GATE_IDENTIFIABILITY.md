# Smoothing density-gate identifiability

## Scope and verdict

This audit asks whether the density gate in the conservative pore-bin
redistribution law is physical and identifiable. It does not optimize the gate
or claim validation.

The logistic audit spans:

- `smoothing_rho_mid = 0.72, 0.76, 0.79, 0.82, 0.86`;
- width `0.0075` (narrow), `0.015` (baseline), and `0.03` (broad);
- `rho0 = 0.65, 0.70, 0.75, 0.80, 0.85`;
- targets 0.85, 0.88, 0.90, and 0.92;
- 75 gate/initial-density conditions and 300 fixed-budget trajectories.

A bounded comparison also evaluates logistic and linear-clipped gate forms at
the baseline center and width. All other initial and mechanism parameters are
identical. Canonical protocol budgets are unchanged.

The central result is unusually sharp: a positive HR crossover occurs if and
only if `rho0 < smoothing_rho_mid` for all 75 logistic conditions. Every
observed crossover first appears at the sampled target 0.90. Gate placement
therefore moves which green densities benefit; it does not destroy the effect.

## Local gate formulation

The model now exposes `smoothing_gate_form` with two validated local forms:

```text
logistic:
    gate = sigmoid((rho_mid - rho) / rho_width)

linear_clipped:
    gate = clip(0.5 + (rho_mid - rho)/(2 rho_width), 0, 1)
```

Both are 0.5 at `rho=rho_mid`, bounded, monotone, and depend only on
instantaneous density and material parameters. The redistribution law remains
schedule-label-free.

## Gate-center and width response at rho = 0.90

Representative baseline-width HR results are:

| rho_mid | rho0=0.65 | 0.70 | 0.75 | 0.80 | 0.85 |
|---:|---:|---:|---:|---:|---:|
| 0.72 | 4.48 | 3.28 | -3.89 | -6.18 | -7.26 |
| 0.76 | 4.56 | 4.87 | 2.45 | -5.21 | -7.22 |
| 0.79 | 4.58 | 5.01 | 4.95 | -1.11 | -6.95 |
| 0.82 | 4.58 | 5.04 | 5.45 | 4.04 | -5.10 |
| 0.86 | 4.60 | 5.07 | 5.56 | 5.87 | 2.71 |

Narrow and broad gates change the transition magnitude but preserve the same
crossover classification. A broad gate leaves more redistribution above its
center and makes the HR transition smoother; a narrow gate approaches a sharp
cutoff.

At 0.90, correlation with `rho0-rho_mid` is strong:

| Width | corr(delta rho, HR_pct) | corr(delta rho, slow redistribution) |
|---|---:|---:|
| narrow | -0.81 | -0.89 |
| baseline | -0.85 | -0.92 |
| broad | -0.90 | -0.95 |

The points do not collapse perfectly because changing `rho0` also changes total
pore volume and instantaneous topology. Still, relative position to the gate is
the dominant organizing coordinate.

## Functional-form comparison

At baseline `rho_mid=0.79`, both forms preserve the same positive/negative
classification. Far from the center they are nearly identical. Close to the
center they differ materially:

| rho0 | Logistic HR_pct | Linear HR_pct | Logistic redistribution | Linear redistribution |
|---:|---:|---:|---:|---:|
| 0.65 | 4.58 | 4.58 | 0.2444 | 0.2444 |
| 0.70 | 5.01 | 5.03 | 0.2135 | 0.2140 |
| 0.75 | 4.95 | 5.50 | 0.1727 | 0.1828 |
| 0.80 | -1.11 | -3.66 | 0.0577 | 0.0282 |
| 0.85 | -6.95 | -7.26 | 0.00225 | 0 |

Thus endpoint processing results can locate a transition region but do not
uniquely identify width or functional form.

## Attainability and other diagnostics

All eligible points at 0.85 and 0.88 are attainable for both comparisons. At
0.90 all HR comparisons and 74/75 TS comparisons are attainable. At 0.92 none
of the 75 conditions supports both paths for either score. Gate changes do not
repair missing late-stage physics.

Every result row includes time to target, grain size, cumulative redistribution,
removable fine-pore fraction, pore mean radius, large-pore fraction, `f_pore`,
isolation, and `E_G` for slow, fast, high-T, and two-step paths. Failed targets
are retained with `NaN` scores.

## Is the density gate identifiable from processing data alone?

Only partially. If multiple known green densities straddle the transition,
the presence or absence of an HR crossover brackets `rho_mid`: in this model
the crossover set directly follows `rho0 < rho_mid`. Processing data alone does
not uniquely identify gate width, functional form, or smoothing rate because
each can change the integrated redistribution near the transition. TS data add
almost no constraint.

A single initial density is especially non-identifying: moving the center,
broadening the gate, or changing the redistribution rate can produce similar
integrated pore shifts and HR values.

## Best interrupted-ramp calibration measurement

Prepare otherwise matched compacts with `rho0` around 0.75, 0.80, and 0.85,
apply the same 0.2 C/min ramp, and quench replicate specimens near 800, 950,
1050, and 1200 C before the final hold. At each interruption measure the full
pore-volume distribution—especially volume-weighted mean radius and removable
fine-pore fraction—using a consistent pore-sensitive method such as calibrated
SAXS combined with microscopy/tomography.

The most constraining observable is the onset and slope of fine-pore depletion
versus temperature for several known `rho0`, not final density alone. The model
predicts strong separation between `rho0<=0.75` and `rho0>=0.80` in the
900--1200 C range for the baseline gate. Simultaneously recording shrinkage
keeps activity/density evolution aligned with the pore measurements.

## Does shifting the gate move or destroy the HR effect?

It moves the effect. At `rho_mid=0.72`, only `rho0=0.65` and 0.70 benefit at
0.90. At 0.86, all five starting densities benefit. The maximum positive HR
remains roughly 4--6% across centers; gate placement changes eligibility for
memory rather than the existence of the underlying redistribution response.

## Is TS insensitive to gate placement?

Yes, within the attainable region. At 0.90, TS is generally about 10.9--11.5%
across centers, widths, and starting densities. Logistic-to-linear changes are
below 0.001 percentage point in the baseline form comparison. One broad/high
gate case fails target and is unscored. This confirms that TS and HR constrain
different parts of the model.

## Should density be replaced by explicit topology?

Probably. The sharp rule `rho0 < rho_mid` shows that density is acting as a
proxy switch. A more physical state would gate smoothing by connected removable
fine-pore population, possibly combined with pore-boundary coverage and
isolation, so the transition emerges from the evolving pore network rather than
an imposed density threshold.

Replacement should follow, not precede, interrupted-ramp measurements. Those
data are needed to define and calibrate the connected fine-pore observable and
to determine whether density retains any independent role.

Exact tables and plots are in `results/smoothing_gate_identifiability/`.
