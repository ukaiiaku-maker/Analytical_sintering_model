# Initial-condition factorial sensitivity

## Scope

This audit separates the initial descriptors that were confounded in the first
processing map. It uses the current `pore_bin_redistribution` mechanism with one
unchanged parameter set and fixed canonical protocol budgets.

The bounded design contains:

- 11 one-at-a-time (OAT) conditions around `rho0=0.75`, `G0=150 nm`, pore
  scale `25 nm`, and log-width `0.65`;
- 16 low/high corner combinations using `rho0={0.65,0.85}`,
  `G0={75,300} nm`, pore scale `{15,40} nm`, and width `{0.45,0.85}`;
- 27 conditions, 108 trajectories, and 135 condition-target rows in total.

The large-pore tail fraction is not independently exposed by the current
initializer. Width and pore scale change the generated distribution, but this
audit does not pretend they are an independent tail control.

## Fixed protocols and scoring

Every condition uses 0.2 and 20 C/min heating ramps, 1350 C isothermal, and
1350→1250 C two-step with switch density 0.85. Budgets are identical across
conditions: 126.75 h slow ramp, 9.1875 h fast ramp, and 96 h for both high-T and
two-step schedules. Targets are 0.80, 0.85, 0.88, 0.90, and 0.92. A metric is
reported only if both paths reach an eligible target above `rho0`.

## Main result at rho = 0.90

### One-at-a-time effects

| Descriptor | Level | HR_pct | TS_pct |
|---|---:|---:|---:|
| rho0 | 0.65 | 4.58 | 11.08 |
| rho0 | 0.70 | 5.01 | 11.14 |
| rho0 | 0.75 | 4.95 | 11.21 |
| rho0 | 0.80 | -1.11 | 11.30 |
| rho0 | 0.85 | -6.95 | 11.46 |
| G0 (nm) | 75 | 6.11 | 11.24 |
| G0 (nm) | 150 | 4.95 | 11.21 |
| G0 (nm) | 300 | 1.20 | 9.11 |
| pore scale (nm) | 15 | 2.51 | 11.79 |
| pore scale (nm) | 25 | 4.95 | 11.21 |
| pore scale (nm) | 40 | 6.40 | 10.58 |
| log-width | 0.45 | 7.46 | 11.45 |
| log-width | 0.65 | 4.95 | 11.21 |
| log-width | 0.85 | 2.92 | 10.93 |

Every OAT condition that has positive HR first crosses at the sampled target
0.90. The requested density grid cannot resolve a finer crossover location.

### Corner-factorial main effects

At 0.90, averaged over the other corner factors:

| Factor change | Mean HR change | Mean TS change |
|---|---:|---:|
| rho0 0.85 → 0.65 | +11.24 points | +0.04 points |
| G0 300 → 75 nm | +4.23 points | +3.62 points |
| pore scale 15 → 40 nm | +3.96 points | +0.88 points |
| width 0.85 → 0.45 | +1.98 points | -0.12 points |

Initial density is therefore the dominant fast-heating control in this model.
Grain size and pore scale are secondary, while distribution width/removable
fine-pore inventory has a smaller but consistent effect. Two-step behavior has
a different sensitivity: `G0` is the clearest control, while `rho0` has almost
no factorial main effect.

## 1. What controls the fast-heating crossover?

The primary control is `rho0` relative to the pre-densification smoothing gate.
In the rho-only sweep, all pore-distribution and grain descriptors are fixed,
yet HR falls from about +5% for `rho0=0.65--0.75` to -1.11% at 0.80 and -6.95%
at 0.85. Slow-ramp cumulative redistribution at 0.90 falls correspondingly:

```text
rho0 0.65: 0.244
rho0 0.70: 0.213
rho0 0.75: 0.173
rho0 0.80: 0.0577
rho0 0.85: 0.00225
```

The fine-pore inventory is identical in this OAT sweep, so insufficient initial
fine pores cannot explain the rho effect. Narrowing the distribution separately
raises the removable fine fraction from 0.486 (width 0.85) through 0.572
(baseline) to 0.691 (width 0.45), and raises HR from 2.92 through 4.95 to 7.46%.
Thus removable inventory is a real secondary control, not the dominant one.

## 2. Is there a meaningful optimum initial state?

Within the tested bounds, the strongest joint result is corner condition
`rho0=0.65`, `G0=75 nm`, pore scale `40 nm`, and width `0.45`:

```text
HR_pct = 10.65%
TS_pct = 10.83%
```

The broader interpretation is more useful than this single corner: positive HR
is favored by green densities at or below about 0.75, small starting grains,
an appreciable absolute pore scale, and a narrow/fine-removable distribution.
The OAT rho sweep peaks weakly at 0.70 rather than at the lowest rho, so the
model does not imply “looser is always better.” The optimum remains a model
prediction, not a calibrated material prescription.

## 3. Does two-step benefit use the same descriptors?

No. Two-step benefit is nearly independent of initial density in the OAT sweep
(11.08--11.46%) and in the corner average. It is more sensitive to initial
grain size: increasing `G0` from 75 to 300 nm lowers the factorial mean TS by
3.62 points. Larger pore scale slightly improves HR but slightly reduces OAT TS,
showing an explicit tradeoff between the two processing objectives. Width has
only a weak TS effect.

Two broad/large-pore corner cases fail the two-step target at 0.90 and remain
unscored. At 0.92 neither canonical comparison is generally attainable.

## 4. Why does the predensified/coarser state fail?

The rho-only experiment isolates the main cause. Raising `rho0` to 0.85 while
holding `G0`, pore scale, width, and removable fine fraction fixed changes HR
from +4.95 to -6.95%. Raising `G0` alone to 300 nm still leaves HR positive at
1.20%. Therefore large `G0` and insufficient initial fine pores are not enough
to explain failure.

Starting at 0.85 lies above the smoothing density gate (`rho_mid=0.79`), so the
slow path accumulates almost no redistribution memory. Initial connectivity
also falls from 0.993 at rho0 0.75 to 0.731 at 0.85, and is a coupled secondary
effect of density in the topology law. The current design cannot vary
connectivity independently of density, but the explicit gate suppression and
the OAT result identify starting above the gate as the primary modeled cause.

## 5. What should be calibrated first?

Initial density itself is directly measurable; the quantity requiring
calibration is the density dependence of the smoothing gate. Interrupted-ramp
experiments across several green densities should measure the removable
fine-pore inventory and pore-size shift before densification. Those data can
locate or replace `smoothing_rho_mid` and its width. Among initial distribution
descriptors, removable fine-pore fraction (or the underlying distribution
width) should be measured first because it is observable and monotonically
shifts HR in the OAT audit.

## 6. What must change before rho > 0.95?

No initial-condition combination repairs the high-density failure. Before
attempting 0.95--0.98, the model needs a late-stage closed/isolated-pore module:
gas-pressure and capillary balance, vacancy transport from isolated pores to
remote sinks, and a topology transition for closed-pore shrinkage. Expanding
the initial-condition search cannot substitute for that missing physics.

## Anti-cheat audit

- All non-initial `Params` fields are asserted identical across all 27 cases.
- Protocol budgets are stored in every result row and tested for uniqueness.
- Unattainable comparisons have `NaN` scores.
- The default remains `pore_bin_redistribution`; `memory_model="none"` remains
  covered by the existing ablation tests.
- Pore nonnegativity, exact `rho=1-sum(phi)`, and schedule-free redistribution
  locality remain regression-tested.

Exact results and plots are in `results/initial_condition_sensitivity/`.
