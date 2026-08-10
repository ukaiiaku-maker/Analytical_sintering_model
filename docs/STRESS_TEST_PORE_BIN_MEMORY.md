# Stress test of observable pore-bin memory

## Verdict

The conservative pore-bin redistribution mode reproduces both required signs
at `rho_target = 0.90` without empirical topology damage:

```text
HR_pct = +4.5990%
TS_pct = +11.5224%
```

All four reference protocols reach target. Unlike empirical topology damage,
the new mode produces a visible schedule-dependent pore distribution. The
heating-rate sign also survives all held-out rates at 0.90, and the two-step
benefit remains confined to the same finite processing window.

This is still not a validation claim. The heating-rate sign is negative at the
earlier target `rho = 0.88`, and the redistribution rate and fine-pore exponent
have not been calibrated to pore-distribution measurements.

## Memory modes

The model now validates and exposes four modes:

- `none`: no empirical or pore-bin memory;
- `empirical_topology_damage`: the previous effective-topology state;
- `pore_bin_redistribution`: the new default and candidate physical memory;
- `combined`: diagnostic only, not used as an accepted result.

The three non-combined modes give the following at `rho = 0.90`:

| Mode | HR_pct | TS_pct | Slow mean r (nm) | Fast mean r (nm) | Slow large fraction | Fast large fraction |
|---|---:|---:|---:|---:|---:|---:|
| None | -6.4503 | 11.5345 | 47.94 | 47.97 | 0.0938 | 0.0940 |
| Empirical damage | 38.0887 | 11.2949 | 47.93 | 47.94 | 0.0938 | 0.0938 |
| Pore-bin redistribution | 4.5990 | 11.5224 | 57.27 | 48.06 | 0.1379 | 0.0942 |

The empirical model changes inferred topology without changing the observable
pore distribution. The new mode depletes the slow-ramp removable fine-pore
fraction to 0.2149 versus 0.3549 fast. Its cumulative adjacent-bin transferred
volume is 0.1810 slow versus 0.00234 fast. This cumulative diagnostic counts
each adjacent-bin crossing, so pore volume transferred through multiple bins
can contribute more than once.

## Local conservative law

For each bin except the largest,

```text
J_i = k_smooth
      * Gaussian(T; T_mid, width)
      * (1 - activity)^m
      * pre_densification_gate(rho)
      * phi_i
      * (r_min/r_i)^q

phi_dot_i     -= J_i
phi_dot_{i+1} += J_i
```

The function receives only instantaneous state, temperature, material
parameters, and renewal kinetics. It has no protocol, schedule, ramp-rate, or
slow/fast input. Tests inspect the function source for schedule-label leakage.

The flux sums to zero to numerical precision and has `rho_dot = 0`.
Redistribution is included as a named power channel in the nonnegative
dissipation partition. After evolution, pore number is recomputed from bin
volume and bin radius, representing fewer pores when volume moves to a larger
bin. Density continues to be calculated only as `1 - sum(phi)`.

## Held-out heating rates at rho = 0.90

Every rate reaches target, and every held-out rate is better than 0.2 C/min:

| Rate (C/min) | G (nm) | HR_pct vs 0.2 | Cumulative transfer | Mean r (nm) | Large fraction |
|---:|---:|---:|---:|---:|---:|
| 0.2 | 473.10 | 0.00 | 0.1810 | 57.27 | 0.1379 |
| 0.5 | 469.66 | 0.73 | 0.0838 | 51.90 | 0.1082 |
| 1 | 462.95 | 2.14 | 0.0443 | 50.00 | 0.1006 |
| 5 | 451.98 | 4.46 | 0.00928 | 48.36 | 0.0951 |
| 10 | 452.47 | 4.36 | 0.00467 | 48.17 | 0.0946 |
| 20 | 451.34 | 4.60 | 0.00234 | 48.06 | 0.0942 |
| 50 | 451.13 | 4.64 | 0.00094 | 48.00 | 0.0940 |

Integrated redistribution decreases smoothly with heating rate. Grain size is
not perfectly monotonic between 5 and 10 C/min, but the small variation occurs
after the benefit has nearly saturated and is physically interpretable as
competition among redistribution, growth, and densification.

## Held-out two-step window

For each of the three memory modes, 18 of 27 two-step/high-T comparisons reach
0.90 and 15 are beneficial. In the pore-bin mode:

- all nine `T2 = 1200 C` paths fail to reach target;
- all nine `T2 = 1250 C` paths reach target and are beneficial;
- `T2 = 1300 C` is beneficial when `T1` is 1350 or 1400 C;
- the three `T1 = T2 = 1300 C` controls give exactly zero benefit.

Thus redistribution neither manufactures universal two-step success nor hides
failed low-temperature holds.

## Target-density sweep

| Mode | Target | HR_pct | TS_pct | Attainment |
|---|---:|---:|---:|---|
| None | 0.88 | -17.15 | 11.02 | all four |
| Empirical | 0.88 | 29.33 | 10.87 | all four |
| Redistribution | 0.88 | -5.39 | 11.03 | all four |
| None | 0.90 | -6.45 | 11.53 | all four |
| Empirical | 0.90 | 38.09 | 11.29 | all four |
| Redistribution | 0.90 | 4.60 | 11.52 | all four |

At 0.92, only the high-T isothermal run reaches target in every mode, so HR and
TS metrics are not reported. No target was lowered to rescue a case.

The negative redistribution HR sign at 0.88 is the main qualification. The
observable pore memory accumulates enough to overcome the original growth
ordering between densities 0.88 and 0.90, not throughout the entire trajectory.

## Parameters and next checks

The new law reuses the broad 1025 C, 180 C-width intermediate-temperature
window and uses `smoothing_rate_s = 2e-5 1/s` and fine-radius exponent 2. These
are named, local material parameters rather than hidden efficiency factors, but
they remain uncalibrated. Before claiming validation:

1. compare interrupted-ramp pore distributions with the predicted slow/fast
   bin shifts;
2. test bin-resolution convergence and replace fixed adjacent-bin jumps with a
   conservative finite-volume radius-space flux;
3. determine whether the sign crossover between 0.88 and 0.90 is observed;
4. infer the smoothing rate and activation window from surface-area or pore-size
   measurements rather than protocol outcomes.

Exact tables and plots are in `results/pore_bin_memory_stress/`.
