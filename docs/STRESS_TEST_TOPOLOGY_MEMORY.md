# Stress test of the topology-memory mechanism

## Scope and verdict

This is a mechanism audit, not a calibration or validation claim. The existing
memory law was not expanded and no new mechanism was added. Across the tested
schedule space, the heating-rate sign is robust and the two-step benefit exists
in a finite processing window. The result is not restricted to the original
0.2/20 C/min pair or the original temperature-window center. However, the
damage magnitude and its topology couplings remain uncalibrated, and the model
does not yet connect accumulated damage to a changed pore-bin distribution.

## Locality and anti-cheat audit

`topology_damage_rate(state, T_C, params, kinetics)` receives no protocol or
schedule object. It uses only instantaneous temperature, density, renewal
activity, accumulated damage, and fixed material parameters:

```text
damage_rate = k_damage
              * Gaussian(T; T_mid, width)
              * (1 - activity)^m
              * pre_densification_gate(rho)
              * (1 - topology_damage)
```

There is no protocol name, ramp-rate value, elapsed-time test, schedule class,
or slow/fast label in the rate law. The state is path dependent because this
local rate is integrated in time. Damage then acts through two named topology
channels: reduced removable pore coverage and increased pore isolation. Source
inspection in the tests rejects schedule-label identifiers in the function.

## Held-out heating rates

All seven ramps reach `rho_target = 0.90`.

| Rate (C/min) | G at target (nm) | Damage | Max damage rate (1/s) | Median E_G | HR_pct vs 0.2 |
|---:|---:|---:|---:|---:|---:|
| 0.2 | 733.25 | 0.8481 | 7.76e-6 | 0.615 | 0.00 |
| 0.5 | 604.44 | 0.5188 | 1.13e-5 | 0.894 | 17.57 |
| 1 | 532.98 | 0.3015 | 1.31e-5 | 0.961 | 27.31 |
| 5 | 466.41 | 0.0679 | 1.48e-5 | 1.176 | 36.39 |
| 10 | 459.17 | 0.0345 | 1.51e-5 | 1.120 | 37.38 |
| 20 | 453.97 | 0.0174 | 1.52e-5 | 1.086 | 38.09 |
| 50 | 452.34 | 0.0071 | 1.53e-5 | 1.061 | 38.31 |

The accumulated damage decreases smoothly with heating rate and the grain-size
benefit saturates above roughly 5--10 C/min. The maximum instantaneous damage
rate increases with ramp rate because slow ramps have already consumed much of
the remaining-undamaged fraction before reaching the rate maximum. Integrated
damage, rather than peak rate, controls the trajectory. This behavior is
physically interpretable and survives every held-out rate tested.

## Held-out two-step window

Of 27 combinations, 18 reach the target and 15 have positive `TS_pct`.

- All nine `T2 = 1200 C` cases fail to reach 0.90 within the fixed time budget.
- All nine `T2 = 1250 C` cases reach target and are beneficial. `TS_pct` spans
  approximately 5.38--16.99%.
- At `T2 = 1300 C`, the three `T1 = 1300 C` cases are isothermal controls and
  give exactly 0%; the six cases with `T1 = 1350` or `1400 C` are beneficial,
  with `TS_pct` approximately 5.78--12.14%.
- Switch densities 0.82, 0.85, and 0.88 preserve the qualitative result when
  the second-step temperature is viable.

The two-step advantage is therefore a finite processing window, not a universal
result. Failed 1200 C paths are retained as failures rather than scored.

## Target-density sweep

At targets 0.88 and 0.90, all four reference protocols reach target and both
percentage metrics remain positive:

| Target | HR_pct | TS_pct |
|---:|---:|---:|
| 0.88 | 29.33 | 10.87 |
| 0.90 | 38.09 | 11.29 |

At 0.92, only the 1350 C isothermal path reaches target. Final densities are
0.9069 (slow), 0.9091 (fast), 0.9200 (high-T), and 0.9031 (two-step). No HR or
TS score is reported at that target.

## Ablations

All eight ablations reach 0.90 for all four reference protocols.

| Ablation | HR_pct | TS_pct | Slow damage | Fast damage |
|---|---:|---:|---:|---:|
| Memory disabled | -6.45 | 11.53 | 0 | 0 |
| Memory enabled | 38.09 | 11.29 | 0.848 | 0.017 |
| Coverage only | 23.62 | 11.39 | 0.844 | 0.017 |
| Isolation only | 12.03 | 11.44 | 0.817 | 0.017 |
| Window center 850 C | 40.19 | 11.50 | 0.908 | 0.023 |
| Window center 1200 C | 29.32 | 11.15 | 0.644 | 0.011 |
| Damage rate / 2 | 27.03 | 11.41 | 0.597 | 0.009 |
| Damage rate x 2 | 42.78 | 11.18 | 0.979 | 0.035 |

Removable coverage is the stronger single channel, though isolation alone also
flips the HR sign. Their combined effect is nonlinear. Shifting the window
center by -175 C or +175 C and changing the damage rate by a factor of two all
preserve positive HR and TS metrics. The result is therefore not narrowly tied
to the original 1025 C center.

## Physical concerns and next step if the result fails validation

The current values are still suspicious in the calibration sense:
`surface_damage_rate_s = 2e-5 1/s` (a 13.9 h intrinsic timescale), a 180 C
Gaussian width, coverage strength 0.65, and isolation strength 0.35 have not
been inferred from experiment. The slow ramp reaches damage 0.85, so the effect
is large. Moreover, pore mean radius and large-pore fraction histories are
nearly identical among heating rates: memory currently changes effective
topology without transferring pore volume among bins. This is the main physical
gap, despite the schedule robustness.

Before any larger search, compare the damage state against measured surface
area, pore-size distribution, or removable fine-pore fraction during interrupted
ramps. If the topology-only memory fails that comparison, the next mechanism to
test should be a conservative surface-diffusion pore-bin redistribution law
that produces the same memory through observable pore coarsening. It should
replace, not stack on top of, the empirical damage coupling during ablation.

Exact tables and plots are in `results/topology_memory_stress/`.
