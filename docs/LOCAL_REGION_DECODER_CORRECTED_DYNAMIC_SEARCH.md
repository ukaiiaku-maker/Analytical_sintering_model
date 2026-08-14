# Decoder-corrected local-region dynamic search

## Outcome

The completed campaign sampled 1,000,000 decoded vectors, retained 20,000 unique Stage-1 dynamic fingerprints, and exactly reconfirmed 1,000 unique fingerprints. Exact runs use a 30-minute maximum step, adaptive density increments, a shared matched first-step state at `rho_switch = 0.88`, and 50 °C boundary discovery followed by 25 °C and 10 °C refinement. Uniform stagnation and 6,000-step criteria reject exhausted or stiff trajectories; targets and time budgets are common to all cases.

The exact trajectory/window filter produced 184 provisional Tier-B cases. The preparation audit rejected 178, primarily because first-step growth exceeded 20%. Six candidates retained a complete finite Chen window, active closed-pore support, exact switch transfer, matched-density high-density benefit, and first-step growth below 20%. No Tier A candidate was found, so the result is conditional Tier B rather than validation.

Best candidate 693168 has:

- exact switch density 0.87999936;
- `G0 = 103.0 nm`, `G1 = 117.1 nm`, first-step growth 13.7%;
- median reduction 89.4%, minimum reduction 88.3%, and `span20 = 0.030` over `rho = 0.95–0.98`;
- second-step success from 925 to 1200 °C;
- density-exhaustion failures through 920 °C and grain-growth failures from 1210 °C;
- finite window width 275 °C.

The very large predicted reduction and broad window are prototype-scale and should not be interpreted quantitatively. The successful candidates include closed fractions at the switch from roughly 0.02 to 0.94, which is a broad and partly suspicious range requiring experimental calibration.

## Numerical corrections

Three corrections were necessary before accepting the map: bounded thermal exponents/analytic grain updates, recording the state that crosses the target, and adaptive integration to prevent first-step overshoot. Recomputing each T2 from a single cloned first-step state enforces identical preparation across the map.
