# Experimental calibration plan

## Objective and evidence discipline

The objective is to test the two dominant mechanism hypotheses independently before fitting the full reduced model. Fast firing should be tested as nucleation-limited onset behavior. Two-step sintering should be tested as PR-prepared closed-pore accommodation memory with a lower density-exhaustion boundary and upper grain-growth boundary.

Calibration should not begin by fitting all model parameters to one shrinkage curve. Each major state or timescale should first be constrained by an observable that was not used to select candidate 693168. The final schedule comparison should then be treated as a blind prediction. Candidate 693168 remains conditional Tier B throughout this process.

## Priority 1 — Matched-density grain distributions

### Measurements

- `G_mean`, `G50`, and `G90` for slow and fast heating paths;
- the same quantities for high-temperature isothermal and two-step schedules;
- density and temperature histories sufficient to interpolate at matched density;
- distribution width and abnormal-growth tail, not mean grain size alone.

### Minimum protocol matrix

Use at least two widely separated heating rates, including the canonical slow/fast pair where experimentally practical, and at least three two-step second temperatures representing lower failure, expected success, and upper growth. Retain a common target density and time budget. Sample multiple densities through the predicted effect interval rather than only the final state.

### Mechanism constrained

This directly tests `R_fast(rho)` and the high-temperature/two-step `G(rho)` separation. The fast-firing hypothesis predicts that the fast path becomes finer after the onset interval. The two-step hypothesis predicts that a successful second step retains the prepared grain distribution while the high-temperature reference coarsens.

### Falsification criterion

The mechanism is weakened if matched-density separation appears only in the mean but not in `G50/G90`, if either path fails to attain the interval, or if the sign and density span are not reproducible.

## Priority 2 — Interrupted-ramp densification onset

### Measurements

- shrinkage/density rate during ramps;
- interrupted density and grain size before, through, and after onset;
- a nucleation proxy if available, such as event statistics, defect-source activation, or a reproducible incubation time;
- replicate ramps to distinguish onset distributions from instrument lag.

### Mechanism constrained

This constrains `tau_nuc`, the low-activity exposure, and the relative timing of exchange and transport. A complementary nucleation-facile material or processing condition is especially discriminating.

### Expected signature

Slow and fast ramps should cross the same material onset regime with different dwell times. Making nucleation facile should substantially reduce the fast-firing grain-size benefit.

### Falsification criterion

The nucleation-limited interpretation is rejected if the same fast-firing separation persists when no incubation/onset interval is detectable or when nucleation is independently shown to be facile.

## Priority 3 — Three-dimensional open/closed pore fraction

### Measurements

- 3D open and closed pore volumes at the initial state, first-step switch, selected second-step interruptions, and target density;
- pore connectivity to the external surface and grain-boundary network;
- uncertainty from segmentation, limited volume, and resolution;
- the fraction of remaining pore volume that is closed, reported separately from absolute porosity.

### Mechanism constrained

This is the primary test of the PR-prepared closed-store hypothesis and the mapping between internal model stores and measurable porosity.

### Expected signature

The high-temperature first step should prepare a schedule-dependent closed inventory that persists into the successful second step and differs from a matched-density comparator.

### Falsification criterion

Candidate-693168's interpretation is rejected if no persistent schedule-dependent closed inventory is observed, or if measured closed fractions are incompatible with the large modeled state even after accounting for the proxy mapping.

## Priority 4 — Pore D50, D90, and large-pore tail

### Measurements

- volume-weighted and number-weighted pore D50 and D90;
- large-pore-tail volume fraction;
- pore number density and surface-area proxy;
- distributions at matched density for first-step, direct-ramp, high-temperature continuation, and two-step paths.

### Mechanism constrained

These measurements constrain conservative PR/surface redistribution and test whether topology memory has an observable pore-distribution signature.

### Expected signature

Slow exposure and first-step preparation should change the pore distribution at matched density rather than only changing an inferred scalar topology state.

### Falsification criterion

The pore-redistribution interpretation is weakened if distributions remain schedule independent within uncertainty while the model requires substantial PR preparation.

## Priority 5 — Connected fine-pore fraction and percolation

### Measurements

- connected fine-pore volume and number fractions;
- pore-bearing grain-boundary coverage;
- connectivity/percolation to free surfaces and along grain boundaries;
- pore-location fractions at grain-boundary segments, triple junctions, and isolated interiors where resolvable.

### Mechanism constrained

These data constrain densification eligibility, removable pore inventory, connected-to-closed transition, and migration drag without relying on density as a topology proxy.

### Expected signature

The prepared first-step state should differ in connected removable inventory and retain part of that difference during the second step.

### Falsification criterion

The topology interpretation is weakened if paths with different predicted behavior have indistinguishable connected fine-pore and percolation states at matched density.

## Priority 6 — Trapped gas or accommodation proxy

### Candidate measurements

- internal pore pressure by a suitable spectroscopy, mechanical inference, or post-quench reconstruction;
- gas content and species;
- pore-shape/volume response during isothermal holds;
- creep or diffusional accommodation response under comparable temperature and pressure;
- a calibrated proxy for accommodation capacity, use, and recovery.

### Mechanism constrained

This is the most direct test of the finite closed-accommodation state that supports high-density shrinkage and creates the lower T2 boundary.

### Expected signature

Low-T2 density exhaustion should coincide with insufficient closed shrinkage or depleted/slow accommodation. Successful T2 conditions should retain enough accommodation to reach the target. The state should remain finite rather than acting as an unlimited sink.

### Falsification criterion

The current candidate interpretation is rejected if low-T2 exhaustion occurs without a measurable accommodation limitation, or if full accommodation remains available while the model requires depletion.

## Priority 7 — Grain-growth mobility across the upper T2 boundary

### Measurements

- isothermal grain-growth mobility over a fine T2 grid bracketing the predicted upper boundary;
- grain-distribution evolution, not only mean size;
- pore/topology state at the start of each mobility measurement;
- sufficient temporal resolution to separate an onset threshold from accumulated hold time.

### Mechanism constrained

This constrains the upper-bound migration/growth activation and its separation from closed-pore densification.

### Expected signature

Grain growth should increase sharply enough across the upper boundary to transform successful densification into growth failure while density remains attainable.

### Falsification criterion

A complete Chen-window interpretation is rejected if no upper migration activation is observed or if every T2 remains below the fixed growth tolerance.

## Priority 8 — Interrupted first/second-step tomography

### Measurements

Acquire 3D states at:

1. the initial condition;
2. several points during first-step heating;
3. immediately at the switch density;
4. early, middle, and late second-step holds for lower-failure, successful, and upper-growth T2;
5. the target or exhausted final state.

Track grain distribution, pore distribution, open/closed fraction, connectivity, pore location, and any accommodation proxy in the same specimens or statistically matched replicates.

### Mechanism constrained

This tests whether first-step topology memory persists long enough to alter second-step migration without directly suppressing densification.

### Falsification criterion

The memory hypothesis is rejected if the prepared state immediately collapses onto the isothermal comparator before the grain trajectories diverge.

## Calibration sequence

1. **Calibrate observables, not outcomes:** determine segmentation/resolution uncertainty and map model pore stores to 3D measurements.
2. **Fit fast-firing timescales independently:** constrain nucleation onset, exchange, and transport using ramp data without two-step outcomes.
3. **Fit pore preparation independently:** constrain PR redistribution and connected-to-closed transition using interrupted first-step pore data.
4. **Fit high-density support:** constrain closed shrinkage and accommodation with low/intermediate T2 holds.
5. **Fit upper migration separately:** constrain grain-growth mobility from the upper-bound grid.
6. **Blind validation attempt:** predict held-out heating rates, switch conditions, and T2 values with all parameters frozen.

The word “validation” should be reserved for the final held-out comparison and only if the observable state trajectories and processing outcomes are both reproduced. Until then, results remain calibration and falsification tests.

## Recommended calibration table

| Model quantity | Measurement | Fit stage | Held-out check | Primary uncertainty |
|---|---|---|---|---|
| `tau_nuc` / onset | interrupted ramp onset | fast kinetics | intermediate heating rates | onset distribution/instrument lag |
| exchange/transport times | relaxation or tracer proxy | fast kinetics | matched-density attainment | mechanism separation |
| closed fraction | 3D tomography | pore preparation | second-step trajectory | segmentation/resolution |
| pore D50/D90/tail | 3D imaging/scattering | PR redistribution | held-out ramps | sampling volume |
| connected fine fraction | connectivity analysis | topology eligibility | held-out first steps | percolation threshold |
| accommodation state | pressure/mechanical proxy | high-density support | low-T2 exhaustion | proxy-to-state mapping |
| growth mobility | isothermal T2 grid | migration | upper-bound location | initial topology control |

## Deliverable decision points

- **Supports fast mechanism:** reproducible matched-density benefit plus measurable onset timing and nucleation-facile weakening.
- **Supports two-step mechanism:** persistent prepared closed/topology state, low-T2 accommodation exhaustion, and upper-T2 migration activation.
- **Rejects candidate 693168:** required closed-store/accommodation state is absent or incompatible with measured magnitude.
- **Ready for quantitative validation attempt:** independently calibrated state trajectories predict held-out schedules without retuning.
