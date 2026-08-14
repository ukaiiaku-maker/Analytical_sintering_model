# Matched-topology and fast-firing mechanism discrimination

## Fixed design

This bounded audit follows the process-window calculation without recalibration. Matched-density histories use G0 = 75, 150, 300 nm, density = 0.80/0.85, slow/fast ramps = 0.2/20 C/min, and a common 1250 C, 96 h follow-up. M0/M1/M2 are represented by the existing density, connectivity, and hybrid-topology smoothing gates; the underlying topology-aware densification physics is unchanged.

The fast-firing audit uses G0 = 50--450 nm, rho0 = 0.65/0.70/0.75, and 0.2/5/20 C/min for density and connectivity gates. Complete state trajectories are evolved during heating. The small two-step initial-condition probe varies only G0, rho0, and gate at T1=1300 C, switch density 0.85, and the seven existing T2 values.

## Matched density: is density sufficient?

All 36 history states reached their requested matched density. Slow and fast histories at equal density retain visibly different states. For example, at G0=150 nm and rho=0.80, density-gate connected coverage is 0.332 after the slow history and 0.436 after the fast history; G1 is 160.5 versus 176.3 nm. The same 1250 C follow-up then gives density gains of 0.1006 and 0.1022.

Across matched G0/density pairs, slow-fast differences in subsequent density gain span:

- density gate: 0.0010--0.0021;
- connectivity gate: 0.00048--0.00153;
- hybrid gate: 0.00028--0.00109.

Thus different histories at the same density do produce measurable but small response differences. Density is not a complete state descriptor. However, topology and G1 change together, and nearest-coverage pairs across sizes still have substantial G1 differences. The separate causal contribution of topology is therefore `NON_IDENTIFIABLE` in this dataset.

Initial size acts primarily through the first-step state, but persistent grain-state memory remains. Approximately matched density/coverage does not eliminate large G1 differences, which strongly alter fractional follow-up growth. The current reduced state cannot construct truly equal-density/equal-topology/equal-G comparisons without synthesizing states, which this audit deliberately avoids.

## Fast firing

The density gate produces positive HR at rho=0.90 for 18/21 size-density combinations. It is positive from 50 through 300 nm across all three initial densities, with HR from about +1.2% to +6.3%, and becomes slightly negative near 450 nm (-0.9% to -1.2%). This is a substantially broader domain than the two-step window.

The connectivity gate is positive for 11/21 combinations and nonmonotonic with size. It is harmful for the smallest particles, positive mainly around 100/150--300 nm depending on rho0, and negative again at 450 nm. Its HR range is -4.19% to +3.80%. Connectivity therefore does not improve cross-size robustness; it over-suppresses slow-ramp smoothing for some nanoscale states.

The initial-density dependence is weak for the density gate and material for the connectivity gate. The small two-step probe finds no 5% success for any G0/rho0/gate combination, so changing initial density does not rescue the missing two-step scale separation.

## Do the two phenomena share a mechanism?

The fixed model supports partially overlapping but distinct origins:

- fast firing is driven by transient ramp competition and accumulated pore redistribution; under the density gate it survives across a broad 50--300 nm domain;
- two-step behavior requires an overlap between continuing densification and suppressed fractional grain growth during T2; that overlap is absent at 5% and appears only narrowly at 10% for the largest size.

The contrast is scientifically useful but not the desired physical picture. The fast-firing mechanism is qualitatively supported. The nanoscale two-step mechanism is not.

## Mechanism assessments

| mechanism statement | assessment |
|---|---|
| density alone is a sufficient state | INCONSISTENT_WITH_SWEEP |
| topology independently identifies follow-up response | NON_IDENTIFIABLE |
| G1/initial-size memory affects follow-up growth | SUPPORTED |
| density-gate fast firing across 50--300 nm | SUPPORTED |
| connectivity gate improves fast-firing robustness | INCONSISTENT_WITH_SWEEP |
| upper-T grain-growth failure | SUPPORTED |
| lower-T densification exhaustion | SUPPORTED |
| nanoscale robust two-step window | INCONSISTENT_WITH_SWEEP |
| late-stage closed-pore mechanism for rho>=0.95 | REQUIRED |

## Measurements that discriminate the state variables

Prepare equal-density specimens with intentionally different pore connectivity, then measure 3D connected GB pore coverage, isolated-pore fraction, full pore-volume distribution, and G before a common T2 hold. The decisive response variables are initial rho-dot, fractional G-dot/G, and density gain over the early second step. To separate size from topology, specimens must be matched in density and G as well as connected coverage; the present trajectories do not naturally supply all three matches.

## Outputs

`results/mechanism_discrimination/` contains matched histories, nearest state pairs, fast-firing metrics and downsampled trajectories, the bounded initial-condition window probe, figures, and the requested model-style robustness matrix.

No calibration or validation is claimed. PR #2 remains draft.
