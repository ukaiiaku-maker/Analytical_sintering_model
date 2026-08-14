# Expanded fixed-model phase-space topology

## Scope and execution

This is mechanism exploration with the existing fixed model, not calibration or validation. No activation energy, prefactor, topology threshold, time budget, or success criterion was changed. PR #2 remains draft.

The main calculated campaign contains 992 first-step groups, 18,848 coarse second-step trajectories, and 8,134 automatically selected 10 C boundary trajectories. It also contains 3,780 heating-rate/peak-temperature trajectories, 60 expanded matched histories, and 234 combined rapid-T1/two-step trajectories. The calculated upper-size and onset extensions add 2,852 and 700 second-step trajectories, respectively. Total two-step trajectories are 30,534, expanded into 1,282,428 explicit target/tolerance classifications. Main-campaign wall time was 2,720 s using four local workers; extensions add about 270 s.

Every second step starts from a simulated first-step state. There are 2,926 trajectory rows with unattainable first steps in the main design. Targets already reached during step 1 are explicitly labeled `INELIGIBLE_TARGET_ALREADY_REACHED` and are not scored as second-step windows.

## Two-step phase-space topology

### Size onset and persistence

At rho target 0.90 and the primary 5% growth criterion:

- no gate produces a window through G0=300 nm;
- nominal-state refinement gives a zero-width onset near 500 nm for hybrid topology and 550 nm for density/connectivity;
- across the broader T1/switch space, isolated zero-width states appear near 450 nm;
- finite windows are repeatable by about 600 nm;
- the window continues through 2000 nm and remains present at diagnostic extensions of 3, 5, and 10 micrometers.

The model therefore predicts a large-particle domain rather than the desired nanoscale domain. It does not show an upper closure within its defensible range. At 3--10 micrometers the reduced model is outside its intended fine-particle interpretation; those points show that the mathematical window persists/asymptotes rather than bounding a physically credible conventional-sintering regime.

At the nominal T1=1300 C, rho switch=0.825 state, the density-gate minimum required tolerance decreases from 23.8% at 325 nm to 9.3% at 400 nm, 6.0% at 475 nm, 5.1% at 500 nm, and 4.8% at 550 nm. Connectivity follows nearly the same descent. Nanoscale values are approximately 50--52%, not near 5%.

### Boundary surfaces

For a representative nominal state:

| G0 | gate | lower density T2 | upper 5% T2 | window |
|---:|---|---:|---:|---:|
| 300 nm | density | 1215 C | 1150 C | none |
| 300 nm | connectivity | 1215 C | 1165 C | none |
| 600 nm | density | 1175 C | 1175 C | 0 C |
| 600 nm | connectivity | 1190 C | 1195 C | 5 C |
| 1000 nm | density | 1165 C | 1185 C | 20 C |
| 1000 nm | connectivity | 1200 C | 1205 C | 5 C |
| 2000 nm | density | 1165 C | 1195 C | 30 C |
| 2000 nm | connectivity | 1215 C | 1225 C | 10 C |

Topology changes boundary location and width but not the fundamental size scaling. Connectivity shifts both bounds upward by roughly 15--50 C and generally narrows the window relative to density/hybrid gates. It does not restore nanoscale eligibility.

At the primary target/tolerance, the main classification table contains 551 density-gate, 351 connectivity-gate, and 182 hybrid-gate successes. Failures span both required regions. Most low-temperature failures are classified kinetic exhaustion, with a smaller time-budget-exhaustion subset; high-temperature failures are grain-growth limited. The explicit reduced topology thresholds rarely become the terminal failure criterion, so the large-size onset is not primarily a connectivity cutoff.

### Density target

Among eligible core states, 5% windows occur at rho target 0.88 and 0.90. Density 0.92 is reached by many schedules but never with <=5% growth. No core schedule reaches rho target 0.94, 0.95, 0.96, or 0.98. This locates the current late-stage deficit between roughly 0.92 and 0.94 rather than merely saying full density fails.

Late-stage/closed-pore physics is therefore `REQUIRED` before full-density Chen-style maps are credible.

### Initial density, T1, and switch history

Changing rho0 from 0.60 to 0.80 shifts widths modestly but does not change the main threshold: 75 nm always fails; 1000 and 2000 nm succeed for all five rho0 values. At 450 nm, only an isolated density-gate rho0=0.80 condition gives a zero-width window.

T1 from 1150--1450 C changes window location but does not rescue 75 nm. At 450 nm, zero-width windows require switch density near 0.875. At 1000/2000 nm, windows occur across all T1 values and switch densities up to about 0.875. Switches 0.90/0.925 are often ineligible for a 0.90 target or unattainable and are never used to manufacture success.

## Why nanoscale particles are penalized

The dominant cause is the existing fractional grain-growth scaling. Clean growth uses approximately dG/dt proportional to 1/G, so fractional growth dG/G scales approximately as 1/G^2. Small G therefore incurs a very large fractional penalty before T2 is hot enough to densify. Densification also contains a transport time proportional to G^2, which eventually penalizes large grains and shifts their lower T2 upward, especially with the connectivity gate. Over the explored range, however, the reduction in fractional growth dominates enough to keep a large-G window open.

Assessment:

- fractional grain-growth scaling as nanoscale penalty: `SUPPORTED`;
- topology as the source of the onset: `INCONSISTENT_WITH_SWEEP`;
- topology as a boundary displacement: `SUPPORTED`;
- density/T1 history as a modulator: `SUPPORTED`;
- density/T1 history as the fundamental cause: `INCONSISTENT_WITH_SWEEP`.

## Matched-history identifiability

Expanded histories use G0=35, 150, 600, 1000, and 2000 nm; matched densities 0.80/0.85; and 0.1/100 C/min histories. Maximally different natural states show large topology, G1, and follow-up differences. As connectivity and G1 are progressively matched, follow-up density differences fall to roughly 0.0003--0.0014 in the best available pairs.

This supports the conclusion that most apparent size/history memory is mediated by the state reached during step 1. It does not prove topology has no independent content, because natural trajectories do not span density, connectivity, pore bins, and G1 independently. Independent topology causality remains `NON_IDENTIFIABLE`.

## Fast-firing response surface

The 0.1 C/min cases are retained, but many cannot reach the requested peak within the fixed time cap. HR is therefore reported against both 0.1 and the established 0.2 C/min reference, only when both paths reach rho=0.90.

At peak 1450 C and rho0=0.70:

- density gate: positive HR from 25--300 nm (about +1.5 to +6.4%), reverses near 450--600 nm, then becomes weakly positive above about 800 nm;
- connectivity gate: negative below about 100 nm, positive around 100--300 nm, negative near 450--600 nm, and weakly positive above about 800 nm;
- hybrid gate: mostly negative, with only narrow/weak positive regions.

For favorable cases the response largely saturates by 20--100 C/min; the calculated optimum is often at 50 or 100 C/min but differences across that high-rate band are small. Extremely high rate does not create a universal reversal, but the sign changes nonmonotonically with G0 and gate.

Peak-temperature surfaces show that rates cannot be interpreted independently of target/actual peak and attainment. High rate creates a different trajectory through redistribution/topology state, but it does not rescue an insufficient peak temperature.

## Rapid T1 heating followed by T2

At G0=75 nm, no 0.2/5/20 C/min T1 history creates a 5% two-step window. At 800 and 2000 nm, changing T1 heating rate can add/remove individual successful T2 points, especially for connectivity, but it does not move nanoscale states into the eligible regime. Rapid preparation is therefore a window modulator, not the missing nanoscale mechanism.

## Answers to the mechanism questions

1. Any 5% window first appears as an isolated zero-width state around 450--550 nm, depending on gate/history.
2. It broadens from roughly 600 nm upward.
3. It does not disappear by 2000 nm and persists mathematically to 10 micrometers; the latter is outside the defensible fine-particle domain.
4. Required tolerance is ~50% at nanoscale sizes, ~20% near 300 nm, ~5% near 500 nm, and ~2--3% above 600--2000 nm.
5. Small G is grain-growth/mixed limited; low T2 is kinetic exhaustion; high T2 is grain-growth limited; some high-switch states are first-step unattainable/ineligible.
6. Topology moves both T2 boundaries and narrows connectivity windows but does not control onset scaling.
7. Initial density modestly shifts width/location but not the size threshold.
8. T1/switch history selects isolated onset states and shifts boundaries but does not rescue nanoscale G0.
9. Nanoscale particles are penalized by dG/G scaling approximately as 1/G^2.
10. The penalty is primarily the grain-growth law coupled to the temperature needed for densification, not topology depletion.
11. Density-gate fast firing is favorable mainly at 25--300 nm and weakly above ~800 nm; other gates have narrower/nonmonotonic domains.
12. Benefit generally saturates over 20--100 C/min and can reverse with size/gate.
13. Rapid T1 heating modulates large-size eligibility but does not create nanoscale eligibility.
14. The phenomena partially overlap through pore redistribution/history, but two-step eligibility is dominated by fractional growth/densification overlap while fast firing is a transient ramp competition.
15. The minimal missing ingredients suggested are a physically revised nanoscale grain-growth suppression/timescale-separation mechanism and explicit late-stage closed-pore densification physics. These should be tested separately after this fixed-model audit.

## Files

`results/expanded_phase_space/` contains trajectory tables, the 1.28-million-row compressed classification table, calculated boundary surfaces, manifest/runtime files, matched-state tiers, fast-firing sign boundaries, upper/onset extensions, and boundary-focused figures. Run the main campaign with `expanded_phase_space_exploration.py`; targeted extensions and analysis have separate reproducible scripts.
