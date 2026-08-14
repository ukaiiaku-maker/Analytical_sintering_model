# Publication-style sintering figures for candidate 693168

## 1. Motivation

The earlier candidate audit established numerical reproducibility but emphasized model diagnostics. This rebuild presents the fixed result in the order used for experimental sintering data: temperature history, density, grain size, grain-size–density trajectories, final outcomes versus second-step temperature, and filled processing maps. Model-internal states follow only after those observables.

## 2. Scientific status

Candidate 693168 remains a **conditional Tier B prototype**. It is not Tier A, calibrated, validated, or paper-ready. Its parameters, local laws, density target (0.98), growth tolerance (20%), and 500 h budgets are frozen. No optimization or parameter search is performed.

The switch is near density 0.88, with initial mean grain size near 103 nm, prepared mean size near 117 nm, and first-step growth near 13.7%. The original-route window is approximately 925–1205 °C. The candidate predicts an approximately 89% high-density grain-size reduction, but this is driven partly by extreme high-temperature reference growth. Roughly 65% of the remaining pore volume is already assigned to the closed store at the switch. These magnitudes are shown, not hidden.

## 3. Experimental-style conventions

All main trajectory figures lead with (T(t)), ρ(t), G(t), and G(ρ). Curves use continuous physical time across the step change. Failure paths remain visible. Filled categorical maps use blue for density exhaustion, green for success, orange for grain-growth failure, purple for mixed failure, and gray for unattainable/ineligible preparation.

Every figure is exported as vector PDF and 600-dpi PNG. Primary panels avoid decoder names and internal completion variables. PR memory, accommodation, and topology closures appear in the mechanism or supplement figures.

## 4. Fast firing

The E0021 material envelope is reconstructed from the archived frozen parameter record at 1, 20, 50, and 100 °C/min. Temperature, density, grain size, G(ρ), and matched-density ratios are shown separately from the two-step result. A bounded heating-rate/peak-temperature/hold map is deterministic presentation mapping, not reoptimization.

Fast firing remains a separate nucleation-limited material-kinetics result: nucleation-facile ablation removes it and PR-off preserves it. Candidate 693168’s topology parameters are not inserted into this material model.

## 5. Two-step trajectories

The main two-step figure compares the 1400 °C reference with lower failure, successful, and upper failure paths. Time is never re-zeroed. The low-temperature isothermal comparator is included in the source history even where it fails to attain the target. Separate figures expose G(ρ), reduction versus density, final density/growth versus T2, and the lower/success/upper triplet.

## 6. Filled Chen maps

Three fixed-candidate maps are generated:

- T1=1200–1550 °C by T2, with 25/10 °C spacing;
- many prepared G1 states from G0, T1, and switch-density routes by T2;
- switch density 0.75–0.94 by T2.

Each map uses exact candidate parameters, cloned first-step states, the common target and budget, and unmodified classifications. Success regions are filled categorical cells rather than sparse success markers. Unattainable and failed points are retained. The T1–T2 plot displays the practical T2<T1 diagonal.

## 7. Mechanism/state figures

Only after the primary trajectories are established do the figures show open, connected, isolated, and closed pore stores; closed accommodation; open/closed removal; PR-prepared memory; and the causal ablation waterfall. The closed-pore fraction and accommodation trajectory remain the key falsification targets.

## 8. Ablation and robustness

The supplement preserves the exact registered ablation definitions and the six-candidate comparison. The initial-state heatmap shows median reduction, window width, and joint acceptance across the bounded ρ0–G0 audit. No failed case is hidden.

## 9. Missing data

The local-region model does not expose calibrated pore D50/D90, gas pressure, absolute PR or densifying energy, or an independent closed-pore Λ/K event model. Pore size is therefore shown only as a clearly labeled volume/number proxy. PR/densification plots use conservative volume-transfer quantities, not invented energies. The frozen archive provides fast-firing histories but not experimental uncertainty.

## 10. Limitations and experimental comparison

The next comparison should measure absolute G(ρ), open/closed porosity, and a 3D connected-pore fraction during interrupted first and second steps. The predicted ~65% closed-store share at density 0.88 and near-closed residual pore population at density 0.98 are stringent falsification targets. High-temperature grain-growth kinetics must also be calibrated before the ~89% reduction magnitude can be interpreted quantitatively.
