# Observable topology gate identifiability audit

## Scope

This is a mechanism prototype and identifiability test, not validation or a new calibration. It keeps the pore-bin redistribution kinetics, thermal protocols, time budgets, initial grain size (150 nm), pore scale (25 nm), and pore width (0.65) fixed. Only `rho0` and `smoothing_gate_mode` vary. PR #2 remains a draft.

The bounded design contains four gates, five initial densities, and four target densities. The canonical comparisons are 0.2 versus 20 C/min and 1350 C isothermal versus 1350 -> 1250 C with switching at rho=0.85. A result is scored only when both paths reach the requested density.

## Gate definitions

- `density`: the previous gate, unchanged: `sigmoid((smoothing_rho_mid-rho)/smoothing_rho_width)`.
- `fine_pore`: a sigmoid of instantaneous removable fine-pore fraction.
- `connectivity`: a sigmoid of instantaneous connected pore-boundary coverage, `f_pore * connectivity`.
- `hybrid_topology`: the product of the fine-pore, connected-coverage, and non-isolation sigmoid gates.
- `none` is available as an ablation (unit gate), but was not added to the requested four-mode comparison.

All topology gates are computed from the instantaneous pore bins and inferred topology. They receive no protocol object, schedule label, ramp rate, or slow/fast flag. The recorded `smoothing_gate_value` makes their action observable and ablatable.

The fixed prototype thresholds are not fitted by gate mode: fine midpoint/width = 0.45/0.08, connected-coverage midpoint/width = 0.45/0.06, and isolation midpoint/width = 0.08/0.03. These are transparent provisional scales, not identified material constants.

## Main results at rho=0.90

| gate | HR_pct, rho0=0.65 | 0.70 | 0.75 | 0.80 | 0.85 | TS range among attained cases |
|---|---:|---:|---:|---:|---:|---:|
| density | +4.58 | +5.01 | +4.95 | -1.11 | -6.95 | 11.08 to 11.46 |
| fine pore | +1.10 | +1.30 | +1.48 | +1.54 | +1.22 | 10.99 to 11.01 (2/5 attained) |
| connectivity | +3.80 | +3.35 | +1.91 | -1.51 | -6.58 | 11.08 to 11.46 |
| hybrid topology | +0.62 | +0.43 | -0.38 | -2.78 | -6.90 | 11.08 to 11.46 |

At rho=0.85, only four initial states are eligible; all eligible comparisons attain the target. At rho=0.88 all 20 conditions attain both comparisons. At rho=0.90, heating-rate comparisons attain in all conditions; density, connectivity, and hybrid modes attain all two-step comparisons, whereas the fine-pore mode attains only 2/5. No mode attains rho=0.92 for both paths, so no 0.92 metric is scored.

## Interpretation

### 1. Can topology reproduce the density-gate crossover without using density?

Yes, connected coverage reproduces it qualitatively. It gives positive HR for rho0=0.65--0.75 and negative HR for rho0=0.80--0.85, close to the density baseline. At rho=0.90, HR correlates with initial connected coverage at 0.989 in connectivity mode, compared with -0.938 versus rho0. The hybrid also retains a crossover, but shifts it downward between rho0=0.70 and 0.75 because multiplying three gates suppresses redistribution more strongly.

Fine-pore inventory alone does not reproduce the crossover. Every initial condition starts with the same normalized removable fine-pore fraction (0.5724), because this audit varies pore volume through rho0 while holding the normalized pore distribution fixed. Consequently its initial gate is 0.8220 for all five states. It predicts a small positive HR benefit everywhere and makes three rho=0.90 two-step targets unattainable. This is a useful falsification: normalized fine fraction alone lacks the extensive/topological information needed to distinguish predensification.

### 2. What controls the crossover now?

Within this deliberately restricted design, connected pore-boundary coverage is the best observable predictor. Initial connected coverage decreases from 0.636 at rho0=0.65 to 0.257 at rho0=0.85. The connectivity gate value correspondingly falls from 0.957 to 0.039, cumulative redistribution falls from 0.224 to 0.005, and HR changes from +3.80% to -6.58%.

This does not prove that coverage is independent of density. In the present topology closure, both `f_pore` and connectivity are inferred partly from density, and fixed normalized pore bins make them strongly collinear with rho0. The result replaces a direct density cutoff with an observable topology pathway, but the current initial-condition family cannot statistically separate the two explanations.

### 3. Why does the predensified case fail?

For connectivity and hybrid gates, rho0=0.85 begins with low connected coverage (0.257), appreciable isolation (0.0585), and a gate value of only 0.039 or 0.021. It therefore accumulates almost no redistribution before densification. The failure is not caused by depleted normalized fine-pore fraction, which is identical across rho0. It is caused by depleted connected/removable boundary topology in this closure.

The fine-pore-only result confirms this distinction: allowing the same fine gate at rho0=0.85 removes the negative HR sign, but it also harms target attainment for the canonical two-step path. Treating a normalized fine fraction as sufficient therefore overstates available removable pore volume.

### 4. Is two-step behavior sensitive to the gate?

For density, connectivity, and hybrid modes, TS at rho=0.90 remains in the narrow 11.08--11.46% range. The gate changes HR much more than TS. This is consistent with the earlier initial-condition audit, which associated the two-step response more strongly with grain-growth state/G0. G0 is held fixed here, so this audit tests gate sensitivity rather than re-estimating the G0 effect.

The fine-pore-only gate is the exception: three two-step cases fail to reach rho=0.90. Those failures are reported, not scored, and show that this gate is not an acceptable standalone replacement under the fixed parameters.

### 5. Does topology reduce non-identifiability?

It improves physical observability but does not yet resolve parameter identifiability. Connected coverage can be measured from 3D pore/GB topology, whereas a fitted density midpoint is only a proxy. However, density and connected coverage remain confounded in this audit because pore shape/distribution is fixed at each rho0 and the closure computes connectivity from density. The hybrid adds thresholds and therefore worsens parameter identifiability unless independent isolation data are supplied.

A decisive calibration must vary topology at matched density: samples with equal rho0 but different fine-pore volume, GB pore coverage, and isolation. If their interrupted-ramp response follows connected coverage rather than density, the topology gate has earned its extra structure.

### 6. Measurements needed

The most useful experiment is an interrupted slow-ramp series on green bodies prepared at matched density but different pore topology, paired with a fast-ramp control. Quench around 800, 950, 1050, and 1200 C and measure:

1. full 3D pore-size distribution and fine-pore *volume per specimen volume*, not only normalized fraction;
2. fraction of grain-boundary area covered by connected pores;
3. connected-cluster and isolated-pore volume fractions;
4. shrinkage/density and grain size at each interruption.

FIB-SEM or diffraction-contrast tomography with segmentation tied to stereological GB coverage would directly constrain the connectivity midpoint and width. Gas adsorption or SAXS can supplement the fine-pore inventory, but cannot alone identify whether pores remain boundary-connected.

## Model decision

`density` remains the default and is exactly preserved. `connectivity` is the strongest observable competitor and should be carried forward as the next experimental hypothesis. `fine_pore` alone fails the crossover discrimination and target-attainment check. `hybrid_topology` is mechanistically plausible but not justified as the default because it adds non-identifiable thresholds and suppresses the benefit too aggressively.

The next scientific step is not a wider parameter search or late-stage physics. It is a matched-density, varied-topology calibration/validation dataset. No result here resolves the rho>0.92 limitation, and no validation is claimed.

## Reproduction

```bash
python3 topology_gate_identifiability.py
python3 -m pytest -q
```

Outputs are under `results/topology_gate_identifiability/`.
