# Visual Inspection Plots V2

## Why regenerate

Version 2 replaces sparse and placeholder-heavy figures with dense solver-history reconstructions, labeled multi-panel plots, continuous two-step physical time, and an automated output audit. No model physics, material kinetics, topology parameters, or success classifications changed.

## Candidates and histories

All five unique E0142 Tier B windows and the best E0021 Tier C comparison are read from the strict-tier tables. Fast paths were rerun at 1, 20, 50, and 100 C/min with frozen parameters. Histories contain at least 1,200 uniformly interpolated time points plus all solver steps and endpoints. The E0142 triplet carries first-step time into the second step without re-zeroing.

## Visual results

E0021 remains the cleaner nucleation-delay fast-firing material and E0142 the stronger Chen overlap. Initial-condition maps explicitly cover rho0=0.60--0.80 and G0=50--300 nm; unattained comparisons remain rejected. E0142 is Tier B, not Tier A, while E0021 remains Tier C for Chen behavior.

The two-step plots separate temperature/density/grain histories from available topology and pore descriptors. The Chen maps use finite filled bands only and retain distinct density, growth, and mixed failures. Fast-firing and two-step ablations are separate figures with explicit labels.

## Missing data and limitations

The separated model does not expose Pcomp,TJ, CTJ, CGBseg, clean-GB fraction, or decomposed stress/power histories. Those panels are omitted rather than fabricated and are documented in `VISUAL_INSPECTION_PLOTS_V2_MISSING_DATA.md`. These figures support inspection, not new validation.
