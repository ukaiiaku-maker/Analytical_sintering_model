# Separate Fast-Firing and Chen Mechanisms

## Scope

This prototype deliberately separates `MaterialKinetics` from `TopologyGrowthClosure`. The former owns serial nucleation/exchange/transport densification and conservative pore redistribution. The latter returns a migration multiplier only. Neither local law receives a schedule label, ramp rate, target density, or protocol identity.

The 256-set fractional screen is mechanism discovery, not calibration. It used three representative peak/hold pairs (1350 C/2 h, 1450 C/8 h, 1550 C/20 h), a 1 C/min reference, and 20/50/100 C/min fast ramps. This is a bounded subset of the requested schedule space.

## Result

Six of 256 material sets produced at least one strict fast-firing result: `G_reference/G_fast >= 1.5` continuously across at least 0.03 density, with both trajectories attaining 0.75--0.92 and no numerical censoring. The topology audit then found complete lower/upper-bracketed practical windows for all six selected materials; both qTJ=0 and qTJ=1 occurred. Thus the prototype has a nonempty coexistence set with unchanged material kinetics.

This does **not** validate the model. Several extreme grain ratios (maximum about 148) reveal stiff, weakly constrained kinetics. The Chen stage is a reduced fixed-state rate/capacity map rather than a full second-step integration. Both facts require experimental constraints and higher-fidelity verification.

## Architectural conclusion

- Fast-firing separation can arise with topology disabled, through time spent in a surface-redistribution regime before nucleation-limited densification activates.
- Chen-style bounds can arise after fixing material parameters and changing migration-only topology parameters.
- The effects are mathematically compatible in this prototype, but the magnitude and parameter identifiability remain unvalidated.
- Spatial/network physics is still required to interpret connected-pore and junction populations rather than treating them as reduced observables.

## Reproduction

```bash
python separated_fast_chen_search.py --materials 256
python -m pytest -q tests/test_separated_fast_chen.py
```

Outputs are in `results/separate_fast_firing_and_chen_mechanisms/`.
