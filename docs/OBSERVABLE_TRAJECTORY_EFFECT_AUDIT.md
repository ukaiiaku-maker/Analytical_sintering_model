# Observable Trajectory Effect Audit

## Scope and decision rule

This is a correction to the interpretation of the production PR/de-sintering
campaign, not a new mechanism or calibration. The audit uses matched-density
interpolation and regards a response as experimentally meaningful only when
the reference/comparison grain-size ratio is at least 1.5 continuously over
`Delta rho >= 0.03`. Ratios below 1.2 are negligible; 1.2--1.5 are weak; and
ratios at least 2 are strong. An isolated target point is never promoted to a
trajectory result.

The committed fast-firing table supplies targets 0.85, 0.88, 0.90, and 0.92.
Its explicit `reference_rate_C_min` field keeps a reachable 0.2 C/min reference
separate from the 1 C/min fallback. The representative histories supply a
continuous diagnostic comparison. Densities above 0.92 were not sampled by
the production campaign, and densities at or above 0.95 are additionally
outside the branch's physically supported open-pore scope.

## Result

The current branch is mechanistically useful but does not yet produce an
experimentally meaningful sintering-trajectory separation.

Across 12,347 attained production fast-firing comparisons, no sampled point
below density 0.92 reaches a grain-size ratio of 1.5. The largest sampled ratio
is 1.301 (23.1% improvement), which remains in the weak effect-size band even
before applying the finite-density-span rule. The representative 1 versus 20
C/min trajectory has maximum/median ratios of 1.012/0.998 in the fully attained
early window (0.75--0.85). Its intermediate window is not jointly attained and
therefore is not scored.

The representative high-temperature/two-step comparison reaches a maximum
ratio of 1.372 below 0.92. It therefore does not meet the 1.5 criterion either.
The fully attained early-window maximum/median are 1.028/1.001. This does not
invalidate the finite Chen-style kinetic windows; it says the selected
representative trajectory does not yet demonstrate a paper-scale grain-size
separation under the new criterion.

## Bounded PR rescue screen

Because the baseline failed, a 24-point OAT/partition screen varied only the
named PR calibration targets: rate prefactor, activation energy, renewal gate,
renewal power, conservative flux partition, and topology power. Targets,
96-hour budgets, schedules, densification physics, and Chen definitions were
unchanged. No candidate achieved a 1.5 ratio over `Delta rho >= 0.03` while
jointly attaining 0.85--0.92. The largest observed ratio was 1.373 for
`k_PR_ref_s=2e-3`, but the full interval was not jointly attained; it is
rejected rather than rescued. Since no trajectory candidate passed the first
screen, no candidate was promoted to an expensive Chen recheck.

## Required questions

1. **Meaningfully different trajectories?** No.
2. **Maximum and median over 0.85--0.92?** The exhaustive discrete production
   table has a global maximum of 1.301. Medians are reported per schedule in
   `density_window_effects.csv`; aggregating unlike schedules into one median
   would conceal reference and attainment differences. The representative
   comparison does not attain the entire interval and has no scored median.
3. **Any fast-firing ratio >= 1.5 below 0.92?** No.
4. **Any representative two-step ratio >= 1.5 below 0.92?** No; maximum 1.372.
5. **Largest differences only above 0.95?** No such conclusion can be drawn:
   the production campaign stops at 0.92.
6. **Are high-density differences supported?** No. Targets 0.95, 0.98, and
   0.99 are explicitly marked unsampled and unsupported without closed-pore
   physics.
7. **Likely experimentally measurable?** The representative separation is
   negligible. Even the production maximum is only weak under the stipulated
   uncertainty bands, so the branch does not support that claim.

## Files and interpretation

`fast_firing_ratio_curves.csv` preserves each schedule and reference choice;
`two_step_ratio_curves.csv` contains the representative matched interpolation;
`density_window_effects.csv` records attainment before scoring;
`high_density_attainment.csv` prevents extrapolation; and the rescue tables
record every accepted/rejected bounded trial. Figures use common, non-magnified
axes and threshold lines. Internal pore-memory diagnostics are retained only
to show that a mechanistic state difference need not create an observable
`G(rho)` separation.

This is a negative/insufficient-mechanism result, not model validation.

