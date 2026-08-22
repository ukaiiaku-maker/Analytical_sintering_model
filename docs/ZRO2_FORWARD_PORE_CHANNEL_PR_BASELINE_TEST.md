# ZrO2 pore-channel / PR baseline test

## Decision

The three literature-guided topology modes are a `diagnostic_negative_result`.
They pass conservation and limiting checks, but none produces a naturally prepared
lower/success/upper two-step topology. No bounded process map was justified or run.

The modes are `surface_coarsening_only_v1`, `PR_pinchoff_v1`, and
`PR_regularization_damage_v1`. They are implemented in a standalone audit wrapper;
accepted forward-model physics is unchanged. The barrier JSON, GB and surface
diffusivities, open renewal law, intrinsic mobility, density identity, schedules,
and success definitions remain fixed. The failed global mobility fit is inactive.
No material parameters were fitted and no physical `Q_closed` was introduced.
`Q_closed_app` remains post-run diagnostic only.

All topology operations are conservative. Surface diffusion, coarsening,
regularization, damaging PR, pinch-off, precursor formation, isolation, and closure
do not directly densify. Only named open and closed shrinkage fluxes alter density.

The fixed-state scan contains 316,800 combinations. The (r^4) surface time makes
the strongest topology response nanoscale; large-particle cases are progressively
slower. This scale sensitivity does not rescue density attainment. Heating-rate
families fail to jointly attain the prescribed matched-density targets, so no
matched-density grain-size ratio is claimed.

This is a bounded physics test, not validation. No validation claim is made.
