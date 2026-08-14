# TJ Pore and Structural-Constraint Interpretation

## Scientific status

The Class-B multihit closure represents compatibility-limited migration of
structurally constrained triple junctions during shear-coupled grain-boundary
motion. Treating a pore-filled TJ as a relaxed accommodation site is a modeling
hypothesis: the source literature motivates structural TJ back stress and
disconnection compatibility, but does not calibrate pore-surface relaxation.

The implementation therefore distinguishes `C_TJ_total`, `C_TJ_pore`,
`C_TJ_structural`, `C_TJ_constraint`, `C_TJ_relaxed`, and `C_TJ_pinned`.
Only `C_TJ_constraint` enters `Lambda_TJ_structural`. Pore drag enters a
separate Class-A resistance and power channel. TJ-assisted densification is
unchanged and cannot be driven by either migration channel.

## Frozen bounded audit

No PR, densification, pore-redistribution, target, schedule, adaptive-boundary,
or 96 h budget parameter was changed. The four modes were tested on both q0 and
q1 frozen bases:

| mode | joint-positive bases | complete canonical Chen bases | beneficial fast cases |
|---|---:|---:|---:|
| current_all_TJ_multihit | 2/4 | 2/4 | 75 |
| pore_relaxed_TJ | 3/4 | 3/4 | 72 |
| pore_pinned_drag | 4/4 | 4/4 | 70 |
| mixed_relaxed_pinned | 4/4 | 4/4 | 64 |

All modes retain harmful and neutral fast-firing regions. The old limiting case
is therefore not uniquely required. Pure relaxation fails only for the
`mech_009_q0` canonical Chen route; pinning or a mixed partition restores that
q0 route. Both q variants remain positive for the `mech_019` family under all
four interpretations.

Across the focused attained points, no occupancy scalar alone is a strong
success predictor. The strongest simple association is lower
`Lambda_TJ/K_TJ` (correlation about -0.43 with the beneficial indicator), which
is consistent with suppressed migration. `C_TJ_pore`, `C_TJ_constraint`,
`C_TJ_pinned`, and pore-drag power have near-zero marginal correlations because
temperature, topology, q, and base mechanism condition their effects.

## Answers

1. The joint result does not require pore-filled TJs to impose structural
   multihit constraints.
2. Pore relaxation remains joint-positive for three frozen bases.
3. Explicit pore drag is more robust in this bounded audit: all four bases
   remain joint-positive, as they do for the mixed partition.
4. q0 is more sensitive in the `mech_009` family; q0 and q1 behave similarly in
   the `mech_019` family.
5. `Lambda_TJ/K_TJ` is the most informative single diagnostic here, but it is
   not an independently calibrated predictor. Conditional histories and the
   separate pore-drag channel are more defensible than raw occupancy.
6. The safest interpretation is that multihit kinetics describe structurally
   constrained shear-coupled TJ migration. Pore-filled TJs can relax part of
   that compatibility requirement and/or impose Class-A capillary drag. The
   model does not claim that all pore-filled TJs are structural constraints.

Location-resolved TJ mobility, pore occupancy, back-stress relaxation, and
interrupted-sintering measurements are needed to calibrate this partition.
