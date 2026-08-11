# Production PR/De-sintering Assessment

## Frozen confirmation

This production confirmation uses one parameter set only:
`early_memory_mode="PR_plus_connected_fine_attrition"` and
`k_PR_ref_s=2e-4`, layered on the four frozen persistent-junction/TJ-multihit
bases. Targets, kinetics, and 96 h budgets are unchanged. The disabled control
is the previously persisted production negative control and remains
joint-negative.

The full preparation calculation produced 6,488 admissible routes and 6,284
unique instantaneous PR-memory states. Adaptive classification retained lower
densification-exhaustion and upper grain-growth failures, separated kinetic
from practical `T2<T1` maps, and did not score censored boundaries.

| base | complete practical 5%/5% | all complete tolerance combinations | beneficial fast cases |
|---|---:|---:|---:|
| mech_009 | 563 | 4,633 | 1,437 |
| mech_019 | 633 | 5,974 | 655 |
| mech_009_q0 | 106 | 2,537 | 1,518 |
| mech_019_q0 | 643 | 7,037 | 605 |

All four frozen candidates are production joint-positive and non-universal.
They retain thousands of harmful and neutral comparisons and roughly 25,900
unattainable target records per base. Thus the result is not caused by lowering
a target, extending selected paths, or making every fast schedule successful.

## Fast-firing pathway

Benefits occur at every attainable target from 0.85 to 0.92 and in all initial
topologies. They are concentrated at 1400–1500 C, 8–20 h holds, and initial
grain sizes below about 150 nm; 225 nm cases are sparse and 300 nm cases do not
form a meaningful beneficial region. Across beneficial cases, median `HR_pct`
is about 3%, slow-reference minus fast PR work is `1.1e6–1.6e6` in model work
units, connected fine-pore fraction is higher by roughly 0.006–0.007, and mean
connected-pore radius is lower by roughly 0.7–1.0 nm.

Across all attainable points, HR has only weak linear correlation with either
PR exposure difference or fine-pore difference because growth regime, topology,
and target strongly stratify the response. The mechanistic criterion is instead
confirmed conditionally: beneficial matched-density cases consistently show
less PR work and more removable connected fine pores. Interrupted-path figures
make the temporal ordering visible.

## Two-step pathway

The full Chen map remains positive. The representative two-step path retains
the prepared persistent `X_J` population, a distinct `Lambda_TJ/K_TJ` history,
and temperature-sensitive multihit completion. PR memory changes preparation
topology, but the finite lower/upper second-step boundaries continue to arise
from densification exhaustion and migration reactivation—not from a hidden
two-step bonus.

## Local PR robustness

The bounded OAT audit used representative, production-positive Chen and
fast-firing routes for each base. Of 64 variants, 28 remain joint-positive.
The counts by group are: activation energy 7/12, PR rate 4/12, renewal gate
midpoint 6/12, renewal exponent 6/12, and flux partition 5/16. This is a finite
robustness neighborhood, not universal robustness. It also identifies the PR
rate and GB/TJ/isolation partition as important experimental calibration
targets. No target, budget, or base kinetics were changed.

## Pore-occupied versus constrained triple junctions

The original multihit closure used pore-connected TJ coverage directly as the
constraint population. That is a useful limiting case, but it is not assumed
to be unique. The new migration-only ablation distinguishes:

- `C_TJ_pore`: pore-connected TJ occupancy;
- `C_TJ_constraint`: structurally constrained junction population;
- `C_TJ_relaxed`: pore occupancy that relaxes/bypasses compatibility;
- `C_TJ_pinned`: occupancy contributing to Class-A pore drag.

`current_all_TJ_multihit` exactly recovers the prior closure. In the bounded
audit, current-all is joint-positive for 2/4 bases, pore-relaxed for 3/4,
mixed for 3/4, and pore-pinned drag for 4/4. Relaxation generally reduces the
number of beneficial cases; explicit pore drag preserves fewer but still
finite benefits and reports its dissipation separately. The q0 and q1 variants
do not give identical mode rankings. Therefore the best current evidence favors
pore occupancy acting at least partly as a drag site, but does not identify a
unique partition between relaxation and structural constraint.

TJ pore densification remains a separate unchanged channel. `P_TJ_multihit`,
`P_TJ_pore_drag`, and TJ-assisted densification power are reported separately;
none of the TJ constraint modes directly removes pore volume.

## Answers

1. **Full Chen map positive?** Yes, for all four frozen bases.
2. **Full fast map positive and non-universal?** Yes, with beneficial, harmful,
   neutral, and unattainable regions retained.
3. **Same parameters?** Yes: the same moderate PR parameters generate both.
4. **Where?** Primarily small initial grains, 1400–1500 C peaks, and longer
   holds, with topology-dependent frequencies.
5. **Observable memory?** Beneficial fast paths retain connected fine pores,
   smaller connected pores, and lower cumulative PR work at matched density.
6. **Two-step origin?** Persistent junction drag plus TJ multihit reactivation
   still defines migration suppression and its upper boundary.
7. **Necessary ingredients?** Conservative early PR memory, separate
   densification/migration channels, persistent junction state, multihit or
   pore-drag migration resistance, and censor-aware adaptive classification.
8. **Measurements?** Location-resolved interrupted-ramp pore distributions,
   TJ occupancy, junction mobility/back-stress proxies, surface-area loss,
   dilatometry, and matched-density grain size across heating rates. These
   would constrain the PR rate, relocation partition, and whether pore-filled
   TJs relax compatibility or act primarily as pinning sites.

This is a production robustness result, not material validation or a fit to a
specific experiment.
