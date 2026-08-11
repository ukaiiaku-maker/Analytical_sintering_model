# Rejected-Case Failure Decomposition

## Scope

Only one bounded-search case ever exceeds a mean-grain ratio of 1.5: combined
case C091. It is not a success. This audit reruns its reference and fast paths
without changing parameters, schedules, targets, or budgets and resolves time,
temperature, mean/median/tail grain sizes, pore D50/D90, pore topology, PR work,
residual stresses, densification rate, and grain-growth rate.

## Why the interval terminates

C091 reaches a maximum `G_mean,reference/G_mean,fast` ratio of 1.662, but the
ratio exceeds 1.5 only from density 0.850 to 0.858 (`Delta rho=0.008`). The
reference path ends at density 0.859 and the fast path at 0.880 under their
unchanged schedules. Thus the high ratio does not collapse within a shared
density interval; it becomes unscorable because the reference path no longer
attains higher density. Its final densification rate remains positive, so this
is fixed-schedule target nonattainment rather than proof of a thermodynamic
densification-exhaustion point. Selectively extending the schedule would
violate the audit rules and was not done.

At density 0.85 the reference/fast pore D90 values are 450/218 nm, connected
fine-pore fractions are 0.024/0.057, large-pore fractions are 0.979/0.935, and
cumulative PR works are approximately 503,000/29,400 model units. The pore
memory is real. It strongly changes the brief grain trajectory, but it does
not supply a jointly attained `Delta rho >= 0.03` interval.

## Failure questions

1. **Does the ratio collapse because the slow/reference path stops attaining
   density?** The ratio does not collapse; scoring stops because the reference
   reaches only 0.859 under the fixed schedule.
2. **Does the fast path later coarsen?** It does coarsen, but there is no
   reference trajectory at the same later densities, so this cannot explain a
   matched-density ratio collapse.
3. **Does large-pore damage affect pores but not growth?** It affects both:
   C091 has large pore and mean/median/tail grain differences. The limitation
   is duration and joint attainment, not complete decoupling from growth.
4. **Does residual stress relax too quickly?** Yes. Large-pore stress falls to
   roughly `8e-13` of its peak on the reference path, while PR/pore memory
   remains. The current stress state cannot carry the early disparity forward.
5. **Is missing late-stage closed-pore persistence the cause?** No evidence
   supports that diagnosis: failure occurs at densities 0.859--0.880, well
   below 0.95.

## Decision

Select **persistent defect/topology-stress memory below density 0.92**. The
late-stage closed/isolated-pore path is deferred because no qualifying or
terminating behavior was observed near density 0.95, and the present model
must not extrapolate there. The next mechanism should store PR-created defect
population and stress/topology history conservatively, with slow relaxation,
then test whether that state changes mean-grain trajectories over a fully
attained finite interval.

