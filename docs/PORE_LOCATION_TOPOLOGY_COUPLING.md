# Pore-location topology coupling

## Status and outcome

This branch starts from PR #2 head and leaves that PR as the aggregate-model
negative control. The new model explicitly resolves pore volume and pore number
by size bin and location: GB segments, triple junctions, and isolated/closed
pores. It is a successful architecture and mechanism-ladder implementation,
but it does **not** produce a robust 5--10% Chen-style nanoscale window.

The strongest nanoscale observation is an isolated 10% success at 300 nm for
the static GBseg-rich ablation. Its sampled window width is zero. No evolving
case produces any success below 450 nm, so this is not validation.

## Architecture

`pore_location_mode="disabled"` delegates directly to the PR #2 solver and is
bit-for-bit identical. `static` enables location-specific eligibility, stress,
pinning, and power without relocation. `evolving` additionally enables:

- conservative GB-segment smoothing `GBseg_i -> GBseg_i+1`;
- conservative `GBseg_i -> TJ_i` relocation;
- conservative `TJ_i -> isolated_i` conversion;
- conservative TJ smoothing `TJ_i -> TJ_i+1`.

Only explicit GB-segment and TJ densification fluxes reduce pore volume.
Isolated pores have no open-pore removal flux. At every step,

```text
rho = 1 - sum_i(phi_GBseg_i + phi_TJ_i + phi_iso_i).
```

Coverage and occupancy use bounded exponential saturation. Intermediate
projected area, line density, reduced GB area, and reduced TJ line scale are
reported. GB-segment and TJ pinning have separate thermal and grain-size
scales. The unsuppressed migration propensity defines clean-GB and drag power;
mobility suppression therefore does not directly alter densification.

Named power channels are `P_GBseg_dens`, `P_TJ_dens`, `P_clean_GB`,
`P_GBseg_drag`, `P_TJ_drag`, `P_GB_to_TJ_relocation`,
`P_TJ_iso_conversion`, `P_iso`, and `P_persistent_drag`. Named stresses are
reported separately for GB-segment pores, TJ pores, clean GBs, isolated pores,
and total activation.

All local flux functions depend only on instantaneous state, temperature, and
material parameters. They contain no protocol, ramp-rate, target, or schedule
labels.

## Level 0: static placement ablation

All cases start at the same aggregate density and aggregate pore-size
distribution and use a 1300 C isothermal path with a fixed 96 h budget.

| location class | target 0.90 | final rho | G final (nm) | time (h) | median E_G |
|---|---|---:|---:|---:|---:|
| GBseg-rich | reached | 0.9000 | 213 | 15.1 | 0.383 |
| mixed GBseg/TJ | reached | 0.9001 | 325 | 48.7 | 0.178 |
| TJ-rich | reached | 0.9000 | 460 | 77.9 | 0.127 |
| isolated-rich | failed | 0.8044 | 954 | 96 | 0.032 |
| clean-GB-rich | failed | 0.7941 | 906 | 96 | 0.037 |

This gives the required physical ordering: GB-segment placement couples most
efficiently to pore removal; TJ placement supplies stress and junction drag but
has lower direct removal efficiency; isolated pores exhaust open-pore
densification; pore-poor clean GB area coarsens without equivalent removal.

## Level 1: conservative evolution and history

At the end of the canonical ramps, the evolving slow path reaches `rho=0.864`
with location fractions GBseg/TJ/isolated = `0.424/0.351/0.226`; the evolving
fast path reaches `rho=0.835` with `0.456/0.361/0.183`. Thus instantaneous
location fluxes create a measurable history difference, particularly more
isolation during the slow exposure.

However, evolving and static final fractions differ only modestly at a matched
protocol. Selective densification and density-driven isolation dominate over
the present relocation rates. Placement memory is observable but weak.

## Bounded parameter screen

The 64-case deterministic fractional screen varies only named quantities:
initial location fractions, location size bias, conservative smoothing,
relocation and isolation rates, TJ densification efficiency, size exponents,
pinning strength, and stress concentration. It uses a fixed reduced map at
100/225 nm and fixed targets/budgets.

No screen case succeeds at 5% or 10% on the reduced map. Thirty-seven cases
are rejected: 31 lose the upper grain-growth boundary and six lose the lower
densification boundary. There is no universal-success case. The two retained
screen cases selected for full maps are negative controls with both boundaries,
not fitted winners.

## Full Chen-style maps

Four cases use the requested `G0=50...600 nm`, `T1=1250/1300/1350 C`,
switch density `0.75/0.80/0.85`, `T2=900...1300 C` in 25 C increments,
`rho_target=0.90`, 5%/10% tolerances, and uniform 96 h step budgets.

| case | 5% successes | 10% successes | finite nanoscale window |
|---|---:|---:|---|
| default evolving | 0 | 0 | none |
| screen_00 evolving | 0 | 0 | none |
| screen_04 evolving | 9 | 27 | none; finite only 450--600 nm |
| static GBseg-rich | 3 | 26 | none; isolated point at 300 nm |

Every selected map retains both densification-exhaustion and grain-growth
failure regions. Failed first steps and targets are retained and never scored.

## Why the nanoscale window still fails

Explicit placement solves an interpretability defect but not the size-scaling
defect. GBseg-rich states densify efficiently, yet their pore-bearing coverage
depletes as those pores are removed. Clean GB migration then dominates. TJ-rich
states retain drag and stress but densify too inefficiently and accumulate too
much growth before reaching the target. Evolving states additionally convert
TJ volume to isolated volume, exhausting the open-pore channel. The current
relocation flux is too weak to replenish GB-segment pore-bearing area, while
increasing it or pinning constants after seeing the map would be parameter
forcing.

The next discriminating physics is not another scalar drag increase. It is an
experimentally constrained rule for whether migrating/bowing boundaries sweep
pores, detach from them, or continually repopulate pore-bearing segments, plus
eventual closed-pore removal for densities above the open-pore regime.

## Reproduction

```bash
python3 pore_location_topology_sensitivity.py
python3 -m pytest -q
python3 -m py_compile *.py
```

All requested CSV tables and diagnostic plots are in
`results/pore_location_topology/`.
