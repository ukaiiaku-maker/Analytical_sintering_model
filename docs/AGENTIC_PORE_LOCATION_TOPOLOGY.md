# Local pore-location action layer: bounded mechanism audit

## Scope and status

This is a mechanism prototype and negative-control audit, not a calibration or
validation claim.  It starts from commit `71e404d` on the preserved
`codex/pore-placement-topology-search` baseline and introduces an explicit,
ablatable action layer.  Densification kinetics, density target (0.90), and
the fixed 96 h first- and second-step budgets were not retuned by case.

The central result is negative: the action layer changes pore placement and
grain growth, but no tested action parameter set creates a finite Chen-style
success window at either 5% or 10% allowed growth.  This rejects the current
closure cleanly and narrows the next scientific question.

## Mechanism architecture

`TopologyAction` records a named source, destination, conservation flag,
density/growth role, nonnegative propensity, driving power, resistance, and
diagnostics.  At every integration step, `score_actions` uses only the
instantaneous pore-location state, temperature, and fixed material/action
parameters.  No protocol name, ramp-rate label, or schedule class enters it.

The constrained allocation has three compatible groups:

1. GB-segment versus TJ pore removal shares one densification capacity;
2. GB smoothing versus GB-to-TJ transfer shares one GB redistribution
   capacity;
3. TJ-to-GB-segment capture versus TJ isolation shares one TJ redistribution
   capacity.

All weights and fluxes are nonnegative.  Smoothing and transfers are exactly
volume conservative; only the two removal actions change density.  Grain
growth remains a clean-GB migration action opposed by GB-segment and TJ drag.
The special `TJ_to_GBseg_capture` action is therefore an observable local
topology transition, not an efficiency multiplier.

The modes are:

- `fixed_flux_baseline`: exact delegated evolving-flux baseline;
- `action_static`: exact delegated static placement ablation;
- `action_evolving_no_capture`: constrained action allocation with capture
  prohibited;
- `action_evolving_capture`: capture competes with isolation;
- `action_evolving_capture_high_TJ_drag`: same action model with doubled,
  explicitly named TJ migration resistance.

## Bounded calculation

The 64-set physical screen varied capture/isolation competition, clean-GB
growth weight, TJ densification participation, capture activation energy, and
TJ drag.  Its reduced map used 75--600 nm grains, two switch densities, five
second-step temperatures, and both growth tolerances.  A set was rejected for
universal success, a missing lower/upper failure boundary, or (after auditing
the results) no success window.  The full map retained the five required
ablations and the initially top-ranked diagnostic set, covering 50--600 nm,
three T1 values, three switch densities, T2 = 900--1300 C in 25 C increments,
and 5%/10% tolerances.  This produced 14,688 classifications in 4,975 s.

## Results

All 64 screen sets have zero 5% and zero 10% successes and are now explicitly
listed in `rejected_action_parameter_sets.csv` as `no_success_window`.  The
pre-rejection ranking selected `action_00` only because all scores tied at
zero; it is not a scientific winner.  Its full map is retained to make that
failure reproducible.

Across every full-map mode and tolerance, `SUCCESS` count is zero.  The maps
still contain the desired physical failure types: densification exhaustion,
grain-growth failure, and mixed failure.  Thus the action layer does not make
all temperatures successful, but it also does not make the lower and upper
bounds overlap to form a window.  There is consequently no nanoscale window
and no large-grain window in this bounded map.

At a representative 1300 C isothermal target run, all five ablations reach
0.90.  Final grain sizes are 387.1 nm (fixed baseline), 365.6 nm (static),
384.6 nm (evolving without capture), 382.5 nm (capture), and 339.0 nm
(capture plus high TJ drag).  Capture alone therefore changes final grain
size by only about 0.5% relative to the no-capture action mode, whereas the
named high-TJ-drag ablation changes it by about 12%.  Even that larger
migration response does not create a 5% or 10% window.

The slow-ramp integrated capture flux is approximately 3.13e-4 pore-volume
fraction, versus 1.84e-4 for the fast ramp.  The corresponding isolation
fluxes are 8.40e-4 and 2.47e-4.  GB-segment removal is about 0.122 (slow) and
0.104 (fast).  Capture is therefore observable and schedule-history
dependent through local state, but is two to three orders smaller than the
dominant removal channel.  It shifts the final slow-ramp location fractions
from `(0.3976, 0.3752, 0.2272)` without capture to
`(0.3993, 0.3756, 0.2251)` with capture: real but too small to reorganize the
window topology.

## Interpretation

The failure is not caused by a missing arbitrary drag multiplier.  Doubling
TJ drag substantially suppresses grain growth while leaving densification
kinetics direct and unchanged, yet the density-attainment boundary still
does not overlap the no-growth boundary.  The tested local capture closure is
also too weak relative to pore removal and smoothing to maintain a sufficiently
large connected, removable population during the second step.

The next mechanism should not be a broader search over these same weights.
The useful next hypothesis is a state-resolved *junction population balance*
or explicit grain-boundary segment-length state: capture and junction loss
must alter the number/density of pinning junctions and their release kinetics,
not merely move a small pore-volume fraction among fixed geometric stores.
That state could couple strongly to migration while remaining separate from
GB diffusion.  Experimental interrupted-ramp measurements should resolve
pore occupancy at boundaries versus triple junctions together with junction
density and segment length; pore-size distribution alone cannot identify
this closure.

## Reproducibility and anti-cheat checks

The generated tables and figures are under `results/pore_location_agentic/`.
Tests verify local scorer inputs and source text, nonnegative propensities,
partitioned weights, exact capture-off ablation, nonnegative pore/number bins,
`rho = 1 - sum(phi)` to 1e-12, conservative redistribution, and mutually
exclusive strict classifications.  Failed first steps are never scored as a
second-step success.  No density target or time budget was changed to rescue
a case.

Key files:

- `action_parameter_screen.csv` and `screen_summary.csv`: all 64 parameter
  definitions and scores;
- `rejected_action_parameter_sets.csv`: explicit rejected-set registry;
- `action_flux_histories.csv`: slow/fast action, power, stress, and topology
  histories;
- `chen_style_action_map.csv` and `window_boundaries_by_action_mode.csv`: full
  classifications and boundaries;
- `selected_case_maps.csv`: the tied, rejected `action_00` diagnostic map;
- `action_history_diagnostics.png`, `chen_style_action_modes.png`, and
  `window_width_vs_G0_by_action_mode.png`: deterministic diagnostics.
