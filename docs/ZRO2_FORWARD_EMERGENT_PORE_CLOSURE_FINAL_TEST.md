# ZrO2 emergent pore-closure final promotion test

## Outcome

`emergent_pore_closure_v1` is a **diagnostic negative result**, not a promoted
production law. The candidate passes its limiting and dimensional checks, but the
actual naturally prepared first-step state has no success points, no finite window,
and no correctly ordered near-window. The inherited reconstructed energy ledger
also contains unresolved budget violations. A process map was therefore not
justified and was not run.

No accepted forward-model physics changed. The barrier JSON, GB and surface
diffusivities, intrinsic mobility, open-pore renewal branch, schedules, density
identity, PR conservation, and classification thresholds are unchanged. The failed
global mobility calibration remains inactive. No physical `Q_closed` was introduced;
`Q_closed_app` is calculated only as a post-run slope.

## Tested structure

The primary candidate multiplies closed inventory, shrinkability/connectivity,
finite accommodation, a pore-size penalty, capillary stress after gas
counterpressure, and a stress-resolved renewal/transport kernel. Radius exponents
3 and 4 were tested. A unit-correct GB-diffusion alternative and a strictly
non-densifying surface-accommodation mode were tested as diagnostics.

The actual exported state contains only (6.59\times10^{-6}) closed pore volume and
accommodation (2.82\times10^{-3}). This inventory–availability product prevents
the density target from being reached. A moderate bracket and a high injected state
produce two low-temperature GB-diffusion success points, but success begins at the
lowest scanned temperature, so there is no lower density-exhaustion boundary. Those
states therefore do not establish a finite Chen window.

All tested temperatures are below the fitted barrier range and are flagged as
extrapolative. This analytical audit is not validation, and no validation claim is
made.
