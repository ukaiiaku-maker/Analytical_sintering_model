# Emergent closure test results

## Analytical gates

All limiting checks pass for renewal, GB-diffusion, and surface-only modes at
(m=3,4): rates vanish with zero inventory or nonpositive stress, decrease with
radius, and fall with gas counterpressure. Surface accommodation alone has zero
density flux. Conservative preparation transfers preserve pore volume,
accommodation stays in [0,1], and the stored-state density identity closes.

All unit-audit entries pass. That result does not make every coefficient physical:
the GB coefficient has the required dimensions but an underived magnitude, and the
shape recovery and size normalizations remain bounded or semi-phenomenological.

## Boundary and promotion gates

The actual selected state produces only density-exhaustion and mixed failures for
both kernels and both exponents. It has zero success points and zero finite windows.
The moderate and injected states yield two GB-diffusion success points at 850 and
875 °C, followed by growth failure. Because success begins at the bottom of the
scan, the required lower boundary is absent. Renewal variants do not succeed.

Consequently no naturally prepared lower/success/upper topology exists, no candidate
is promoted, no broad process map is run, and the strict finite-window count is zero.
The requested process-map panels are failure-mode diagnostics, not success maps.

## Ablations and ledger

Removing closed shrinkage or using surface accommodation alone eliminates density
gain. Removing preparation or precursor support reduces gain. High gas pressure
strongly suppresses shrinkage, while infinite accommodation increases it. The
GB-diffusion alternative exhausts the bracketed closed inventory and thereby exposes
its overactive/no-lower-boundary behavior.

The selected-path energy ledger retains named open/closed work, PR, smoothing,
coarsening, growth, drag, gas, other, and residual channels. Reconstructed histories
still contain budget violations, so ledger consistency does not pass. Rates are not
silently rescaled to hide this result.

Every candidate temperature lies below the fitted barrier range. These findings are
not validation, and no validation claim is made.
