# ZrO2 forward distributional pore-population model

## Scope

The previous archive-only audit was not sufficient to test the distributional hypothesis. This branch therefore implements actual forward evolution of topology-labelled pore populations. It is a bounded model-form test, not calibration or validation. Candidate 693168 is a response-form target only and supplies no ZrO2 parameter values.

The state has 16 logarithmic radius bins and four explicit stores: `open_connected`, `precursor`, `isolated`, and `closed`. Each bin retains pore volume, number proxy, radius, surface-area proxy, topology age, and origin. Closed bins additionally retain shrinkability, accommodation capacity/availability/use/recovery, gas pressure, and driving stress.

The exact identity is

`rho = 1 - sum_alpha sum_i phi_alpha,i`.

Regularization, coarsening, PR pinch-off, precursor formation, isolation, and closure are conservative. Only named open and closed shrinkage remove pore volume. Surface accommodation alone cannot densify.

## Representations

- `lognormal`: conservative bin evolution followed by single-lognormal moment projection.
- `bimodal`: conservative bin evolution followed by two-mode projection.
- `discrete_bin`: full bin population retained without projection.

Initial distributions are generated from the prescribed lognormal, bimodal, or controlled discrete families. Synthetic states are tagged `synthetic_forward_state`; the candidate-like state is tagged `candidate_response_target_only`. Reconstructed positive distributions are never labeled measured.

## Fixed physical inputs

The barrier JSON and material laws remain unchanged:

- `D_GB = 0.056 exp[-380000/(RT)] m2/s`
- `D_s = 0.10 exp[-380000/(RT)] m2/s`
- accepted intrinsic mobility branch unchanged and failed global fit inactive
- existing renewal and GB-diffusion closed kernels, tested with `m=3,4`

No physical `Q_closed` was introduced. `Q_closed_app` remains a post-run diagnostic only. Geometry and gate constants are classified in `distribution_parameter_registry.csv`.

No validation claim is made.
