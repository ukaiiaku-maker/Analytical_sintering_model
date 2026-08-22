# Distributional evolution laws

## Surface time and conservative evolution

The surface coefficient is

`B_s(T) = D_s(T) gamma_s Omega^(4/3)/(k_B T)`,

and `tau_s,i = C_s r_i^4/B_s`. The geometry coefficient `C_s=1` is geometry-derived, not an activation-energy fit.

Regularization moves volume toward the next smaller bin. Damaging coarsening moves it toward the next larger bin. Exponential bin operators are used so transfers remain nonnegative, conservative, and dependent on physical time. Neither branch directly changes density.

## PR and topology transfer

`I_PR = lambda_seg/(2 pi r F_GB) - 1` and `P_pinch = logistic(I_PR/w_I)`. Pinched volume is partitioned conservatively as 0.50 precursor, 0.30 isolated, and 0.20 closed. Precursor and isolated stores may subsequently transfer to closed bins without changing total pore volume.

## Named shrinkage

Open shrinkage uses the fixed stress-dependent barrier renewal construction binwise. Closed shrinkage uses the existing renewal-limited or GB-diffusion-limited candidate:

`sigma_c,i = max(2 gamma_s/r_i - P_gas,i, 0)`.

The closed flux contains closed inventory, shrinkability, finite accommodation, the selected `r^-m` geometry, and the fixed physical kernel. It vanishes for zero closed inventory or nonpositive stress. Surface diffusion may recover accommodation but has no standalone density flux.

## Distributional Zener coupling

The distributional metric is `sum_i phi_open,i/r_i`; the mean-radius comparator is reported separately. Zener and migration modifiers enter grain growth only and never alter open or closed density rates.

All local constitutive functions omit schedule and protocol labels. No validation claim is made.
