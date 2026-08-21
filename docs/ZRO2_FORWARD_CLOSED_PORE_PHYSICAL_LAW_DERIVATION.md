# ZrO2 closed-pore physical-law derivation

## Bin state and force

Each candidate uses closed inventory, number/radius, shrinkability, geometry,
trapped-gas pressure, and available/used/recovered accommodation. The density
identity remains

\[
\rho=1-\sum_i(\phi_{o,i}+\phi_{iso,i}+\phi_{c,i}),
\]

and the bin force is

\[
\sigma_{c,i}=\max(C_{geom,i}2\gamma_s/r_{c,i}-P_{gas,i},0).
\]

Only named open or closed shrinkage removes pore volume. PR, precursor formation,
closure transfer, surface smoothing, pinning, and migration do not directly change
density.

## Candidate laws

`renewal_limited_closed_shrinkage` reuses the fitted stress-resolved barrier,
GB transport time, surface-exchange time, and cycle activity. It tests radius
exponents 3 and 4. Its event-scale normalization remains semi-phenomenological.

`GB_diffusion_closed_shrinkage` uses

\[
\dot\rho_c=C_{GB,c}\phi_c\chi A_c
 D_{GB}\Omega\sigma_c(k_BT)^{-1}r_c^{-m}.
\]

Dimensional analysis requires (C_{GB,c}) to have units m\(^{m-2}\). The audit
uses (r_{ref}^{m-2}) only to bracket rates. That factor is a dimensional geometry
calibration target, not an activation energy or fitted ZrO2 property.

`surface_diffusion_accommodation_only` changes pore shape/accommodation through a
surface-diffusion scale and has identically zero density rate.

`gas_limited_closed_shrinkage` applies explicit idealized counterpressure to the GB
transport law. It stops when gas pressure reaches capillary pressure.

`empirical_reduced_closure` is retained only as a labeled diagnostic comparator.
Its `Q_closed_emp` is not a physical ZrO2 material property and is not compared as
though it were the earlier reduced-model `Q_closed`.

## Preparation and growth

The proposed PR flux scales with (D_s r^{-4}), low renewal activity, topology,
and available open inventory. Transfers are conservative. A physical coefficient
and topology mapping remain unresolved. Intrinsic growth stays
\(\dot G_{int}=M_{GB}\gamma_{GB}/G\); observed growth is multiplied by pore, Zener,
junction, preparation, and accommodation factors. The Zener scale is
\(R_Z=C_Zr_p/f_v\), preferably evaluated from bin moments. These modifiers alter
migration only.

All candidate temperatures (850–1350 °C) are below the lowest fitted barrier slice
(1557 °C). Barrier-shape parameters are therefore clamped to the nearest fitted
slice while Arrhenius temperature dependence remains active. This extrapolation is
an important unresolved limitation.

This derivation is **not validation**. No validation claim is made.
