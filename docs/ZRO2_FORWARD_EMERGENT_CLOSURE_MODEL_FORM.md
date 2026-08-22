# Emergent pore-closure model form

For closed bin (i), the candidate state is
\(\{\phi_{c,i},N_{c,i},r_{c,i},\chi_{c,i},C_{geom,i},P_{gas,i},A_{c,i}\}\),
with available, used, and recovered accommodation recorded separately. Density is

\[
\rho=1-\sum_i(\phi_{open,i}+\phi_{iso,i}+\phi_{closed,i}).
\]

Only named open and closed shrinkage fluxes alter this identity. PR/coarsening,
precursor formation, closure transfer, surface smoothing, and migration are
non-densifying.

The driving stress and renewal kernel are

\[
\sigma_{c,i}=\max(C_{geom,i}2\gamma_s/r_{c,i}-P_{gas,i},0),
\]
\[
r_{nuc,i}=\nu_0e^{-G^*(\sigma_{c,i},T)/(k_BT)},\qquad
\tau_{sink,i}=C_c\frac{k_BT}{\sigma_{c,i}\Omega}\frac{\ell_{c,i}^2}{D_{GB}},
\]
\[
\dot\rho_{c,i}=\phi_{c,i}\chi_{c,i}A_{c,i}
\left(\frac{r_{ref}}{r_{c,i}}\right)^m
\frac{1}{\tau_{sink,i}}\frac{r_{nuc,i}\tau_{sink,i}}
{1+r_{nuc,i}\tau_{sink,i}}.
\]

The GB-diffusion diagnostic uses

\[
\dot\rho_{c,i}=C_{GB,c}\phi_{c,i}\chi_{c,i}A_{c,i}D_{GB}
\frac{\Omega\sigma_{c,i}}{k_BT}r_{c,i}^{-m}.
\]

Dimensional closure requires (C_{GB,c}) in m\(^{m-2}\). The audit uses
\(r_{ref}^{m-2}\); this magnitude is semi-phenomenological and is not a material
activation energy.

Surface diffusion supplies only bounded shape/accommodation recovery. Closed
shrinkage consumes accommodation, while recovery cannot raise it above its maximum.
Gas pressure is an explicit counterforce. The strict GB-area-loss equality remains
diagnostic only; the energy ledger never silently forces all released energy into
densification.

This model form is not validation. No validation claim is made.
