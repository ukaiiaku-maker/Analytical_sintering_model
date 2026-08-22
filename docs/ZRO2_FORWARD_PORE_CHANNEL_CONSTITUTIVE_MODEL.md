# Pore-channel constitutive model

The dimensionally audited surface coefficient is

\[
B_s(T)=\frac{D_s(T)\gamma_s\Omega^{4/3}}{k_BT},\qquad [B_s]=\mathrm{m^4/s},
\]

with \(\tau_s=C_sr^4/B_s\). The geometry product (C_s=1) is a
semi-phenomenological baseline, not a material activation barrier.

Surface coarsening transfers open-pore volume between size classes. PR pinch-off
uses

\[
I_{PR}=\frac{\lambda_{seg}}{2\pi rF_{GB}}-1,\qquad
P_{pinch}=\{1+\exp[-I_{PR}/w_I]\}^{-1}.
\]

The bounded (F_{GB}) surrogate is geometry-derived uncertainty; (w_I=0.15) and
the conservative precursor/isolated/closed partition are bounded uncertainties.
They were not fitted.

The regularization/damage mode compares a moderate-activity, not-overwide
regularization flux with a low-activity, wide-tail damaging flux. Regularization
narrows the width and preserves connectivity. Damage broadens the tail, reduces
connected fine pores, and transfers volume toward precursor, isolated, and closed
stores. Every transfer closes exactly in pore volume.

Pore-size pinning enters only through migration. (R_Z=4r/(3f_v)),
(P_Z\propto f_v/r), and a bounded `Gamma_Zener` multiplies intrinsic growth. It
does not modify open or closed density fluxes. The broader named energy ledger is
diagnostic; no equality forces released interfacial power into densification, and
the former GB-area-only equality is not a main closure.

This constitutive audit is not validation. No validation claim is made.
