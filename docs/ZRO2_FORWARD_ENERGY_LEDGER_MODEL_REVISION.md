# ZrO2 forward-model energy-ledger revision

## Scope and decision

The previous equality

\[
P_{dens}=\zeta\eta_A\gamma_{GB}C_{GB}\dot G/G^2,
\qquad \zeta\eta_A=1,
\]

is too narrow to serve as the fundamental thermodynamic closure. It counts only
modeled grain-boundary-area loss, omits pore and external surface energies and gas
work, and feeds the grain-growth law back into activation stress. It is retained as
the `strict_GB_area_loss_balance` diagnostic ablation so that earlier figure
trajectories remain reproducible. This branch does not change those trajectories.

## Revised diagnostic ledger

The stored-energy inventory contains external surface, grain-boundary, pore-surface,
and gas proxy terms. Their signed time derivatives form the available-power ledger.
Expenditures are named separately: open densification, closed densification, surface
smoothing, conservative PR/coarsening, GB migration, pore/Zener/junction drag, gas
work, and an explicit unresolved channel. `Pi_dens`, `Pi_total`, excess power, and
budget violations are reported. Violations are not repaired by rescaling a rate.

The exported histories lack pore-bin arrays, so historical pore area is reconstructed
from total porosity and D50, and external surface uses a stated compact-scale proxy.
Those terms are geometry diagnostics, not new material properties.

## Stress interpretation

The primary physical candidate is local capillary force minus gas counterpressure.
A bounded-power variant checks that the resulting work is consistent with the full
ledger. Neither has been activated in the accepted integrator. Densification remains
local pore-volume removal coupled to activation/renewal kinetics; power ratios are
diagnostics or possible constraints, not a hidden efficiency multiplier.

The reconstructed histories contain frequent budget violations, including severe
violations when reconstructed available power approaches zero. This rejects the
claim that the current strict equality is a sufficient general ledger; it does not
by itself identify a unique replacement closure.

This is a diagnostic model revision and **not validation**. No validation claim is
made.
