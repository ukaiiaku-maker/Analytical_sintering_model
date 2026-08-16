# Binwise closed-pore evolution addendum

Closed pores are now represented by bin-resolved number proxy, evolving radius, shrinkable fraction, geometry factor, gas pressure, and finite accommodation used/recovered. The density identity remains `rho = 1 - sum(phi_open + phi_iso + phi_closed)`. Closure and PR preparation are conservative; only named open and closed removal change density.

Local driving stress is `max(C_geom 2 gamma_s/r - P_gas, 0)`. The renewal mode uses the unchanged ZrO2 barrier and GB diffusivity in a serial nucleation/exchange/transport cycle, with `r^-4` default and `r^-3` sensitivity. The GB alternative uses dimensional transport with its prefactor labeled semi-phenomenological. Gas pressure uses a transparent ideal-compression proxy when enabled. Surface diffusion changes accommodation only.

Accommodation consumption is proportional to closed volume removed. Recovery and shape relaxation are local state/temperature laws. Infinite accommodation is not used in the comparison.

The identical 950 °C ramps formed very little late closed inventory and showed no measurable closed density contribution for any physical mapping. Injected GB/gas states densified universally and erased the lower boundary. Naturally prepared GB/gas states restored closed removal but produced no strict success; renewal retained the lower boundary but remained below strict success. Surface accommodation alone produced no closed density gain. No mode passed the lower-success-upper acceptance gate, so no Chen map was run.

Post-run `Q_closed_app` values are saved only as apparent diagnostics and are not inputs or ZrO2 properties. The outcome indicates incomplete first-step state preparation and unconstrained geometry/gas/accommodation mapping, not a need to invent another mechanism. No validation is claimed.
