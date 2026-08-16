# Resolved-rule forward-model implementation

The new `mechanism_mode="resolved_rules"` is separate from the earlier baseline. It keeps the ZrO2 barrier JSON, GB and surface diffusivities, energy balance, schedules, and targets unchanged.

Densification uses a serial nucleation/exchange/transport cycle and a connected-removable eligibility factor. Density changes only through named open- and closed-pore shrinkage. PR adjacent-bin coarsening and relocation conserve pore volume while building bounded preparation memory. A named conservative transition moves PR-prepared open/precursor porosity into the closed store. Late open-path eligibility exhausts at closure, and further density requires closed shrinkage with finite accommodation.

Grain growth is intrinsic Arrhenius capillarity migration multiplied by pore/Zener, closed-accommodation, PR-memory, and event-completion activity factors. Mobility never enters the density flux. No Dong–Chen low-temperature mobility transition or schedule label is used. This implementation is conditional and is not validated.
