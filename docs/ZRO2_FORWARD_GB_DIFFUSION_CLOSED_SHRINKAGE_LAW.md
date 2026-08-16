# GB-diffusion closed-pore shrinkage candidate

The candidate uses `phi_closed A F_activity (P_cap-P_gas) D_GB Omega/(k_B T) r^-m` with `m=3` or `4`. For `m=3`, `C_closed_GB` is dimensionless; for `m=4`, it has units of length. The nominal values are geometric estimates, and factors 0.03–100 are explicitly labeled prefactor uncertainty.

`D_GB` retains the transferred ZrO2/ScSZ form (`D0=0.056 m²/s`, `Q=380 kJ/mol`). Capillary pressure is `2 gamma_s/r`; the gas-free GB submode uses that pressure directly and the gas candidate applies a bounded proxy. The law removes only `phi_closed` and reports inventory, radius, accommodation, activity, pressure, transport, prefactor, and density rate.

Nominal GB shrinkage restored strong closed removal. In the injected candidate state it exhausted the available closed inventory and erased the low-T density boundary. In naturally prepared forward maps it yielded no strict success and no finite window. It is neither calibrated nor validated.
