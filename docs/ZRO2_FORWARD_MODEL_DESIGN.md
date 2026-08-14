# ZrO2 forward model design

The explicit state contains density, grain size, fixed logarithmic pore-radius bins, open/isolated/closed pore volumes, and a bounded closed-pore accommodation proxy. Density is always the identity `1 - sum(phi_open + phi_iso + phi_closed)`.

Triple-line geometry uses `L_TJ/V=C_TJ/G^2`, `A_GB/V=C_GB/G`, `rho_TL_area=(C_TJ/C_GB)/G`, and event strain `b*rho_TL_area`. Effective stress solves `(1-rho)*sigma*edot=P_surf` on bounded stress limits and records excess surface power.

Connected pore removal is apportioned by `phi_i/r_i^4`. Surface redistribution transfers connected volume to the next radius bin, conserving volume and connectivity. Separate instantaneous density gates move volume into isolated and closed reservoirs. Only a named closed-reservoir shrinkage flux removes closed pore volume.

Grain growth uses the high-temperature Arrhenius mobility and a smooth combination of distribution-aware Zener pinning and mobile-pore drag. Thermal path objects are kept outside constitutive modules.
