# ZrO2 forward parameters

Diffusivities use `D0 exp(-Q/RT)` with GB `(0.056 m2/s, 380 kJ/mol)` and surface `(0.10 m2/s, 380 kJ/mol)`. These originate from 10 mol% Sc2O3 cubic zirconia and carry a composition-transfer caveat. Mobility uses only the 4.2 eV high-temperature activation energy. Its provisional absolute prefactor is not calibrated while the required inputs are absent.

Initial relative density is `2.88/5.95`, initial grain size is 10.20 nm, and connected pore volume is initialized as a lognormal distribution centered at 24.5 nm. Initial isolated and closed volumes are zero. Geometry defaults are `C_TJ=10`, `C_GB=2`. A single CS-calibrated site-density multiplier of 39.5 and `M0=5.8e-3 m4/(J s)` are used unchanged for every other thermal history.

The barrier is evaluated by PCHIP within its fitted 1557–2052 °C range. Because all target schedules lie partly or wholly below that range, parameters are conservatively clamped to the nearest fitted slice below 1557 °C and every such evaluation is flagged. This substantial extrapolation limitation prevents a validation claim.
