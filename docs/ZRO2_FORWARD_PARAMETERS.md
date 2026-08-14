# ZrO2 forward parameters

Diffusivities use `D0 exp(-Q/RT)` with GB `(0.056 m2/s, 380 kJ/mol)` and surface `(0.10 m2/s, 380 kJ/mol)`. These originate from 10 mol% Sc2O3 cubic zirconia and carry a composition-transfer caveat. Mobility uses only the 4.2 eV high-temperature activation energy. Its provisional absolute prefactor is not calibrated while the required inputs are absent.

Initial relative density is `2.88/5.95`, initial grain size is 10.20 nm, and connected pore volume is initialized as a lognormal distribution centered at 24.5 nm. Initial isolated and closed volumes are zero. Geometry defaults are `C_TJ=10`, `C_GB=2`.
