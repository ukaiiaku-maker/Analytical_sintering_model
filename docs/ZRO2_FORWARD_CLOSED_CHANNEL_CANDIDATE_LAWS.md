# Candidate physical closed-channel laws

This registry is analysis only, not validation, and installs no production law. Candidate mappings include GB-diffusion-controlled closed-pore removal, surface-diffusion shape accommodation, renewal-limited serial shrinkage, gas/accommodation-limited shrinkage, and an explicitly empirical reduced closure.

GB-diffusion and renewal-limited laws are recommended for a next implementation comparison because they provide dimensional transport and state requirements. Surface diffusion may relax shape but must not count as density without a removal path. Gas-limited closure requires trapped-pressure state. An empirical `k0_closed exp(-Q_closed_emp/RT)` law is permissible only if labeled empirical; `Q_closed_emp` would not be a ZrO2 material property.
