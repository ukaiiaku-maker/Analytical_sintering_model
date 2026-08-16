# Defined-law port into the ZrO2 forward model

`mechanism_mode="defined_laws_port"` assembles the already-defined serial renewal, bounded capillary/Onsager stress, conservative PR topology transfer, named open and closed removal, finite closed accommodation, and intrinsic-growth-times-migration laws. It does not add a mechanism.

The renewal channel evaluates the fitted `G*(sigma,T)`, `r_nuc`, GB-diffusion sink time, `Lambda`, and `Lambda/(1+Lambda)` once. Eligibility is explicit, and open density change is tied to removable open pore volume. Surface work is used in the stress/power accounting; it does not directly set density. PR bin crossing and open-to-precursor/closed transfers are conservative.

`closed_mapping_mode` provides source-traceable mappings: current baseline, defined-law port, reduced candidate-law transfer, GB diffusion, renewal-limited, gas accommodation, and empirical diagnostic. The reduced and empirical modes are not physical ZrO2 calibration.

The named discussion DOCX was unavailable locally; its controlling rules were taken from the supplied continuation brief. The renewal/Onsager equations were traced to the available Python and MATLAB source scripts. No validation is claimed.
