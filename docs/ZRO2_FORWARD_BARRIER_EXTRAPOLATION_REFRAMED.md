# Barrier extrapolation, reframed after conditioning

Initialization and barrier extrapolation are now separated. Starting every comparison at the observed 950 °C state removes the need for the model to generate the missing green-to-950 evolution before evaluating post-950 kinetics.

For the conditioned trajectory-calibration parameters, the JSON nearest-slice clamp gives final densities 0.9950 and 0.9592 for 5 and 50 °C/min. The fixed-low-temperature-slope diagnostic gives 0.9957 and 0.9621; the mathematical continuation gives 0.9949 and 0.9591. These modest changes show that the conditioned post-950 discrepancy is not dominated by the clamp in the way the full-process onset was. The generic anchored barrier gives much higher density but is non-JSON diagnostic evidence only.

Barrier extrapolation below 1557 °C remains a primary uncertainty. No barrier was refitted, the generic anchor is not final evidence, the low-temperature Dong–Chen transition is absent, and no validation claim is made.
