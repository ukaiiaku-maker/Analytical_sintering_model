# Fixed-parameter ZrO2 processing-window prediction

No data fitting or parameter optimization was performed, and no model physics was changed. The ZrO2 barrier, GB diffusivity, surface diffusivity, PR/topology laws, and defined open/closed shrinkage laws remain those at source commit `2249f02`. The failed global mobility fit was not used; intrinsic mobility is the current default provisional branch and never enters density flux.

The campaign evaluated 5,020 unique ZrO2 initial/preparation states and promoted 24 states for exact cloned-switch-state scans. Candidate 693168 is not parameterized as ZrO2 and was used only as a mechanism-response reference. Its state values were not copied. ZrO2 independently searched `rho0`, `G0`, pore state, `rho_switch`, and the resulting `G1`.

This is a fixed-physics forward-mapping campaign, not calibration or validation. The model is not validated.
