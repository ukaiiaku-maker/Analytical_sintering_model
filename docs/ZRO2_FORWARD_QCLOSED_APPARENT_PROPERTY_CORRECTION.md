# Correction: apparent closed-channel coordinate

This correction is diagnostic and not validation. `Q_closed` is not an independently defined physical ZrO2 input. Prior `Q_closed_eff` values are renamed `Q_closed_app` or `Q_closed_proxy`: apparent finite-difference slopes of a composite model rate. They cannot be used as direct material-property pass/fail tests.

The trusted physical inputs remain the fitted stress-resolved nucleation function `G*(sigma,T)`, GB diffusivity, surface diffusivity, and the uncertain high-temperature growth branch. The corrected statement is: the forward model's apparent closed-channel coordinate does not map cleanly to the reduced successful coordinate; the physical closed-channel law has not yet been defined.

The apparent slopes vary with normalization. Raw, inventory-normalized, and inventory/accommodation-normalized rates give about 126–135 kJ/mol in supported states. Dividing additionally by the implemented closed activity removes essentially the full apparent thermal slope, giving approximately zero. This demonstrates emergence from activity coupling rather than a standalone material barrier.
