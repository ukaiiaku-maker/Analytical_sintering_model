# Mobility assumptions

`M_GB(T)=M0 exp(-Q_M/RT)` is used only in intrinsic grain growth. The tested modes are `fixed_highT_literature`, `CS_endpoint_calibrated`, `CS_curve_regularized`, `bounded_uncertainty_factor`, and `activation_energy_envelope`. Calibration-labeled modes inherit the one global CS calibration; they are not method-specific fits.

Mobility does not enter nucleation, exchange, transport, renewal activity, effective stress, PR transfer, open shrinkage, closed shrinkage, or density identity. There is no low-temperature transition and no schedule-dependent mobility. Activation-energy cases preserve the 1500 °C mobility anchor. All results are uncertainty diagnostics, not validation.
