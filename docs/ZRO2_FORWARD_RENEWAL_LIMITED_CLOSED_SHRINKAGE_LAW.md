# Renewal-limited closed-pore shrinkage candidate

This candidate uses a serial cycle: `tau_cycle_closed = tau_nuc_closed + tau_exchange_closed + tau_transport_closed`, followed by `rho_dot_closed = phi_closed A_closed eps_closed/tau_cycle_closed`. Local stress is the nonnegative capillary-minus-gas pressure scaled by `C_sigma_closed`. The existing fitted stress/temperature barrier supplies `Gstar_closed`; GB diffusion supplies the transport time. Exchange is explicitly represented by a fast 1 s diagnostic value because no independently constrained closed-pore exchange law is available.

The implementation reports stress, barrier, nucleation rate/time, exchange time, transport time, cycle time, Lambda, and the closed rate. Its nominal candidate-state maximum closed contribution was 0.00320 and it produced no strict success. This is an interpretable candidate, not a tuned or validated law.
