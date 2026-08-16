# ZrO2 open/closed rate-handoff audit

This focused implementation audit is conditional and is not validation. The barrier, GB and surface diffusivities, intrinsic-mobility audit, schedules, time budgets, and strict targets were held fixed. Eleven designs covered the unchanged resolved default, diagnostic open recovery, closed-rate factors 3–100, closed-inventory factors 3–30, a local availability-balanced handoff, and an externally injected candidate-like state.

The unchanged default reproduced the prior 50 °C/min result (`rho=0.900504`, `G=0.258471 µm`). Diagnostic open recovery reached `rho=0.971813` and retained a 76% matched-density smaller-grain sign, showing that open-channel rate recovery is necessary. It did not create a strict Chen success or finite window, so it is not sufficient.

Closed-rate and inventory factors did not recover fast density because the fast ramp creates too little closed inventory for that channel to carry meaningful density. Even 100x closed rate contributes only about 1.7e-6 integrated fast-ramp density. The balanced law correctly leaves open eligibility essentially unchanged when closed availability is negligible, but consequently cannot repair the missing rate by itself.

No forward-eligible mode, and not even candidate-state injection, generated a strict success or finite bracketed window. The next defect is therefore broader than switch-state magnitude alone: useful closed-state persistence/rate and the open-rate formulation must be reconciled while upper-bound intrinsic growth remains controlled.
