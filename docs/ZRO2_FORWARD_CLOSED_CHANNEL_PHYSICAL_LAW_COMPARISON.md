# ZrO2 closed-channel physical-law comparison

This branch implements a schedule-independent registry of closed-pore laws while preserving the fitted barrier, the 380 kJ/mol GB and surface diffusivities, and the inherited high-temperature mobility law. The default remains `resolved_proxy_current`. `Q_closed_app` is not an input.

The PDF-conditioned ramps share `rho=0.66`, `G=50 nm`, and the same open-pore population at 950 °C. Nominal final densities for 5/50 °C min⁻¹ were 0.937011/0.900504 for the inherited proxy, 0.939140/0.900671 for GB shrinkage, 0.936878/0.900504 for renewal-limited shrinkage, and 0.939140/0.900671 for the bounded-gas candidate. Surface accommodation alone and the nominal empirical diagnostic gave essentially the no-new-removal result. Matched-density median `G_5/G_50` was about 0.60 for all nominal laws.

The injected candidate-like state (`rho=0.88`, `G=117 nm`, closed fraction 0.649, `A_closed=0.152`, PR memory 1) is a state-sufficiency diagnostic, not a forward prediction. Nominal GB and gas candidates removed all available closed inventory and produced some strict T2 successes, but they also produced universal density attainment across the candidate-state scan and lost the lower boundary. Renewal-limited shrinkage was much smaller (maximum inventory-bounded closed contribution 0.00320) and produced no strict success.

In naturally prepared forward states, every eligible nominal law produced zero strict mini-map successes and zero finite windows. GB and gas laws restored substantial closed removal (up to the available inventory) but did not create the required lower-success-upper topology. Thus closed shrinkage can be restored dimensionally, but the present nominal mappings do not recover a strict Chen window.

Candidate 693168 remains a conditional comparator: its closed fraction (~0.649), accommodation (~0.152), substantial closed density contribution (~0.244 in the inherited reduced layer), and conditional interval describe a qualitative state/rate effect. Its parameters were not copied.

The prefactor envelopes are uncertainty audits, not calibration. The empirical closure is diagnostic only; its `Q_closed_emp` is not a ZrO2 material property. These results are not validation.
