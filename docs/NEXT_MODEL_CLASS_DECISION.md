# Next Model-Class Decision

No new model is implemented here. This is a proposal requiring user approval.

## Options

### 1. Minimal spatial pore-network / percolation model — recommended if data confirm the effect

Pros: represents connected transport paths, pinch-off, pore clusters, local
matrix densification, and interacting defect-rich regions. It directly
addresses the failure of weighted independent cohorts and scalar topology
gates. Cons: introduces graph topology, exchange rules, and new calibration
parameters.

### 2. Stochastic local-region ensemble

Pros: cheaper and can represent distributions of density, pore topology,
stress, and grain size. Cons: remains mean-field unless regions exchange pore
volume, stress, or boundary motion; the independent-cohort audit shows that
missing interactions may be decisive.

### 3. Phase-field or graph-based reduced model

Pros: directly follows pore pinch-off, clustering, grain topology, and boundary
migration. Cons: substantially more expensive and harder to calibrate; likely
premature without discriminating data.

### 4. Experimental calibration first — immediate recommendation

Pros: determines whether the required 1.5 ratio is real, its density window,
and which topology observables carry memory. It prevents another unconstrained
search. Cons: requires data not yet encoded.

## Recommendation

Do not add another scalar closure. Move to a minimal spatial/network or
interacting stochastic local-region model only if experimental data confirm a
large `G(rho)` separation that the reduced model must reproduce. The immediate
next step is experimental calibration and uncertainty encoding. If confirmed,
the first new model should be a small graph of connected pore/GB regions with
explicit pinch-off, local transport, and inter-region exchange—not a full
phase-field calculation.

