# Equation-to-code traceability

## How to use the trace

The authoritative machine table is results/equation_functional_form_audit/equation_to_source_function.csv. Every row contains an equation ID, source file, source function, resolved line number, code excerpt, result family, and evidence role. Line numbers are regenerated from stable source anchors by equation_functional_form_audit.py.

## Final fast-firing path

| IDs | Source function | Role |
|---|---|---|
| FF-01–FF-10 | separated_fast_chen_model.material_rates | stress, serial times, activity, density, PR ablation, growth |
| FF-11 | separated_fast_chen_model.topology_growth_factor | migration-only diagnostic closure |
| FF-12 | separated_fast_chen_model.run | conservative adjacent-bin update |
| MET-01 | observable_trajectory_effect_audit.matched_curve | joint-support interpolation |
| MET-02 | relative_material_property_window_search.fast_metric | ratio/span/pass rule |

The exact fast-firing result uses FF-01–FF-10, MET-01, and MET-02. FF-09/FF-12 remain necessary to state the PR-off ablation, not to claim PR control.

## PR and pore-location controls

PR-01–PR-09 trace to pr_desintering_memory_model.local_competition. PL-01–PL-06 trace to pore_location_topology_model. AL-01–AL-03 trace to pore_location_agentic_model. These equations document conservative topology evolution and action allocation, but their final role is negative-control or diagnostic.

## Persistent junction and multihit

TJ-01–TJ-06 and TJ-08 trace to agentic_mechanism_model.local_mechanism. TJ-07 traces to agentic_mechanism_model.poisson_completion. The exact Poisson equation belongs to this earlier mechanism family.

The decoder-corrected local model uses LR-09, a sigmoid completion proxy. A paper must not cite TJ-07 as the candidate-693168 closed-accommodation law.

## Candidate 693168

| IDs | Source function | Role |
|---|---|---|
| LR-01, LR-07, LR-08, LR-11, LR-12 | interacting_local_region_model.advance | conservation, state updates, bounded accommodation, growth step |
| LR-02–LR-06, LR-09–LR-10 | interacting_local_region_model.local_fluxes | local activity, open/closed shrinkage, PR, migration |
| MET-03–MET-04 | audit_candidate_693168_closed_accommodation.classify | growth fraction and point class |
| MET-05 | adaptive_T2_boundary_search.status | complete boundary topology |
| MET-06 | audit_candidate_693168_closed_accommodation.score_histories | high-density trajectory reduction |

The candidate audit reports accommodation fraction as P_comp_closed for plotting compatibility, but explicitly stores Lambda_closed and K_closed as unavailable. That diagnostic naming does not create a Poisson law.

## Property attribution

PROP-01–PROP-03 trace the additive/multiplicative design and Latin-hypercube generation. PROP-04 and PROP-05 are dimensionless diagnostic groups. PROP-06 and PROP-07 are screening-only response surfaces. PROP-08 is the exact classification composed from exact fast and exact two-step flags.

## Missing path

No single file named local_region_decoder_corrected_dynamic_search.py exists. The current trace uses:

- interacting_local_region_decoder.py for the complete decoded parameter vector;
- interacting_local_region_model.py for state evolution;
- local_region_decoder_corrected_postprocess.py for corrected postprocessing;
- audit_candidate_693168_closed_accommodation.py for exact candidate reproduction and classification.

No archived result directory was restored.
