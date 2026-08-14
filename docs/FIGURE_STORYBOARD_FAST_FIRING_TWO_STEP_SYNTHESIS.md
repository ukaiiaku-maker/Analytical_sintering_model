# Figure storyboard: fast-firing and two-step synthesis

## Story arc

The figure sequence should move from mechanism separation, to exact overlap, to the two individual property windows, to the Chen boundaries, and finally to family robustness, screening limitations, and experimental falsification. The eight final synthesis figures are the primary sequence. Candidate-693168 presentation figures supply detailed SI support.

## Main synthesis figures

### Figure 1 — `mechanism_chain_fast_vs_twostep`

- **Scientific question:** Can fast firing and two-step behavior coexist without sharing one causal scalar?
- **Source table:** conceptual synthesis constrained by `mechanism_attribution_summary.csv`, exact ablation tables, and `final_property_window_summary.csv`.
- **Main message:** Fast firing is a nucleation-onset timing effect; two-step behavior is a PR-prepared closed-accommodation effect; serial attainment and bounded states are shared constraints.
- **Caption draft:** Fast firing and two-step sintering are compatible but arise from different dominant pathways. Nucleation waiting creates slow-ramp growth exposure, whereas PR-prepared closed topology creates lower exhaustion and upper growth boundaries during the second step. The same perturbation vector is evaluated in separate model layers.
- **Caveat:** This schematic does not imply a new dynamically coupled model or a universal causal hierarchy for all Tier-B candidates.

### Figure 2 — `relative_property_phase_map_exact`

- **Scientific question:** Is there a finite exact relative-property region in which both behaviors survive?
- **Source table:** `material_property_window_exact_promotions.csv` merged with `material_property_window_scorecard.csv` for dimensionless coordinates.
- **Main message:** Seventy-three exact perturbations pass both mechanisms, distinct from 485 fast-only, 119 two-step-only, and 1,226 neither cases.
- **Caption draft:** Exact-promoted perturbations classified in nucleation-dominance and closed-shrinkage/growth coordinates. The overlap is finite and heterogeneous rather than a unique parameter point.
- **Caveat:** The axes are model-derived dimensionless groups. Their observed ranges are coverage limited, and the fast and two-step states come from separate existing model layers.

### Figure 3 — `fast_firing_property_window`

- **Scientific question:** What causal and material-property evidence supports nucleation-limited fast firing?
- **Source tables:** `fast_firing_ratio_curves.csv`, `fast_firing_OAT_window.csv`, and exact promotion ablation fields.
- **Main message:** The base ratio remains above 1.5 over an attained interval; `Delta Q_nuc` has a finite 0 to +50 kJ/mol survival window; nucleation-facile ablation weakens the effect while PR-off often preserves it.
- **Caption draft:** Exact fast-firing evidence identifies a finite nucleation-onset window. Lowering the barrier makes nucleation too facile, whereas sufficiently increasing it removes attainable separation. PR is not the primary causal channel in this envelope.
- **Caveat:** The trajectory panel is the representative exact base pass. Exact OAT perturbations are summarized by their trajectory metrics rather than plotted as full stored density curves.

### Figure 4 — `two_step_property_window`

- **Scientific question:** Which closed-store and growth properties preserve candidate-693168 two-step behavior?
- **Source tables:** `dense_candidate_693168_matched_density_curves.csv`, `two_step_OAT_window.csv`, and exact dimensionless groups.
- **Main message:** Candidate 693168 retains a large matched-density grain-size separation. `Delta Q_closed = -25…+100 kJ/mol` survives, the PR prefactor threshold is 0.3x, and a finite window depends on closed/growth selectivity.
- **Caption draft:** Two-step behavior is retained when first-step PR preparation produces a closed store and closed shrinkage remains selective against growth. The lower boundary is lost below the closed-barrier window, while inadequate PR preparation destroys the complete response.
- **Caveat:** The strong absolute reduction is not quantitatively calibrated; the closed-store fraction and accommodation trajectory remain internal model states requiring measurement.

### Figure 5 — `chen_window_mechanism_boundaries`

- **Scientific question:** Does candidate 693168 exhibit a complete window with physically distinct lower and upper failures?
- **Source table:** `candidate_693168_T2_classification_fine.csv` and `candidate_693168_Chen_boundaries_fine.csv`.
- **Main message:** The fine exact map contains density-exhaustion failures below the 925 °C boundary, success through 1205 °C, and grain-growth failures above it.
- **Caption draft:** Candidate 693168 produces a finite second-step window. Closed shrinkage/accommodation controls density exhaustion at low T2, whereas thermally activated grain growth controls the upper boundary.
- **Caveat:** The classification uses the fixed Tier-B target, time budget, and 20% growth tolerance; it is not a quantitative fit to a specific experimental map.

### Figure 6 — `six_TierB_family_property_summary`

- **Scientific question:** Is the conditional two-step response unique to candidate 693168?
- **Source table:** `six_TierB_family_mechanism_summary.csv`.
- **Main message:** Six exact Tier-B base candidates span different closed fractions, preparation growth, window widths, reductions, and destructive-ablation sets.
- **Caption draft:** The exact Tier-B family demonstrates multiple conditional topologies. Candidate 693168 is the strongest separation comparator, while other candidates provide low-closed-fraction and intermediate-topology tests.
- **Caveat:** Material-property robustness across the family is a reduced transfer from exact 693168 OAT, not six independent exact material searches.

### Figure 7 — `surrogate_vs_exact_warning`

- **Scientific question:** Can the 50k dimensionless screen be used as final evidence?
- **Source table:** `surrogate_vs_exact_comparison.csv`.
- **Main message:** The screen predicts 19,880 both-pass rows, but the exact promoted union contains only 73, a 272.3-fold discrepancy.
- **Caption draft:** Surrogate feasibility metrics substantially overpredict exact overlap. All final classifications and property-window claims therefore use exact-promoted trajectories only.
- **Caveat:** Exact evaluation was performed on promoted subsets, so the comparison is a warning about screen calibration rather than an exhaustive exact labeling of all 50,655 rows.

### Figure 8 — `experimental_falsification_targets`

- **Scientific question:** Which measurements can discriminate the proposed mechanisms?
- **Source table:** `experimental_falsification_targets.csv`.
- **Main message:** The model makes observable predictions for ramp onset, grain distributions, 3D pore state, connectivity, closed accommodation, and upper-bound growth mobility.
- **Caption draft:** Experimental measurements prioritized by the mechanism they constrain and the signature that would falsify the current attribution. Closed-pore fraction and accommodation trajectory are the primary unresolved calibration targets.
- **Caveat:** Several model stores currently map only to proxies; the measurement-to-state conversion must be defined before quantitative fitting.

## Candidate-693168 presentation figures for SI

### SI Figure A — `candidate_693168_final_panel`

- **Scientific question:** Does the strongest candidate reproduce, attain, and remain bounded under the complete audit?
- **Source tables:** candidate reproduction summary, matched-density curves, fine T2 classifications, ablation summary, and timestep audit.
- **Main message:** Candidate 693168 combines a large attained trajectory separation with a complete finite Chen window and timestep stability.
- **Caption draft:** Consolidated exact audit of candidate 693168 showing absolute trajectories, window boundaries, closed-store evolution, and robustness checks.
- **Caveat:** Conditional Tier B; the absolute magnitude and closed-store trajectory are not calibrated.

### SI Figure B — `candidate_693168_mechanism_schematic`

- **Scientific question:** How do PR preparation, closed transition, shrinkage, accommodation, and growth connect?
- **Source tables:** exact ablation summary and closed-accommodation history.
- **Main message:** The first step prepares a closed store; finite accommodation and thermal shrinkage control the lower response; growth controls the upper response.
- **Caption draft:** Causal mechanism schematic derived from exact destructive ablations for candidate 693168.
- **Caveat:** The diagram represents this candidate, not every Tier-B topology.

### SI Figure C — `candidate_693168_physical_time_histories`

- **Scientific question:** Are lower failure, success, and upper failure distinguishable in continuous physical time?
- **Source table:** `dense_candidate_693168_T2_scan_histories.csv`.
- **Main message:** Temperature, density, grain size, and state trajectories separate the three regimes without resetting second-step time.
- **Caption draft:** Continuous physical-time histories for representative lower failure, success, and upper growth failure.
- **Caveat:** The common time budget is part of the classification and must be retained in experimental comparisons.

### SI Figure D — `candidate_693168_fine_Chen_filled_window` and classification map

- **Scientific question:** Are success points contiguous and bracketed?
- **Source tables:** candidate fine T2 classification and boundary tables.
- **Main message:** A filled 925–1205 °C band lies between retained density-exhaustion and grain-growth failures.
- **Caption draft:** Fine exact Chen-style classification and filled success window for candidate 693168.
- **Caveat:** The displayed band uses the fixed Tier-B growth tolerance; changing the acceptance rule would change its width.

### SI Figure E — `candidate_693168_PR_damage_closed_transition_history`

- **Scientific question:** Does first-step preparation leave a persistent state that supports the second step?
- **Source table:** `candidate_693168_closed_accommodation_history.csv` and dense candidate histories.
- **Main message:** PR memory, closed transition, closed fraction, and accommodation evolve coherently across the switch.
- **Caption draft:** Evolution of the PR-prepared closed-store state through first-step preparation and second-step densification.
- **Caveat:** These are conservative model-store variables and proxies, not directly measured pore fractions or energies.

### SI Figure F — `candidate_693168_ablation_waterfall`

- **Scientific question:** Which channels are necessary for the joint trajectory/window response?
- **Source table:** `candidate_693168_ablation_summary_final.csv`.
- **Main message:** Removing PR damage, closed transition, closed shrinkage, or finite accommodation destroys the joint result.
- **Caption draft:** Exact ablation waterfall identifying the causal closed-store chain in candidate 693168.
- **Caveat:** Ablation survival does not establish that a channel is universally irrelevant; it defines necessity only within this candidate and tested protocol.

### SI Figure G — six-candidate ablation summary

- **Scientific question:** Do all Tier-B candidates depend on the same secondary channels?
- **Source tables:** six-candidate comparison and family mechanism summary.
- **Main message:** The family shares closed-support requirements but differs in junction, drag, stress, and topology sensitivities.
- **Caption draft:** Candidate-specific destructive-ablation patterns across the six exact Tier-B base cases.
- **Caveat:** The candidate family is qualitative and not independently calibrated.

## Recommended manuscript ordering

Use synthesis Figures 1, 3, 4, 5, and 8 in the main text if space permits. Place the exact phase map and surrogate warning together to establish evidence discipline. Use the six-candidate family summary to prevent over-identification with candidate 693168. Move detailed physical-time histories, state trajectories, timestep checks, decoder provenance, and full ablations to the SI.
