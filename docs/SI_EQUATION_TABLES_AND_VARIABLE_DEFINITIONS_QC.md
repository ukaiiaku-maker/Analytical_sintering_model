# SI equation tables and variable definitions

The complete row-level registry is in results/equation_functional_form_audit/. This document provides the paper-facing subset. Units marked “code-normalized” require calibration before dimensional interpretation.

## Table S1. Variables and units

| Symbol | Definition | Units/status |
|---|---|---|
| \(T,T_C\) | absolute and Celsius temperature | K, °C |
| \(\rho\) | relative density reconstructed from pore volume | dimensionless |
| \(G\) | mean/local grain size | m in fast model; nm-like local code unit |
| \(r_i,r_0\) | pore-bin and reference radii | m |
| \(\phi_i\) | pore-volume fraction in bin \(i\) | dimensionless |
| \(\phi^{GB},\phi^{TJ},\phi^{iso},\phi^{closed}\) | location-resolved pore stores | dimensionless |
| \(N_i,N^{closed}\) | pore-number proxy | code-normalized |
| \(\sigma_{\rm loc}\) | fast-model activation stress | Pa |
| \(\sigma_{\rm res}\) | local residual-stress state | normalized proxy |
| \(\tau_{\rm nuc},\tau_{\rm ex},\tau_{\rm tr},\tau_{\rm cyc}\) | serial waiting times | s |
| \(a\) | renewal activity | dimensionless |
| \(\dot\rho,\dot G\) | density and grain-size rates | s\(^{-1}\); grain unit s\(^{-1}\) |
| \(E_G\) | density gain per logarithmic grain growth | dimensionless ratio |
| \(C_{GB},C_{TJ}\) | GB pore coverage and TJ occupancy | dimensionless |
| \(X_J\) | persistent junction population | bounded dimensionless state |
| \(\Lambda_{TJ},K_{TJ}\) | expected and required TJ event counts | dimensionless counts |
| \(P_{\rm comp}\) | completion probability | dimensionless |
| \(A_{\rm closed},A_{\rm cap}\) | available/capacity accommodation | dimensionless |
| \(R_{\rm fast}\) | matched-density grain-size ratio | dimensionless |
| \(g_2\) | second-step grain-growth fraction | dimensionless |
| \(\Theta_{\rm nuc}\) | nucleation dominance ratio | dimensionless |
| \(S_{\rm closed/growth}\) | closed-shrinkage/growth selectivity | dimensionless |

## Table S2. Material kinetic parameters

| Parameter | Meaning | Unit | Final role |
|---|---|---|---|
| \(Q_{\rm nuc},\nu_0,v^*\) | nucleation barrier, attempt rate, activation volume | J mol\(^{-1}\), s\(^{-1}\), m\(^3\) | fast firing |
| \(Q_{\rm ex},\tau_{\rm ex,0}\) | exchange barrier and time prefactor | J mol\(^{-1}\), s | serial cycle |
| \(Q_{\rm tr},\tau_{\rm tr,0}\) | transport barrier and \(G^2\)-scaled prefactor | J mol\(^{-1}\), s m\(^{-2}\) | serial cycle |
| \(Q_{\rm GB},D_{\rm GB,0}\) | intrinsic growth barrier/prefactor | J mol\(^{-1}\), m\(^2\) s\(^{-1}\) | growth |
| \(Q_s,D_{s,0},k_{\rm PR}\) | surface redistribution parameters | J mol\(^{-1}\), m\(^2\) s\(^{-1}\), registered scale | PR ablation |
| \(\epsilon_{\rm event},\zeta_\eta\) | density gain and scale ratio | dimensionless | densification |
| \(Q_{\rm closed},k_{\rm closed}\) | local closed shrinkage barrier/prefactor | J mol\(^{-1}\), code s\(^{-1}\) | candidate 693168 |
| \(Q_g,k_g\) | local growth barrier/prefactor | J mol\(^{-1}\), code grain\(^2\) s\(^{-1}\) | upper boundary |

## Table S3. Pore-state variables

| Store/state | Source | Update type | Density role |
|---|---|---|---|
| \(\phi^{GB}\) | initial connected pore inventory | PR, sweep, recapture, open loss | open shrinkage can remove |
| \(\phi^{TJ}\) | initial/PR TJ inventory | PR input, sweep, closed transition | transition conservative |
| \(\phi^{iso}\) | isolated inventory | PR, detachment, closed transition | not directly removed by open flux |
| \(\phi^{closed}\) | initial/PR/transition closed inventory | conservative input, closed loss | closed shrinkage can remove |
| \(C_{\rm rem}\) | GB store divided by total local pore | reconstructed | open eligibility |
| \(A_{\rm closed}\) | finite accommodation state | recovery minus closed loss | multiplies closed shrinkage |
| PR damage memory | normalized cumulative PR state | production/relaxation | indirect topology role |
| \(X_J\) | persistent junction state | production/relaxation | migration only |

## Table S4. PR/de-sintering flux terms

| ID | Term | Functional form | Conservation |
|---|---|---|---|
| FF-09 | separated PR propensity | surface Arrhenius × \((1-a)^2\) × fine fraction | yes |
| PR-01–PR-05 | production PR source | radius, activity, thermal, topology gates | source definition |
| PR-06 | smoothing/GB-to-TJ/TJ-to-isolated | nonnegative shares sum to one | yes |
| PR-07–PR-08 | PR/density competition | normalized hazards; removal scaled by density share | PR itself does not densify |
| LR-06–LR-07 | local PR and five-way partition | local activity gate and decoded partitions | yes |
| LR-08 | isolated/TJ-to-closed transition | Arrhenius × source × \(\rho^3\) | yes |

## Table S5. Migration activity factors

| Class | Registry ID | Factor | Evidence role |
|---|---|---|---|
| Intrinsic | FF-10 | \(D_{\rm GB}\gamma_{\rm GB}/G\) | final fast/growth |
| Class A drag | PL-06 | \((1+R_{GB}+R_{TJ})^{-1}\) | negative control |
| Persistent drag | TJ-01–TJ-03 | \((1+R_J)^{-1}\) contribution | prior Chen diagnostic |
| Class B multihit | TJ-04–TJ-07 | exact Poisson completion | prior Chen diagnostic |
| Coupled | TJ-08 | \(P_{\rm comp}/(1+R_J+R_{pore})\) | prior Chen diagnostic |
| Local candidate | LR-09–LR-10 | sigmoid completion divided by drag | conditional final two-step |
| Separated topology | FF-11 | structural plus pore drag proxy | diagnostic |

## Table S6. Closed-pore accommodation variables

| Variable/parameter | Implemented meaning | Limitation |
|---|---|---|
| \(\phi^{closed}\) | modeled closed pore-volume store | not directly a stereological measurement |
| \(N^{closed}\) | closed pore-number proxy | code-normalized |
| \(\xi_r=\phi^{closed}/N^{closed}\) | closed radius proxy | not calibrated diameter |
| \(A_{\rm closed}\) | currently available accommodation | phenomenological state |
| \(A_{\rm cap}\) | maximum accommodation | requires calibration |
| \(\tau_A\) | recovery time | requires independent transient data |
| \(r_g\) | gas/sintering-pressure ratio proxy | no explicit gas equation |
| \(P_{\rm comp,closed}\) | audit alias \(A_{\rm closed}/A_{\rm cap}\) | not Poisson completion |
| \(\Lambda_{\rm closed},K_{\rm closed}\) | not implemented | intentionally unavailable |

## Table S7. Classification criteria

| Label | Density condition | Growth condition | Additional rule |
|---|---|---|---|
| SUCCESS | target attained | \(g_2\le g_{\rm tol}\) | first step eligible |
| GRAIN_GROWTH_FAILURE | target attained | \(g_2>g_{\rm tol}\) | — |
| DENSIFICATION_EXHAUSTION_FAILURE | target not attained | \(g_2\le g_{\rm tol}\) | — |
| MIXED_FAILURE | target not attained | \(g_2>g_{\rm tol}\) | — |
| UNATTAINABLE_FIRST_STEP | switch not attained | not scored | excluded |
| INELIGIBLE_TARGET_ALREADY_REACHED | target reached in step 1 | not scored | excluded |
| COMPLETE_WINDOW | lower and upper brackets | finite success band | \(T_2<T_1\), width ≥25 °C |
| Fast pass | both paths attained | \(R_{\rm fast}\ge1.5\) | continuous \(\Delta\rho\ge0.03\) |

## Table S8. Mechanism outcomes and tested conditions

| Mechanism family | Tested conditions/count | Outcome | Final role |
|---|---:|---|---|
| Aggregate topology/PR | historical bounded maps | negative fast control | superseded architecture |
| Pore placement and action allocation | historical bounded maps | insufficient | negative controls |
| Persistent junction + TJ multihit | q0/q1 dynamic maps | finite prior Chen windows | diagnostic prior family |
| Production PR memory | production maps | observable memory, weak fast trajectory | negative control |
| Heterogeneity/defect/connected-sink | bounded campaigns | short, unattained, or censored effects | negative controls |
| Separated nucleation material | 1,000 exact fast promotions | controlling fast-firing mechanism | final |
| Decoder-corrected local region | six exact Tier-B base candidates | conditional two-step family | final conditional |
| Candidate 693168 | 101 fine T2 points, 12 ablations, 30 initial-condition rows | complete finite window; strong trajectory | conditional Tier B |
| Relative property campaign | 50,655 screen; 1,000 + 1,000 exact promotions | 1,903 exact: 485/119/73/1,226 | final attribution |

## Source and evidence note

The CSV registry contains the complete variable and parameter lists, source excerpts, and evidence tags. Surrogate equations PROP-06 and PROP-07 are screening-only. Exact classifications use PROP-08. Candidate 693168 remains conditional Tier B, and its closed/accommodation state is the principal calibration target.

## Table S9. Complete equation-to-source and evidence map

| Equation ID | Equation name | Source file | Source function | Evidence role |
|---|---|---|---|---|
| FF-01 | Capillary/pore activation stress | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-02 | Stress-assisted nucleation time | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-03 | Exchange time | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-04 | Transport time | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-05 | Serial cycle time | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-06 | Renewal activity | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-07 | Connected-fine eligibility | separated_fast_chen_model.py | material_rates | implemented_proxy |
| FF-08 | Serial-cycle densification | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-09 | Low-activity surface redistribution | separated_fast_chen_model.py | material_rates | causal_ablation |
| FF-10 | Intrinsic grain growth | separated_fast_chen_model.py | material_rates | final_evidence |
| FF-11 | Migration-only topology factor | separated_fast_chen_model.py | topology_growth_factor | diagnostic_only |
| FF-12 | Conservative adjacent-bin PR update | separated_fast_chen_model.py | run | causal_ablation |
| PR-01 | Fine-bin logistic weight | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-02 | Low-renewal gate | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-03 | Relative PR thermal factor | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-04 | Connected fine topology gate | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-05 | PR source by bin | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-06 | Conservative PR partitions | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-07 | Densification/PR competition | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-08 | Competition-scaled removal | pr_desintering_memory_model.py | local_competition | negative_control |
| PR-09 | PR/densifying work proxies | pr_desintering_memory_model.py | local_competition | diagnostic_only |
| PL-01 | Three-store pore conservation | pore_location_topology_model.py | run | negative_control |
| PL-02 | GB coverage | pore_location_topology_model.py | topology_diagnostics | negative_control |
| PL-03 | TJ occupancy | pore_location_topology_model.py | topology_diagnostics | negative_control |
| PL-04 | Location densification eligibility | pore_location_topology_model.py | instantaneous | negative_control |
| PL-05 | Location density rate | pore_location_topology_model.py | instantaneous | negative_control |
| PL-06 | Class-A pore pinning | pore_location_topology_model.py | instantaneous | negative_control |
| AL-01 | Nonnegative action propensity | pore_location_agentic_model.py | score_actions | negative_control |
| AL-02 | Normalized action allocation | pore_location_agentic_model.py | score_actions | negative_control |
| AL-03 | Constrained pair allocation | pore_location_agentic_model.py | allocated_fluxes | negative_control |
| TJ-01 | Persistent-junction production | agentic_mechanism_model.py | local_mechanism | negative_control |
| TJ-02 | Junction relaxation | agentic_mechanism_model.py | local_mechanism | negative_control |
| TJ-03 | Persistent drag | agentic_mechanism_model.py | local_mechanism | negative_control |
| TJ-04 | Required TJ packet count | agentic_mechanism_model.py | local_mechanism | negative_control |
| TJ-05 | Pore/structural constraint split | agentic_mechanism_model.py | local_mechanism | negative_control |
| TJ-06 | Expected TJ events | agentic_mechanism_model.py | local_mechanism | negative_control |
| TJ-07 | Poisson completion | agentic_mechanism_model.py | poisson_completion | negative_control |
| TJ-08 | Coupled TJ mobility | agentic_mechanism_model.py | local_mechanism | negative_control |
| LR-01 | Four-store density identity | interacting_local_region_model.py | advance | final_evidence |
| LR-02 | Local activity proxy | interacting_local_region_model.py | local_fluxes | implemented_proxy |
| LR-03 | Open-pore shrinkage | interacting_local_region_model.py | local_fluxes | final_evidence |
| LR-04 | Closed-radius proxy | interacting_local_region_model.py | local_fluxes | implemented_proxy |
| LR-05 | Closed-pore shrinkage | interacting_local_region_model.py | local_fluxes | final_evidence |
| LR-06 | Local PR damage | interacting_local_region_model.py | local_fluxes | final_evidence |
| LR-07 | Conservative PR store partition | interacting_local_region_model.py | advance | final_evidence |
| LR-08 | Conservative closed transition | interacting_local_region_model.py | local_fluxes | final_evidence |
| LR-09 | Local TJ completion proxy | interacting_local_region_model.py | local_fluxes | implemented_proxy |
| LR-10 | Local migration/growth | interacting_local_region_model.py | local_fluxes | final_evidence |
| LR-11 | Closed accommodation update | interacting_local_region_model.py | advance | final_evidence |
| LR-12 | Inverse-size exact growth step | interacting_local_region_model.py | advance | final_evidence |
| MET-01 | Matched-density interpolation | observable_trajectory_effect_audit.py | matched_curve | final_evidence_metric |
| MET-02 | Fast-firing ratio and rule | relative_material_property_window_search.py | fast_metric | final_evidence_metric |
| MET-03 | Second-step growth fraction | audit_candidate_693168_closed_accommodation.py | classify | final_evidence_metric |
| MET-04 | Chen point classification | audit_candidate_693168_closed_accommodation.py | classify | final_evidence_metric |
| MET-05 | Complete Chen window | adaptive_T2_boundary_search.py | status | final_evidence_metric |
| MET-06 | Two-step matched-density reduction | audit_candidate_693168_closed_accommodation.py | score_histories | final_evidence_metric |
| PROP-01 | Activation-energy perturbation | relative_material_property_window_search.py | apply_fast | final_evidence_metric |
| PROP-02 | Prefactor perturbation | relative_material_property_window_search.py | design | final_evidence_metric |
| PROP-03 | Latin-hypercube design | relative_material_property_window_search.py | design | screening_only |
| PROP-04 | Nucleation dominance | mechanism_dimensionless_groups.py | fast_groups | diagnostic_only |
| PROP-05 | Closed/growth selectivity | mechanism_dimensionless_groups.py | two_step_groups | diagnostic_only |
| PROP-06 | Surrogate fast score | relative_material_property_window_search.py | screen | screening_only |
| PROP-07 | Surrogate two-step score | relative_material_property_window_search.py | screen | screening_only |
| PROP-08 | Exact joint classification | relative_material_property_window_search.py | main | final_evidence_metric |
| NC-01 | Nonnegative dissipation partition | topology_constrained_sintering.py | solve_dissipation_partition | superseded |
| NC-02 | Aggregate smoothing redistribution | topology_constrained_sintering.py | surface_smoothing_redistribution | superseded |
| NC-03 | Density efficiency | topology_constrained_sintering.py | run | diagnostic_only |

## Manuscript non-claims after QC

- Candidate 693168 is **conditional Tier B, not validation** and not a calibrated Tier-A material model.
- The large attained high-temperature/two-step grain-size separation is **not inherently unphysical**; its magnitude remains an experimental-scale prediction requiring calibration.
- The modeled closed-pore/accommodation trajectory is the primary calibration and falsification target.
- Surrogate and screening-only equations are not final evidence; exact-promoted rows control classifications.
- No hidden closed-pore Lambda/K law was implemented. Closed accommodation is an implemented bounded proxy, not a derived closed-pore Poisson or gas-transport law.
