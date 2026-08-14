# Supporting Information structure for the sintering mechanism model

## SI overview

The Supporting Information should preserve the evidence hierarchy and the progression from negative controls to exact-promoted conditional candidates. The main text can present the mechanistic conclusion and compact property windows; the SI should show how alternatives failed, how candidate decoding and exact promotion were controlled, and which results remain calibration-sensitive.

## SI Section S1 — Model hierarchy

### S1.1 State variables and observables

Document density, mean/distributional grain size, pore-bin volume and number, connected/open/closed pore inventories, pore-location fractions, junction state, stress/work variables, and accommodation stores. Distinguish physical observables from internal proxies.

### S1.2 Serial densification kinetics

Define nucleation, exchange, and transport times and their serial activity partition. Introduce `Theta_nuc`, the serial time fractions, low-activity exposure, and matched-density efficiency metrics.

### S1.3 Migration and growth layer

Document the migration-only closures, pore and junction drag, optional multihit constraints, residual stress, and the separation between densification and grain growth.

### S1.4 Pore redistribution and high-density support

Document conservative pore-bin redistribution, PR/surface evolution, connected-to-closed transition, closed-pore shrinkage, accommodation capacity/recovery, and exact pore-volume conservation.

### S1.5 Separate fast-firing and two-step evaluation layers

State explicitly that the final material perturbation vector was evaluated in the existing fast-firing and two-step layers without adding a hidden shared density channel. Joint pass means both exact evaluations pass, not that a new coupled model was fitted.

Suggested tables: state-variable definitions, flux registry, parameter provenance, observable/proxy mapping, and conservation identities.

## SI Section S2 — Negative-control progression

Organize the negative controls chronologically and state only the conclusions established in the existing reports.

| Stage | Negative control or limitation | Established lesson | SI source |
|---|---|---|---|
| Lambda-only optimization | `r_nuc * tau_sink` can rise in coarse grains while densification efficiency falls | renewal activity alone is insufficient | architecture and mechanism-option reports |
| Empirical topology damage | flips the fast-heating sign but changes inferred topology more than observable pore distribution | schedule memory must be linked to observable state | topology-memory stress report |
| Conservative pore-bin redistribution | creates schedule-dependent pore distributions and a density-window crossover | pore redistribution can supply observable ramp memory | pore-bin and density-window reports |
| Density smoothing gate | crossover follows `rho0 - smoothing_rho_mid` | the density gate is a useful proxy but is not uniquely identifiable | smoothing-gate identifiability report |
| Observable topology gates | tests fine-pore/connectivity replacements | topology gates reduce direct density dependence but remain conditional | topology-gate identifiability report |
| Fixed-model Chen maps | practical windows appear mainly at coarse grain size | baseline growth scaling does not create the desired nanoscale window | expanded phase-space report |
| Junction-limited migration | shifts the size onset but does not create a robust strict nanoscale window | migration suppression is necessary but not sufficient | nanoscale growth suppression report |
| Pore/junction pinning and placement | topology moves boundaries but does not alone determine the size onset | pore location must be treated as state-resolved, not a scalar multiplier | pinning and pore-placement reports |
| Separated-mechanism studies | fast firing survives PR-off and fails under nucleation-facile ablation | fast firing is primarily nucleation limited | separated production reports |
| Closed-store local-region candidates | finite two-step windows emerge with closed support and upper growth activation | closed accommodation is the leading two-step hypothesis | decoder-corrected and 693168 audit reports |

Do not present failed models as obsolete data. Present them as falsification steps that narrowed the admissible mechanism structure.

## SI Section S3 — Decoder correction and exact promotion

### S3.1 Decoder audit

Explain the sample-to-parameter mapping correction and fingerprint checks. Provide the deterministic seed, candidate reconstruction procedure, and parameter-vector hash where available.

### S3.2 Search funnel

Report the established counts: 1,000,000 decoded vectors, 20,000 retained Stage-1 dynamic fingerprints, 1,000 exact reconfirmations, 184 provisional Tier-B cases, and six candidates surviving the preparation audit. No Tier-A candidate was found.

### S3.3 Exact trajectory rules

Document common density targets and time budgets, cloned first-step state transfer, adaptive integration, stagnation rejection, matched-density interpolation support, and the requirement for both lower and upper Chen boundaries.

### S3.4 Rejection provenance

Summarize preparation-growth, attainment, censoring, missing-boundary, and physical-state rejection counts. Retain a reason for every rejection.

Primary sources: `docs/LOCAL_REGION_DECODER_CORRECTED_DYNAMIC_SEARCH.md`, `docs/LOCAL_REGION_DECODER_AUDIT.md`, and `results/local_region_decoder_corrected_dynamic_search/`.

## SI Section S4 — Candidate 693168 audit

### S4.1 Parameter reconstruction and timestep convergence

Report the deterministic fingerprint and 30-, 15-, and 5-minute maximum-step comparison. State the small boundary shifts and preservation of the switch state.

### S4.2 Absolute trajectories

Show high-temperature and two-step `G(rho)` on absolute and logarithmic axes. Over `rho = 0.95–0.98`, report high-temperature grain sizes of approximately 858–1123 nm and two-step sizes near 118 nm, with both paths attained.

### S4.3 Closed-pore/accommodation trajectory

Show closed fraction, open/closed shrinkage contributions, accommodation capacity/availability/use, and PR/closed-transition history. State that model-store fractions are not stereological porosity.

### S4.4 Fine Chen map

Show density-exhaustion failures, the 925–1205 °C success band, and grain-growth failures. Retain the fixed target and growth tolerance.

### S4.5 Causal ablations

Show that no PR damage, no closed transition, no closed shrinkage, and infinite accommodation destroy the joint result. Report retained non-primary channels without turning survival into proof of irrelevance.

Primary sources: `docs/CANDIDATE_693168_CLOSED_ACCOMMODATION_AUDIT.md`, its final figure manifest, and `results/audit_candidate_693168_closed_accommodation/`.

## SI Section S5 — Six-candidate Tier-B family

Present candidates 693168, 822940, 581668, 295003, 366094, and 85161 in a common table. Include `G0`, `G1`, preparation growth, median reduction, window width, success range, closed fraction at switch/target, closed-shrinkage share, robustness count, and destructive ablations.

Use the family to show that exact Tier-B behavior is not unique to one topology. Preserve the distinctions: 693168 is the strongest separation comparator, 822940 provides a low-closed-fraction case, and 581668 provides an intermediate topology. Do not imply that all six share an identical causal hierarchy or an exact material-property robustness audit.

Primary table: `results/final_mechanism_synthesis_and_property_windows/source_tables/six_TierB_family_mechanism_summary.csv`.

## SI Section S6 — Relative material-property window

### S6.1 Perturbation and promotion design

Report 50,655 screened rows, including 50,000 Latin-hypercube points and bounded OAT/pairwise/diagnostic rows. Clearly label screening metrics as surrogate/dimensionless feasibility evidence.

### S6.2 Exact classifications

Report 1,000 exact fast promotions, 1,000 exact two-step promotions, and 1,903 unique exact cases. Final counts are 485 fast-only, 119 two-step-only, 73 both-pass, and 1,226 neither.

### S6.3 Local OAT windows

Report:

- `Delta Q_nuc = 0…+50 kJ/mol` survives; `-25` and `+75 kJ/mol` fail.
- `Delta Q_closed = -25…+100 kJ/mol` survives; the lower boundary disappears at `-50 kJ/mol`.
- `Delta Q_growth` remains viable throughout the tested `±100 kJ/mol`; its limit was not found.
- PR prefactor threshold: `0.3x` base.
- Growth prefactor threshold: `0.1x` base to retain the upper boundary.

### S6.4 Joint dimensionless envelope

Present exact observed ranges with their coverage-limited caveat. Do not convert the scatter envelope into independent rectangular parameter tolerances.

### S6.5 Surrogate warning

Show the 19,880 surrogate both predictions beside the 73 exact both-pass cases. State that the 272.3-fold discrepancy prevents surrogate rows from supporting final mechanistic claims.

Primary sources: the final synthesis source tables and `docs/MATERIAL_PROPERTY_WINDOW_MECHANISM_ATTRIBUTION.md`.

## SI Section S7 — Figure and table inventory

### Main-text candidates

1. Mechanism chain: fast firing versus two step.
2. Exact relative-property phase map.
3. Fast-firing property window.
4. Two-step property window.
5. Fine Chen boundary mechanism map.
6. Six-candidate Tier-B family summary.
7. Surrogate-versus-exact warning.
8. Experimental falsification targets.

### Candidate-693168 presentation figures

- final multi-panel result;
- mechanism schematic;
- physical-time histories;
- fine Chen filled/classification maps;
- closed-accommodation and PR/transition histories;
- ablation waterfall;
- six-candidate comparison.

### Core machine-readable tables

- `final_property_window_summary.csv`;
- `exact_behavior_classification_counts.csv`;
- `fast_firing_OAT_window.csv`;
- `two_step_OAT_window.csv`;
- `dimensionless_group_thresholds_final.csv`;
- `six_TierB_family_mechanism_summary.csv`;
- `surrogate_vs_exact_comparison.csv`;
- `experimental_falsification_targets.csv`;
- candidate-693168 fine T2 classification, boundary, ablation, and history tables.

## SI Section S8 — Falsification measurements

Organize the experimental section around measurements that distinguish mechanisms rather than merely fit shrinkage curves:

1. matched-density `G_mean/G50/G90` under multiple ramps and two-step schedules;
2. interrupted-ramp densification onset;
3. exchange and transport relaxation;
4. 3D open/closed pore fraction;
5. pore D50/D90 and large-pore tail;
6. connected fine-pore fraction and percolation;
7. trapped-gas or accommodation proxy;
8. grain-growth mobility across the upper T2 boundary;
9. interrupted first/second-step tomography.

End the SI with a calibration table mapping each model state to its measurement, uncertainty, fitted/not-fitted status, and falsification criterion. Preserve the explicit statement that the current work is conditional mechanism attribution, not validation.
