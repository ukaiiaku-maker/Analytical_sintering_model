# Material-property window mechanism attribution

## Status and scope

This is a **conditional mechanism-attribution audit, not validation or paper-ready calibration**. The topology and geometric parameters of candidate 693168 were frozen. A common relative material-property perturbation vector was evaluated in two existing, separate model layers: the E0021 fast-firing model and the candidate-693168 two-step model. Therefore, `both_pass` means that the same perturbation vector passed both exact evaluations; it does not imply a newly coupled state trajectory.

The campaign screened 50,655 designs (50,000 Latin-hypercube points plus bounded OAT, pairwise, and diagnostic rows). It exactly evaluated 1,000 fast-firing promotions and 1,000 two-step promotions. Their union contains 1,903 points because 97 points were promoted by both rankings. Exact classification gave 485 fast-only, 119 two-step-only, 73 both-pass, and 1,226 neither cases. The screen had predicted 19,880 both-pass cases, so the reduced screen is useful for promotion but is not itself mechanistic evidence.

## Main attribution

Fast firing is controlled primarily by a nucleation-limited waiting interval relative to exchange, transport, and non-densifying grain growth. The full base path reaches a maximum matched-density ratio of 1.796 over a density span of 0.17. Among exact fast promotions, 55.8% pass the full rule, only 19.3% still pass the nucleation-facile ablation, and 71.7% pass with PR disabled. Thus PR is not the causal fast-firing channel in this envelope.

Candidate 693168's two-step response is controlled by a different chain: PR prepares a closed-pore store, closed shrinkage/accommodation sustains high-density densification at the lower second-step temperatures, and thermally activated growth supplies the upper boundary. The exact base reproduction gives a reduction of 0.906 over a 0.03 high-density span and a 250 °C bracketted Chen window. Its large grain-size separation is not flagged as an artifact because both paths attain the interval and the stored states remain bounded.

The shared contributor is high-density/matched-density attainment with bounded states. Nucleation waiting is specific to fast firing. PR-prepared closed topology, finite accommodation, and closed shrinkage versus migration are specific to the two-step result.

## Relative-property windows

OAT evidence is the cleanest local attribution:

- Fast firing survives `Q_nuc` shifts of 0 to +50 kJ/mol around the base. It fails at -25 kJ/mol and at +75 kJ/mol. This is a finite, nonmonotonic onset window, not permission to quote ±100 kJ/mol as an allowable uncertainty.
- Fast firing survives `Q_growth` shifts from -100 through 0 kJ/mol and fails at +25 kJ/mol. `Q_exchange` survives through +50 kJ/mol but fails at +75 kJ/mol. The full tested `Q_transport` OAT range (±75 kJ/mol) passes.
- Two-step behavior survives `Q_closed` shifts from -25 through +100 kJ/mol. At -50 kJ/mol and below the apparent density response loses the required lower boundary, so it is not a Chen window.
- Two-step behavior survives the full tested `Q_growth` activation-energy OAT range (±100 kJ/mol), but a growth prefactor of 0.03× loses the upper boundary; 0.1× through 30× retain both boundaries.
- The PR prefactor must be at least 0.3× in the tested OAT. At 0.03× and 0.1× the full two-step rule fails. The full tested PR activation-energy OAT range (±100 kJ/mol) passes, indicating prefactor/topology production is more identifiable than its barrier in this parameterization.
- The closed-shrinkage prefactor passes from 0.03× through 30×, although reduction and window width move. This broad local result does not remove the independent concern about closed-accommodation capacity and recovery calibration.

Across the 73 joint exact cases, the observed, coverage-limited ranges are:

| Group | Minimum | Median | Maximum |
|---|---:|---:|---:|
| `Q_nuc - Q_growth` (kJ/mol) | -50.0 | -50.0 | 96.0 |
| `Q_nuc - Q_transport` (kJ/mol) | 15.0 | 90.0 | 201.8 |
| `Q_closed - Q_growth` (kJ/mol) | -251.8 | -226.8 | -126.8 |
| `log10(k_closed/k_growth)` | -1.523 | 0.000 | 1.498 |
| `log10(k_PR/k_growth)` | -1.477 | 0.000 | 1.477 |
| `Theta_nuc` | 23.4 | 1.07e3 | 2.96e5 |
| `S_closed_growth` | 0.116 | 3.85 | 126 |

These are the ranges sampled among promoted exact successes, not universal necessary-and-sufficient bounds.

## Six-candidate family

All six exact Tier-B candidates are preserved in the fixed registry. The family comparison is deliberately a **reduced transfer from exact 693168 OAT**, not six new exact searches. Under that proxy, 693168 is the most robust comparator. Candidate 581668 retains a moderate transferred margin but remains calibration-sensitive; 822940 has nine destructive-ablation robustness counts yet a narrower property margin. Candidates 295003, 366094, and 85161 have weaker transferred margins. Candidate 693168 remains representative of the strongest conditional closed-store mechanism, not representative of a quantitatively calibrated material class.

## Calibration and falsification needs

The most important measurements are independent nucleation onset times during ramps, exchange/transport relaxation times, closed-pore fraction at the first-step switch, pore-size/location distributions, closed-pore shrinkage rate versus temperature, accommodation capacity/recovery, and grain-growth mobility over the same T2 interval. Interrupted first-step microscopy plus in-situ or interrupted high-density pore evolution would distinguish PR-prepared storage from a fitted state variable.

The principal remaining concern is the closed-store/accommodation trajectory. The base candidate reaches closed fractions near unity at the target and draws almost all high-density densification from closed shrinkage. That behavior is a falsifiable conditional hypothesis, not validation.

Primary data are in `results/relative_material_property_window_attribution/source_tables/`; exact checkpoints and the compressed full screen are retained with the run state.
