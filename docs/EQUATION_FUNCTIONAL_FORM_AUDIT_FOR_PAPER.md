# Equation and functional-form audit for the paper

## Scope and result

This is a static audit of the implemented equations that generated the final synthesis. It executes no trajectory, changes no parameter, and introduces no physics. The generated registry contains 67 equations, 54 variable definitions, and 46 parameter definitions. Twenty-seven requested source files were reviewed; 12 directly furnish registered equation rows, while the remaining search and plotting files provide orchestration, classification use, or result provenance. The single requested historical path not present as one script is the local-region decoder-corrected dynamic search; its equations are recoverable from interacting_local_region_decoder.py, interacting_local_region_model.py, local_region_decoder_corrected_postprocess.py, and the candidate-693168 audit.

The audit separates four evidence classes:

1. Final-evidence closures: the serial material law, candidate-693168 local closed-store law, matched-density metrics, and exact classification logic.
2. Causal ablations: the PR-off and nucleation-facile modifications used to attribute fast firing.
3. Diagnostic and negative-control closures: aggregate partitioning, pore placement, action allocation, persistent-junction drag, and TJ multihit.
4. Screening-only equations: the Latin-hypercube surrogate response scores. They are not final evidence.

Machine-readable provenance is in results/equation_functional_form_audit/.

## A. Serial nucleation, exchange, and transport

Registry equations FF-01 through FF-10 contain the exact fast-material expressions. Stress-assisted nucleation, exchange, and \(G^2\)-scaled transport are summed in a serial cycle. Renewal activity is the fraction of that cycle not spent waiting for nucleation. Connected density and fine-pore fractions provide the explicit geometry eligibility, and density gain per completed cycle gives \(\dot\rho\).

The controlling final interpretation is timing-based: the slow path accumulates more grain growth during a low-activity interval, whereas the fast path crosses it quickly and still activates soon enough to attain the scored density interval. The intrinsic growth law is independent of the density law except through the shared instantaneous state.

The no-nucleation-limitation ablation makes the nucleation time \(10^{-6}\) of the exchange-plus-transport time. The exchange-limited and transport-only variants alter the serial structure explicitly and are retained as causal alternatives, not final material equations.

## B. Fast-firing observable

MET-01 and MET-02 restrict interpolation to the density interval attained by both paths. The controlling criterion is \(G_{\rm ref}/G_{\rm fast}\ge1.5\) continuously over at least 0.03 density. Censored and unattained intervals are excluded.

The exact base reaches 1.796 over a 0.17 span. Among 1,000 exact fast promotions, 55.8% pass the full rule, 19.3% pass after nucleation is made facile, and 71.7% pass with PR disabled. The final fast-firing role is therefore assigned to nucleation-limited onset, not PR.

## C. PR/de-sintering and conservative pore memory

FF-09 and FF-12 are the simple conservative surface redistribution used in the separated fast model. PR-01 through PR-09 describe the more explicit production competition closure:

- a logistic fine-radius weight;
- a low-renewal gate;
- a reference-temperature Arrhenius factor;
- connected fine-pore topology eligibility;
- conservative smoothing, GB-to-TJ, and TJ-to-isolated partitions;
- a nonnegative densification/PR competition partition;
- separate densifying and non-densifying work proxies.

PR does not directly remove pore volume. In final attribution it is a preparation channel for the candidate-693168 closed state, not the controlling fast-firing channel.

## D. Pore location and action allocation

PL-01 through PL-06 trace the GB-segment, TJ, and isolated stores, the density identity, exponential coverage laws, location-resolved densification eligibility, and Class-A pinning. AL-01 through AL-03 trace the nonnegative local action scores and the three constrained action pairs:

- GB-segment versus TJ shrinkage;
- GB smoothing versus GB-to-TJ relocation;
- TJ-to-GB capture versus TJ isolation.

These layers are scientifically useful negative controls. They made topology observable and conservative but were insufficient as the final joint mechanism.

## E. Persistent junction and TJ multihit

TJ-01 through TJ-08 document junction-state production, Arrhenius relaxation, persistent drag, required packet count, expected event count, exact Poisson completion, pore/structural TJ separation, and coupled mobility. The q_TJ=0 and q_TJ=1 forms remain explicit.

The local-region model does not reuse exact Poisson completion. Its LR-09 equation is a sigmoid normal-approximation-style proxy. The two forms must not be conflated.

Persistent junction plus TJ multihit produced earlier complete Chen windows. In candidate 693168, removing TJ multihit, attached-pore drag, or persistent junction still retains a complete window. Those closures are not its primary cause.

## F. Decoder-corrected local-region and closed accommodation

LR-01 through LR-12 are the decisive two-step equations. Each local region contains GB-segment, TJ, isolated, and closed pore stores. Open shrinkage removes GB-connected volume. PR and sweep terms redistribute pore volume conservatively. A conservative isolated/TJ-to-closed transition prepares the closed store. Closed shrinkage is proportional to closed volume, a relative-temperature factor, available accommodation, a gas-ratio factor, and a code-normalized radius proxy.

Available accommodation recovers toward a finite capacity and is consumed by closed-pore removal. This state is an implemented bounded proxy. There is no derived closed-pore Poisson event model and no exact closed \(\Lambda/K\) equation. Candidate audit columns named P_comp_closed are the normalized accommodation fraction, while Lambda_closed and K_closed remain unavailable.

Candidate 693168 has 13.71% first-step growth and 64.94% of remaining pore volume in the closed store at the switch. Its fine exact map contains 101 T2 points and a complete 925–1205 °C success band; timestep-controlled summaries give 930–1190/1200 °C depending on the maximum step. The mechanism is conditional Tier B because the large closed-store/accommodation trajectory is not experimentally calibrated.

## G. Chen classification and adaptive boundaries

MET-03 through MET-06 document growth fraction, mutually exclusive point classes, complete-window requirements, and matched-density two-step reduction. Complete means both lower exhaustion and upper growth brackets, a finite success band, practical \(T_2<T_1\), and width at least 25 °C. Failed first steps, target-already-reached states, and censored bounds are not promoted.

The adaptive search begins with a coarse grid, extends downward if density is already attained at the lowest point, extends upward until an upper failure is found, and uses 10 °C refinement inside changing 25 °C brackets.

## H. Relative-property perturbations and exact promotion

PROP-01 through PROP-08 document additive activation-energy shifts, multiplicative prefactors, the Latin-hypercube construction, dimensionless summaries, screen surrogates, and exact classification. The screen evaluated 50,655 rows and predicted 19,880 both-pass cases. The final exact union contains 1,903 rows: 485 fast-only, 119 two-step-only, 73 both-pass, and 1,226 neither. The 272.3-fold both-pass discrepancy prevents screen equations from being evidence.

The exact local windows are:

- \(\Delta Q_{\rm nuc}=0\) through \(+50\) kJ mol\(^{-1}\) survives; \(-25\) and \(+75\) fail.
- \(\Delta Q_{\rm closed}=-25\) through \(+100\) kJ mol\(^{-1}\) survives; at \(-50\), the lower boundary is absent.
- The tested \(\Delta Q_{\rm growth}=\pm100\) kJ mol\(^{-1}\) remains viable; its limit was not found.
- The tested PR prefactor threshold is 0.3 times base.
- The tested growth prefactor threshold for retaining an upper boundary is 0.1 times base.

These are observed model windows, not universal constants.

## Mechanism-to-result conclusion

The final fast-firing equations are FF-01 through FF-10 plus MET-01 and MET-02. The final conditional two-step equations are LR-01, LR-03, LR-05 through LR-12, and MET-03 through MET-06. PROP-01, PROP-02, and PROP-08 define the exact property attribution. All other families remain visible as ablations, diagnostics, negative controls, or superseded architecture.

No validation is claimed. No model physics changed in this audit.
