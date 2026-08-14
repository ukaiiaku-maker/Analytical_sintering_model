# Methods equations: fast firing and two-step sintering

## Equation hierarchy

This compact hierarchy links the paper Methods equations to the machine registry.

| Paper role | Registry IDs | Implemented content | Evidence status |
|---|---|---|---|
| Fast activation stress | FF-01 | capillary stress plus large-pore concentration | final |
| Serial kinetic times | FF-02–FF-06 | nucleation, exchange, transport, cycle, activity | final |
| Fast density law | FF-07–FF-08 | connected/fine eligibility and event density gain | final; gate is a proxy |
| Fast PR ablation | FF-09, FF-12 | non-densifying conservative redistribution | causal ablation |
| Intrinsic growth | FF-10 | inverse-grain-size Arrhenius growth | final |
| Fast trajectory criterion | MET-01–MET-02 | matched density, ratio, continuous span | final |
| Closed-store conservation | LR-01, LR-07–LR-08 | four stores and conservative transfers | final conditional two-step |
| Open/closed shrinkage | LR-03–LR-05 | explicit density channels | final conditional two-step |
| Local migration/growth | LR-09–LR-10, LR-12 | completion proxy, drag, inverse-size step | final conditional two-step |
| Finite accommodation | LR-11 | bounded recovery/consumption state | final conditional two-step; proxy |
| Chen classification | MET-03–MET-06 | growth fraction, classes, brackets, reduction | final |
| Property attribution | PROP-01–PROP-08 | perturbation, screening, exact promotion | exact rows final; screen not final |

## Fast-firing functional chain

\[
\sigma_{\rm loc}\rightarrow\tau_{\rm nuc},\qquad
(\tau_{\rm nuc},\tau_{\rm ex},\tau_{\rm tr})\rightarrow\tau_{\rm cyc},a,
\]
\[
(\rho,\phi_i,r_i,a)\rightarrow\dot\rho,\dot G,\qquad
[\rho(t),G(t)]\rightarrow R_{\rm fast}(\rho).
\]

The density equation contains no schedule label. Schedule history enters only by integration of instantaneous state and temperature. Nucleation-facile is an explicit alternative equation; PR-off sets only the conservative PR flux to zero.

## Two-step functional chain

\[
(T,\phi^{GB},a)\rightarrow J_{\rm PR}
\rightarrow(\phi^{GB},\phi^{TJ},\phi^{iso},\phi^{closed}),
\]
\[
(\phi^{closed},A_{\rm closed},T)\rightarrow\dot\rho_{\rm closed},\qquad
(T,G,\mathcal T)\rightarrow\dot G,
\]
\[
[\rho_2,G_2;G_1]\rightarrow
\{\text{exhaustion, success, growth failure, mixed}\}.
\]

PR and the closed transition move volume conservatively. Open and closed shrinkage are the only local density-removal channels. The second-step lower boundary is produced by failure of closed shrinkage/accommodation to attain target within the common budget; the upper boundary is produced by excessive grain growth.

## Migration closure classes

Class A uses \(\Gamma=(1+D)^{-1}\). Class B uses exact Poisson completion in agentic_mechanism_model. The decoder-corrected local model instead uses a sigmoid TJ-completion proxy and \(\Gamma=P_{TJ}/(1+D)\). Class C exchange-limited closure was registered historically but is not a final candidate equation. Coupled persistent-junction plus multihit equations are prior Chen mechanisms, not the primary candidate-693168 cause.

## Conservation statements

1. Adjacent-bin smoothing has equal source loss and destination gain.
2. GB-to-TJ, TJ-to-isolated, detachment/recapture, and open-to-closed transitions are conservative.
3. PR work does not directly reduce pore volume.
4. Density is reconstructed from total pore volume.
5. Migration-only closures modify \(\dot G\), not \(\dot\rho\), at a shared state.
6. All successful classifications require target attainment and bounded growth.

## Exact evidence

The property campaign screened 50,655 rows, exactly promoted 1,000 fast and 1,000 two-step rows, and produced 1,903 unique exact cases. Counts are 485 fast-only, 119 two-step-only, 73 both-pass, and 1,226 neither. These exact rows, not the 19,880 surrogate joint predictions, control the final classification.
