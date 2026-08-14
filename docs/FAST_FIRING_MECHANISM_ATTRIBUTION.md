# Fast-firing mechanism attribution

## Result

The current fast-firing effect is a nucleation-limited timing effect. A slow ramp spends longer in a low-activity interval before densifying events activate, allowing more non-densifying growth per density increment. A fast ramp crosses that interval rapidly and subsequently densifies once exchange and transport can complete events.

This remains a **model attribution, not validation**. The exact base result has `R_fast = 1.796` and a 0.17 density span above the 1.5 ratio threshold. Across 1,000 exact fast promotions, 558 pass the full rule. The nucleation-facile ablation passes only 193/1,000, while PR-off passes 717/1,000. PR-off survival is expected under the current interpretation and is not a rejection.

## Local material window

For frozen topology and OAT changes around E0021:

- `Q_nuc`: pass at 0, +25, and +50 kJ/mol; fail at -25 and +75 kJ/mol. The largest contiguous tested shift from the base is +50 kJ/mol.
- nucleation prefactor: pass from 0.1× through 3×; fail at 0.03× and 10×. Too slow loses attainable separation, while too facile removes the waiting interval.
- `Q_exchange`: pass through +50 kJ/mol and fail at +75 kJ/mol; all negative tested shifts pass.
- `Q_transport`: the full ±75 kJ/mol tested range passes locally.
- `Q_growth`: shifts from -100 to 0 kJ/mol pass; +25 kJ/mol fails.
- PR barrier and prefactor changes do not change the base fast-firing metric in this separated layer.

The leading exact rank associations for `R_fast` include the nucleation prefactor. `Theta_nuc`, the nucleation serial-time fraction, and low-activity exposure organize the response, but no single dimensionless group is sufficient across the promoted subset.

## Causal interpretation

Necessary in the present model:

1. A finite nucleation waiting regime.
2. A fast/slow difference in exposure to that regime.
3. Exchange and transport fast enough that both trajectories attain the scored density interval.
4. A competing non-densifying growth channel during slow exposure.

Not necessary in the present fast-firing layer:

- PR redistribution. It may coexist with the effect, but its removal does not automatically destroy it.
- Candidate-693168 closed-pore accommodation. That belongs to the separate two-step model layer.

## Measurements

Measure densification-event onset and pore/grain evolution during interrupted ramps at several heating rates. Independent estimates of nucleation activation, exchange, and transport times are more discriminating than a single fitted shrinkage curve. The defining falsification is whether a nucleation-facile material retains the same matched-density fast-firing separation; the model predicts that it should not.
