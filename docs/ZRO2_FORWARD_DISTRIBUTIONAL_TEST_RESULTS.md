# Distributional model test results

## Limiting and fixed-state tests

All conservation and limiting tests passed. The full prescribed fixed-state grid contains 79,200 points spanning 700–1500 °C, six D50 values, five widths, four tail weights, four segment ratios, and five activities. It records regularization, damage, pinch-off, distribution-shape changes, Zener changes, and conservative topology production.

## Synthetic boundary test

The bounded forward test evaluated 2,016 trajectories: eight synthetic states, three representations, two physical closed kernels, `m=3,4`, and T2 from 850–1350 °C in 25 °C increments. It produced 36 strict success points but zero finite Chen windows.

Every success occurred at 850 and 875 °C and therefore had an upper boundary but no lower density-exhaustion boundary. The high-useful-closed, gas-limited, and candidate-response-like states were the successful families. Candidate-like states remain diagnostic and are not ZrO2 predictions.

The remaining points comprised 156 density-exhaustion failures, 625 mixed failures, and 1,199 grain-growth failures. Because the synthetic topology gate failed, the broad initial-state/process grid was not justified or run; failure-mode tables were generated instead.

## Heating-rate response

The bounded campaign ran 525 paths across three representations, five peak temperatures, seven heating rates, and five holds. Of 525 matched-density rows, 399 pairs attained the requested density. The maximum slow/reference-to-fast grain-size ratio was about 1.26 at rho=0.75 for the lognormal family. However, only 12 of the 399 attained comparisons had a ratio above one; 387 had the wrong sign, with the fast ramp producing more grain growth at matched density. The median ratio was 0.152. No density span reached ratios 1.5 or 2.0. Thus the tested distribution evolution does not reproduce the required heating-rate response.

## Ablations

For the best low-temperature synthetic success:

- removing closed shrinkage or imposing high gas pressure destroyed density attainment;
- disabling precursor-to-closed transfer reduced density but did not destroy success;
- disabling PR pinch-off reduced density but did not destroy success;
- regularization and damaging-coarsening ablations had negligible effect;
- distributional-Zener removal and mean-radius substitution had no effect at the low-temperature best case;
- infinite accommodation increased density and did not restore a lower boundary;
- removing the energy-ledger gate retained success.

The intended coupled distributional mechanism is therefore not causally supported. No validation claim is made.
