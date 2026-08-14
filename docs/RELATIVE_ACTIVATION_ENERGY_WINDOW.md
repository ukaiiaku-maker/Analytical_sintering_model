# Relative activation-energy window

## Interpretation

The audit perturbs barriers around the current base values while topology remains fixed. Reported exact envelopes are coverage-limited observed ranges, not universal constants. OAT results give local identifiability; Latin-hypercube successes show conditional combinations.

## OAT windows

| Parameter shift | Fast firing | Two step | Interpretation |
|---|---|---|---|
| `Delta Q_nuc` | 0 to +50 kJ/mol pass; -25 and +75 fail | insensitive in separated layer | finite nucleation-onset window |
| `Delta Q_exchange` | -75 to +50 pass; +75 fails | insensitive | completion cannot become too slow |
| `Delta Q_transport` | full ±75 range passes | insensitive | weak locally in this envelope |
| `Delta Q_growth` | -100 to 0 pass; +25 fails | full ±100 range passes | fast timing and two-step upper boundary constrain different aspects |
| `Delta Q_PR` | full ±100 range passes | full ±100 range passes | barrier poorly identifiable locally; prefactor is decisive |
| `Delta Q_closed` | insensitive | -25 to +100 pass; -50 loses lower boundary | lower boundary constrains closed kinetics |

The largest allowable local `Q_nuc` shift is +50 kJ/mol from base before the next tested point fails. The two-step local `Q_closed` range is -25 to +100 kJ/mol. The two-step `Q_growth` range spans the full tested ±100 kJ/mol, so its barrier limit was not found; however, reducing its prefactor to 0.03× removes the upper boundary.

## Joint exact envelope

Seventy-three exact perturbations pass both separate layers. Their observed relative ranges include:

- `Q_nuc-Q_growth`: -50 to +96 kJ/mol;
- `Q_nuc-Q_transport`: +15 to +202 kJ/mol;
- `Q_nuc-Q_PR`: +51 to +251 kJ/mol;
- `Q_closed-Q_growth`: -252 to -127 kJ/mol;
- `Q_PR-Q_closed`: -74 to +212 kJ/mol;
- `log10(k_closed/k_growth)`: -1.52 to +1.50;
- `log10(k_PR/k_growth)`: -1.48 to +1.48.

These ranges are conditional combinations from a promoted subset. They must not be read as independent rectangular tolerances.

## Dimensionless thresholds

Joint successes occupy `Theta_nuc = 23.4–2.96e5` and `S_closed_growth = 0.116–126`, with medians 1.07e3 and 3.85. The breadth confirms that the outcomes depend on combinations and histories, not a unique scalar threshold.

The complete numerical tables are `dimensionless_thresholds.csv` and `relative_activation_energy_window.csv`.
