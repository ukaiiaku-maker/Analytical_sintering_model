# Figure source data package

Built from frozen evidence at `codex/figure-source-data-package@579367dddf34b5d6b167e4cb0f0b81157f6b1ad7`. No simulation, search, parameter
tuning, classification change, or model-physics import is performed by the builder.

The package contains 52 figure/source tables plus citation, panel,
provenance, dictionary, missing-channel, and QC metadata. Final synthesis claims use
only `evidence_level=exact`. Surrogate rows are explicitly labeled
`surrogate_screen` and are promotion evidence only.

Candidate 693168 remains conditional Tier B, not validation. Its bounded
closed-pore accommodation trajectory is the primary calibration target. Large
high-temperature/two-step grain-size separation is not inherently unphysical, but
the magnitude remains uncalibrated. No closed-pore Poisson/gas-transport law and no
hidden closed-pore Lambda/K law are claimed.

## Loading in Python

```python
import pandas as pd
df = pd.read_csv("results/figure_source_data_package/figure_02_exact_property_phase_map/exact_property_phase_map.csv")
dense = pd.read_csv("results/figure_source_data_package/candidate_693168/dense_time_histories.csv")
```

Use `all_columns_dictionary.csv` for units/labels and `plotting_recommendations.csv`
for panel mappings. No journal aesthetics are hard-coded.
