# Chen map plot data dictionary

- `UNATTAINABLE_FIRST_STEP`: requested switch state was not reached within the common 500 h budget.
- `TARGET_REACHED_DURING_FIRST_STEP`: density 0.98 was reached before a valid second step.
- `DENSIFICATION_EXHAUSTION_FAILURE`: target missed and second-step growth <=20%.
- `SUCCESS`: target reached and second-step growth <=20%.
- `GRAIN_GROWTH_FAILURE`: target reached but growth >20%.
- `MIXED_FAILURE`: target missed and growth >20%.
- `NUMERICAL_CENSOR`: nonfinite integration output.

All tiles retain failures; no schedule-specific parameters or filtering are used.
