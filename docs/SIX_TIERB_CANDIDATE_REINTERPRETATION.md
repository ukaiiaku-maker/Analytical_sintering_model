# Reinterpretation of the six exact Tier-B candidates

The family contains candidates 693168, 822940, 581668, 295003, 366094, and 85161. All have exact high-density attainment, complete finite Chen windows, and recorded causal ablations. None is classified as an artifact solely because of reduction magnitude.

- **693168:** calibration-sensitive Tier B; strongest large experimental-scale separation, broad window, stable timestep audit, but high closed-store support and 13.7% first-step growth.
- **822940:** plausible Tier B; moderate reduction, low switch closed fraction, broad initial-condition robustness, and useful comparison against closed-dominated solutions.
- **581668:** calibration-sensitive Tier B; moderate reduction and five robustness passes, with a near-saturated target closed store.
- **295003:** calibration-sensitive Tier B; low preparation growth and broad window, but near-saturated target closed store and no retained bounded robustness passes in the archived neighborhood.
- **366094:** calibration-sensitive Tier B; lowest median reduction, low closed fraction, nearly 20% preparation growth, and limited robustness coverage. It is a useful open-shrinkage comparator.
- **85161:** calibration-sensitive Tier B; high switch closed fraction, near-saturated target store, and a narrower window.

Closed fraction alone does not rank reduction or window width. The useful family diversity lies in its different closed/open shrinkage partitions, preparation growth, and causal migration-side ablations. Quantitative values are in `tierB_candidate_reinterpretation.csv` and `tierB_candidate_plausibility_scorecard.csv`.

No candidate is Tier A or validated. Independent timestep sweeps are available only for 693168; the other five are exact production reconfirmations but retain an explicit timestep-evidence gap.
