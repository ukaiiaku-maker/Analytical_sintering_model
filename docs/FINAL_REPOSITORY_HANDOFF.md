# Final repository handoff

## Current final scientific status

Fast firing and two-step sintering are compatible in the current reduced
framework, but the exact evidence attributes them to different dominant
mechanisms. Fast firing is controlled primarily by nucleation-limited material
kinetics: the fast ramp crosses the low-activity interval with less
pre-densification grain growth. Two-step behavior in the current best candidate
is controlled primarily by PR-prepared, bounded closed-pore accommodation
memory, which creates a density-exhaustion lower boundary and a grain-growth
upper boundary.

Candidate 693168 is **conditional Tier B**, not validation and not a paper-ready
calibration. No model physics, topology parameters, material parameters,
candidate IDs, density targets, classifications, time budgets, or final
evidence rules were changed for the final synthesis, equation audit, source-data
package, or this handoff.

## Exact final counts

| Quantity | Count | Evidence status |
|---|---:|---|
| Screened rows | 50,655 | surrogate screening/promotion |
| Unique exact cases | 1,903 | exact promoted union |
| Exact fast-only | 485 | exact final classification |
| Exact two-step-only | 119 | exact final classification |
| Exact both-pass | 73 | exact final classification |
| Exact neither | 1,226 | exact final classification |
| Surrogate both | 19,880 | screening only; not final evidence |

Candidate 693168 status: conditional Tier B; not validation.

## Branch provenance

The verified local/remote branch heads, purposes, and controlling reports are
listed in [FINAL_GITHUB_BRANCH_INDEX.md](FINAL_GITHUB_BRANCH_INDEX.md).

## Key reports

- [Final mechanism synthesis](FINAL_FAST_FIRING_AND_TWO_STEP_MECHANISM_SYNTHESIS.md)
- [Final synthesis caption drafts](FINAL_MECHANISM_SYNTHESIS_CAPTIONS.md)
- [Publication-style candidate figures](PUBLICATION_STYLE_SINTERING_FIGURES_693168.md)
- [Candidate 693168 closed-accommodation audit](CANDIDATE_693168_CLOSED_ACCOMMODATION_AUDIT.md)
- [Tier-B experimental-plausibility reframe](TIERB_EXPERIMENTAL_PLAUSIBILITY_REFRAME.md)
- [Relative material-property window](RELATIVE_MATERIAL_PROPERTY_WINDOW_REFRAMED.md)
- [Equation functional-form audit](EQUATION_FUNCTIONAL_FORM_AUDIT_FOR_PAPER.md)
- [QC Methods text](METHODS_TEXT_WITH_EQUATIONS_FOR_PAPER_QC.md)
- [QC SI equation tables](SI_EQUATION_TABLES_AND_VARIABLE_DEFINITIONS_QC.md)
- [Figure source-data package README](FIGURE_SOURCE_DATA_PACKAGE_README.md)

## Key source-data packages

- Figure source-data package: `results/figure_source_data_package/`
- Compressed figure package: `results/figure_source_data_package_20260814.zip`
- Equation audit: `results/equation_functional_form_audit/`
- Final synthesis tables: `results/final_mechanism_synthesis_and_property_windows/source_tables/`
- Candidate 693168 tables: `results/figure_source_data_package/candidate_693168/`
- Fast-firing tables: `results/figure_source_data_package/fast_firing/`
- Complete Chen-map tables: `results/figure_source_data_package/chen_maps/`
- Six Tier-B family tables: `results/figure_source_data_package/figure_06_six_tierB_family/`

The package provides complete source data for 33 figure references through 52
scientific datasets and 61 CSV source/metadata tables. The archive is about
8.3 MB. No journal aesthetics are encoded.

## Final evidence versus screening

Exact-promoted rows control all final pass/fail counts. The surrogate screen is
only a promotion tool; its 19,880 `both` rows overpredict the 73 exact both-pass
cases and are not final evidence. Candidate 693168 remains conditional Tier B.
Its large high-temperature/two-step grain-size separation is not automatically
an artifact, but its magnitude is not experimentally calibrated.

The primary calibration and falsification targets are the closed-pore fraction
and bounded accommodation trajectory during interrupted first and second steps.
The accommodation law is an implemented bounded proxy, not a derived
closed-pore Poisson or gas-transport law. No hidden closed-pore Lambda/K law is
claimed.

## Regenerating key outputs

Use the project virtual environment where dependencies are installed:

```bash
python3 final_mechanism_synthesis_and_property_windows.py
python3 final_mechanism_synthesis_plots.py
python3 equation_audit_qc.py
python3 build_figure_source_data_package.py
python3 audit_figure_source_data_package.py
python3 -m pytest -q -m "not requires_archived_results"
python3 -m py_compile *.py
```

The first two commands reproduce existing final synthesis outputs and may be
computationally more expensive. The figure-source builder itself performs no
simulation and reads frozen compact evidence. The publication/handoff task did
not rerun the synthesis or any search.

## Archive policy and known limitations

- No validation claim is made.
- No paper-ready calibration is claimed.
- The closed-pore/accommodation state is not experimentally calibrated.
- Historical result folders are not restored by default; compact archives are
  used where appropriate.
- Tests requiring archived historical fixtures are marked
  `requires_archived_results`; current tests use
  `-m "not requires_archived_results"`.
- Candidate histories do not expose TJ-specific Lambda/K, `P_comp_TJ`, `X_J`,
  separately named residual stress, trapped-gas pressure, or absolute-energy
  channels. Missing values are documented rather than invented.
- Two existing repository files exceed 50 MB: the historical backup ZIP and a
  previously tracked local-region history table. Neither is part of this
  handoff commit.
- The 862 intentional historical-result deletions remain unstaged and are not
  mixed into publication or science commits.
