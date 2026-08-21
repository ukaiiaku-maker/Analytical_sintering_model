# Final ZrO2 summary figures

This branch regenerates final heating-rate and two-step summary figures only. It changes no model physics, barrier data, diffusivities, closed-pore laws, PR/coarsening laws, schedules, targets, or success definitions. No material parameters were fitted, no mobility optimization was used, and no new broad process-state search was run.

The selected cases come from the existing fixed-parameter processing-window result tables at source commit `255446c`. P001 provides the strong fast-firing heating-rate response. P014 provides the finite two-step window and the low-temperature failure, central strict success, and high-temperature failure trajectories. Selected paths were rerun only to reconstruct dense histories with the same stored parameters and selected first-step state; that reconstruction was not a search.

Candidate 693168 remains a prior mechanism reference only and is not a ZrO2 parameterization. These figures are publication-quality candidates, but the model is not validated.

The main deliverables are in `results/zro2_forward_final_summary_figures/figures_main`; supplements and figure-specific source CSVs are in the adjacent `figures_supplement` and `source_tables` directories. `final_figure_inventory.csv` links each rendered figure to its source table, and `final_figure_qc_report.csv` records the automated quality checks.
