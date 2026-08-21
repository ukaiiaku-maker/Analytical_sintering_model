# Final figure source data

Every plotted main or supplemental figure is linked to one CSV by `results/zro2_forward_final_summary_figures/final_figure_inventory.csv`. The figure-specific tables in `source_tables/` are direct assemblies of the selected reconstructed histories and existing classification results. `source_file_resolution_manifest.csv` records how requested process-search filenames map to the available processing-window result filenames.

`final_heating_rate_histories.csv` and its panel/source derivative contain P001 at 0.2, 1, and 100 °C min⁻¹. `final_twostep_histories.csv` and its derivative contain P014's identical first-step state followed by the 850, 925, 1025, and 1500 °C paths with continuous physical time. `final_chen_map_source.csv` preserves the existing classification points and boundary fields. `missing_history_fields.csv` explicitly records unavailable requested history fields; no missing field was fabricated.

The extraction changes no physics and performs no fitting, mobility optimization, or broad search. Dense selected histories use the same stored parameters. Candidate 693168 is not used as a ZrO2 parameterization. The model is not validated.
