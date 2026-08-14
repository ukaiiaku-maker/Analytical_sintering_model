# Visual Inspection Candidate Plots

This package is a plotting-only assembly built from the strict-tier tables and frozen E0021/E0142 parameters. It selects all five unique E0142 Tier B windows, the widest/smallest-G1 E0142 representative, and the widest E0021 Tier C comparison. Missing representative histories were rerun only for those frozen paths.

The intended inspection is visual rather than inferential: fast-firing panels show time, temperature, density, grain size, serial kinetic times, activity, and ablations; Chen panels show distinct failure classes, filled success bands, q0/q1 comparisons, and a lower/success/upper triplet.

E0142 is labeled Tier B, never Tier A. E0021 remains Tier C despite its cleaner fast-firing timing signature. Unavailable decomposed power, stress, and topology-coverage fields are listed in `visual_inspection_missing_data.md` and were not fabricated.
