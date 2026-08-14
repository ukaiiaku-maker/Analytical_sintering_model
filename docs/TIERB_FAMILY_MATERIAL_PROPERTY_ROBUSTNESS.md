# Tier-B family material-property robustness

## Scope

All six exact Tier-B candidates are included: 693168, 822940, 581668, 295003, 366094, and 85161. This is a **reduced family transfer**, not an exact OAT or Latin-hypercube campaign for every candidate, and it is **not validation**. Exact material perturbations were run for the frozen 693168 topology; those survival margins were transferred according to each candidate's base reduction and Chen-window margin.

## Comparison

| Candidate | Base reduction | Window (°C) | Closed fraction at switch | Reduced conclusion |
|---:|---:|---:|---:|---|
| 693168 | 0.881 | 260 | 0.649 | robust comparator within this audit |
| 822940 | 0.319 | 125 | 0.022 | narrower, calibration-sensitive transfer |
| 581668 | 0.339 | 265 | 0.460 | moderate transfer, calibration-sensitive |
| 295003 | 0.423 | 290 | 0.621 | calibration-sensitive despite broad base window |
| 366094 | 0.231 | 215 | 0.147 | weak transferred margin |
| 85161 | 0.395 | 60 | 0.942 | narrow-window, weak transferred margin |

Candidate 693168 remains the most useful representative of the strong closed-store response, but not of a quantitatively calibrated class. Candidate 822940 is the cleaner lower-closed-fraction comparator and has nine destructive-ablation robustness checks, although the property-margin proxy is narrower. Candidate 581668 is an important intermediate topology.

## What is robust and what is not

The existence of six distinct Tier-B topologies shows that a finite qualitative window is not unique to 693168. It does not establish material-property robustness. The reduced transfer predicts only 693168 as a robust comparator under the chosen margin rule; all other candidates remain calibration-sensitive.

An exact family OAT was deliberately not run in this campaign. The `tierB_family_exact_OAT.csv` table is empty and records that limitation rather than silently converting a proxy into exact evidence. A future bounded exact check should prioritize 822940 and 581668, not repeat another topology search.
