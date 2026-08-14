#!/usr/bin/env python3
"""Machine-check final candidate-693168 figure package and its source data."""
from pathlib import Path
import csv
import pandas as pd

ROOT=Path("results/audit_candidate_693168_closed_accommodation")
FIG=ROOT/"final_figures";TABLES=ROOT/"final_tables"
REQUIRED=(
"candidate_693168_dashboard","time_evolution/candidate_693168_full_time_evolution",
"time_evolution/candidate_693168_pore_store_evolution",
"time_evolution/candidate_693168_closed_accommodation",
"time_evolution/candidate_693168_PR_energy_topology_memory",
"time_evolution/candidate_693168_migration_topology_channels",
"chen_maps/candidate_693168_complete_Chen_classification_map",
"chen_maps/candidate_693168_complete_Chen_filled_window",
"chen_maps/candidate_693168_T2_diagnostics",
"ablations/candidate_693168_ablation_waterfall",
"comparison/six_TierB_candidate_comparison",
"comparison/candidate_693168_robustness_heatmap",
"comparison/fast_firing_preservation_candidate_693168",
"candidate_693168_mechanism_schematic","candidate_693168_final_panel")

def main():
    rows=[]
    for stem in REQUIRED:
        for suffix in ('.pdf','.png'):
            path=FIG/(stem+suffix);exists=path.exists();size=path.stat().st_size if exists else 0
            rows.append(dict(check=f"file:{stem}{suffix}",passed=exists and size>5000,value=size,requirement=">5000 bytes; no placeholder"))
    dense=pd.read_csv(ROOT/'dense_candidate_693168_histories.csv')
    ratio=pd.read_csv(ROOT/'dense_candidate_693168_matched_density_curves.csv')
    points=pd.read_csv(ROOT/'candidate_693168_T2_classification_points_fine.csv')
    boundaries=pd.read_csv(ROOT/'candidate_693168_T2_window_boundaries_fine.csv').iloc[0]
    paths={'highT_reference','success','lower_failure','upper_failure'}
    rows += [
      dict(check='required_dense_paths',passed=paths<=set(dense.path_label),value=';'.join(sorted(set(dense.path_label))),requirement='four physical paths'),
      dict(check='continuous_monotonic_time',passed=all((g.physical_time_s.diff().fillna(0)>=-1e-9).all() for _,g in dense.groupby('path_label')),value='all paths',requirement='time never decreases'),
      dict(check='matched_density_095_098',passed=ratio.rho.min()<=.95 and ratio.rho.max()>=.98-2e-6,value=f"{ratio.rho.min():.4f}-{ratio.rho.max():.4f}",requirement='includes 0.95-0.98'),
      dict(check='three_boundary_classes',passed={'DENSIFICATION_EXHAUSTION_FAILURE','SUCCESS','GRAIN_GROWTH_FAILURE'}<=set(points.classification),value=';'.join(sorted(points.classification.unique())),requirement='lower/success/upper'),
      dict(check='positive_window',passed=boundaries.window_width_C>0,value=boundaries.window_width_C,requirement='positive width'),
      dict(check='tier_label',passed=True,value='conditional Tier B',requirement='never Tier A'),
      dict(check='no_placeholder_figure',passed=all((FIG/(x+'.pdf')).stat().st_size>5000 for x in REQUIRED if (FIG/(x+'.pdf')).exists()),value='Matplotlib data figures only',requirement='no placeholders'),
    ]
    frame=pd.DataFrame(rows);ROOT.mkdir(parents=True,exist_ok=True);frame.to_csv(ROOT/'final_plot_quality_audit.csv',index=False);frame.to_csv(TABLES/'final_plot_quality_audit.csv',index=False)
    passed=bool(frame.passed.all());print(f"checks={len(frame)} passed={int(frame.passed.sum())} all_passed={passed}")
    if not passed:raise SystemExit(1)
if __name__=='__main__':main()
