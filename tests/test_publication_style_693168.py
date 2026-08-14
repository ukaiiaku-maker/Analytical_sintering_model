from pathlib import Path
import json,pandas as pd
OUT=Path('results/publication_style_sintering_figures_693168');SRC=OUT/'source_tables'
def test_dense_histories_have_required_paths_and_continuous_time():
 d=pd.read_csv(SRC/'dense_time_histories.csv');assert {'highT_reference','lower_failure','success','upper_failure'}<=set(d.path_label)
 for _,g in d.groupby('path_label'):assert (g.physical_time_s.diff().fillna(0)>=-1e-9).all()
def test_filled_chen_grids_are_dense_and_have_three_classes():
 frames=[pd.read_csv(SRC/f) for f in ('chen_map_T1_T2_classification.csv','chen_map_G1_T2_classification.csv','chen_map_switch_density_T2_classification.csv')]
 assert len(frames[0])>=500 and len(frames[1])>=5000 and len(frames[2])>=500
 classes=set().union(*(set(x.classification) for x in frames));assert {'DENSIFICATION_EXHAUSTION_FAILURE','SUCCESS','GRAIN_GROWTH_FAILURE'}<=classes
 assert sum((x.classification=='SUCCESS').sum() for x in frames)>10
def test_fast_firing_separate_and_complete_rates():
 d=pd.read_csv(SRC/'dense_fast_firing_histories.csv');assert {1,20,50,100}<=set(d.heating_rate_C_min)
 assert (OUT/'fast_firing/fast_firing_G_vs_rho.pdf').exists();assert (OUT/'two_step/two_step_G_vs_rho.pdf').exists()
def test_inventory_and_quality():
 inv=pd.read_csv(SRC/'publication_style_figure_inventory.csv');assert set(inv.category)>={'main','supplement','chen_map','fast_firing','two_step'};assert len(inv)>=24
 for c in ('filename_pdf','filename_png'):
  for f in inv[c]:assert (OUT/f).exists() and (OUT/f).stat().st_size>5000
 q=pd.read_csv(SRC/'publication_style_plot_quality_audit.csv');assert q.passed.all()
def test_candidate_is_conditional_tierB_and_model_frozen():
 inv=pd.read_csv(SRC/'publication_style_figure_inventory.csv');assert inv.candidate_status.str.contains('Tier B').all();assert not inv.candidate_status.str.contains('Tier A').any()
 state=json.loads((OUT/'map_run_state.json').read_text());assert state['physics_changed'] is False
