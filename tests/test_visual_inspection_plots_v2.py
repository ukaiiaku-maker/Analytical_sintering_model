from pathlib import Path
import json
import pandas as pd

ROOT=Path("results/visual_inspection_candidate_plots_v2")

def test_candidates_and_tiers():
    d=pd.read_csv(ROOT/"tables/selected_candidates_for_visualization.csv")
    assert {"E0021","E0142"}<=set(d.material_id)
    assert set(d[d.material_id=="E0021"].tier)=={"Tier_C"}
    assert set(d[d.material_id=="E0142"].tier)=={"Tier_B"}

def test_dense_histories_and_continuous_two_step_time():
    idx=pd.read_csv(ROOT/"tables/dense_fast_histories_index.csv");assert idx.n_points.min()>=1000
    d=pd.read_csv(ROOT/"histories/dense_two_step_histories.csv")
    for _,g in d.groupby("role"):
        assert (g.physical_time_s.diff().fillna(0)>=0).all()
        first=g[g.stage=="first_step"].physical_time_s.max();second=g[g.stage=="second_step"].physical_time_s.min();assert second>=first

def test_robustness_grid():
    d=pd.read_csv(ROOT/"tables/fast_firing_initial_condition_robustness.csv")
    assert d.rho0.nunique()==5 and d.G0_nm.nunique()==6
    assert {20,50,100}==set(d.fast_rate_C_min)

def test_inventory_and_quality_audit():
    inv=pd.read_csv(ROOT/"tables/figure_inventory.csv");audit=pd.read_csv(ROOT/"tables/plot_quality_audit.csv")
    assert len(inv)==18 and audit.passed.all()
    for _,r in inv.iterrows():
        assert (ROOT/"figures"/r.filename_pdf).exists();assert (ROOT/"figures"/r.filename_png).exists()
    assert json.loads((ROOT/"quality_audit_summary.json").read_text())["two_step_time_monotonic"]

def test_filled_windows_and_ablation_labels():
    b=pd.read_csv(ROOT/"tables/E0142_TierB_Chen_window_boundaries_v2.csv");assert (b.window_width_C>0).all()
    a=pd.read_csv(ROOT/"tables/fast_firing_ablation_labels.csv")
    for c in ("PR_redistribution_enabled","nucleation_limitation_enabled","topology_enabled","growth_before_activation_enabled"):assert c in a
    assert Path("docs/VISUAL_INSPECTION_PLOTS_V2_MISSING_DATA.md").exists()
