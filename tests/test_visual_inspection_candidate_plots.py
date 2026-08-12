from pathlib import Path
import pandas as pd

ROOT=Path("results/visual_inspection_candidate_plots")

def test_selected_candidates_include_both_materials_and_correct_tiers():
    d=pd.read_csv(ROOT/"selected_candidates_for_visualization.csv")
    assert {"E0021","E0142"}<=set(d.material_id)
    assert set(d[d.material_id=="E0021"].tier)=={"Tier_C"}
    assert set(d[d.material_id=="E0142"].tier)=={"Tier_B"}

def test_inventory_files_exist():
    d=pd.read_csv(ROOT/"figure_inventory.csv")
    assert len(d)>=18
    for _,r in d.iterrows():
        stem=ROOT/"figures"/r.filename_pdf.removesuffix('.pdf')
        assert stem.with_suffix('.pdf').exists()
        assert stem.with_suffix('.png').exists()

def test_filled_windows_have_positive_width_and_no_isolated_points():
    d=pd.read_csv(ROOT/"tables/selected_Chen_window_boundaries.csv")
    assert (d.window_width_C>0).all()
    assert (d.T2_last_success_C>d.T2_first_success_C).all()

def test_ratio_thresholds_and_missing_data_are_explicit():
    src=Path("visual_inspection_candidate_plots.py").read_text()
    assert "(1.2,1.5,2)" in src
    missing=(ROOT/"visual_inspection_missing_data.md").read_text()
    assert "not reconstructed or inferred" in missing

def test_visualization_does_not_modify_physics():
    src=Path("visual_inspection_candidate_plots.py").read_text()
    assert "replace(mat,ablation_mode" in src
    assert "MaterialKinetics(" not in src
