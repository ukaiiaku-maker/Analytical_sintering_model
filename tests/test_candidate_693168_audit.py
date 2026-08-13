from pathlib import Path
import numpy as np
import pandas as pd
import audit_candidate_693168_closed_accommodation as audit
import interacting_local_region_model as model

ROOT=Path("results/audit_candidate_693168_closed_accommodation")

def test_candidate_parameter_vector_loaded_exactly():
    _,decoded=audit.candidate_parameters()
    assert len(decoded)==47
    assert audit.decoder.fingerprint(decoded)=="b46304a8d7c7055a"

def test_exact_reproduction_switch_and_G1():
    d=pd.read_csv(ROOT/"candidate_693168_reproduction_summary.csv")
    exact=d[d.run_label=="original_exact"].iloc[0]
    assert abs(exact.exact_switch_density-.87999936)<1e-4
    assert abs(exact.G1_nm-117.066668)/117.066668<.01

def test_tighter_timestep_comparison_persisted():
    d=pd.read_csv(ROOT/"candidate_693168_reproduction_summary.csv")
    assert {30.,15.,5.}<=set(d.max_timestep_min)

def test_dense_histories_and_continuous_time():
    d=pd.read_csv(ROOT/"dense_candidate_693168_histories.csv")
    assert {"highT_reference","success","lower_failure","upper_failure"}<=set(d.path_label)
    for _,g in d.groupby("path_label"):
        assert (g.physical_time_s.diff().fillna(0)>=-1e-9).all()

def test_matched_density_only_both_attained():
    d=pd.read_csv(ROOT/"dense_candidate_693168_matched_density_curves.csv")
    assert d.rho.min()<=.95 and d.rho.max()>=.98-2e-6
    assert d.both_paths_attained.all()
    assert d.highT_interpolation_supported.all() and d.two_step_interpolation_supported.all()

def test_T2_map_has_cloned_state_boundary_topology():
    d=pd.read_csv(ROOT/"candidate_693168_T2_classification_points_fine.csv")
    assert {"DENSIFICATION_EXHAUSTION_FAILURE","SUCCESS","GRAIN_GROWTH_FAILURE"}<=set(d.classification)
    b=pd.read_csv(ROOT/"candidate_693168_T2_window_boundaries_fine.csv").iloc[0]
    assert b.complete and b.window_width_C>0

def test_closed_flux_zero_without_closed_volume():
    p,_=audit.candidate_parameters();s=model.initial(p["N_regions"],p=p);s.phi_closed[:]=0
    assert np.all(model.local_fluxes(s,1200,p)["rho_dot_closed"]==0)

def test_open_shrinkage_does_not_remove_closed_or_isolated():
    p,_=audit.candidate_parameters();p.update(k_closed=0.,k_PR=0.,closed_transition=0.,detachment=0.,recapture=0.,k_sweep_damaged=0.,k_sweep_connected=0.)
    s=model.initial(p["N_regions"],p=p);before=(s.phi_closed.copy(),s.phi_iso.copy())
    model.advance(s,1400,p,60,model.network_adjacency(p["N_regions"],p))
    assert np.allclose(before[0],s.phi_closed) and np.allclose(before[1],s.phi_iso)

def test_finite_accommodation_nonnegative_and_bounded():
    d=pd.read_csv(ROOT/"dense_candidate_693168_histories.csv")
    assert (d.closed_accommodation_available>=0).all()
    assert (d.closed_accommodation_available<=d.closed_accommodation_capacity+1e-9).all()

def test_infinite_accommodation_labeled_as_ablation():
    d=pd.read_csv(ROOT/"candidate_693168_ablation_audit.csv")
    assert (d.ablation=="infinite_closed_accommodation").sum()==1

def test_required_final_figures_nonempty_and_tierB():
    inventory=pd.read_csv(ROOT/"final_figure_inventory.csv")
    assert len(inventory)>=15
    for col in ("filename_pdf","filename_png"):
        for name in inventory[col]:assert (ROOT/name).exists() and (ROOT/name).stat().st_size>5000
    quality=pd.read_csv(ROOT/"final_plot_quality_audit.csv")
    assert quality.passed.all()

def test_model_files_reported_frozen():
    import json
    state=json.loads((ROOT/"audit_run_state.json").read_text())
    assert state["frozen_model_files_changed"] is False
