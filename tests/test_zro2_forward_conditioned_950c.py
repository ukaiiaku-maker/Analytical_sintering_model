from dataclasses import fields
import inspect
from pathlib import Path

import numpy as np
import pandas as pd

from zro2_forward.conditioned_950c import make_pdf_conditioned_initial_state
from zro2_forward.integrator import ForwardModel, ModelParameters
from run_zro2_forward_conditioned_chen_map import chen_success, finite_window


ROOT = Path(__file__).parents[1]
OUT = ROOT / "results/zro2_forward_pdf_conditioned_950C_comparison"


def states():
    return {method: make_pdf_conditioned_initial_state() for method in ("CS", "LMS", "HMS")}


def test_conditioned_initial_state_is_common_and_not_curve_specific():
    ss = states()
    for state in ss.values():
        assert np.isclose(state.T_K, 1223.15)
        assert np.isclose(state.rho, .66)
        assert np.isclose(state.G_m, 50e-9)
    reference = ss["CS"].pores
    for method in ("LMS", "HMS"):
        candidate = ss[method].pores
        for name in ("radii_m", "phi_open", "phi_iso", "phi_closed", "number_open"):
            assert np.array_equal(getattr(reference, name), getattr(candidate, name))
    source = inspect.getsource(make_pdf_conditioned_initial_state)
    assert not any(label in source for label in ('"CS"', '"LMS"', '"HMS"', '"TSS"'))


def test_separate_tables_and_canonical_start_are_persisted():
    full = pd.read_csv(OUT / "full_process_dense_histories.csv")
    conditioned = pd.read_csv(OUT / "pdf_conditioned_dense_histories.csv")
    assert full.T_C.min() < 950
    starts = conditioned.sort_values("t_s").groupby("case").first()
    assert starts.T_C.ge(950).all()
    canonical = pd.read_csv(OUT / "pdf_conditioned_initial_states.csv").query("state_id == 'nominal'").iloc[0]
    assert np.isclose(canonical.T_start_C, 950)
    assert np.isclose(canonical.rho_start, .66)
    assert np.isclose(canonical.G_start_nm, 50)
    run_source = (ROOT / "run_zro2_forward_pdf_conditioned_comparison.py").read_text()
    assert run_source.count("make_pdf_conditioned_initial_state()") >= 2


def test_pre950_mismatch_is_reframed_and_missing_data_are_not_overlaid():
    reinterpretation = pd.read_csv(OUT / "baseline_reinterpretation.csv")
    wording = " ".join(reinterpretation.corrected_interpretation.fillna("").astype(str)).lower()
    assert "initialization/pre-950 c prediction mismatch" in wording
    plot_source = (ROOT / "plot_zro2_forward_conditioned_results.py").read_text()
    assert "q=g[g.T_C<=950]" in plot_source
    assert "common target state" in plot_source
    observed = pd.read_csv(OUT / "target_curves_observed_interval_950C.csv")
    assert observed.loc[observed.T_C.notna(), "T_C"].ge(950).all()


def test_no_method_labels_enter_local_constitutive_law():
    source = inspect.getsource(ForwardModel.rates)
    assert not any(label in source for label in ("CS", "LMS", "HMS", "TSS"))


def test_barrier_mode_is_recorded_for_conditioned_diagnostics():
    modes = pd.read_csv(OUT / "pdf_conditioned_barrier_mode_comparison.csv")
    expected = {"nearest_slice_clamp", "pchip_extrapolate", "fixed_lowT_slope", "generic_anchor_barrier"}
    assert set(modes.barrier_mode) == expected
    chen = pd.read_csv(OUT / "pdf_conditioned_chen_classification_points.csv")
    assert chen.barrier_mode.notna().all()


def test_chen_success_and_window_guardrails_are_joint_and_bracketed():
    assert chen_success(.976, .29)
    assert not chen_success(.95, .29)
    assert not chen_success(.976, .60)
    assert finite_window([1100, 1200], True, True)
    assert not finite_window([1100, 1200], False, True)
    assert not finite_window([1100], True, True)
    table = pd.read_csv(OUT / "pdf_conditioned_chen_classification_points.csv")
    assert (table.chen_success == (table.density_ok & table.strict_target_ok)).all()
    windows = pd.read_csv(OUT / "pdf_conditioned_chen_window_boundaries.csv")
    assert (~windows.finite_window | (windows.lower_boundary_present & windows.upper_boundary_present)).all()


def test_microwave_multiplier_is_disabled_by_default():
    names = {field.name for field in fields(ModelParameters)}
    assert not any("microwave" in name.lower() for name in names)
