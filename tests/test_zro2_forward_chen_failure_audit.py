import inspect
import re
from pathlib import Path
import numpy as np
import pandas as pd
from zro2_forward.integrator import ForwardModel
from zro2_forward.pore_population import transfer_fluxes
from zro2_forward.densification import kinetic_state,density_rate,connectivity
from zro2_forward.grain_growth import growth_state
from zro2_forward.sensitivity_audit import chen_success,finite_window

ROOT=Path(__file__).parents[1];OUT=ROOT/"results/zro2_forward_950C_sensitivity_chen_failure_audit"

def test_strict_success_is_joint_and_finite_window_is_bracketed_nonisolated():
    assert chen_success(.976,.29);assert not chen_success(.95,.29);assert not chen_success(.976,.5)
    assert finite_window([1100,1200],True,True);assert not finite_window([1100],True,True);assert not finite_window([1100,1200],False,True)

def test_persisted_windows_require_boundaries_and_reject_zero_width():
    b=pd.read_csv(OUT/"chen_boundary_ordering_table.csv");assert (~b.finite_window_present|(b.lower_boundary_present&b.upper_boundary_present&b.success_band_present&~b.zero_width_or_isolated)).all()
    assert b.loc[b.T_lower_density_C.isna()|b.T_upper_growth_C.isna(),"gap_C"].isna().all()

def test_relaxed_diagnostics_never_overwrite_strict():
    r=pd.read_csv(OUT/"chen_target_relaxation_diagnostics.csv");assert not r.overwrites_strict.any();strict=r.query("evidence_type=='strict'");assert len(strict)==1 and strict.success_count.iloc[0]==0 and strict.finite_window_count.iloc[0]==0;assert r.query("evidence_type!='strict'").evidence_type.eq("relaxed_diagnostic").all()

def test_no_method_or_schedule_labels_in_local_constitutive_functions():
    forbidden=("cs","lms","hms","tss","fast","slow","protocol","target","schedule","ramp_rate")
    for fn in (ForwardModel.rates,transfer_fluxes,kinetic_state,density_rate,connectivity,growth_state):
        src=inspect.getsource(fn).lower();assert not any(re.search(rf"\b{word}\b",src) for word in forbidden)

def test_every_conditioned_run_records_barrier_and_extrapolation():
    for name in ("common_state_fast_rate_summary.csv","common_state_matched_density_curves.csv","chen_failure_decomposition_full.csv","chen_failure_OAT_summary.csv"):
        x=pd.read_csv(OUT/name);assert x.barrier_mode.notna().all();assert x.barrier_extrapolated.notna().all()

def test_failure_and_triage_are_nonvalidation_diagnostics():
    x=pd.read_csv(OUT/"chen_failure_decomposition_full.csv");assert len(x)==195 and not x.chen_success.any();t=pd.read_csv(OUT/"candidate_triage_950C_sensitivity.csv");assert not t.accepted_for_future_calibration.any()
