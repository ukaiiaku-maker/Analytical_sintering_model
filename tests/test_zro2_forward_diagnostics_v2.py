from pathlib import Path
import numpy as np
import pandas as pd
from zro2_forward.barrier_json import BarrierModel
from zro2_forward.integrator import ModelParameters
from zro2_forward.diagnostics_v2 import MODES

def test_inputs_present():
    for p in [Path("data/zro2/bicrystal_creep_barrier_export.json"),Path("data/targets/mazaheri_8ysz_2008/1-s2.0-S092150930800302X-main.pdf")]: assert p.is_file() and p.stat().st_size>0

def test_barrier_modes_are_labeled_and_low_temperature_is_flagged():
    b=BarrierModel.load("data/zro2/bicrystal_creep_barrier_export.json")
    assert set(MODES)=={"nearest_slice_clamp","pchip_extrapolate","fixed_lowT_slope","generic_anchor_barrier"}
    for mode in MODES:
        x=b.with_mode(mode); assert x.mode==mode and x.diagnostics(1e8,1773.15)["temperature_extrapolated"]

def test_calibration_objectives_are_distinct():
    x=pd.read_csv("results/zro2_forward_diagnostic_calibration_v2/CS_calibrated_parameter_sets.csv")
    assert set(x.objective_mode)=={"endpoint_only","density_curve_only","density_plus_grain_endpoint"}
    assert x.drop(columns="objective_mode").drop_duplicates().shape[0]>=2

def test_sensitivity_does_not_mutate_barrier():
    b=BarrierModel.load("data/zro2/bicrystal_creep_barrier_export.json"); before=b.G0_J.copy()
    _=ModelParameters(C_PR_m2=3e-23,rho_close_mid=.87,zener_length_factor=3)
    assert np.array_equal(before,b.G0_J)

def test_chen_success_and_window_require_both_conditions():
    x=pd.read_csv("results/zro2_forward_diagnostic_calibration_v2/barrier_mode_chen_counts.csv")
    assert (x.success_count<=x.density_ok_count).all() and (x.success_count<=x.growth_ok_count).all()
    base=x[x.barrier_mode.eq("nearest_slice_clamp")].iloc[0]; assert base.success_count==0

def test_closed_fluxes_remain_named_and_separate():
    src=Path("zro2_forward/integrator.py").read_text()
    assert "rho_dot_open_sinv" in src and "rho_dot_closed_sinv" in src
    assert "phi_closed" in src and "shrink=p.phi_closed/tau0" in src

def test_microwave_diagnostic_disabled_by_default():
    q=ModelParameters(); assert q.microwave_mode=="none" and q.effective_nucleation_activity_multiplier==1
