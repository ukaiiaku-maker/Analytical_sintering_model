from dataclasses import replace
import inspect
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd

from zro2_forward.conditioned_950c import make_pdf_conditioned_initial_state
from zro2_forward.pore_population import transfer_fluxes
from zro2_forward.integrator import ModelState
from zro2_forward.sensitivity_audit import calibrated_model

ROOT=Path(__file__).parents[1];OUT=ROOT/"results/zro2_forward_950C_sensitivity_chen_failure_audit"

def test_common_paths_clone_identical_arrays_and_record_all_inputs():
    d=pd.read_csv(OUT/"common_state_factorial_design.csv");required={"rho_start","G_start_nm","pore_D50_nm","pore_log_width","phi_iso_fraction","phi_closed_fraction"};assert required<=set(d)
    r=d.iloc[0];a=make_pdf_conditioned_initial_state(r.rho_start,r.G_start_nm,pore_D50_nm=r.pore_D50_nm,pore_log_width=r.pore_log_width,phi_iso_fraction=r.phi_iso_fraction,phi_closed_fraction=r.phi_closed_fraction);b=make_pdf_conditioned_initial_state(r.rho_start,r.G_start_nm,pore_D50_nm=r.pore_D50_nm,pore_log_width=r.pore_log_width,phi_iso_fraction=r.phi_iso_fraction,phi_closed_fraction=r.phi_closed_fraction)
    for name in ("radii_m","phi_open","phi_iso","phi_closed","number_open"):assert np.array_equal(getattr(a.pores,name),getattr(b.pores,name))

def test_primary_closed_zero_and_method_specific_initialization_absent():
    d=pd.read_csv(OUT/"common_state_factorial_design.csv");assert d.query("state_class=='primary_common_state'").phi_closed_fraction.eq(0).all();assert d.query("phi_closed_fraction>0").state_class.eq("diagnostic_closed_initial_state").all()
    src=(ROOT/"run_zro2_forward_950C_common_state_sensitivity.py").read_text();assert not any(f'"{x}"' in src for x in ("CS","LMS","HMS","TSS"))

def test_PR_is_conservative_and_named_transitions_are_separate():
    s=make_pdf_conditioned_initial_state();p=s.pores;po,pi,pc,_=transfer_fluxes(p,.85,1e-14,.1,.5)
    assert abs((po+pi+pc).sum())<1e-18
    source=inspect.getsource(transfer_fluxes);assert "isolation_rate" in source and "closure_rate" in source

def test_doubling_radius_gives_sixteenfold_closed_removal_time():
    _,q=calibrated_model();r=25e-9;activity=.2
    tau1=q.closed_tau0_s*(r/25e-9)**4/activity;tau2=q.closed_tau0_s*(2*r/25e-9)**4/activity
    assert np.isclose(tau2/tau1,16,rtol=.02)

def test_open_and_closed_shrinkage_touch_only_named_stores():
    m,_=calibrated_model({"C_PR_m2":0.,"C_iso_m2":0.,"C_close_m2":0.});s=make_pdf_conditioned_initial_state(phi_closed_fraction=.1);od,ii,cc,closed,_,diag=m.rates(s,1300)
    assert np.allclose(ii,0);assert closed>0;assert np.all(cc<=0);assert diag["rho_dot_open_sinv"]>=0 and diag["rho_dot_closed_sinv"]>=0

def test_barrier_status_and_full_process_separation():
    s=pd.read_csv(OUT/"common_state_fast_rate_summary.csv");assert s.barrier_mode.notna().all() and s.barrier_extrapolated.all();assert not (OUT/"full_process_dense_histories.csv").exists();assert (ROOT/"results/zro2_forward_pdf_conditioned_950C_comparison/full_process_dense_histories.csv").exists()

def test_legacy_mechanism_search_files_are_unmodified():
    names=subprocess.check_output(["git","diff","--name-only","b60f2bf"],cwd=ROOT,text=True).splitlines();allowed=("zro2_forward/","run_zro2_forward_","plot_zro2_forward_","tests/test_zro2_forward_","docs/ZRO2_FORWARD_","results/zro2_forward_950C_");assert all(n.startswith(allowed) for n in names)
