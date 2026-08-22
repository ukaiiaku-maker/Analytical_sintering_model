from __future__ import annotations

import hashlib
import inspect
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd

import promote_zro2_emergent_pore_closure_final_test as p

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/zro2_forward_emergent_pore_closure_final_test"


def test_01_barrier_hash_unchanged():
    assert hashlib.sha256(p.BARRIER_PATH.read_bytes()).hexdigest()=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"


def test_02_03_diffusivity_laws_unchanged():
    T=1473.15
    assert np.isclose(p.MAT.D_GB(T),.056*np.exp(-380000/(p.R_GAS*T)))
    assert np.isclose(p.MAT.D_s(T),.10*np.exp(-380000/(p.R_GAS*T)))


def test_04_failed_global_mobility_fit_inactive():
    assert p.MAT.mobility_prefactor_status=="calibrated once to conventional-sintering final grain size"
    assert p.MAT.M0_m4_J_s==5.8e-3


def test_05_06_no_physical_qclosed_and_apparent_is_postrun():
    assert "q_closed" not in {x.lower() for x in inspect.signature(p.emergent_pore_closure_v1).parameters}
    q=pd.read_csv(OUT/"closed_law_parameter_classification.csv")
    row=q[q.term=="Q_closed_app"].iloc[0]
    assert row.classification=="empirical diagnostic" and "post-run" in row.mapping


def test_07_08_closed_shrinkage_requires_inventory_and_stress():
    for kernel in ("renewal","GB_diffusion","surface_only"):
        assert p.emergent_pore_closure_v1(1473.15,25e-9,0,1,.5,.25,3,kernel=kernel)["rho_dot_closed_sinv"]==0
        assert p.emergent_pore_closure_v1(1473.15,25e-9,.1,1,.5,1.1,3,kernel=kernel)["rho_dot_closed_sinv"]==0


def test_09_gas_pressure_reduces_stress():
    a=p.emergent_pore_closure_v1(1473.15,25e-9,.1,1,.5,0,3)
    b=p.emergent_pore_closure_v1(1473.15,25e-9,.1,1,.5,.9,3)
    assert b["sigma_c_Pa"]<a["sigma_c_Pa"] and b["rho_dot_closed_sinv"]<a["rho_dot_closed_sinv"]


def test_10_radius_law_decreases_rate():
    for m in (3,4):
      for kernel in ("renewal","GB_diffusion"):
        a=p.emergent_pore_closure_v1(1473.15,10e-9,.1,1,.5,.25,m,kernel=kernel)
        b=p.emergent_pore_closure_v1(1473.15,100e-9,.1,1,.5,.25,m,kernel=kernel)
        assert b["rho_dot_closed_sinv"]<a["rho_dot_closed_sinv"]


def test_11_surface_accommodation_alone_does_not_densify():
    z=p.emergent_pore_closure_v1(1473.15,25e-9,.1,1,.2,.25,4,kernel="surface_only")
    assert z["A_dot_shape_recovery_sinv"]>0 and z["rho_dot_closed_sinv"]==0


def test_12_preparation_transfers_conserve_volume():
    x=np.array([.2,.1,.05]); y=p.conservative_transfer(x,.001,10,0,2)
    assert np.isclose(x.sum(),y.sum()) and (y>=0).all()


def test_13_only_named_shrinkage_changes_density():
    reg=pd.read_csv(OUT/"PR_preparation_flux_registry.csv")
    assert not reg.changes_density.any()
    h=pd.read_csv(OUT/"final_emergent_closure_ablation_histories.csv")
    assert not h[h.ablation=="surface_accommodation_only"].changes_density.any()


def test_14_density_identity_at_every_stored_boundary_state():
    q=pd.read_csv(OUT/"final_boundary_preservation_test.csv")
    assert np.max(np.abs(q.density_identity_residual))<1e-12


def test_15_accommodation_bounded():
    h=pd.read_csv(OUT/"final_emergent_closure_ablation_histories.csv")
    assert h.A_closed.between(0,1).all()


def test_16_energy_ledger_channels_named():
    q=pd.read_csv(OUT/"final_energy_ledger_channel_registry.csv")
    required={"P_open_dens","P_closed_dens","P_PR","P_surface_smooth","P_pore_coarsen","P_GB_growth","P_drag","P_gas","P_other","P_residual"}
    assert required==set(q.channel) and q.named.all()


def test_17_strict_GB_balance_diagnostic_only():
    h=pd.read_csv(OUT/"final_energy_ledger_selected_paths.csv")
    assert (h.strict_GB_area_balance_status=="diagnostic_only").all()


def test_18_any_empirical_closure_labeled_diagnostic():
    q=pd.read_csv(OUT/"closed_law_parameter_classification.csv")
    empirical=q[q.classification=="empirical diagnostic"]
    assert len(empirical)>0 and empirical.term.str.contains("app").all()


def test_19_finite_window_requires_lower_and_upper_boundaries():
    q=pd.read_csv(OUT/"final_boundary_topology_summary.csv")
    assert ((~q.strict_finite_window)|(q.lower_boundary&q.upper_boundary&(q.success_points>1))).all()


def test_20_T2_paths_clone_identical_first_step_state():
    q=pd.read_csv(OUT/"final_boundary_preservation_test.csv")
    for _,z in q.groupby("state_id"):
        for c in ("rho","G_nm","phi_closed","A","chi","r_nm"): assert z[c].nunique()==1


def test_21_local_law_contains_no_processing_labels():
    tokens=set(re.findall(r"[a-z_]+",inspect.getsource(p.emergent_pore_closure_v1).lower()))
    forbidden={"cs","lms","hms","tss","fast","slow","protocol","schedule","ramp_rate","target"}
    assert not tokens.intersection(forbidden)


def test_22_reports_explicitly_disclaim_validation():
    names=["ZRO2_FORWARD_EMERGENT_PORE_CLOSURE_FINAL_TEST.md","ZRO2_FORWARD_EMERGENT_CLOSURE_MODEL_FORM.md",
           "ZRO2_FORWARD_EMERGENT_CLOSURE_PARAMETER_MAPPING.md","ZRO2_FORWARD_EMERGENT_CLOSURE_TEST_RESULTS.md",
           "ZRO2_FORWARD_EMERGENT_CLOSURE_NEXT_DECISION.md"]
    for name in names:
        text=(ROOT/"docs"/name).read_text().lower(); assert "not validation" in text or "no validation claim" in text


def test_figures_have_sources_and_no_success_process_map():
    q=pd.read_csv(OUT/"figure_inventory.csv")
    assert len(q)==8 and q.pdf_nonempty.all() and q.png_nonempty.all() and not q.placeholder.any()
    assert not q.success_colored_process_map.any()
    for source in q.source_table: assert (OUT/source).is_file()


def test_23_staged_payload_excludes_forbidden_files_when_present():
    names=subprocess.run(["git","diff","--cached","--name-only"],cwd=ROOT,text=True,capture_output=True,check=True).stdout.splitlines()
    allowed_prefix="results/zro2_forward_emergent_pore_closure_final_test/"
    for name in names:
        assert not name.endswith(".DS_Store")
        assert not (name.lower().endswith(".pdf") and not name.startswith(allowed_prefix))
        if name.startswith("results/"): assert name.startswith(allowed_prefix)
