from __future__ import annotations
import hashlib, inspect, json, subprocess
from pathlib import Path
import numpy as np
import pandas as pd

import build_zro2_distributional_pore_population_model as m

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results/zro2_forward_distributional_pore_population_model"
SHA="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"

def test_01_barrier_hash_unchanged():assert hashlib.sha256((ROOT/"data/zro2/bicrystal_creep_barrier_export.json").read_bytes()).hexdigest()==SHA
def test_02_gb_diffusivity_unchanged():assert (m.MAT.D_GB0_m2_s,m.MAT.Q_GB_J_mol)==(.056,380000.)
def test_03_surface_diffusivity_unchanged():assert (m.MAT.D_s0_m2_s,m.MAT.Q_s_J_mol)==(.10,380000.)
def test_04_failed_global_mobility_inactive():assert not json.loads((OUT/"run_state.json").read_text())["failed_global_mobility_fit_active"]
def test_05_no_physical_qclosed():assert "Q_closed" not in m.PopulationParameters.__dataclass_fields__
def test_06_qclosed_app_diagnostic_only():
    q=pd.read_csv(OUT/"distribution_parameter_registry.csv");r=q[q.parameter.eq("Q_closed_app")].iloc[0];assert r.classification=="empirical diagnostic"
def test_07_density_identity():
    s=m.initial_state();assert abs(s.rho-(1-s.phi.sum()))<1e-14
def test_08_transfers_conserve():
    q=pd.read_csv(OUT/"distribution_conservation_tests.csv");assert q[q.test.str.contains("conserved")].passed.all()
def test_09_surface_alone_no_density():
    q=pd.read_csv(OUT/"distribution_conservation_tests.csv");assert q[q.test.eq("surface_PR_no_density")].passed.all()
def test_10_pr_alone_no_density():
    q=pd.read_csv(OUT/"distribution_conservation_tests.csv");assert q[(q.mode=="PR_only")&q.test.eq("surface_PR_no_density")].passed.all()
def test_11_zero_closed_inventory():
    q=pd.read_csv(OUT/"distribution_conservation_tests.csv");assert q[q.test.eq("zero_closed_inventory")].passed.all()
def test_12_nonpositive_stress_stops_closed():
    q=pd.read_csv(OUT/"distribution_conservation_tests.csv");assert q[q.test.eq("nonpositive_stress_stops_closed_shrinkage")].passed.all()
def test_13_gas_reduces_stress():
    q=pd.read_csv(OUT/"distribution_conservation_tests.csv");assert q[q.test.eq("gas_reduces_closed_stress")].passed.all()
def test_14_r4_surface_time():
    q=pd.read_csv(OUT/"distribution_conservation_tests.csv");r=q[q.test.eq("r4_time_increases")].iloc[0];assert r.passed and r.residual_or_metric>1
def test_15_small_pores_pin_more():
    a=m.initial_state(D50_nm=10);b=m.initial_state(D50_nm=100);assert m.metrics(a)["Zener_pinning_metric_minv"]>m.metrics(b)["Zener_pinning_metric_minv"]
def test_16_mean_and_distributional_zener_reported():
    q=pd.read_csv(OUT/"source_tables/distributional_Zener_metric_source.csv");assert {"mean_radius_Zener","distributional_Zener"}<=set(q)
def test_17_migration_does_not_change_density_rate():
    s=m.initial_state();a=m.rates(s,1473.15,m.PopulationParameters());b=m.rates(s,1473.15,m.PopulationParameters(distributional_zener=False));assert a["rho_dot_open"]==b["rho_dot_open"] and a["rho_dot_closed"]==b["rho_dot_closed"]
def test_18_t2_paths_clone_state():
    q=pd.read_csv(OUT/"synthetic_distribution_boundary_test.csv");assert q.groupby(["state_id","representation","closed_kernel","m"]).initial_state_fingerprint.nunique().eq(1).all()
def test_19_finite_window_requires_boundaries():
    q=pd.read_csv(OUT/"distribution_window_boundaries.csv");assert not q.strict_finite_window.any();assert not ((q.strict_finite_window)&~(q.lower_boundary&q.upper_boundary)).any()
def test_20_local_laws_have_no_schedule_labels():
    source="\n".join(inspect.getsource(x) for x in (m.local_terms,m.rates,m._closed_unit_rate));forbidden=("CS","LMS","HMS","TSS","fast","slow","protocol","schedule","ramp_rate","target");assert not any(x in source for x in forbidden)
def test_21_reconstructed_positive_not_measured():
    q=pd.read_csv(OUT/"initial_distribution_families.csv");assert q.provenance.eq("synthetic_forward_state").all()
def test_22_candidate_is_response_target_only():
    q=pd.read_csv(OUT/"synthetic_distribution_boundary_test.csv");z=q[q.state_id.eq("candidate_response_target_like")];assert len(z) and z.candidate_response_target_only.all()
def test_23_figures_have_sources_no_placeholders():
    q=pd.read_csv(OUT/"figure_inventory.csv");assert (~q.placeholder).all();assert all(Path(x).exists() for x in q.pdf);assert all(Path(x).exists() for x in q.png);assert all(Path(x).exists() for x in q.source_table)
def test_24_reports_nonvalidation():
    names=["ZRO2_FORWARD_DISTRIBUTIONAL_PORE_MODEL.md","ZRO2_FORWARD_DISTRIBUTIONAL_EVOLUTION_LAWS.md","ZRO2_FORWARD_DISTRIBUTIONAL_TEST_RESULTS.md","ZRO2_FORWARD_DISTRIBUTIONAL_CHEN_TOPOLOGY.md","ZRO2_FORWARD_DISTRIBUTIONAL_NEXT_DECISION.md"]
    assert all("No validation claim" in (ROOT/"docs"/x).read_text() for x in names)
def test_25_staged_scope_excludes_unrelated_files():
    staged=subprocess.run(["git","diff","--cached","--name-only"],cwd=ROOT,text=True,capture_output=True,check=True).stdout.splitlines()
    allowed=("build_zro2_distributional_pore_population_model.py","run_zro2_distributional_pore_population_tests.py","plot_zro2_distributional_pore_population_results.py","tests/test_zro2_forward_distributional_pore_population_model.py","docs/ZRO2_FORWARD_DISTRIBUTIONAL_","results/zro2_forward_distributional_pore_population_model/")
    assert all(any(x==a or x.startswith(a) for a in allowed) for x in staged)
