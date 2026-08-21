from pathlib import Path
import hashlib,inspect,subprocess
import numpy as np,pandas as pd
from zro2_forward.material_zro2 import MaterialParameters
from zro2_forward.resolved_rules import ResolvedRuleModel
from zro2_forward.pore_population import initial_population
import run_zro2_forward_processing_window_prediction as run
ROOT=Path(__file__).parents[1];OUT=ROOT/"results/zro2_forward_processing_window_prediction_figures"
def test_01_barrier_hash_unchanged():assert hashlib.sha256((ROOT/"data/zro2/bicrystal_creep_barrier_export.json").read_bytes()).hexdigest()=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
def test_02_gb_diffusivity_unchanged():m=MaterialParameters();assert (m.D_GB0_m2_s,m.Q_GB_J_mol)==(.056,380000.)
def test_03_surface_diffusivity_unchanged():m=MaterialParameters();assert (m.D_s0_m2_s,m.Q_s_J_mol)==(.10,380000.)
def test_04_closed_laws_unchanged():assert not subprocess.check_output(["git","diff","--name-only","2249f02","--","zro2_forward/closed_channel_laws.py","zro2_forward/closed_pore_evolution.py"],cwd=ROOT,text=True).strip()
def test_05_failed_mobility_not_active():s=pd.read_json(OUT/"run_state.json",typ="series");assert not s.failed_global_mobility_active and run.params().M0_factor==1 and run.params().Q_M_J_mol_override is None
def test_06_no_optimizer_or_fit_called():s=(ROOT/"run_zro2_forward_processing_window_prediction.py").read_text().lower();assert not any(x in s for x in ("scipy.optimize","least_squares","differential_evolution","curve_fit"))
def test_07_candidate_values_not_copied():s=pd.read_csv(OUT/"prior_solution_manifest.csv").iloc[0];assert not s.state_values_copied and not s.use_as_ZrO2_parameter and pd.read_csv(OUT/"initial_state_grid.csv").initial_state_id.str.startswith("S").all()
def test_08_candidate_prior_only():s=pd.read_csv(OUT/"prior_solution_manifest.csv").iloc[0];assert s.response_form_target_only and not s.parameterized_as_ZrO2
def test_09_search_axes_present():x=pd.read_csv(OUT/"first_step_processing_state_screen.csv");assert {"rho0","G0_nm","rho_switch","G1_nm"}<=set(x)
def test_10_identical_cloned_first_state_across_T2():x=pd.read_csv(OUT/"twostep_second_step_classification_points.csv");assert x.groupby("case_id").switch_state_hash.nunique().eq(1).all()
def test_11_continuous_second_step_time():x=pd.read_csv(OUT/"twostep_second_step_classification_points.csv");assert (x.physical_start_time_s>0).all() and (x.physical_end_time_s>x.physical_start_time_s).all()
def test_12_named_density_fluxes_only():src=inspect.getsource(ResolvedRuleModel.rates);assert "rho_dot_open_sinv" in src and "rho_dot_closed_sinv" in src and "rho_dot_total_sinv" in src
def test_13_pr_transfer_conserves_volume():
 p=initial_population(rho0=.66);f,c=run.ResolvedRuleModel.__mro__[0].__module__,None
 from zro2_forward.resolved_rules import conservative_adjacent_PR
 flux,_=conservative_adjacent_PR(p.phi_open,np.ones_like(p.phi_open)*1e-3);assert abs(flux.sum())<1e-14
def test_14_surface_accommodation_nondensifying():s=(ROOT/"zro2_forward/closed_channel_laws.py").read_text();assert "surface_diffusion_accommodation_only" in s and "A_dot_closed_sinv" in s
def test_15_no_schedule_labels_in_local_rate_law():
 src=inspect.getsource(ResolvedRuleModel.rates).lower();assert not any(x in src for x in ('"cs"','"lms"','"hms"','"tss"','protocol','ramp_rate','schedule'))
def test_16_finite_window_has_both_boundaries():x=pd.read_csv(OUT/"twostep_window_boundaries.csv");f=x[x.finite_window];assert (f.lower_boundary_present&f.upper_boundary_present&f.window_width_C.gt(0)&f.success_count.ge(2)).all()
def test_17_figure_qc_all_pass():assert pd.read_csv(OUT/"figure_qc_report.csv").qc_pass.all()
def test_18_every_figure_has_source():
 x=pd.read_csv(OUT/"figure_inventory.csv");assert x.source_table.notna().all() and all((OUT/p).exists() for p in x.source_table)
def test_19_reports_nonvalidation():
 for p in ROOT.glob("docs/ZRO2_FORWARD_*PROCESSING_WINDOW*.md"):assert "not validation" in p.read_text().lower() or "not validated" in p.read_text().lower()
def test_20_no_forbidden_staged_files():
 x=subprocess.check_output(["git","diff","--cached","--name-only"],cwd=ROOT,text=True).splitlines();assert not any(".DS_Store" in p or "mechanism_search" in p or (p.lower().endswith(".pdf") and not p.startswith("results/zro2_forward_processing_window_prediction_figures/")) for p in x)
