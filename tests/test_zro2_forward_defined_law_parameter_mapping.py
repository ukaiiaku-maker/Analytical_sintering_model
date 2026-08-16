from dataclasses import replace
from pathlib import Path
import hashlib,inspect,json,subprocess
import numpy as np,pandas as pd
from zro2_forward.conditioned_950c import BARRIER,make_pdf_conditioned_initial_state
from zro2_forward.grain_growth import growth_state
from zro2_forward.material_zro2 import MaterialParameters
from zro2_forward.resolved_rules import ResolvedRuleModel,ResolvedRuleParameters,resolved_initial_state,conservative_adjacent_PR
from zro2_forward.closed_channel_laws import closed_channel_rates
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results/zro2_forward_defined_law_parameter_mapping";SHA="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
def state():
 s=resolved_initial_state(make_pdf_conditioned_initial_state(rho=.88,G_nm=117));p=s.pores.copy();p.phi_closed=.65*p.phi_open;p.phi_open*=.35;return replace(s,rho=1-p.total,pores=p,A_closed=.152,PR_memory=1.)
def model(**kw):return ResolvedRuleModel(parameters=replace(ResolvedRuleParameters(),mechanism_mode="defined_laws_port",closed_mapping_mode="defined_laws_port",**kw))
def test_01_barrier_hash_unchanged():assert hashlib.sha256((ROOT/BARRIER).read_bytes()).hexdigest()==SHA
def test_02_gb_diffusivity_unchanged():assert (MaterialParameters().D_GB0_m2_s,MaterialParameters().Q_GB_J_mol)==(.056,380000.)
def test_03_surface_diffusivity_unchanged():assert (MaterialParameters().D_s0_m2_s,MaterialParameters().Q_s_J_mol)==(.10,380000.)
def test_04_growth_mobility_unchanged():assert (MaterialParameters().M0_m4_J_s,MaterialParameters().Q_M_J_mol)==(5.8e-3,4.2*96485.33212)
def test_05_q_closed_not_physical_input():
 x=pd.read_csv(OUT/"defined_law_parameter_mapping.csv");q=x[x.parameter.eq("Q_closed_app")].iloc[0];assert q.parameter_class=="empirical_diagnostic" and not q.physical_ZrO2_input
def test_06_every_law_has_source():
 x=pd.read_csv(OUT/"defined_law_registry.csv");assert len(x)==12 and x.source_file.notna().all() and x.current_forward_function.notna().all()
def test_07_nonphysical_parameters_labeled():
 x=pd.read_csv(OUT/"defined_law_parameter_mapping.csv");bad=x[~x.physical_ZrO2_input & x.requires_calibration];assert bad.parameter_class.isin(["global_calibration","bounded_uncertainty","reduced_phenomenological","empirical_diagnostic","missing_physical_mapping"]).all()
def test_08_pr_transfer_conserves_volume():
 p=state().pores;f,_=conservative_adjacent_PR(p.phi_open,np.ones_like(p.phi_open));assert abs(f.sum())<1e-15
def test_09_pr_transfer_no_direct_density():
 p=state().pores;f,_=conservative_adjacent_PR(p.phi_open,np.ones_like(p.phi_open));assert np.isclose(1-(p.total+f.sum()*.01),1-p.total)
def test_10_open_shrinkage_only_open_store():
 s=replace(state(),PR_memory=0.);od,ii,_,_,_,d=model(no_PR_damage=True,no_closed_transition=True).rates(s,1473.15);assert np.isclose(od.sum()+d["rho_dot_open_sinv"],0) and np.all(ii==0)
def test_11_closed_shrinkage_only_closed_store():
 s=state();s.pores.phi_closed[:]=0;s=replace(s,rho=1-s.pores.total);assert model().rates(s,1473.15)[3]==0
def test_12_isolated_store_not_open_removed():
 s=replace(state(),PR_memory=0.);before=s.pores.phi_iso.copy();od,ii,_,_,_,d=model(no_PR_damage=True,no_closed_transition=True).rates(s,1473.15);assert d["rho_dot_open_sinv"]>=0 and np.all(ii==0) and np.array_equal(before,s.pores.phi_iso)
def test_13_density_identity_every_step():
 s=state();m=model()
 for _ in range(20):s,_=m.step(s,1473.15,.01);assert np.isclose(s.rho,1-s.pores.total)
def test_14_migration_modifiers_do_not_change_density_rate():
 s=state();a=model(M0_factor=.1).rates(s,1473.15)[5];b=model(M0_factor=10.,no_pore_drag=True).rates(s,1473.15)[5];assert np.isclose(a["rho_dot_total_sinv"],b["rho_dot_total_sinv"])
def test_15_smaller_pores_pin_more():
 m=MaterialParameters();phi=np.array([.05]);small=growth_state(1e-6,np.array([10e-9]),phi,1473.15,m);large=growth_state(1e-6,np.array([100e-9]),phi,1473.15,m);assert small["P_Z_Pa"]>large["P_Z_Pa"] and small["Gamma_growth"]<large["Gamma_growth"]
def test_16_finite_window_requires_boundaries():
 x=pd.read_csv(OUT/"mini_map_window_boundaries.csv");assert len(x)==0 or ((~x.finite_window)|(x.lower_boundary_present&x.upper_boundary_present)).all()
def test_17_injection_diagnostic_tag():
 x=pd.read_csv(OUT/"candidate_state_T2_scan_by_mode.csv");assert x.candidate_state_injected.all() and x.diagnostic_only.all()
def test_18_empirical_tag_nonvalidated():
 x=pd.read_csv(OUT/"candidate693168_defined_law_comparison.csv").query("model=='empirical_rate_scale_diagnostic'");assert not x.validated.any();p=pd.read_csv(OUT/"defined_law_parameter_mapping.csv");assert (p.parameter_class=="empirical_diagnostic").any()
def test_19_local_law_has_no_path_labels():
 sig=str(inspect.signature(closed_channel_rates));assert not any(v in sig for v in ("schedule","protocol","ramp_rate","target","success"))
def test_20_no_old_search_files_modified():
 names=subprocess.check_output(["git","diff","--name-only","83eab35f95a90e0d65db05d9cddee8dced8cf55b"],cwd=ROOT,text=True).splitlines();assert not any("mechanism_search" in n or n.endswith("_search.py") for n in names)
