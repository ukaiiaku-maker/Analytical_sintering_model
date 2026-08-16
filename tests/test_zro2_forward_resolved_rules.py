import hashlib,inspect,re
from pathlib import Path
import numpy as np,pandas as pd
from zro2_forward.conditioned_950c import make_pdf_conditioned_initial_state
from zro2_forward.resolved_rules import ResolvedRuleModel,ResolvedRuleParameters,resolved_initial_state,conservative_adjacent_PR,ABLATIONS
from run_zro2_forward_resolved_rules import classify,window_rows

ROOT=Path(__file__).parents[1];OUT=ROOT/'results/zro2_forward_resolved_rules'
def state():return resolved_initial_state(make_pdf_conditioned_initial_state())
def test_fixed_inputs_and_locality():
 assert hashlib.sha256((ROOT/'data/zro2/bicrystal_creep_barrier_export.json').read_bytes()).hexdigest()=='fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37'
 src=inspect.getsource(ResolvedRuleModel.rates).lower();forbidden=('cs','lms','hms','tss','fast','slow','protocol','schedule','target','ramp_rate');assert not any(re.search(rf'\b{x}\b',src) for x in forbidden)
def test_PR_is_conservative_and_not_density_flux():
 phi=np.array([.1,.08,.04]);rate=np.array([.2,.1,.05]);flux,cross=conservative_adjacent_PR(phi,rate);assert abs(flux.sum())<1e-15
 src=inspect.getsource(conservative_adjacent_PR);assert 'rho' not in src
def test_named_stores_and_density_identity_in_histories():
 h=pd.read_csv(OUT/'resolved_rule_smoke_histories.csv');assert np.allclose(h.rho,1-h.open_fraction-h.isolated_fraction-h.closed_fraction,atol=2e-6);assert (h[['open_fraction','isolated_fraction','closed_fraction']]>=0).all().all()
def test_open_closed_isolated_rules_and_bounded_accommodation():
 src=inspect.getsource(ResolvedRuleModel.rates);assert 'p.phi_closed/tau_closed' in src and 'p.phi_open' in src and 'p.phi_iso' in src;h=pd.read_csv(OUT/'resolved_rule_smoke_histories.csv');assert h.A_closed.between(0,1).all()
def test_migration_modifiers_and_M0_do_not_change_shared_state_density_rate():
 s=state();a=ResolvedRuleModel(parameters=ResolvedRuleParameters(M0_factor=.1));b=ResolvedRuleModel(parameters=ResolvedRuleParameters(M0_factor=10,no_pore_drag=True));ra=a.rates(s,1500);rb=b.rates(s,1500);assert np.isclose(ra[5]['rho_dot_open_sinv'],rb[5]['rho_dot_open_sinv']);assert np.isclose(ra[5]['rho_dot_closed_sinv'],rb[5]['rho_dot_closed_sinv']);assert rb[4]['G_dot_intrinsic_m_s']/ra[4]['G_dot_intrinsic_m_s']>90
def test_strict_chen_and_window_rules():
 class F:rho=.976;G_m=.29e-6
 class X:rho=.8
 assert classify(F(),X(),True)=='SUCCESS';F.rho=.95;assert classify(F(),X(),True)!='SUCCESS'
 x=pd.DataFrame([{'T1_C':1400,'T2_C':1100,'switch_density':.8,'hold_h':20,'classification':'DENSIFICATION_EXHAUSTION_FAILURE','density_ok':False,'grain_ok':True},{'T1_C':1400,'T2_C':1200,'switch_density':.8,'hold_h':20,'classification':'SUCCESS','density_ok':True,'grain_ok':True},{'T1_C':1400,'T2_C':1300,'switch_density':.8,'hold_h':20,'classification':'GRAIN_GROWTH_FAILURE','density_ok':True,'grain_ok':False}]);assert not window_rows(x).finite_window.iloc[0]
def test_ablation_inventory_and_prior_baseline_separate():
 a=pd.read_csv(OUT/'resolved_rule_ablation_summary.csv');assert set(a.ablation)==set(ABLATIONS);assert set(ABLATIONS[:4])==set(a[a.controlling_expected].ablation);assert (ROOT/'results/zro2_forward_required_chen_physics_gap_analysis').exists()
