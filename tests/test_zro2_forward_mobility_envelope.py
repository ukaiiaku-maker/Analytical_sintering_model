import hashlib,inspect,re
from pathlib import Path
import numpy as np,pandas as pd
from zro2_forward.material_zro2 import MaterialParameters
from zro2_forward.resolved_rules import ResolvedRuleModel,ResolvedRuleParameters,resolved_initial_state
from zro2_forward.conditioned_950c import make_pdf_conditioned_initial_state

ROOT=Path(__file__).parents[1];OUT=ROOT/'results/zro2_forward_mobility_envelope_audit'
def test_barrier_and_diffusivity_source_unchanged():
 assert hashlib.sha256((ROOT/'data/zro2/bicrystal_creep_barrier_export.json').read_bytes()).hexdigest()=='fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37';m=MaterialParameters();assert m.D_GB0_m2_s==.056 and m.Q_GB_J_mol==380000 and m.D_s0_m2_s==.10 and m.Q_s_J_mol==380000
def test_mobility_changes_growth_only_at_shared_state():
 s=resolved_initial_state(make_pdf_conditioned_initial_state());a=ResolvedRuleModel(parameters=ResolvedRuleParameters(M0_factor=.03));b=ResolvedRuleModel(parameters=ResolvedRuleParameters(M0_factor=30));ra=a.rates(s,1500);rb=b.rates(s,1500);assert np.isclose(ra[5]['rho_dot_open_sinv'],rb[5]['rho_dot_open_sinv']);assert np.isclose(ra[5]['activity'],rb[5]['activity']);assert rb[5]['G_dot_intrinsic_m_s']/ra[5]['G_dot_intrinsic_m_s']>900
def test_no_low_temperature_transition_or_schedule_labels():
 src=inspect.getsource(ResolvedRuleModel).lower();assert 'dong' not in src and 'mobility_transition' not in src;assert not any(re.search(rf'\b{x}\b',src) for x in ('cs','lms','hms','tss','fast','slow','schedule','protocol'))
def test_common_initial_state_and_strict_outputs_preserved():
 s=[resolved_initial_state(make_pdf_conditioned_initial_state()) for _ in range(4)];assert all(np.array_equal(s[0].pores.phi_open,x.pores.phi_open) for x in s[1:]);p=pd.read_csv(OUT/'mobility_chen_classification_points.csv');assert (p.strict_success==(p.density_ok&p.grain_ok)).all()
def test_windows_require_boundaries_and_are_uncertainty_not_validation():
 b=pd.read_csv(OUT/'mobility_chen_window_boundaries.csv');assert (~b.finite_window|(b.lower_boundary_present&b.upper_boundary_present&b.window_width_C.gt(0))).all();assert b.non_validation_flag.all()
def test_baseline_not_overwritten_and_parameters_recorded_everywhere():
 assert (ROOT/'results/zro2_forward_required_chen_physics_gap_analysis/boundary_gap_strict_summary.csv').exists()
 for name in ['mobility_envelope_design.csv','mobility_CS_conditioned_summary.csv','mobility_fast_rate_summary.csv','mobility_chen_classification_points.csv','mobility_chen_window_boundaries.csv','mobility_boundary_gap_summary.csv','mobility_pathway_consistency.csv']:
  x=pd.read_csv(OUT/name);assert {'gb_mobility_mode','M0_factor','Q_M_kJ_mol'}<=set(x)
