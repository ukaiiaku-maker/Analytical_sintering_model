from dataclasses import replace
from pathlib import Path
import hashlib,inspect,json,subprocess
import numpy as np,pandas as pd
import run_zro2_forward_closed_channel_physical_law_comparison as base
from zro2_forward.conditioned_950c import BARRIER
from zro2_forward.material_zro2 import MaterialParameters
from zro2_forward.resolved_rules import ResolvedRuleModel,conservative_adjacent_PR
from zro2_forward.closed_channel_laws import closed_channel_rates
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/zro2_forward_port_defined_closed_laws';SHA='fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37'
def state(mode='mechanistic_renewal_closed'):
 s=base.initial();p=s.pores.copy();shape=p.phi_open/p.phi_open.sum();p.phi_closed=.65*p.total*shape;p.phi_open=.35*p.total*shape;s=replace(s,rho=1-p.total,pores=p,A_closed=.152,PR_memory=1.);q=base.params(mechanism_mode='defined_closed_laws_port',closed_mapping_mode=mode);return s,ResolvedRuleModel(parameters=q)
def test_01_barrier_unchanged():assert hashlib.sha256((ROOT/BARRIER).read_bytes()).hexdigest()==SHA
def test_02_gb_diffusivity_unchanged():assert (MaterialParameters().D_GB0_m2_s,MaterialParameters().Q_GB_J_mol)==(.056,380000.)
def test_03_surface_diffusivity_unchanged():assert (MaterialParameters().D_s0_m2_s,MaterialParameters().Q_s_J_mol)==(.10,380000.)
def test_04_mobility_unchanged():assert (MaterialParameters().M0_m4_J_s,MaterialParameters().Q_M_J_mol)==(5.8e-3,4.2*96485.33212)
def test_05_no_physical_q_closed():
 x=pd.read_csv(OUT/'defined_closed_parameter_mapping.csv');assert not ((x.parameter_name.str.contains('Q_closed'))&x.parameter_class.isin(['fixed_ZrO2_input','literature_input'])).any()
def test_06_laws_registered():
 x=pd.read_csv(OUT/'defined_closed_law_registry.csv');required={'renewal_limited_closed_shrinkage','GB_diffusion_closed_shrinkage','surface_diffusion_accommodation_only','gas_accommodation_pressure_limit','candidate693168_reduced_closed_channel'};assert required<=set(x.law_id)
def test_07_phenomenological_labeled():
 x=pd.read_csv(OUT/'phenomenological_to_mechanistic_mapping_status.csv');assert x.mapping_status.notna().all() and (x.mapping_status=='phenomenological_unmapped').any()
def test_08_empirical_labeled():
 x=pd.read_csv(OUT/'defined_closed_parameter_mapping.csv');q=x[x.parameter_class.eq('empirical_diagnostic')];assert len(q)>0 and q.notes.str.contains('diagnostic').all()
def test_09_pr_transfer_conservative():
 p=state()[0].pores;f,_=conservative_adjacent_PR(p.phi_open,np.ones_like(p.phi_open));assert abs(f.sum())<1e-15
def test_10_closure_preparation_conservative():
 s,m=state();o,i,c,r,_,d=m.rates(s,1473.15);assert abs(o.sum()+i.sum()+c.sum()+d['rho_dot_open_sinv']+r)<1e-10
def test_11_open_shrinkage_only_open():
 s,m=state();s=replace(s,PR_memory=0.);o,i,c,r,_,d=m.rates(s,1473.15);assert np.all(np.array(json.loads(d['rho_dot_open_i_json']))>=0)
def test_12_closed_shrinkage_only_closed():
 s,m=state();s.pores.phi_closed[:]=0;s=replace(s,rho=1-s.pores.total);assert m.rates(s,1473.15)[3]==0
def test_13_surface_accommodation_non_densifying():
 s,m=state('mechanistic_surface_accommodation');assert m.rates(s,1473.15)[3]==0
def test_14_density_identity_stored_times():
 x=pd.read_csv(OUT/'closed_law_bin_histories.csv.gz',nrows=100)
 for r in x.itertuples():assert np.isclose(r.rho,1-sum(json.loads(r.phi_open_i))-sum(json.loads(r.phi_iso_i))-sum(json.loads(r.phi_closed_i)),atol=2e-8)
def test_15_accommodation_bounded():
 x=pd.read_csv(OUT/'closed_law_bin_histories.csv.gz',nrows=500)
 for r in x.itertuples():
  a=np.array(json.loads(r.A_closed_available_i));am=np.array(json.loads(r.A_closed_max_i));assert np.all(a>=0) and np.all(a<=am+1e-12)
def test_16_gas_pressure_nonnegative():
 s,m=state('mechanistic_gas_accommodation');d=m.rates(s,1473.15)[5];assert np.all(np.array(json.loads(d['P_gas_i_json']))>=0)
def test_17_injection_tagged():
 x=pd.read_csv(OUT/'candidate_state_T2_scan_by_mode.csv');assert x.candidate_state_injected.all() and x.diagnostic_only.all()
def test_18_windows_require_boundaries():
 x=pd.read_csv(OUT/'mini_map_window_boundaries.csv');assert len(x)==0 or ((~x.finite_window)|(x.lower_boundary_present&x.upper_boundary_present)).all()
def test_19_local_laws_have_no_path_labels():
 sig=str(inspect.signature(closed_channel_rates));assert not any(v in sig for v in ('schedule','protocol','ramp_rate','target','success'))
def test_20_no_old_search_modified():
 names=subprocess.check_output(['git','diff','--name-only','986aabb'],cwd=ROOT,text=True).splitlines();assert not any('mechanism_search' in n or n.endswith('_search.py') for n in names)
