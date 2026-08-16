from dataclasses import replace
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import run_zro2_forward_closed_channel_physical_law_comparison as base
from zro2_forward.conditioned_950c import BARRIER
from zro2_forward.closed_pore_evolution import ClosedPoreEvolution
from zro2_forward.resolved_rules import ResolvedRuleModel
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/zro2_forward_closed_pore_evolution_addendum';SHA='fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37'
def injected(mode='mechanistic_renewal_limited',mexp=4):
 s=base.initial();p=s.pores.copy();shape=p.phi_open/p.phi_open.sum();p.phi_closed=.65*p.total*shape;p.phi_open=.35*p.total*shape;s=replace(s,rho=1-p.total,pores=p,A_closed=.152,PR_memory=1.);q=base.params(mechanism_mode='defined_laws_port',closed_mapping_mode=mode,closed_radius_exponent=mexp);return s,ResolvedRuleModel(parameters=q)
def test_fixed_inputs():
 assert hashlib.sha256((ROOT/BARRIER).read_bytes()).hexdigest()==SHA;s,m=injected();assert (m.material.Q_GB_J_mol,m.material.Q_s_J_mol)==(380000.,380000.)
def test_bin_state_and_density_identity():
 s,m=injected();n,d=m.step(s,1473.15,.01);assert isinstance(n.closed_pores,ClosedPoreEvolution) and np.isclose(n.rho,1-n.pores.total)
def test_number_radius_and_bounded_factors():
 s,m=injected();n,d=m.step(s,1473.15,.01);c=n.closed_pores;assert np.all(c.N_closed>=0) and np.all(c.radii_m(n.pores.phi_closed)>0) and np.all((c.chi_shrink>=0)&(c.chi_shrink<=1)) and np.all((c.C_geom>=0)&(c.C_geom<=1))
def test_finite_accommodation_consumed():
 s,m=injected();n,d=m.step(s,1473.15,1.);assert np.all(n.closed_pores.A_used>=0) and np.all(n.closed_pores.A_available()<=n.closed_pores.A_max)
def test_gas_pressure_opposes_capillarity():
 s,m=injected('mechanistic_gas_accommodation');_,d=m.step(s,1473.15,.01);pg=np.array(json.loads(d['P_gas_i_json']));sig=np.array(json.loads(d['sigma_closed_i_json']));assert np.all(pg>=0) and np.all(sig>=0) and d['gas_model']=='ideal_compression_proxy'
def test_surface_accommodation_no_density_removal():
 s,m=injected('mechanistic_surface_accommodation');_,_,_,r,_,d=m.rates(s,1473.15);assert r==0 and sum(json.loads(d['rho_dot_closed_i_json']))==0
def test_conservative_transfers_and_named_removal():
 s,m=injected();o,i,c,r,_,d=m.rates(s,1473.15);assert abs(o.sum()+i.sum()+c.sum()+d['rho_dot_open_sinv']+r)<1e-10
def test_all_required_bin_diagnostics_saved():
 x=pd.read_csv(OUT/'closed_pore_bin_histories.csv',nrows=1);required=['N_closed_i_json','r_closed_i_json','chi_shrink_i_json','C_geom_i_json','P_gas_i_json','sigma_closed_i_json','Gstar_closed_i_json','r_nuc_closed_i_json','tau_transport_closed_i_json','tau_cycle_closed_i_json','A_closed_i_json','A_closed_used_i_json','A_closed_recovered_i_json','rho_dot_closed_i_json','rho_dot_open_i_json'];assert set(required)<=set(x.columns)
def test_q_closed_app_diagnostic_only():
 x=pd.read_csv(OUT/'apparent_closed_rate_slopes.csv');assert x.diagnostic_only.all() and not x.material_input.any()
def test_acceptance_and_no_map():
 x=pd.read_csv(OUT/'boundary_topology_acceptance.csv');assert not x.plausible_topology.any();state=json.loads((OUT/'run_state.json').read_text());assert not state['mini_map_run'] and not state['validation']
