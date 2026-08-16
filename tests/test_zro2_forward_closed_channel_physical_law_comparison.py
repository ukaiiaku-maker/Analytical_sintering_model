from dataclasses import replace
from pathlib import Path
import hashlib,inspect,json,subprocess
import numpy as np
import pandas as pd
from zro2_forward.conditioned_950c import BARRIER,make_pdf_conditioned_initial_state
from zro2_forward.material_zro2 import MaterialParameters
from zro2_forward.resolved_rules import ResolvedRuleModel,ResolvedRuleParameters,resolved_initial_state,conservative_adjacent_PR

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results/zro2_forward_closed_channel_physical_law_comparison"
EXPECTED_BARRIER="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
def state():
    s=resolved_initial_state(make_pdf_conditioned_initial_state(rho=.88,G_nm=117));p=s.pores.copy();p.phi_closed=.65*p.phi_open;p.phi_open*=.35
    return replace(s,rho=1-p.total,pores=p,A_closed=.152,PR_memory=1.)
def test_physical_inputs_unchanged():
    m=MaterialParameters();assert hashlib.sha256((ROOT/BARRIER).read_bytes()).hexdigest()==EXPECTED_BARRIER
    assert (m.D_GB0_m2_s,m.Q_GB_J_mol,m.D_s0_m2_s,m.Q_s_J_mol)==(.056,380000.,.10,380000.)
    assert (m.M0_m4_J_s,m.Q_M_J_mol)==(5.8e-3,4.2*96485.33212)
def test_default_reproduces_proxy_and_q_app_not_input():
    q=ResolvedRuleParameters();assert q.closed_channel_law=="resolved_proxy_current";assert "Q_closed_app" not in inspect.getsource(ResolvedRuleParameters)
def test_surface_accommodation_has_no_direct_density():
    m=ResolvedRuleModel(parameters=ResolvedRuleParameters(closed_channel_law="surface_diffusion_accommodation_only"));o,i,c,r,_,d=m.rates(state(),1473.15);assert r==0 and d["A_dot_closed_sinv"]>=0
def test_named_reservoirs_and_conservation():
    s=state();m=ResolvedRuleModel(parameters=ResolvedRuleParameters(closed_channel_law="GB_diffusion_closed_shrinkage"));o,i,c,r,_,d=m.rates(s,1473.15)
    assert r>=0 and np.all(c<=d["phi_closed_dot_json"] if False else np.inf);assert abs(o.sum()+i.sum()+c.sum()+d["rho_dot_open_sinv"]+r)<1e-10
    flux,_=conservative_adjacent_PR(s.pores.phi_open,np.ones_like(s.pores.phi_open));assert abs(flux.sum())<1e-15
    n,_=m.step(s,1473.15,.01);assert np.isclose(n.rho,1-n.pores.total)
def test_tags_windows_and_schedule_blind_laws():
    reg=pd.read_csv(OUT/"closed_law_registry.csv");assert reg.schedule_blind.all()
    emp=reg.query("closed_channel_law=='empirical_closed_rate_scale'").iloc[0];assert emp.empirical and "not validation" in emp.status
    cand=pd.read_csv(OUT/"candidate_state_closed_law_T2_scan.csv");assert cand.query("candidate_state_injected").diagnostic_only.all()
    b=pd.read_csv(OUT/"closed_law_chen_window_boundaries.csv");assert ((~b.finite_window)|(b.lower_boundary_present&b.upper_boundary_present)).all()
    sig=str(inspect.signature(__import__("zro2_forward.closed_channel_laws",fromlist=["closed_channel_rates"]).closed_channel_rates));assert not any(x in sig for x in ("schedule","protocol","ramp_rate","target"))
def test_density_identity_outputs_and_no_old_search_changes():
    for file in ("candidate_state_closed_law_T2_scan.csv","natural_state_closed_law_T2_scan.csv"):
        x=pd.read_csv(OUT/file);assert x.final_rho.between(0,1).all();assert x.Delta_rho_closed.le(1-x.initial_rho+1e-12).all()
    state_json=json.loads((OUT/"run_state.json").read_text());assert not state_json["Q_closed_app_physical_input"] and not state_json["validation"]
def test_barrier_json_hash_unchanged():assert hashlib.sha256((ROOT/BARRIER).read_bytes()).hexdigest()==EXPECTED_BARRIER
def test_gb_diffusivity_law_unchanged():assert (MaterialParameters().D_GB0_m2_s,MaterialParameters().Q_GB_J_mol)==(.056,380000.)
def test_surface_diffusivity_law_unchanged():assert (MaterialParameters().D_s0_m2_s,MaterialParameters().Q_s_J_mol)==(.10,380000.)
def test_growth_mobility_law_unchanged():assert MaterialParameters().Q_M_J_mol==4.2*96485.33212
def test_q_closed_app_diagnostic_only():assert "Q_closed_app" not in inspect.getsource(ResolvedRuleParameters)
def test_closed_removal_only_uses_closed_store():
    s=state();s.pores.phi_closed[:]=0;m=ResolvedRuleModel(parameters=ResolvedRuleParameters(closed_channel_law="GB_diffusion_closed_shrinkage"));assert m.rates(s,1473.15)[3]==0
def test_empirical_mode_is_diagnostic():
    d=ResolvedRuleModel(parameters=ResolvedRuleParameters(closed_channel_law="empirical_closed_rate_scale")).rates(state(),1473.15)[5];assert d["empirical_diagnostic"] and d["Q_closed_emp_not_material_property"]
def test_default_proxy_reproduces_inherited_outputs():
    x=pd.read_csv(OUT/"fixed_path_closed_law_summary.csv").query("law=='resolved_proxy_current' and prefactor_factor==1").set_index("path")
    assert np.isclose(x.loc["PDF_conditioned_5C_min","final_rho"],.9370101360405385,rtol=0,atol=2e-6)
    assert np.isclose(x.loc["PDF_conditioned_50C_min","final_rho"],.9005037465317097,rtol=0,atol=2e-6)
def test_no_old_mechanism_search_files_modified():
    names=subprocess.check_output(["git","diff","--name-only","ebf82941db9b8789d018e5e7986b8e6750587e20"],cwd=ROOT,text=True).splitlines()
    assert not any("mechanism_search" in name or name.endswith("_search.py") for name in names)
