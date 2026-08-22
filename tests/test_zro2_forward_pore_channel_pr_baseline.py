from __future__ import annotations
import hashlib,inspect,re,subprocess
from pathlib import Path
import numpy as np,pandas as pd
import test_zro2_pore_channel_pr_baseline as p

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results/zro2_forward_pore_channel_pr_baseline_test"

def test_01_barrier_hash(): assert hashlib.sha256(p.BARRIER_PATH.read_bytes()).hexdigest()=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
def test_02_03_diffusivities():
 T=1473.15; R=8.31446261815324
 assert np.isclose(p.MAT.D_GB(T),.056*np.exp(-380000/(R*T))) and np.isclose(p.MAT.D_s(T),.10*np.exp(-380000/(R*T)))
def test_04_failed_mobility_inactive(): assert p.MAT.M0_m4_J_s==5.8e-3 and "calibrated once" in p.MAT.mobility_prefactor_status
def test_05_06_no_physical_qclosed():
 q=pd.read_csv(OUT/"pore_channel_parameter_registry.csv"); x=q[q.parameter=="Q_closed_app"].iloc[0]
 assert x.classification=="empirical diagnostic" and not (q.classification.eq("physical")&q.parameter.str.contains("closed",case=False)).any()
def test_07_08_topology_does_not_densify():
 z=p.pore_channel_terms(1473.15,25e-9,10,.5,.1,.8,.6)
 assert z["rho_dot_topology"]==0
 for m in p.MODES:
  s=p.conservative_step({"open":.2,"precursor":.03,"isolated":.04,"closed":.05},m,z,100); assert abs(s["conservation_residual"])<1e-14
def test_09_regularization_damage_conserve():
 z=p.pore_channel_terms(1473.15,25e-9,10,.5,.1,.8,.6); s=p.conservative_step({"open":.2,"precursor":.03,"isolated":.04,"closed":.05},"PR_regularization_damage_v1",z,1e5)
 assert abs(s["conservation_residual"])<1e-14 and all(s[k]>=0 for k in ("open","precursor","isolated","closed"))
def test_10_density_identity():
 b=pd.read_csv(OUT/"pore_channel_boundary_preservation_test.csv"); assert abs(b.density_identity_residual).max()<1e-10
def test_11_12_closed_limits():
 for k in ("renewal","GB_diffusion"):
  assert p.emergent_pore_closure_v1(1473.15,25e-9,0,1,.5,.25,3,kernel=k)["rho_dot_closed_sinv"]==0
  assert p.emergent_pore_closure_v1(1473.15,25e-9,.05,1,.5,1.1,3,kernel=k)["rho_dot_closed_sinv"]==0
def test_13_gas_reduces_stress():
 a=p.emergent_pore_closure_v1(1473.15,25e-9,.05,1,.5,0,3); b=p.emergent_pore_closure_v1(1473.15,25e-9,.05,1,.5,.9,3); assert b["sigma_c_Pa"]<a["sigma_c_Pa"]
def test_14_r4_time():
 a=p.pore_channel_terms(1473.15,10e-9,10,.5,.1,.8,.6); b=p.pore_channel_terms(1473.15,100e-9,10,.5,.1,.8,.6); assert np.isclose(b["tau_s_s"]/a["tau_s_s"],1e4)
def test_15_smaller_pores_pin_more():
 a=p.pore_channel_terms(1473.15,10e-9,10,.5,.1,.8,.6); b=p.pore_channel_terms(1473.15,100e-9,10,.5,.1,.8,.6); assert a["P_Z_Pa"]>b["P_Z_Pa"]
def test_16_migration_does_not_change_density():
 h=pd.read_csv(OUT/"pore_channel_heating_rate_histories.csv"); assert {"Gamma_migration","rho_dot_total"}.issubset(h.columns)
 assert pd.read_csv(OUT/"pore_channel_energy_ledger_registry.csv").forced_equality.eq(False).all()
def test_17_cloned_T2_states():
 q=pd.read_csv(OUT/"pore_channel_boundary_preservation_test.csv")
 for _,z in q.groupby(["state_id","mode"]):
  for c in ("rho","G_nm","open","precursor","isolated","closed","r_nm","W","conn","A"): assert z[c].nunique()==1
def test_18_window_requires_boundaries():
 q=pd.read_csv(OUT/"pore_channel_window_boundaries.csv"); assert ((~q.strict_finite_window)|(q.lower_boundary&q.upper_boundary&(q.success_points>1))).all()
def test_19_no_local_processing_labels():
 source=(inspect.getsource(p.pore_channel_terms)+inspect.getsource(p.conservative_step)).lower(); tokens=set(re.findall(r"[a-z_]+",source)); forbidden={"cs","lms","hms","tss","fast","slow","protocol","schedule","ramp_rate","target"}; assert not tokens&forbidden
def test_20_empirical_fallback_diagnostic_only():
 q=pd.read_csv(OUT/"pore_channel_parameter_registry.csv"); assert q[q.classification.str.contains("empirical")].mapping.str.contains("post-run").all()
def test_21_figures_sourced_no_placeholders():
 q=pd.read_csv(OUT/"figure_inventory.csv"); assert len(q)==9 and q.pdf_nonempty.all() and q.png_nonempty.all() and not q.placeholder.any() and not q.success_colored_map.any()
 for x in q.source_table: assert (OUT/x).is_file()
def test_22_reports_nonvalidation():
 names=["ZRO2_FORWARD_PORE_CHANNEL_PR_BASELINE_TEST.md","ZRO2_FORWARD_PORE_CHANNEL_CONSTITUTIVE_MODEL.md","ZRO2_FORWARD_SURFACE_DIFFUSION_PR_LITERATURE_MAPPING.md","ZRO2_FORWARD_PORE_CHANNEL_TEST_RESULTS.md","ZRO2_FORWARD_PORE_CHANNEL_NEXT_DECISION.md"]
 for n in names: assert "not validation" in (ROOT/"docs"/n).read_text().lower() or "no validation claim" in (ROOT/"docs"/n).read_text().lower()
def test_23_staged_scope():
 names=subprocess.run(["git","diff","--cached","--name-only"],cwd=ROOT,text=True,capture_output=True,check=True).stdout.splitlines(); prefix="results/zro2_forward_pore_channel_pr_baseline_test/"
 for n in names:
  assert not n.endswith(".DS_Store") and not (n.lower().endswith(".pdf") and not n.startswith(prefix))
  if n.startswith("results/"): assert n.startswith(prefix)
