from pathlib import Path
import hashlib,inspect,subprocess,sys,re
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from run_zro2_forward_resolved_rules import params,initial
from zro2_forward.resolved_rules import ResolvedRuleModel,conservative_adjacent_PR
from zro2_forward.conditioned_950c import run_path
from zro2_forward.schedules import RampNoHold

OUT=ROOT/"results/zro2_forward_open_closed_rate_handoff_audit"
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def test_resolved_default_reproduction():
 f,_=run_path(ResolvedRuleModel(parameters=params(open_closed_handoff_mode="resolved_default")),RampNoHold(50,1500,start_C=950),initial(),60,60,"default")
 assert abs(f.rho-.9005037465)<2e-6 and abs(f.G_m*1e6-.258471024)<2e-6
def test_barrier_and_diffusivities_unchanged():
 assert sha("data/zro2/bicrystal_creep_barrier_export.json")=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
 assert sha("zro2_forward/material_zro2.py")=="35d6897595a8e8d94a53951c87d05564e984e479a2321cee330c860b789d5b1d"
def test_mobility_audit_not_modified():
 x=pd.read_csv(ROOT/"results/zro2_forward_mobility_envelope_audit/mobility_boundary_gap_summary.csv");assert x.boundary_gap_C.max()==-100 and x.strict_Chen_success_count.sum()==0
def test_PR_conservation():
 phi=np.array([.1,.2,.3]);f,c=conservative_adjacent_PR(phi,np.ones(3));assert abs(f.sum())<1e-15 and np.allclose(f,[-.1,-.1,.2])
def test_named_channel_and_transfer_identities():
 m=ResolvedRuleModel(parameters=params(open_closed_handoff_mode="balanced_handoff"));s=initial();od,ii,cc,closed,_,d=m.rates(s,1500+273.15)
 assert abs(d["rho_dot_total_sinv"]-d["rho_dot_open_sinv"]-d["rho_dot_closed_sinv"])<1e-15
 assert closed==0 and d["open_eligibility_eff"]==d["open_eligibility_base"]
 assert abs(od.sum()+ii.sum()+cc.sum()+d["rho_dot_open_sinv"]+d["rho_dot_closed_sinv"])<1e-10
def test_conservative_transfers_do_not_change_density():
 m=ResolvedRuleModel(parameters=params());s=initial();od,ii,cc,_,_,d=m.rates(s,1500+273.15)
 transfer_sum=od.sum()+ii.sum()+cc.sum()+d["rho_dot_open_sinv"]+d["rho_dot_closed_sinv"]
 assert abs(transfer_sum)<1e-10
def test_open_shrinkage_requires_open_inventory():
 s=inspect.getsource(ResolvedRuleModel.rates).replace(" ","")
 assert "open_shrink=-removal_weights(p)*rho_open" in s
def test_closed_shrinkage_requires_closed_inventory():
 m=ResolvedRuleModel(parameters=params());_,_,_,closed,_,_=m.rates(initial(),1500+273.15);assert closed==0
def test_balanced_handoff_requires_closed_availability():
 m=ResolvedRuleModel(parameters=params(open_closed_handoff_mode="balanced_handoff"));_,_,_,_,_,d=m.rates(initial(),1500+273.15)
 assert d["closed_availability"]==0 and d["open_eligibility_eff"]==d["open_eligibility_base"]
def test_open_and_closed_flux_integrals_sum():
 x=pd.read_csv(OUT/"fast_rate_flux_integrals.csv");assert (x.Delta_rho_open+x.Delta_rho_closed-x.Delta_rho_total).abs().max()<1e-10
def test_candidate_injection_is_diagnostic_only():
 x=pd.read_csv(OUT/"candidate_state_injection_diagnostic.csv");assert x.diagnostic_only.all() and not x.forward_prediction_eligible.any();assert (x.injected_closed_fraction==.65).all()
def test_ablation_and_window_rules():
 a=pd.read_csv(OUT/"handoff_ablation_summary.csv");assert not a.parent_window_present.any() and a.ablation_result.eq("not_interpretable_parent_has_no_window").all()
 b=pd.read_csv(OUT/"chen_handoff_window_boundaries.csv");assert not b.finite_window.any();assert ((~b.finite_window)|(b.lower_boundary_present&b.upper_boundary_present)).all()
def test_no_schedule_labels_in_handoff_local_law():
 s=inspect.getsource(ResolvedRuleModel.rates).lower()
 for token in ["cs","lms","hms","tss","fast","slow","protocol","schedule","ramp_rate","target"]:assert re.search(rf"\b{token}\b",s) is None
def test_old_mechanism_search_files_not_modified():
 names=subprocess.check_output(["git","diff","--name-only","008abd3de98d1540df0abd05ab602d0e9a90f602"],cwd=ROOT,text=True).splitlines()
 assert not any("mechanism_search" in n for n in names)
def test_reports_disclaim_validation():
 for p in ROOT.glob("docs/ZRO2_FORWARD_*HANDOFF*.md"):assert "not validation" in p.read_text().lower()
if __name__=="__main__":
 for n,v in sorted(globals().copy().items()):
  if n.startswith("test_") and callable(v):v();print("PASS",n)
