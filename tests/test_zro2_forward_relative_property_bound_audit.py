from pathlib import Path
import hashlib,inspect,re,sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));OUT=ROOT/"results/zro2_forward_relative_property_bound_audit"
from zro2_forward.barrier_json import BarrierModel
from zro2_forward.resolved_rules import ResolvedRuleModel
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def test_barrier_hash_unchanged():assert sha("data/zro2/bicrystal_creep_barrier_export.json")=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
def test_GB_diffusivity_law_unchanged():
 s=(ROOT/"zro2_forward/material_zro2.py").read_text();assert "D_GB0_m2_s: float = 0.056" in s and "Q_GB_J_mol: float = 380000.0" in s
def test_surface_diffusivity_law_unchanged():
 s=(ROOT/"zro2_forward/material_zro2.py").read_text();assert "D_s0_m2_s: float = 0.10" in s and "Q_s_J_mol: float = 380000.0" in s
def test_high_temperature_mobility_law_unchanged():assert sha("zro2_forward/material_zro2.py")=="35d6897595a8e8d94a53951c87d05564e984e479a2321cee330c860b789d5b1d"
def test_source_windows_loaded_from_synthesis_source():
 x=pd.read_csv(OUT/"source_property_windows.csv");assert x.source.str.contains("material_property_window_scorecard|FINAL_MECHANISM_SYNTHESIS").all();assert len(x)>=13
def test_effective_Gstar_uses_stress_and_temperature():
 x=pd.read_csv(OUT/"effective_barrier_values_vs_T.csv");r=x.iloc[len(x)//3];b=BarrierModel.load(ROOT/"data/zro2/bicrystal_creep_barrier_export.json");assert abs(b.Gstar(r.sigma_eff_Pa,r.T_C+273.15)/1.602176634e-19-r.Gstar_eV)<1e-9;assert x.sigma_eff_Pa.nunique()>1
def test_Q_closed_proxy_status_explicit():
 x=pd.read_csv(OUT/"effective_rate_values_vs_T.csv");assert x.Q_closed_status.isin(["proxy/unidentified","effective finite-difference proxy"]).all();assert x.Q_closed_status.eq("proxy/unidentified").any()
def test_audit_does_not_modify_physics():assert sha("zro2_forward/resolved_rules.py")=="930a8eb5a85723d2da1367d54cfa407eaf3222494927c6dac4477c80c9009a74"
def test_no_schedule_labels_in_local_law():
 s=inspect.getsource(ResolvedRuleModel.rates).lower()
 for t in ["lms","hms","tss","fast","slow","protocol","schedule","ramp_rate","target"]:assert re.search(rf"\b{t}\b",s) is None
def test_reports_disclaim_validation():
 for p in ROOT.glob("docs/ZRO2_FORWARD_*PROPERTY*.md"):assert "not validation" in p.read_text().lower()
if __name__=="__main__":
 for n,v in sorted(globals().copy().items()):
  if n.startswith("test_") and callable(v):v();print("PASS",n)
