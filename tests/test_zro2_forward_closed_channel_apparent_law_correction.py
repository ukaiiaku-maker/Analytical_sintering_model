from pathlib import Path
import hashlib,json,subprocess,sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"results/zro2_forward_closed_channel_apparent_law_correction"
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def test_barrier_unchanged():assert sha("data/zro2/bicrystal_creep_barrier_export.json")=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
def test_GB_diffusivity_unchanged():
 s=(ROOT/"zro2_forward/material_zro2.py").read_text();assert "D_GB0_m2_s: float = 0.056" in s and "Q_GB_J_mol: float = 380000.0" in s
def test_surface_diffusivity_unchanged():
 s=(ROOT/"zro2_forward/material_zro2.py").read_text();assert "D_s0_m2_s: float = 0.10" in s and "Q_s_J_mol: float = 380000.0" in s
def test_mobility_unchanged():assert sha("zro2_forward/material_zro2.py")=="35d6897595a8e8d94a53951c87d05564e984e479a2321cee330c860b789d5b1d"
def test_reports_state_not_physical_input():
 paths=list(ROOT.glob("docs/ZRO2_FORWARD_*CLOSED_CHANNEL*.md"))+list(ROOT.glob("docs/ZRO2_FORWARD_QCLOSED*.md"));text="\n".join(p.read_text() for p in paths);assert "not an independently defined physical" in text and "not validation" in text
def test_corrected_tables_use_apparent_terms():
 for p in OUT.glob("*.csv"):
  cols=pd.read_csv(p,nrows=0).columns
  if p.name!="closed_channel_terminology_correction.csv":assert "Q_closed" not in cols and "Q_closed_eff" not in cols
def test_apparent_slopes_diagnostic_only():
 x=pd.read_csv(OUT/"apparent_closed_channel_slope_derivation.csv");assert x.diagnostic_only.all() and x.not_material_property.all()
def test_no_apparent_quantity_is_physical_pass_fail():
 x=pd.read_csv(OUT/"corrected_statewise_property_window_classification.csv");assert x.closed_channel_classification.isin(["apparent_proxy_not_physical_property","insufficient_state_support"]).all()
def test_candidate_laws_registry_only():
 x=pd.read_csv(OUT/"closed_channel_candidate_law_registry.csv");assert x.validation_status.str.contains("not installed").all();assert sha("zro2_forward/resolved_rules.py")=="930a8eb5a85723d2da1367d54cfa407eaf3222494927c6dac4477c80c9009a74"
def test_required_factors_analysis_only():assert pd.read_csv(OUT/"required_closed_rate_magnitude.csv").interpretation.str.contains("not applied").all()
def test_candidate_is_conditional_not_validated():assert pd.read_csv(OUT/"closed_channel_candidate693168_comparison.csv").query("model=='candidate_693168'").status.str.contains("conditional comparator, not validation").all()
def test_all_reports_disclaim_validation():
 for p in ROOT.glob("docs/ZRO2_FORWARD_*CLOSED_CHANNEL*.md"):assert "not validation" in p.read_text().lower()
def test_no_old_mechanism_search_modified():
 names=subprocess.check_output(["git","diff","--name-only","d4fd1b3128ccda910e8955a1c955d06c1fb6d349"],cwd=ROOT,text=True).splitlines();assert not any("mechanism_search" in n for n in names)
def test_no_broad_parameter_search():
 for p in [ROOT/"correct_zro2_forward_closed_channel_property_interpretation.py",ROOT/"analyze_zro2_forward_closed_channel_physical_mapping.py"]:assert "mechanism search" not in p.read_text().lower() and "parameter_grid" not in p.read_text()
if __name__=="__main__":
 for n,v in sorted(globals().copy().items()):
  if n.startswith("test_") and callable(v):v();print("PASS",n)
