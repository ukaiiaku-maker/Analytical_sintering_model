from pathlib import Path
import hashlib
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/zro2_forward_resolved_rule_rate_balance_audit"
def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def test_model_physics_files_unchanged():
    assert sha("zro2_forward/resolved_rules.py")=="ef0415c41f4041f4294eeebd0ad9c7b8622983867f25e53257370bbb8ad998e4"
def test_barrier_hash_unchanged():
    assert sha("data/zro2/bicrystal_creep_barrier_export.json")=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"
def test_diffusivity_parameters_unchanged():
    assert sha("zro2_forward/material_zro2.py")=="35d6897595a8e8d94a53951c87d05564e984e479a2321cee330c860b789d5b1d"
def test_mobility_audit_preserved():
    x=pd.read_csv(ROOT/"results/zro2_forward_mobility_envelope_audit/mobility_boundary_gap_summary.csv")
    assert x.strict_Chen_success_count.sum()==0 and x.finite_strict_Chen_window_count.sum()==0
    assert x.boundary_gap_C.max()==-100
def test_ablation_interpretation_requires_parent_window():
    x=pd.read_csv(OUT/"resolved_rule_ablation_summary_reinterpreted.csv")
    assert not x.parent_window_present.any() and not x.ablation_interpretable.any()
    assert x.reason_not_interpretable.eq("not_interpretable_parent_has_no_window").all()
def test_density_identity_and_flux_contribution_sum():
    x=pd.read_csv(OUT/"trajectory_integral_summary.csv").query("source_layer=='resolved'")
    assert ((x.Delta_rho_open+x.Delta_rho_closed-x.Delta_rho_total_flux).abs()<1e-10).all()
    h=pd.read_csv(OUT/"density_flux_comparison.csv").query("source_layer=='resolved'")
    assert (h.rho+h.phi_open+h.phi_precursor+h.phi_closed-1).abs().max()<1e-10
    # Recorded-rate trapezoids approximate adaptive-step state increments.
    assert (x.Delta_rho_total_flux-x.Delta_rho_state).abs().max()<7e-3
def test_no_schedule_labels_in_local_law():
    s=(ROOT/"zro2_forward/resolved_rules.py").read_text()
    assert all(token not in s for token in ["5C", "50C", "two_step", "schedule_label"])
def test_reports_explicitly_disclaim_validation():
    for p in ROOT.glob("docs/ZRO2_FORWARD_*RATE_BALANCE*.md"):
        assert "not validation" in p.read_text().lower()
    for name in ["ZRO2_FORWARD_RESOLVED_RULE_FAILURE_DECOMPOSITION.md","ZRO2_FORWARD_RESOLVED_RULE_ABLATION_REINTERPRETATION.md","ZRO2_FORWARD_NEXT_IMPLEMENTATION_DECISION.md"]:
        assert "not validation" in (ROOT/"docs"/name).read_text().lower()

if __name__=="__main__":
    for name,value in sorted(globals().copy().items()):
        if name.startswith("test_") and callable(value): value(); print(f"PASS {name}")
