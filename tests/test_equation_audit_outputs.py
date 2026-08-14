"""Integrity tests for the documentation-only equation audit."""
from pathlib import Path
import subprocess

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "equation_functional_form_audit"
PROTECTED_MODEL_SOURCES = {
    "separated_fast_chen_model.py",
    "pr_desintering_memory_model.py",
    "pore_location_topology_model.py",
    "pore_location_agentic_model.py",
    "agentic_mechanism_model.py",
    "interacting_local_region_model.py",
    "mechanism_dimensionless_groups.py",
    "relative_material_property_window_search.py",
    "adaptive_T2_boundary_search.py",
    "audit_candidate_693168_closed_accommodation.py",
}


def test_registry_is_complete_and_traceable():
    path = OUT / "equation_registry.csv"
    assert path.exists()
    rows = pd.read_csv(path, keep_default_na=False)
    assert len(rows) >= 30
    assert rows.equation_id.is_unique
    assert rows.source_file.ne("").all()
    assert rows.source_function.ne("").all()
    assert rows.source_line.ne("").all()
    assert rows.code_excerpt.ne("").all()


def test_required_tables_and_methods_exist():
    for name in (
        "variable_definitions.csv",
        "parameter_definitions.csv",
        "equation_to_source_function.csv",
        "branch_result_to_equation_map.csv",
        "diagnostic_only_vs_final_equations.csv",
        "missing_equation_sources.csv",
    ):
        assert (OUT / name).exists()
    assert len(pd.read_csv(OUT / "variable_definitions.csv")) >= 30
    assert (ROOT / "docs" / "METHODS_TEXT_WITH_EQUATIONS_FOR_PAPER.md").exists()


def test_final_claims_use_exact_promoted_counts_and_tier_label():
    text = (ROOT / "docs" / "EQUATION_FUNCTIONAL_FORM_AUDIT_FOR_PAPER.md").read_text()
    for token in ("1,903", "485", "119", "73", "1,226", "50,655", "19,880"):
        assert token in text
    combined = "\n".join(
        (ROOT / "docs" / name).read_text()
        for name in (
            "EQUATION_FUNCTIONAL_FORM_AUDIT_FOR_PAPER.md",
            "METHODS_TEXT_WITH_EQUATIONS_FOR_PAPER.md",
            "EQUATION_LIMITATIONS_AND_APPROXIMATIONS.md",
        )
    )
    assert "conditional Tier B" in combined
    assert "Tier A" in combined


def test_surrogate_equations_are_not_final_evidence():
    rows = pd.read_csv(OUT / "equation_registry.csv", keep_default_na=False)
    surrogate = rows[rows.equation_id.isin(["PROP-06", "PROP-07"])]
    assert len(surrogate) == 2
    assert set(surrogate.evidence_role) == {"screening_only"}
    assert not surrogate.evidence_role.str.startswith("final").any()


def test_audit_branch_does_not_change_protected_model_sources():
    def names(args):
        result = subprocess.run(
            ["git", "diff", "--name-only", *args],
            cwd=ROOT, text=True, capture_output=True, check=True,
        )
        return set(result.stdout.splitlines())

    changed = names([]) | names(["--cached"])
    check_base = subprocess.run(
        ["git", "cat-file", "-e", "8ff2262^{commit}"],
        cwd=ROOT, capture_output=True,
    )
    if check_base.returncode == 0:
        changed |= names(["8ff2262...HEAD"])
    assert not (changed & PROTECTED_MODEL_SOURCES)
