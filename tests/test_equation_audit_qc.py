"""Tests for manuscript integration of the equation registry."""
from pathlib import Path
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "equation_functional_form_audit"
METHODS = ROOT / "docs" / "METHODS_TEXT_WITH_EQUATIONS_FOR_PAPER_QC.md"
SI = ROOT / "docs" / "SI_EQUATION_TABLES_AND_VARIABLE_DEFINITIONS_QC.md"
REPORT = ROOT / "docs" / "EQUATION_AUDIT_QC_REPORT.md"


def registry():
    return pd.read_csv(OUT / "equation_registry.csv", keep_default_na=False)


def tagged_ids(text):
    return set(re.findall(r"Eq\. ([A-Z]+-\d+)", text))


def test_registry_has_at_least_67_complete_rows():
    rows = registry()
    assert len(rows) >= 67
    assert rows.equation_id.is_unique
    required = (
        "equation_id", "source_file", "source_function", "evidence_role",
        "changes_density", "changes_migration",
        "conservatively_redistributes_pore_volume", "implementation_note",
    )
    for field in required:
        assert field in rows
        assert rows[field].astype(str).str.strip().ne("").all()


def test_every_tagged_methods_equation_is_registered():
    ids = set(registry().equation_id)
    methods = METHODS.read_text()
    used = tagged_ids(methods)
    assert len(re.findall(r"\\tag\{\d+\}", methods)) == 25
    assert used
    assert used <= ids


def test_every_final_equation_appears_in_methods_or_si():
    rows = registry()
    final_ids = set(rows[rows.evidence_role.str.startswith("final_evidence")].equation_id)
    methods_ids = tagged_ids(METHODS.read_text())
    si_ids = set(re.findall(r"\| ([A-Z]+-\d+) \|", SI.read_text()))
    assert final_ids <= methods_ids | si_ids


def test_screening_only_is_never_final_evidence():
    rows = registry()
    screening = rows[rows.evidence_role.eq("screening_only")]
    assert len(screening) > 0
    assert not screening.evidence_role.str.startswith("final").any()


def test_closed_proxy_and_candidate_nonclaims_are_explicit():
    text = METHODS.read_text() + "\n" + SI.read_text() + "\n" + REPORT.read_text()
    assert "implemented bounded proxy" in text
    assert "not a derived closed-pore Poisson or gas-transport law" in text
    assert "No hidden closed-pore Lambda/K law was implemented" in text
    assert "conditional Tier B, not validation" in text
    assert "not inherently unphysical" in text
    assert "primary calibration and falsification target" in text


def test_machine_qc_checks_all_pass():
    checks = pd.read_csv(OUT / "equation_audit_qc_checks.csv")
    assert len(checks) >= 10
    assert checks.passed.astype(bool).all()
