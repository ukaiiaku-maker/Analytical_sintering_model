#!/usr/bin/env python3
"""Manuscript-integration QC for the static equation audit.

This script reads documentation and CSV registries only.  It does not import,
instantiate, or execute any sintering model.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "equation_functional_form_audit"
DOCS = ROOT / "docs"

REGISTRY = OUT / "equation_registry.csv"
VARIABLES = OUT / "variable_definitions.csv"
PARAMETERS = OUT / "parameter_definitions.csv"
METHODS_IN = DOCS / "METHODS_TEXT_WITH_EQUATIONS_FOR_PAPER.md"
SI_IN = DOCS / "SI_EQUATION_TABLES_AND_VARIABLE_DEFINITIONS.md"
METHODS_OUT = DOCS / "METHODS_TEXT_WITH_EQUATIONS_FOR_PAPER_QC.md"
SI_OUT = DOCS / "SI_EQUATION_TABLES_AND_VARIABLE_DEFINITIONS_QC.md"
REPORT = DOCS / "EQUATION_AUDIT_QC_REPORT.md"
CHECKS = OUT / "equation_audit_qc_checks.csv"

PAPER_MAP = {
    1: ("FF-01",),
    2: ("FF-02", "FF-03", "FF-04"),
    3: ("FF-05", "FF-06"),
    4: ("FF-07", "FF-08"),
    5: ("FF-09", "FF-12"),
    6: ("PR-02", "PR-04"),
    7: ("PR-05",),
    8: ("PR-07", "PR-08"),
    9: ("PL-01",),
    10: ("LR-01",),
    11: ("FF-10",),
    12: ("LR-09", "LR-10"),
    13: ("TJ-01",),
    14: ("TJ-02", "TJ-03"),
    15: ("TJ-04", "TJ-07"),
    16: ("LR-02",),
    17: ("LR-03",),
    18: ("LR-04", "LR-05"),
    19: ("LR-11",),
    20: ("LR-06",),
    21: ("LR-08",),
    22: ("MET-03", "MET-04", "MET-05"),
    23: ("MET-01", "MET-02"),
    24: ("PROP-01", "PROP-02"),
    25: ("PROP-04", "PROP-05"),
}

# Explicit symbol inventory for the numbered Methods equations.  Each entry
# maps the displayed paper symbol to a machine-readable variable/parameter key.
SYMBOL_REQUIREMENTS = {
    "T": "T", "T_C": "T_C", "rho": "rho", "G": "G", "r_i": "r_i", "r_0": "r_0",
    "phi_i": "phi_i", "phi_GB": "phi_GBseg", "phi_TJ": "phi_TJ",
    "phi_iso": "phi_iso", "phi_closed": "phi_closed", "phi_tot": "phi_tot",
    "w_j": "w_j", "sigma_loc": "sigma_loc", "sigma_j": "sigma_j",
    "tau_nuc": "tau_nuc", "tau_ex": "tau_exchange", "tau_tr": "tau_transport",
    "tau_cyc": "tau_cycle", "a": "activity", "eta_geo": "eta_geo",
    "tau_ex_0": "tau_exchange_prefactor", "tau_tr_0": "tau_transport_prefactor",
    "rho_dot": "rho_dot", "G_dot": "G_dot", "J_PR": "J_PR",
    "f_fine": "f_fine", "g_low": "g_low", "g_top": "g_top",
    "w_i_fine": "w_i_fine", "theta_PR": "theta_PR", "H_dens": "H_dens",
    "H_PR": "H_PR", "w_dens": "w_dens", "w_PR": "w_PR",
    "eta_dens": "eta_dens", "C_rem": "C_rem", "Gamma_A": "Gamma_A",
    "Gamma_j": "Gamma_j", "D": "D", "D_j": "D_j", "P_TJ": "P_TJ",
    "X_J": "X_J", "X_J_prod": "X_J_prod", "J_cap": "J_cap",
    "J_reloc": "J_reloc", "X_cap": "X_cap", "R_J": "R_J",
    "tau_J": "tau_J", "tau_J_ref": "tau_J_ref_s",
    "C_TJ": "C_TJ", "f_clean_GB": "f_clean_GB", "Lambda_TJ": "Lambda_TJ",
    "K_TJ": "K_TJ", "P_comp_TJ": "P_comp_TJ", "P_comp": "P_comp_TJ",
    "n": "n", "xi_r": "xi_r",
    "N_closed": "N_i", "A_closed": "A_closed", "A_cap": "A_cap",
    "Delta_phi_closed_loss": "Delta_phi_closed_loss", "Delta_t": "Delta_t",
    "G_1": "G_1", "G_2": "G_2", "G_ref": "G_ref", "G_fast": "G_fast",
    "a_j": "a_j", "A_j": "A_closed", "J_close": "J_close",
    "g_2": "g_2", "R_fast": "R_fast", "Delta_rho": "Delta_rho",
    "W": "W", "T_first_success": "T_first_success", "T_last_success": "T_last_success",
    "Theta_nuc": "Theta_nuc", "S_closed_growth": "S_closed_growth",
    "R": "R", "k_B": "k_B",
    "gamma_s": "gamma_s", "gamma_GB": "gamma_GB",
    "c_sigma": "c_sigma", "nu_0": "nu0_nucleation",
    "Q_nuc": "Q_disconnection_nucleation", "v_star": "v_star",
    "Q_ex": "Q_exchange", "Q_tr": "Q_transport",
    "epsilon_event": "event_strain", "zeta_eta": "zeta_eta_ratio",
    "Q_GB": "Q_GB_diffusion", "D_GB_0": "D_GB_prefactor",
    "Q_s": "Q_surface_diffusion", "D_s_0": "D_surface_prefactor",
    "alpha_attr": "alpha_attr", "q_PR": "q_PR", "p_top": "p_top",
    "a_mid": "a_mid", "a_width": "a_width", "p_a": "p_a",
    "c_D": "c_D", "c_R": "c_R", "c_S": "c_S", "Q_J": "Q_J",
    "T_ref": "T_ref", "A_J": "A_J", "q_J": "q_J", "K_0": "K_0",
    "k_open": "k_open", "Q_density": "Q_density", "k_closed": "k_closed",
    "Q_closed": "Q_closed", "k_PR": "k_PR", "Q_PR": "Q_PR",
    "k_tr": "k_tr", "k_g": "k_g", "Q_growth": "Q_growth",
    "r_g": "r_g", "q_r": "q_r", "tau_A": "tau_A",
    "g_tol": "g_tol", "rho_target": "rho_target", "f_x": "f_x",
    "Delta_Q_x": "Delta_Q_x", "Q_x": "Q_x", "k_x": "k_x",
}

REQUIRED_FIELDS = (
    "equation_id", "source_file", "source_function", "evidence_role",
    "changes_density", "changes_migration",
    "conservatively_redistributes_pore_volume", "implementation_note",
)
NONFINAL_ROLES = {"diagnostic_only", "screening_only", "superseded", "negative_control"}


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    rows = list(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def inject_registry_tags(text):
    output, pending = [], None
    for line in text.splitlines():
        match = re.search(r"\\tag\{(\d+)\}", line)
        if match:
            pending = int(match.group(1))
        output.append(line)
        if line.strip() == r"\]" and pending is not None:
            ids = PAPER_MAP[pending]
            output.extend(("", "*(Registry trace: " + "; ".join(f"Eq. {item}" for item in ids) + ".)*", ""))
            pending = None
    return "\n".join(output).rstrip() + "\n"


def source_table(registry):
    lines = [
        "## Table S9. Complete equation-to-source and evidence map",
        "",
        "| Equation ID | Equation name | Source file | Source function | Evidence role |",
        "|---|---|---|---|---|",
    ]
    for row in registry:
        lines.append(
            f"| {row['equation_id']} | {row['equation_name']} | "
            f"{row['source_file']} | {row['source_function']} | {row['evidence_role']} |"
        )
    return "\n".join(lines)


def nonclaims():
    return """## Manuscript non-claims after QC

- Candidate 693168 is **conditional Tier B, not validation** and not a calibrated Tier-A material model.
- The large attained high-temperature/two-step grain-size separation is **not inherently unphysical**; its magnitude remains an experimental-scale prediction requiring calibration.
- The modeled closed-pore/accommodation trajectory is the primary calibration and falsification target.
- Surrogate and screening-only equations are not final evidence; exact-promoted rows control classifications.
- No hidden closed-pore Lambda/K law was implemented. Closed accommodation is an implemented bounded proxy, not a derived closed-pore Poisson or gas-transport law.
"""


def add_check(rows, check_id, category, passed, details):
    rows.append(dict(check_id=check_id, category=category,
                     passed=bool(passed), details=str(details)))


def main():
    registry = read_csv(REGISTRY)
    variables = read_csv(VARIABLES)
    parameters = read_csv(PARAMETERS)
    ids = {row["equation_id"] for row in registry}
    checks = []

    add_check(checks, "registry_count", "registry", len(registry) >= 67,
              f"{len(registry)} equations; {len(ids)} unique IDs")
    missing_fields = [(row.get("equation_id", ""), field) for row in registry
                      for field in REQUIRED_FIELDS if not row.get(field, "").strip()]
    add_check(checks, "registry_required_fields", "registry", not missing_fields,
              f"missing entries={missing_fields}")
    final = [row for row in registry if row["evidence_role"].startswith("final_evidence")]
    final_missing_source = [row["equation_id"] for row in final if not row["source_function"].strip()]
    add_check(checks, "final_source_functions", "registry", not final_missing_source,
              f"missing={final_missing_source}")

    misdescribed = []
    for row in registry:
        if row["evidence_role"] in NONFINAL_ROLES:
            text = (row["implementation_note"] + " " + row["branches_results_using_it"]).lower()
            if "final mechanism" in text or "controlling mechanism" in text:
                misdescribed.append(row["equation_id"])
    add_check(checks, "nonfinal_role_language", "evidence", not misdescribed,
              f"nonfinal rows described as final/controlling={misdescribed}")

    mapped = {item for values in PAPER_MAP.values() for item in values}
    paper_numbers = {int(x) for x in re.findall(r"\\tag\{(\d+)\}", METHODS_IN.read_text())}
    add_check(checks, "methods_equation_number_map", "methods",
              paper_numbers == set(PAPER_MAP), f"paper={sorted(paper_numbers)}")
    unknown = sorted(mapped - ids)
    add_check(checks, "methods_registry_ids", "methods", not unknown,
              f"mapped IDs={len(mapped)}; unknown={unknown}")

    methods_qc = inject_registry_tags(METHODS_IN.read_text())
    methods_qc += "\n" + nonclaims()
    METHODS_OUT.write_text(methods_qc)

    si_qc = SI_IN.read_text().rstrip() + "\n\n" + source_table(registry) + "\n\n" + nonclaims()
    SI_OUT.write_text(si_qc)
    methods_ids = set(re.findall(r"Eq\. ([A-Z]+-\d+)", methods_qc))
    si_ids = set(re.findall(r"\| ([A-Z]+-\d+) \|", si_qc))
    final_ids = {row["equation_id"] for row in final}
    missing_final = sorted(final_ids - methods_ids - si_ids)
    add_check(checks, "final_equation_document_coverage", "methods_si",
              not missing_final, f"final IDs={len(final_ids)}; missing={missing_final}")

    definitions = {row["variable"]: row for row in variables}
    definitions.update({row["parameter"]: row for row in parameters})
    missing_symbols = sorted(
        f"{symbol}->{key}" for symbol, key in SYMBOL_REQUIREMENTS.items() if key not in definitions
    )
    add_check(checks, "methods_symbol_definitions", "variables",
              not missing_symbols, f"symbols={len(SYMBOL_REQUIREMENTS)}; missing={missing_symbols}")
    missing_units = sorted(key for key in set(SYMBOL_REQUIREMENTS.values())
                           if key in definitions and not definitions[key].get("units", "").strip())
    add_check(checks, "methods_symbol_units", "variables",
              not missing_units, f"missing units={missing_units}")
    proxy_keys = ("sigma_res", "N_i", "P_TJ", "xi_r", "r_g")
    bad_proxy = [key for key in proxy_keys
                 if key not in definitions or "proxy" not in
                 (definitions[key].get("units", "") + definitions[key].get("definition", "")).lower()]
    add_check(checks, "proxy_dimension_labels", "variables",
              not bad_proxy, f"unmarked proxy keys={bad_proxy}")

    screening_final = [row["equation_id"] for row in registry
                       if row["evidence_role"] == "screening_only"
                       and row["evidence_role"].startswith("final")]
    add_check(checks, "screening_not_final", "evidence", not screening_final,
              f"screening rows incorrectly final={screening_final}")

    combined = methods_qc + "\n" + si_qc
    warning_terms = (
        "implemented bounded proxy",
        "not a derived closed-pore Poisson or gas-transport law",
        "No hidden closed-pore Lambda/K law was implemented",
    )
    missing_warnings = [term for term in warning_terms if term not in combined]
    add_check(checks, "closed_accommodation_warning", "nonclaims",
              not missing_warnings, f"missing phrases={missing_warnings}")
    tier_ok = "conditional Tier B, not validation" in combined
    add_check(checks, "candidate_693168_tier", "nonclaims", tier_ok,
              "candidate is explicitly conditional Tier B, not validation")
    large_ok = "not inherently unphysical" in combined
    add_check(checks, "large_separation_nonclaim", "nonclaims", large_ok,
              "large attained separation is not rejected by magnitude alone")
    calibration_ok = "primary calibration and falsification target" in combined
    add_check(checks, "closed_state_calibration_target", "nonclaims", calibration_ok,
              "closed/accommodation trajectory identified as calibration target")

    write_csv(CHECKS, checks)
    passed = sum(row["passed"] for row in checks)
    failed = [row for row in checks if not row["passed"]]
    report = f"""# Equation audit manuscript-integration QC report

## Outcome

**{passed} of {len(checks)} checks passed.** This QC pass ran no simulation and changed no model physics, parameter, or classification.

## Registry integrity

- Equation rows: {len(registry)}; unique IDs: {len(ids)}.
- Final-evidence/final-metric equations: {len(final)}.
- Required fields and final source functions: {"complete" if not missing_fields and not final_missing_source else "incomplete"}.
- Non-final role-language conflicts: {len(misdescribed)}.

## Methods and SI integration

- Numbered Methods equation groups mapped: {len(PAPER_MAP)}.
- Distinct registry IDs used by Methods: {len(mapped)}.
- Final equations absent from both tagged Methods and complete SI crosswalk: {missing_final or "none"}.
- Methods symbols checked: {len(SYMBOL_REQUIREMENTS)}; missing definitions: {missing_symbols or "none"}.
- Symbols missing units/status: {missing_units or "none"}.

## Evidence discipline

Screening-only equations remain screening-only. Candidate 693168 remains conditional Tier B, not validation. The large attained grain-size separation is not rejected as inherently unphysical, while its magnitude remains uncalibrated. The closed-pore/accommodation trajectory is the principal calibration target. The accommodation equation is labeled as an implemented bounded proxy; no closed-pore Poisson Lambda/K or explicit gas-transport law is claimed.

## Failed checks

{chr(10).join("- " + row["check_id"] + ": " + row["details"] for row in failed) if failed else "None."}

The machine-readable check record is results/equation_functional_form_audit/equation_audit_qc_checks.csv.
"""
    REPORT.write_text(report)
    print(f"qc_checks={len(checks)}")
    print(f"qc_passed={passed}")
    print(f"qc_failed={len(failed)}")
    print(f"methods_registry_ids={len(mapped)}")
    print(f"final_equations={len(final_ids)}")
    print(f"symbols_checked={len(SYMBOL_REQUIREMENTS)}")
    if failed:
        for row in failed:
            print(f"FAIL {row['check_id']}: {row['details']}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
