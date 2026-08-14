#!/usr/bin/env python3
"""Audit the manuscript figure-source-data package without running the model."""

from __future__ import annotations

import json
import subprocess
import zipfile
from datetime import date
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "figure_source_data_package"
QC = OUT / "qc"
QC.mkdir(parents=True, exist_ok=True)
CHECKS: list[dict[str, object]] = []


def check(check_id: str, category: str, passed: bool, details: str) -> None:
    CHECKS.append(
        {"check_id": check_id, "category": category, "passed": bool(passed), "details": details}
    )


def read(relative: str) -> pd.DataFrame:
    return pd.read_csv(OUT / relative, low_memory=False)


def main() -> None:
    citation = read("figure_citation_map.csv")
    manifest = read("figure_data_manifest.csv")
    dictionary = read("all_columns_dictionary.csv")

    citation_files = [ROOT / value for value in citation.source_data_file]
    check(
        "citation_sources_exist",
        "coverage",
        all(path.exists() and path.stat().st_size > 0 for path in citation_files),
        f"{len(citation_files)} figure references checked",
    )
    check(
        "manifest_files_exist",
        "coverage",
        all((ROOT / value).exists() and (ROOT / value).stat().st_size > 0 for value in manifest.relative_path),
        f"{len(manifest)} exported scientific datasets checked",
    )

    required_meta = {
        "figure_id",
        "panel_id",
        "candidate_id",
        "evidence_level",
        "model_layer",
        "source_file",
        "source_table",
        "source_branch_or_commit",
        "units_reference",
        "human_readable_label_reference",
        "derivation_notes",
    }
    missing_meta = []
    empty_tables = []
    for row in manifest.itertuples():
        frame = pd.read_csv(ROOT / row.relative_path, nrows=5, low_memory=False)
        if frame.empty:
            empty_tables.append(row.relative_path)
        missing = required_meta - set(frame.columns)
        if missing:
            missing_meta.append(f"{row.relative_path}:{sorted(missing)}")
    check("scientific_tables_nonempty", "schema", not empty_tables, f"empty={empty_tables}")
    check("scientific_table_metadata", "schema", not missing_meta, f"missing={missing_meta[:5]}")
    allowed_evidence = {"exact", "surrogate_screen", "diagnostic", "negative_control", "schematic"}
    levels = set(manifest.evidence_level.dropna().astype(str))
    check("evidence_levels_explicit", "evidence", levels <= allowed_evidence, f"levels={sorted(levels)}")

    exact_map = read("figure_02_exact_property_phase_map/exact_property_phase_map.csv")
    counts = read("figure_02_exact_property_phase_map/exact_behavior_classification_counts.csv")
    expected_counts = {
        "screened_rows": 50655,
        "unique_exact_cases": 1903,
        "fast_only": 485,
        "two_step_only": 119,
        "both_pass": 73,
        "neither": 1226,
        "surrogate_both": 19880,
    }
    observed = dict(zip(counts.behavior_class_or_total, counts["count"]))
    check("exact_counts_reproduced", "evidence", all(int(observed.get(k, -1)) == v for k, v in expected_counts.items()), f"observed={observed}")
    check("exact_union_rows", "evidence", len(exact_map) == 1903, f"rows={len(exact_map)}")
    check("surrogate_not_final", "evidence", not ((manifest.evidence_level.eq("surrogate_screen")) & manifest.notes.astype(str).str.contains("final evidence", case=False, na=False)).any(), "screening tables are explicitly non-final")

    fast_curves = read("figure_03_fast_firing_property_window/fast_firing_matched_density_curves.csv")
    fast_ablation = read("fast_firing/clean_fast_firing_T_rho_G_time.csv")
    required_ablations = {"full", "PR-off", "nucleation-facile"}
    check("fast_firing_ablations_present", "fast_firing", required_ablations <= set(fast_ablation.ablation), f"ablations={sorted(set(fast_ablation.ablation))}")
    check("fast_firing_materials_present", "fast_firing", {"E0021", "E0142"} <= set(fast_ablation.material_id), f"materials={sorted(set(fast_ablation.material_id))}")
    check("fast_matched_density_attained", "fast_firing", fast_curves.both_paths_attained.astype(bool).all(), f"rho_range={fast_curves.rho.min():.3f}-{fast_curves.rho.max():.3f}")

    two = read("figure_04_two_step_property_window/two_step_matched_density_curves_693168.csv")
    supported = two[two.interpolation_support_flag.astype(bool)]
    success = supported[supported.path_pair_id.str.contains("two_step_success")]
    covers = not success.empty and success.rho.min() <= 0.95 and success.rho.max() >= 0.98
    check("candidate_693168_density_interval", "candidate_693168", covers, f"success_rho_range={success.rho.min() if not success.empty else None}-{success.rho.max() if not success.empty else None}")

    candidate_hist = read("candidate_693168/dense_time_histories.csv")
    time_ok = True
    discontinuities = []
    for path, group in candidate_hist.groupby("path_label"):
        times = group.physical_time_s.dropna()
        if not times.is_monotonic_increasing:
            time_ok = False
            discontinuities.append(path)
        stages = group.stage.dropna().astype(str)
        if "step2" in set(stages) and group.loc[stages.eq("step2"), "physical_time_s"].min() < group.loc[~stages.eq("step2"), "physical_time_s"].max():
            time_ok = False
            discontinuities.append(f"{path}:stage reset")
    check("continuous_physical_time", "candidate_693168", time_ok, f"issues={discontinuities}")

    fine = read("figure_05_chen_window_boundaries/candidate_693168_fine_T2_classification.csv")
    required_classes = {"DENSIFICATION_EXHAUSTION_FAILURE", "SUCCESS", "GRAIN_GROWTH_FAILURE"}
    check("fine_T2_three_regions", "Chen_map", required_classes <= set(fine.classification), f"classes={sorted(set(fine.classification))}")
    check("fine_T2_not_success_only", "Chen_map", (~fine.success_band_flag.astype(bool)).any(), f"points={len(fine)} failures={(~fine.success_band_flag.astype(bool)).sum()}")

    chen_files = [
        "chen_maps/chen_map_T1_T2_classification.csv",
        "chen_maps/chen_map_G1_T2_classification.csv",
        "chen_maps/chen_map_switch_density_T2_classification.csv",
    ]
    chen_complete = True
    chen_details = []
    for filename in chen_files:
        frame = read(filename)
        classes = set(frame.classification)
        ok = required_classes <= classes and (~frame.success_flag.astype(bool)).any()
        chen_complete &= ok
        chen_details.append(f"{filename}:{len(frame)}:{sorted(classes)}")
    check("complete_Chen_maps", "Chen_map", chen_complete, "; ".join(chen_details))

    family = read("figure_06_six_tierB_family/six_TierB_candidate_summary.csv")
    expected_ids = {693168, 822940, 581668, 295003, 366094, 85161}
    observed_ids = set(pd.to_numeric(family.candidate_id, errors="coerce").dropna().astype(int))
    check("six_TierB_candidates", "family", observed_ids == expected_ids, f"ids={sorted(observed_ids)}")

    surrogate = read("figure_07_surrogate_vs_exact/surrogate_vs_exact_counts.csv")
    values = dict(zip(surrogate.metric, surrogate.value))
    check("surrogate_exact_warning_counts", "evidence", int(values.get("surrogate_both", -1)) == 19880 and int(values.get("exact_both", -1)) == 73, f"values={values}")

    dict_columns = set(dictionary.column)
    science_columns = set()
    for row in manifest.itertuples():
        science_columns.update(pd.read_csv(ROOT / row.relative_path, nrows=0).columns)
    missing_defs = sorted(science_columns - dict_columns)
    unit_gaps = dictionary.units.isna() | dictionary.units.astype(str).str.strip().eq("")
    label_gaps = dictionary.human_readable_label.isna() | dictionary.human_readable_label.astype(str).str.strip().eq("")
    check("all_columns_defined", "dictionary", not missing_defs, f"missing={missing_defs[:20]}")
    check("all_units_present", "dictionary", not unit_gaps.any(), f"gaps={int(unit_gaps.sum())}")
    check("all_labels_present", "dictionary", not label_gaps.any(), f"gaps={int(label_gaps.sum())}")

    missing_channels = read("candidate_693168/missing_candidate_693168_channels.csv")
    check("missing_channels_honest", "guardrails", missing_channels.handling.str.contains("no value invented", case=False).all(), f"channels={len(missing_channels)}")
    readme = (ROOT / "docs" / "FIGURE_SOURCE_DATA_PACKAGE_README.md").read_text()
    guardrail_text = " ".join([
        readme,
        (ROOT / "docs" / "FIGURE_SOURCE_DATA_MISSING_ITEMS.md").read_text(),
    ]).lower()
    check("conditional_TierB_nonclaim", "guardrails", "conditional tier b" in guardrail_text and "not validation" in guardrail_text, "candidate status stated")
    check(
        "closed_proxy_warning",
        "guardrails",
        "bounded proxy" in guardrail_text
        and "hidden closed-pore lambda/k law" in guardrail_text
        and "not a derived poisson or gas-transport law" in guardrail_text,
        "proxy, non-derived-law, and no-hidden-law warnings present",
    )

    porcelain = subprocess.run(["git", "status", "--porcelain=v1"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.splitlines()
    model_changes = []
    allowed_new = {
        "build_figure_source_data_package.py",
        "audit_figure_source_data_package.py",
        "tests/test_figure_source_data_package.py",
        "conftest.py",
        "tests/test_final_mechanism_synthesis.py",
    }
    for line in porcelain:
        path = line[3:]
        if path.startswith("results/") or path.startswith("docs/FIGURE_SOURCE_DATA") or path in allowed_new:
            continue
        if path.endswith(".py"):
            model_changes.append(path)
    check("no_model_physics_files_changed", "git_hygiene", not model_changes, f"unexpected_python_changes={model_changes}")
    staged = subprocess.run(["git", "diff", "--cached", "--name-status"], cwd=ROOT, check=True, capture_output=True, text=True).stdout.splitlines()
    staged_deletions = [line for line in staged if line.startswith("D\tresults/")]
    check("no_archived_result_deletions_staged", "git_hygiene", not staged_deletions, f"staged_result_deletions={len(staged_deletions)}")

    archive = ROOT / "results" / f"figure_source_data_package_{date.today():%Y%m%d}.zip"
    check("package_archive_exists", "archive", archive.exists() and archive.stat().st_size > 0, f"path={archive.relative_to(ROOT)} bytes={archive.stat().st_size if archive.exists() else 0}")

    checks = pd.DataFrame(CHECKS)
    checks.to_csv(QC / "figure_source_data_audit.csv", index=False)
    summary = {
        "passed": bool(checks.passed.all()),
        "checks_total": len(checks),
        "checks_passed": int(checks.passed.sum()),
        "checks_failed": int((~checks.passed).sum()),
        "figure_references": len(citation),
        "figure_references_complete": int(sum(path.exists() and path.stat().st_size > 0 for path in citation_files)),
        "scientific_tables": len(manifest),
        "exact_counts": expected_counts,
        "candidate_693168": {
            "dense_histories": (OUT / "candidate_693168/dense_time_histories.csv").exists(),
            "fine_Chen_map": (OUT / "candidate_693168/fine_T2_scan.csv").exists(),
            "matched_density_curves": (OUT / "candidate_693168/matched_density_ratio_curves.csv").exists(),
            "ablations": (OUT / "candidate_693168/ablation_summary.csv").exists(),
        },
        "fast_firing": {"full": "full" in set(fast_ablation.ablation), "PR_off": "PR-off" in set(fast_ablation.ablation), "nucleation_facile": "nucleation-facile" in set(fast_ablation.ablation)},
        "Chen_maps": {"T1_T2": (OUT / chen_files[0]).exists(), "G1_T2": (OUT / chen_files[1]).exists(), "switch_density_T2": (OUT / chen_files[2]).exists()},
        "archive": str(archive.relative_to(ROOT)),
        "archive_bytes": archive.stat().st_size if archive.exists() else 0,
    }
    (QC / "figure_source_data_package_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    # Rebuild the deliverable after QC so the audit CSV/JSON are included.  This
    # remains a pure packaging operation over already-exported files.
    archive.unlink(missing_ok=True)
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for path in sorted(OUT.rglob("*")):
            if path.is_file():
                bundle.write(path, path.relative_to(ROOT / "results"))
        for name in [
            "FIGURE_SOURCE_DATA_PACKAGE_README.md",
            "FIGURE_SOURCE_DATA_DICTIONARY.md",
            "FIGURE_SOURCE_DATA_PROVENANCE.md",
            "FIGURE_SOURCE_DATA_MISSING_ITEMS.md",
        ]:
            path = ROOT / "docs" / name
            bundle.write(path, Path("figure_source_data_package") / "docs" / name)
    print(f"checks={len(checks)}")
    print(f"passed={int(checks.passed.sum())}")
    print(f"failed={int((~checks.passed).sum())}")
    print(f"figure_references={len(citation)}")
    print(f"scientific_tables={len(manifest)}")
    print(f"archive_bytes={archive.stat().st_size}")
    if not checks.passed.all():
        print(checks.loc[~checks.passed].to_string(index=False))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
