import csv
import math
from dataclasses import fields
from pathlib import Path

import numpy as np

import initial_condition_factorial_map as audit
import topology_constrained_sintering as model


RESULTS = Path(__file__).parents[1] / "results" / "initial_condition_sensitivity"


def test_design_is_bounded_and_oat_changes_one_descriptor():
    oat, factorial = audit.oat_design(), audit.factorial_design()
    assert len(oat) == 11
    assert len(factorial) == 16
    baseline = audit.BASELINE
    descriptors = ("rho0", "G0_nm", "pore_scale_nm", "log_width")
    for point in oat[1:]:
        changed = [name for name in descriptors if getattr(point, name) != getattr(baseline, name)]
        assert changed == [point.varied_descriptor]


def test_all_noninitial_mechanism_parameters_are_shared():
    base = model.Params(memory_model="pore_bin_redistribution")
    audit.assert_shared_mechanism_parameters(audit.all_design_points(), base)
    allowed = audit.INITIAL_FIELDS
    reference = {field.name: getattr(base, field.name) for field in fields(base) if field.name not in allowed}
    for point in audit.all_design_points():
        params = point.params(base)
        assert {field.name: getattr(params, field.name) for field in fields(params) if field.name not in allowed} == reference


def test_scores_require_both_paths_and_fixed_budgets():
    rows = []
    for filename in ("one_at_a_time_results.csv", "factorial_results.csv"):
        with (RESULTS / filename).open() as stream:
            rows.extend(csv.DictReader(stream))
    budgets = {(row["slow_time_budget_h"], row["fast_time_budget_h"], row["high_time_budget_h"], row["two_step_time_budget_h"]) for row in rows}
    assert len(budgets) == 1
    for row in rows:
        if not math.isnan(float(row["HR_pct"])):
            assert row["eligible_target"] == "True" and row["slow_reached"] == "True" and row["fast_reached"] == "True"
        if not math.isnan(float(row["TS_pct"])):
            assert row["eligible_target"] == "True" and row["high_reached"] == "True" and row["two_step_reached"] == "True"


def test_unattainable_rows_are_not_scored_for_both_comparisons():
    with (RESULTS / "unattainable_cases.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    for row in rows:
        assert row["eligible_target"] == "True"
        assert row["HR_scored"] == "False" or row["TS_scored"] == "False"
        if row["HR_scored"] == "False":
            assert math.isnan(float(row["HR_pct"]))
        if row["TS_scored"] == "False":
            assert math.isnan(float(row["TS_pct"]))


def test_factorial_corner_pore_states_remain_physical():
    base = model.Params(memory_model="pore_bin_redistribution", t_max_s=1800)
    for point in (audit.factorial_design()[0], audit.factorial_design()[-1]):
        result = model.run(point.params(base), model.Iso(1300, 1800))
        assert np.all(result["pore_phi"] >= 0)
        assert np.all(result["pore_N"] >= 0)
        assert np.allclose(result["rho"], 1.0 - result["pore_phi"].sum(axis=1), atol=1e-12)
