import ast
import csv
import inspect
import math
from pathlib import Path

import numpy as np
import pytest

import smoothing_gate_identifiability as audit
import topology_constrained_sintering as model


RESULTS = Path(__file__).parents[1] / "results" / "smoothing_gate_identifiability"


def test_gate_forms_are_bounded_monotone_and_centered():
    for form in ("logistic", "linear_clipped"):
        params = model.Params(smoothing_gate_form=form, smoothing_rho_mid=.79, smoothing_rho_width=.015)
        values = [model.smoothing_density_gate(rho, params) for rho in (.70, .79, .88)]
        assert 0 <= values[2] <= values[1] <= values[0] <= 1
        assert np.isclose(values[1], .5)
    with pytest.raises(ValueError):
        model.smoothing_density_gate(.79, model.Params(smoothing_gate_form="protocol_specific"))


def test_gate_and_redistribution_laws_have_no_schedule_inputs():
    for function in (model.smoothing_density_gate, model.surface_smoothing_redistribution):
        source = inspect.getsource(function)
        names = {node.id for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Name)}
        assert not names & {"protocol", "schedule", "ramp_rate", "heating_rate", "slow", "fast"}


def test_gate_grid_is_complete_and_scores_only_attainable_pairs():
    with (RESULTS / "gate_sensitivity_results.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 5 * 3 * 5 * 4
    assert {float(row["smoothing_rho_mid"]) for row in rows} == set(audit.RHO_MIDS)
    assert {row["rho_width_label"] for row in rows} == set(audit.WIDTHS)
    assert {float(row["rho0"]) for row in rows} == set(audit.RHO0_VALUES)
    budgets = {(row["slow_time_budget_h"], row["fast_time_budget_h"], row["high_time_budget_h"], row["two_step_time_budget_h"]) for row in rows}
    assert len(budgets) == 1
    for row in rows:
        if not math.isnan(float(row["HR_pct"])):
            assert row["eligible_target"] == "True" and row["slow_reached"] == "True" and row["fast_reached"] == "True"
        if not math.isnan(float(row["TS_pct"])):
            assert row["eligible_target"] == "True" and row["high_reached"] == "True" and row["two_step_reached"] == "True"


def test_observed_crossover_tracks_initial_density_relative_to_gate():
    with (RESULTS / "gate_sensitivity_results.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    one_per_condition = {}
    for row in rows:
        one_per_condition.setdefault(row["condition_id"], row)
    for row in one_per_condition.values():
        has_crossover = not math.isnan(float(row["HR_crossover_density"]))
        assert has_crossover == (float(row["rho0_minus_rho_mid"]) < 0)


def test_gate_variants_preserve_pore_conservation():
    base = audit.base_params()
    parameter_sets = [
        audit.condition_params(base, .65, .72, audit.WIDTHS["narrow"]),
        audit.condition_params(base, .85, .86, audit.WIDTHS["broad"]),
        audit.condition_params(base, .75, .79, audit.WIDTHS["baseline"], "linear_clipped"),
    ]
    audit.assert_only_gate_and_rho_vary(parameter_sets, base)
    for params in parameter_sets:
        result = model.run(replace_time(params), model.Iso(1300, 1800))
        assert np.all(result["pore_phi"] >= 0)
        assert np.all(result["pore_N"] >= 0)
        assert np.allclose(result["rho"], 1 - result["pore_phi"].sum(axis=1), atol=1e-12)


def replace_time(params):
    from dataclasses import replace
    return replace(params, t_max_s=1800)
