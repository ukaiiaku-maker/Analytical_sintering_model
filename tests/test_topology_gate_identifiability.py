import ast
import csv
import inspect
import math
from dataclasses import replace
from pathlib import Path

import numpy as np

import topology_constrained_sintering as model
import topology_gate_identifiability as audit

RESULTS = Path(__file__).parents[1] / "results" / "topology_gate_identifiability"


def test_density_mode_preserves_original_gate_exactly():
    params = audit.base_params()
    for rho in np.linspace(.65, .92, 15):
        state = model.initial_state(replace(params, rho0=float(rho)))
        assert model.smoothing_topology_gate(state, params) == model.smoothing_density_gate(state.rho, params)


def test_topology_gates_are_bounded_and_schedule_free():
    for mode in (*audit.GATE_MODES, "none"):
        params = replace(audit.base_params(), smoothing_gate_mode=mode)
        value = model.smoothing_topology_gate(model.initial_state(params), params)
        assert 0 <= value <= 1
    for function in (model.smoothing_topology_gate, model.surface_smoothing_redistribution):
        source = inspect.getsource(function)
        names = {node.id for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Name)}
        assert not names & {"protocol", "schedule", "ramp_rate", "heating_rate", "slow", "fast"}


def test_gate_mode_table_is_complete_and_never_scores_failure():
    with (RESULTS / "gate_mode_results.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(audit.GATE_MODES) * len(audit.RHO0_VALUES) * len(audit.TARGETS)
    assert {row["smoothing_gate_mode"] for row in rows} == set(audit.GATE_MODES)
    for row in rows:
        if not math.isnan(float(row["HR_pct"])):
            assert row["eligible_target"] == "True"
            assert row["slow_reached"] == row["fast_reached"] == "True"
        if not math.isnan(float(row["TS_pct"])):
            assert row["eligible_target"] == "True"
            assert row["high_reached"] == row["two_step_reached"] == "True"


def test_topology_modes_preserve_exact_pore_conservation():
    for mode in audit.GATE_MODES:
        params = replace(audit.condition_params(audit.base_params(), mode, .75), t_max_s=1800)
        result = model.run(params, model.Iso(1300, 1800))
        assert np.all(result["pore_phi"] >= 0)
        assert np.all(result["pore_N"] >= 0)
        assert np.allclose(result["rho"], 1-result["pore_phi"].sum(axis=1), atol=1e-12)
