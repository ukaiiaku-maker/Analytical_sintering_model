import csv
import math
from pathlib import Path

import numpy as np

import density_window_processing_map as processing_map
import topology_constrained_sintering as model


RESULTS = Path(__file__).parents[1] / "results" / "density_window_processing_map"


def test_unattainable_targets_are_reported_without_scores():
    with (RESULTS / "unattainable_cases.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    for row in rows:
        assert row["eligible_target"] == "True"
        assert row["reached_target"] == "False"
        assert math.isnan(float(row["G_at_target_nm"]))
        assert math.isnan(float(row["time_to_target_h"]))


def test_density_resolved_metrics_require_both_paths_to_reach():
    with (RESULTS / "heating_rate_density_curves.csv").open() as stream:
        for row in csv.DictReader(stream):
            if not math.isnan(float(row["HR_pct_vs_0p2"])):
                assert row["eligible_target"] == "True"
                assert row["slow_reached"] == "True"
                assert row["current_reached"] == "True"
    with (RESULTS / "two_step_density_grid.csv").open() as stream:
        for row in csv.DictReader(stream):
            if not math.isnan(float(row["TS_pct"])):
                assert row["eligible_target"] == "True"
                assert row["high_reached"] == "True"
                assert row["two_step_reached"] == "True"


def test_reported_target_attainment_is_not_extrapolated():
    with (RESULTS / "all_protocol_targets.csv").open() as stream:
        for row in csv.DictReader(stream):
            if row["reached_target"] == "True":
                assert float(row["final_density"]) >= float(row["target_density"]) - 1e-12


def test_ensemble_states_preserve_pore_invariants():
    for initial in processing_map.INITIAL_CLASSES:
        params = initial.params()
        result = model.run(params, model.Iso(1300, 3600))
        assert np.all(result["pore_phi"] >= 0)
        assert np.all(result["pore_N"] >= 0)
        assert np.allclose(result["rho"], 1.0 - result["pore_phi"].sum(axis=1), atol=1e-12)
        assert params.memory_model == "pore_bin_redistribution"


def test_failed_sample_is_not_scored():
    result = model.run(model.Params(t_max_s=1), model.Iso(25, 1))
    sampled = processing_map.sample_at_density(result, 0.98, 0.75)
    assert not sampled["reached_target"]
    assert math.isnan(sampled["G_at_target_nm"])
    assert math.isnan(sampled["E_G"])
