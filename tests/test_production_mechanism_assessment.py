import csv
from pathlib import Path

import numpy as np

import agentic_mechanism_model as model
import production_mechanism_assessment as production
import topology_constrained_sintering as aggregate


ROOT = Path(__file__).parents[1] / "results" / "production_mechanism_assessment"


def rows(name):
    with (ROOT / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_frozen_candidate_set_and_common_numerical_budget():
    candidates = production.frozen_mechanisms()
    assert set(candidates) == {"mech_009", "mech_019", "mech_009_q0", "mech_019_q0"}
    assert production.BUDGET == 96 * 3600
    assert all(p.action.location.base.dt_max_s == production.NUMERICAL_DT_MAX_S for p in candidates.values())


def test_tier_rules_are_exclusive_and_require_complete_practical_window():
    base = dict(map_type="practical", boundary_status="COMPLETE_WINDOW", window_width_C=25,
                growth_tolerance=.05, first_step_growth_fraction=.05)
    assert production.tier(base) == "Tier_A"
    assert production.tier({**base, "growth_tolerance": .10, "first_step_growth_fraction": .10}) == "Tier_B"
    assert production.tier({**base, "boundary_status": "UPPER_BOUND_RIGHT_CENSORED"}) == "Tier_C"
    assert production.tier({**base, "window_width_C": 24.9}) == "Tier_C"


def test_persisted_successes_are_practical_and_fully_bracketed():
    successful = rows("successful_practical_windows.csv")
    assert successful
    for row in successful:
        assert row["tier"] in {"Tier_A", "Tier_B"}
        assert row["boundary_status"] == "COMPLETE_WINDOW"
        assert row["lower_bracketed"] == row["upper_bracketed"] == "True"
        assert float(row["T_last_success_C"]) < float(row["T1_C"])


def test_failed_fast_targets_are_not_scored_and_negative_result_is_preserved():
    summary = rows("production_fast_firing_summary.csv")
    assert summary and not any(r["response_class"] == "beneficial" for r in summary)
    for row in summary:
        if row["response_class"] == "unattainable":
            assert not row["HR_pct_median"]
    score = rows("joint_mechanism_scorecard.csv")
    assert all(r["positive_fast_firing_count"] == "0" and r["joint_positive"] == "False" for r in score)


def test_fast_smoke_run_preserves_exact_pore_conservation_and_nonnegative_bins():
    p0 = production.frozen_mechanisms()["mech_009"]
    p = production.fast_params(p0, 75, .70, production.TOPOLOGIES["baseline"])
    h = model.run(p, production.FastSchedule(20, 1350, 0))
    phi = h["phi_GBseg"] + h["phi_TJ"] + h["phi_iso"]
    assert np.all(phi >= 0)
    assert np.all(h["N_GBseg"] >= 0) and np.all(h["N_TJ"] >= 0) and np.all(h["N_iso"] >= 0)
    assert np.max(np.abs(h["rho"] - (1 - phi.sum(axis=1)))) < 1e-12


def test_production_mechanisms_keep_migration_separate_from_densification():
    candidates = production.frozen_mechanisms()
    state = model.initial_state(candidates["mech_009"])
    rates = [model.local_mechanism(state, 1200, p)["rho_dot"] for p in candidates.values()]
    assert np.allclose(rates, rates[0], rtol=0, atol=0)
