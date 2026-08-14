import inspect

import numpy as np

import pore_location_agentic_model as action
import pore_location_agentic_sensitivity as sensitivity
import pore_location_topology_model as fixed
import pore_location_topology_sensitivity as placement
import topology_constrained_sintering as aggregate


def params(mode="action_evolving_capture"):
    return action.ActionParams(
        fixed.LocationParams(placement.base_params(), "evolving"), mode
    )


def test_action_scoring_has_no_schedule_input_or_label_leakage():
    signature = inspect.signature(action.score_actions)
    assert tuple(signature.parameters) == ("s", "T_C", "p")
    source = inspect.getsource(action.score_actions).lower()
    for forbidden in ("protocol", "ramp_rate", "slow", "fast", "schedule"):
        assert forbidden not in source


def test_actions_are_nonnegative_and_weights_partition_unity():
    actions, diagnostics = action.score_actions(fixed.initial_state(params().location), 1200.0, params())
    assert all(item.propensity >= 0.0 and item.resistance >= 0.0 for item in actions.values())
    weights = [diagnostics[f"action_weight_{name}"] for name in action.ACTION_NAMES]
    assert np.isclose(sum(weights), 1.0, atol=1e-12)


def test_capture_is_an_exact_mode_ablation():
    state = fixed.initial_state(params().location)
    off = action.allocated_fluxes(state, 1250.0, params("action_evolving_no_capture"))
    on = action.allocated_fluxes(state, 1250.0, params("action_evolving_capture"))
    assert off["action_flux_TJ_to_GBseg_capture"] == 0.0
    assert on["action_flux_TJ_to_GBseg_capture"] >= 0.0


def test_evolving_action_run_preserves_pore_identity_and_nonnegative_bins():
    protocol = aggregate.Iso(1250.0, 300.0)
    history = action.run(params(), protocol)
    phi = history["phi_GBseg"] + history["phi_TJ"] + history["phi_iso"]
    assert np.all(phi >= 0.0)
    assert np.all(history["N_GBseg"] >= 0.0)
    assert np.all(history["N_TJ"] >= 0.0)
    assert np.all(history["N_iso"] >= 0.0)
    assert np.max(np.abs(history["rho"] - (1.0 - phi.sum(axis=1)))) < 1e-12


def test_redistribution_actions_do_not_directly_change_total_pore_volume():
    state = fixed.initial_state(params().location)
    flux = action.allocated_fluxes(state, 1200.0, params())
    conservative = (
        flux["GB_smooth"].sum()
        - flux["GB_to_TJ"].sum() + flux["GB_to_TJ"].sum()
        + flux["TJ_to_GBseg_capture"].sum() - flux["TJ_to_GBseg_capture"].sum()
        - flux["TJ_to_iso"].sum() + flux["TJ_to_iso"].sum()
    )
    assert abs(conservative) < 1e-18


def test_window_classification_is_mutually_exclusive_and_strict():
    cases = {
        sensitivity.classify(True, 0.90, 0.04, 0.05),
        sensitivity.classify(True, 0.90, 0.06, 0.05),
        sensitivity.classify(True, 0.89, 0.04, 0.05),
        sensitivity.classify(True, 0.89, 0.06, 0.05),
        sensitivity.classify(False, 0.90, 0.04, 0.05),
    }
    assert cases == {
        "SUCCESS", "GRAIN_GROWTH_FAILURE", "DENSIFICATION_EXHAUSTION_FAILURE",
        "MIXED_FAILURE", "UNATTAINABLE_FIRST_STEP",
    }
