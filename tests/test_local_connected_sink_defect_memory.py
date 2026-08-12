import inspect
import numpy as np
import pandas as pd

import local_connected_sink_defect_model as local
import local_connected_sink_defect_search as screen
import persistent_defect_topology_stress_model as persistent
import residual_stress_memory_model as residual
import pr_desintering_memory_model as memory


def base(mode="disabled_local_mixture"):
    row={**screen.design()[0],"mode":mode};return screen.params(row)


def test_disabled_mode_recovers_nonpersistent_local_baseline():
    p=base();protocol=screen.parent.prior.FastSchedule(20,1350,0);_,states=local.run(p,protocol)
    for name,spec in local.normalized_specs(p):
        expected=memory.run(local.cohort_param(p,spec),protocol)
        assert np.array_equal(states[name]["rho"],expected["rho"])


def test_weights_normalized_and_global_density_is_weighted_local_density():
    p=base("static_defect_mixture");specs=local.normalized_specs(p);w=np.array([s.weight for _,s in specs]);assert np.all(w>=0) and np.isclose(w.sum(),1)
    agg,states=local.run(p,screen.parent.prior.FastSchedule(20,1350,0));assert np.isclose(agg["rho"][0],local.global_density_identity(p,states,0))


def test_isolated_pores_not_removed_and_stores_remain_nonnegative():
    p=base("static_defect_mixture");_,states=local.run(p,screen.parent.prior.FastSchedule(20,1350,0))
    for h in states.values():
        assert np.sum(h["phi_iso"][-1])>=np.sum(h["phi_iso"][0])-1e-12
        assert all(np.min(h[k])>=0 for k in ("phi_GBseg","phi_TJ","phi_iso"))
        assert np.allclose(h["rho"],1-np.sum(h["phi_GBseg"]+h["phi_TJ"]+h["phi_iso"],axis=1))


def test_stress_and_defect_states_do_not_directly_contain_density():
    assert "rho" not in residual.ResidualStressState.__dataclass_fields__;assert "rho" not in persistent.PersistentState.__dataclass_fields__


def test_local_functions_have_no_schedule_label_leakage():
    forbidden=("slow","fast","ramp_rate","rho_target","schedule_class")
    for fn in local.LOCAL_FUNCTIONS:
        src=inspect.getsource(fn).lower();assert not any(x in src for x in forbidden)


def test_short_span_or_censored_high_ratios_are_rejected():
    q=pd.read_csv("results/local_connected_sink_defect_memory/rejected_cases.csv");high=q[q.max_ratio>=1.5]
    assert len(high) and (~high.meaningful).all();assert high.rejection_reason.isin(["short_span_below_0.03","reference_nonattainment","numerical_censored"]).all()


def test_chen_not_scored_without_trajectory_gate():
    m=pd.read_csv("results/local_connected_sink_defect_memory/meaningful_trajectory_cases.csv");c=pd.read_csv("results/local_connected_sink_defect_memory/Chen_preservation_summary.csv")
    assert len(m)==0 and len(c)==0
