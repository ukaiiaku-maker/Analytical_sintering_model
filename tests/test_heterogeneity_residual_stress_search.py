import inspect
import numpy as np

import heterogeneous_initial_state_model as hetero
import residual_stress_memory_model as residual
import joint_heterogeneity_residual_stress_search as search
import observable_trajectory_effect_audit as effect
import production_mechanism_assessment as prior
import production_pr_desintering_assessment as production
import pr_desintering_memory_model as memory


def test_disabled_layers_recover_current_baseline():
    p=search.base_params("mech_009_q0",75,.70,"baseline")
    hp=hetero.HeterogeneousParams(p,initial_microstructure_mode="baseline_narrow",G0_mean_nm=75)
    protocol=prior.FastSchedule(20,1400,2)
    a=hetero.run(hp,protocol);b=memory.run(p,protocol)
    assert np.allclose(a["rho"],b["rho"]) and np.allclose(a["G_mean"],b["G"])
    c=residual.run(p,protocol,residual.ResidualStressParams(mode="disabled"))
    assert np.array_equal(c["rho"],b["rho"])


def test_cohort_aggregation_preserves_density_and_nonnegative_stores():
    p=search.base_params("mech_009",75,.70,"baseline")
    hp=hetero.HeterogeneousParams(p,initial_microstructure_mode="large_pore_tail",G0_mean_nm=75,large_pore_tail_fraction=.1)
    items=[]
    for spec,cp in hetero.cohort_params(hp):
        h=memory.run(cp,prior.FastSchedule(20,1350,0));items.append((spec,h))
        assert all(np.min(h[k])>=0 for k in ("phi_GBseg","phi_TJ","phi_iso"))
        assert np.allclose(h["rho"],1-np.sum(h["phi_GBseg"]+h["phi_TJ"]+h["phi_iso"],axis=1))
    a=hetero.aggregate_histories(items);w=np.array([x[0].weight for x in items]);w/=w.sum()
    assert np.isclose(a["rho"][0],sum(w[i]*items[i][1]["rho"][0] for i in range(len(items))))


def test_residual_flux_is_conservative_and_stress_has_no_direct_pore_source():
    p=search.base_params("mech_009",75,.70,"baseline");s=memory.initial_state(p);d=memory.local_competition(s,1000,p)
    rp=residual.ResidualStressParams(mode="initial_only",sigma_res_scale=1,stress_sign="tensile");rs=residual.initial_state(rp,1)
    out=residual.local_residual_coupling(d,rs,1000,rp)
    assert np.isclose(np.sum(out["GB_smooth"]),np.sum(d["GB_smooth"]),atol=1e-14)
    assert out["rho_dot"]==np.sum(out["GBseg_remove"]+out["TJ_remove"])


def test_local_closures_have_no_schedule_leakage():
    forbidden=("slow","fast","ramp_rate","rho_target","schedule_class")
    for fn in hetero.LOCAL_FUNCTIONS+residual.LOCAL_FUNCTIONS:
        src=inspect.getsource(fn).lower();assert not any(word in src for word in forbidden)


def test_observable_rule_and_high_density_flag_are_preserved():
    rho=np.arange(.85,.921,.001);ratio=np.where(rho<=.88,1.6,1.)
    assert effect.longest_span(rho,ratio,1.5)>=.03-1e-12
    import pandas as pd
    q=pd.DataFrame({"rho":np.arange(.95,.991,.001),"ratio":1.6})
    assert effect.classify(q)=="unsupported_high_density"


def test_chen_window_requires_two_boundaries_and_T2_below_T1():
    good=dict(complete_practical_window=True,lower_bracketed=True,upper_bracketed=True,T2_C=1250,T1_C=1350)
    assert search.chen_window_valid(good)
    for key in ("lower_bracketed","upper_bracketed"):bad={**good,key:False};assert not search.chen_window_valid(bad)
    assert not search.chen_window_valid({**good,"T2_C":1400})


def test_all_persisted_rejections_have_reasons():
    import pandas as pd
    q=pd.read_csv("results/heterogeneity_residual_stress_search/rejected_cases.csv")
    assert len(q) and q.rejection_reason.notna().all() and (q.rejection_reason.str.len()>0).all()
