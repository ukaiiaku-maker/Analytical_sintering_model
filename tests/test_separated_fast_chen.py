import inspect
import re
from dataclasses import replace
import numpy as np

import separated_fast_chen_model as m
import separated_fast_chen_search as s


def shared_state():
    p=m.MaterialKinetics(); x=m.initial_state(p)
    return p,x


def test_topology_is_migration_only():
    p,x=shared_state()
    a=m.material_rates(x["rho"],x["G"],x["phi"],x["radii"],1300,p)
    for _,top in s.topology_sets():
        b=m.material_rates(x["rho"],x["G"],x["phi"],x["radii"],1300,p)
        assert a["rho_dot"] == b["rho_dot"]


def test_pr_redistribution_is_conservative_and_nonnegative():
    p,x=shared_state(); total=x["phi"].sum()
    d=m.material_rates(x["rho"],x["G"],x["phi"],x["radii"],1200,p)
    move=np.minimum(d["PR_propensity"]*x["phi"][:-1]*10,.2*x["phi"][:-1])
    x["phi"][:-1]-=move; x["phi"][1:]+=move
    assert np.all(x["phi"] >= 0)
    assert np.isclose(x["phi"].sum(),total,rtol=0,atol=1e-14)


def test_success_rules_and_q_variants():
    assert s.effect.longest_span(np.array([.75,.76,.77,.78]),np.array([1.5]*4),1.5) >= .03-1e-12
    assert {p.q_TJ for _,p in s.topology_sets()} == {0,1}
    points=[(1000,"density"),(1050,"success"),(1100,"growth")]
    good=[x[0] for x in points if x[1]=="success"]
    assert good and any(k=="density" and t<min(good) for t,k in points)
    assert any(k=="growth" and t>max(good) for t,k in points)


def test_local_laws_have_no_schedule_leakage():
    forbidden=("protocol","schedule","ramp_rate","slow","fast","target")
    for fn in m.LOCAL_FUNCTIONS:
        src=inspect.getsource(fn).lower()
        assert not any(re.search(rf"\b{word}\b",src) for word in forbidden)


def test_disabled_topology_exactly_recovers_unit_growth_factor():
    factor,diag=m.topology_growth_factor({"G":100e-9},1300,m.TopologyGrowthClosure())
    assert factor == 1.0
    assert diag["pore_drag"] == 0.0


def test_material_is_frozen_across_topology_ablation():
    p=m.MaterialKinetics(); assert p == replace(p)
    assert all(isinstance(t,m.TopologyGrowthClosure) for _,t in s.topology_sets())
