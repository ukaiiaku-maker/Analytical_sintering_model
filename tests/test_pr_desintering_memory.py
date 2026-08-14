import csv
import inspect
from dataclasses import replace
from pathlib import Path
import numpy as np

import agentic_mechanism_model as baseline
import joint_pr_desintering_search as search
import pr_desintering_memory_model as memory
import production_mechanism_assessment as production


ROOT=Path(__file__).parents[1]/"results"/"pr_desintering_fast_firing_memory"


def candidate(mode="PR_plus_connected_fine_attrition"):
    p=search.frozen_base()["mech_009"]
    return memory.PRMemoryParams(p,mode,k_PR_ref_s=2e-4)


def test_local_law_has_no_schedule_or_protocol_leakage():
    assert tuple(inspect.signature(memory.local_competition).parameters)==("s","T_C","p")
    source=inspect.getsource(memory.local_competition).lower()
    for word in ("protocol","schedule","ramp_rate","slow","fast","rho_target"):
        assert word not in source


def test_disabled_mode_exactly_recovers_production_negative_control_solver():
    p=candidate("disabled");protocol=production.FastSchedule(20,1350,2)
    wrapped=memory.run(p,protocol);direct=baseline.run(p.base,protocol)
    for key in ("rho","G","phi_GBseg","phi_TJ","phi_iso","X_J"):
        assert np.array_equal(wrapped[key],direct[key])


def test_PR_flux_is_conservative_and_does_not_directly_change_density():
    p=candidate();s=memory.initial_state(p);d=memory.local_competition(s,1050,p)
    pr_gb=d["GB_smooth"]-baseline.local_mechanism(s.base,1050,p.base)["GB_smooth"]-d["GB_to_TJ"]+baseline.local_mechanism(s.base,1050,p.base)["GB_to_TJ"]
    pr_tj=d["GB_to_TJ"]-baseline.local_mechanism(s.base,1050,p.base)["GB_to_TJ"]-d["TJ_to_iso"]+baseline.local_mechanism(s.base,1050,p.base)["TJ_to_iso"]
    pr_iso=d["TJ_to_iso"]-baseline.local_mechanism(s.base,1050,p.base)["TJ_to_iso"]
    assert abs(float(np.sum(pr_gb+pr_tj+pr_iso)))<1e-18


def test_enabled_run_preserves_pores_and_nonnegative_bins():
    p=candidate();h=memory.run(p,production.FastSchedule(20,1400,2))
    phi=h["phi_GBseg"]+h["phi_TJ"]+h["phi_iso"]
    assert np.all(phi>=0) and np.all(h["N_GBseg"]>=0) and np.all(h["N_TJ"]>=0) and np.all(h["N_iso"]>=0)
    assert np.max(np.abs(h["rho"]-(1-phi.sum(axis=1))))<1e-12
    assert np.all(np.diff(h["cumulative_PR_desintering_work"])>=-1e-12)


def test_slow_path_accumulates_more_PR_memory_than_fast_path():
    p0=search.frozen_base()["mech_009"]
    p=memory.PRMemoryParams(production.fast_params(p0,75,.70,production.TOPOLOGIES["baseline"]),"PR_plus_connected_fine_attrition",k_PR_ref_s=2e-4)
    paths={rate:memory.run(p,production.FastSchedule(rate,1400,20)) for rate in (1,20)}
    idx={r:int(np.flatnonzero(h["rho"]>=.85)[0]) for r,h in paths.items()}
    assert paths[1]["cumulative_PR_desintering_work"][idx[1]]>paths[20]["cumulative_PR_desintering_work"][idx[20]]
    assert paths[20]["connected_fine_pore_fraction"][idx[20]]>paths[1]["connected_fine_pore_fraction"][idx[1]]


def test_persisted_joint_candidates_are_attainable_nonuniversal_and_keep_chen_windows():
    with (ROOT/"full_joint_scorecard.csv").open(newline="") as stream:score=list(csv.DictReader(stream))
    assert score and all(r["joint_positive"]=="True" and int(r["complete_chen_count"])>0 and int(r["beneficial_fast_count"])>0 and r["universal_fast_benefit"]=="False" for r in score)
    with (ROOT/"reduced_joint_screen.csv").open(newline="") as stream:reduced=list(csv.DictReader(stream))
    disabled=[r for r in reduced if r["variant"]=="disabled"]
    assert len(disabled)==4 and all(int(r["beneficial_fast_count"])==0 for r in disabled)
