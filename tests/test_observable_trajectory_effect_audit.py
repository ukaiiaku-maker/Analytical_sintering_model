import numpy as np
import pandas as pd

import observable_trajectory_effect_audit as audit


def frame(rho, grain):
    return pd.DataFrame({"rho": rho, "G_nm": grain})


def test_matched_interpolation_uses_jointly_attained_density_only():
    ref=frame([.80,.85,.90],[100,120,150]);cmp=frame([.82,.86,.88],[100,105,110])
    q=audit.matched_curve(ref,cmp,step=.01)
    assert q.rho.min()>=.82 and q.rho.max()<=.88+1e-12
    assert np.allclose(q.ratio,q.G_reference_nm/q.G_comparison_nm)


def test_unattainable_window_is_not_scored():
    q=audit.matched_curve(frame([.80,.88],[100,120]),frame([.80,.88],[100,100]),step=.01)
    row=audit.window_row("x","fast",q,"intermediate",.85,.92)
    assert not row["both_paths_attained"] and row["trajectory_class"]=="unattainable"
    assert "median_ratio" not in row


def test_effect_tiers_are_mutually_exclusive():
    assert [audit.tier(x) for x in (1.1,1.3,1.6,2.1)]==["negligible","weak","meaningful","strong"]


def test_meaningful_requires_ratio_and_finite_span():
    rho=np.arange(.85,.921,.001)
    short=pd.DataFrame({"rho":rho,"ratio":np.where((rho>=.88)&(rho<=.89),1.6,1.0)})
    long=pd.DataFrame({"rho":rho,"ratio":np.where((rho>=.86)&(rho<=.90),1.6,1.0)})
    assert audit.classify(short)!="trajectory_meaningful"
    assert audit.classify(long)=="trajectory_meaningful"


def test_high_density_only_effect_is_separate():
    rho=np.arange(.95,.991,.001);q=pd.DataFrame({"rho":rho,"ratio":np.full(len(rho),1.6)})
    assert audit.classify(q)=="unsupported_high_density"


def test_rescue_design_keeps_targets_budgets_and_physics_out_of_parameters():
    forbidden={"rho_target","budget_h","t_max_s","schedule","rate"}
    for _,params in audit.rescue_design(): assert forbidden.isdisjoint(params)


def test_audit_figure_generator_runs(tmp_path,monkeypatch):
    monkeypatch.setattr(audit,"OUT",tmp_path)
    rho=np.linspace(.75,.92,20)
    base=pd.DataFrame({"rho":rho,"G_nm":100+200*(rho-.75),
        "connected_fine_pore_fraction":.3-.5*(rho-.75),
        "connected_mean_radius_nm":30+20*(rho-.75),
        "cumulative_PR_desintering_work":np.linspace(0,1,len(rho))})
    slow=base.copy();fast=base.copy();fast.G_nm*=.98
    high=base.copy();high.G_nm*=1.2;two=base.copy()
    hr=audit.matched_curve(slow,fast);ts=audit.matched_curve(high,two)
    windows=[audit.window_row("x",kind,c,"early",.75,.85) for kind,c in (("fast_firing",hr),("two_step",ts))]
    attainment=pd.DataFrame({"rho_target":[.90,.95],"comparison_attained":[True,False],"n_cases":[10,0]})
    audit.plots(pd.DataFrame(),hr,ts,windows,attainment,slow,fast,high,two)
    assert (tmp_path/"rho_G_trajectories_representative.pdf").is_file()
    assert (tmp_path/"negative_control_internal_vs_observable.png").is_file()
