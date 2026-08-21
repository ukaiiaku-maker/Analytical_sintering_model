from pathlib import Path
import hashlib, subprocess
import pandas as pd
from zro2_forward.material_zro2 import MaterialParameters

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/zro2_forward_final_summary_figures"
BASE=ROOT/"results/zro2_forward_processing_window_prediction_figures"
SHA="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"

def test_01_barrier_hash_unchanged():
    assert hashlib.sha256((ROOT/"data/zro2/bicrystal_creep_barrier_export.json").read_bytes()).hexdigest()==SHA

def test_02_gb_diffusivity_unchanged():
    m=MaterialParameters();assert (m.D_GB0_m2_s,m.Q_GB_J_mol)==(.056,380000.)

def test_03_surface_diffusivity_unchanged():
    m=MaterialParameters();assert (m.D_s0_m2_s,m.Q_s_J_mol)==(.10,380000.)

def test_04_closed_pore_laws_unchanged():
    changed=subprocess.check_output(["git","diff","--name-only","255446c","--","zro2_forward/closed_channel_laws.py","zro2_forward/closed_pore_evolution.py"],cwd=ROOT,text=True)
    assert not changed.strip()

def test_05_no_optimizer_or_fitting_call():
    text="\n".join((ROOT/name).read_text() for name in ("build_zro2_forward_final_summary_histories.py","plot_zro2_forward_final_summary_figures.py"))
    assert not any(token in text for token in ("differential_evolution(","least_squares(","curve_fit(","minimize(",".optimize("))

def test_06_no_broad_search_call():
    text=(ROOT/"build_zro2_forward_final_summary_histories.py").read_text()
    assert not any(token in text for token in ("ProcessPoolExecutor(","screen_one(","second_scan(","run_grid(","parameter_grid("))

def test_07_selected_cases_trace_to_existing_results():
    selected=set(pd.read_csv(OUT/"final_case_selection.csv").case_id)
    prior=set(pd.read_csv(BASE/"selected_representative_cases.csv").case_id)
    assert selected=={"P001","P014"} and selected<=prior
    h=pd.read_csv(BASE/"heating_rate_endpoint_summary.csv")
    assert {(1500.,r) for r in (.2,1.,100.)}<={tuple(x) for x in h[h.case_id.eq("P001")][["T_peak_C","rate_C_min"]].to_numpy()}

def test_08_heating_figure_has_six_panels():
    inv=pd.read_csv(OUT/"final_figure_inventory.csv").set_index("figure_id")
    assert inv.loc["final_fig1_heating_rate_response","panel_count"]==6

def test_09_twostep_figure_has_six_panels():
    inv=pd.read_csv(OUT/"final_figure_inventory.csv").set_index("figure_id")
    assert inv.loc["final_fig2_twostep_vs_isothermal_response","panel_count"]==6

def test_10_twostep_time_is_continuous():
    h=pd.read_csv(OUT/"final_twostep_histories.csv")
    for _,g in h.groupby("path_label"):
        assert (g.physical_time_s.diff().dropna()>=0).all()
        switch=g[g.path_type.ne("first_step")].physical_time_s.min()
        assert switch>0 and switch>=g[g.path_type.eq("first_step")].physical_time_s.max()

def test_11_finite_windows_have_both_boundaries():
    m=pd.read_csv(OUT/"final_chen_map_source.csv")
    finite=m[m.finite_window.fillna(False)]
    assert len(finite)>0 and finite.lower_boundary_C.notna().all() and finite.upper_boundary_C.notna().all()
    assert (finite.lower_boundary_C<finite.upper_boundary_C).all()

def test_12_required_figure_qc_passes():
    q=pd.read_csv(OUT/"final_figure_qc_report.csv")
    required=q[q.figure_file.str.contains(r"final_fig[123]_")]
    assert len(required)==3 and required.pass_qc.all()

def test_13_every_main_figure_has_source_csv():
    inv=pd.read_csv(OUT/"final_figure_inventory.csv")
    main=inv[inv.figure_group.eq("figures_main")]
    assert len(main)==3 and all((OUT/p).is_file() for p in main.source_table)

def test_14_no_forbidden_files_or_deletions_staged():
    staged=subprocess.check_output(["git","diff","--cached","--name-only","--diff-filter=ACDMRTUXB"],cwd=ROOT,text=True).splitlines()
    assert not any(p.endswith(".DS_Store") or (p.endswith(".pdf") and not p.startswith("results/zro2_forward_final_summary_figures/")) or (p.startswith("results/") and "zro2_forward_final_summary_figures/" not in p) for p in staged)
    deleted=subprocess.check_output(["git","diff","--cached","--name-only","--diff-filter=D"],cwd=ROOT,text=True).strip()
    assert not deleted

def test_15_reports_state_nonvalidation():
    docs=[ROOT/"docs/ZRO2_FORWARD_FINAL_SUMMARY_FIGURES.md",ROOT/"docs/ZRO2_FORWARD_FINAL_HEATING_RATE_FIGURE_CAPTION.md",ROOT/"docs/ZRO2_FORWARD_FINAL_TWOSTEP_FIGURE_CAPTION.md",ROOT/"docs/ZRO2_FORWARD_FINAL_FIGURE_SOURCE_DATA.md"]
    assert all("not validated" in p.read_text().lower() for p in docs)
