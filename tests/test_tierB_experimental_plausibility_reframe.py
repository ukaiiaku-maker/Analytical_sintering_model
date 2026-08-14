from pathlib import Path
import pandas as pd


ROOT=Path("results/reframe_tierB_experimental_plausibility")
DOCS=Path("docs")
IDS={693168,822940,581668,295003,366094,85161}


def test_all_six_candidates_remain_explicit_tierB_interpretations():
    d=pd.read_csv(ROOT/"tierB_candidate_reinterpretation.csv")
    assert set(d.candidate_id)==IDS
    assert d.interpretation.str.contains("Tier B").all()
    assert not d.artifact_flag.any()


def test_large_ratio_does_not_create_artifact_flag():
    d=pd.read_csv(ROOT/"tierB_candidate_reinterpretation.csv").set_index("candidate_id")
    c=d.loc[693168]
    assert c.median_reduction > .5
    assert bool(c.high_density_interval_attained)
    assert bool(c.plausible_large_separation_flag)
    assert not bool(c.artifact_flag)


def test_artifact_logic_uses_explicit_conditions_not_ratio_magnitude():
    source=Path("reframe_tierB_experimental_plausibility.py").read_text()
    artifact_block=source.split("artifact = bool(",1)[1].split(")",1)[0]
    assert "median_reduction" not in artifact_block
    assert "numerical_flags" in artifact_block and "ablation_causal" in artifact_block


def test_captions_avoid_unsupported_implausibility_language():
    paths=[DOCS/"PUBLICATION_STYLE_FIGURE_CAPTIONS_693168.md",DOCS/"CANDIDATE_693168_FINAL_CAPTION_DRAFTS.md",DOCS/"UPDATED_CAPTION_LANGUAGE_TIERB.md"]
    text="\n".join(p.read_text().lower() for p in paths)
    assert "physically implausible" not in text
    assert "not because the ratio is large" in text


def test_reports_retain_calibration_and_nonvalidation_limits():
    paths=[DOCS/"TIERB_EXPERIMENTAL_PLAUSIBILITY_REFRAME.md",DOCS/"CANDIDATE_693168_REVISED_INTERPRETATION.md",DOCS/"SIX_TIERB_CANDIDATE_REINTERPRETATION.md"]
    for path in paths:
        text=path.read_text().lower()
        assert "calibrat" in text
        assert "validat" in text


def test_material_window_does_not_modify_topology_parameters():
    d=pd.read_csv(ROOT/"relative_material_property_window_reframed.csv")
    assert len(d)==512
    assert not d.topology_parameters_modified.any()
    source=Path("material_property_window_search.py").read_text()
    assert "interacting_local_region_model" not in source
    assert "topology_parameters_modified=False" in source


def test_required_figure_pairs_exist_and_are_nonempty():
    names=["tierB_candidate_family_summary","large_separation_not_artifact_G_rho",
           "experimental_plausibility_map_highT_vs_twostep_G","closed_fraction_vs_reduction_TierB",
           "TierB_window_width_vs_reduction","relative_material_property_window_summary"]
    for name in names:
        for ext in ("pdf","png"):
            p=ROOT/"figures"/f"{name}.{ext}"
            assert p.exists() and p.stat().st_size>10_000
