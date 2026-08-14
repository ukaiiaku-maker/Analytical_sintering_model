from zro2_forward.targets_mazaheri_8ysz import final_targets, required_inputs
def test_target_rows_and_inputs():
    x=final_targets(); assert set(x.method)=={"CS","LMS","HMS","TSS"}
    assert "reported final state only" in open("data/targets/mazaheri_8ysz_2008/schedules.csv").read()
    assert all(required_inputs().values())
