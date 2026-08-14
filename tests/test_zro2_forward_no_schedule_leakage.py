from pathlib import Path
def test_local_laws_have_no_path_labels():
    forbidden=("protocol","schedule","fast","slow","two_step","target","hms","lms","tss")
    for name in ("densification.py","energy_balance.py","grain_growth.py","pore_population.py"):
        text=(Path("zro2_forward")/name).read_text().lower()
        assert not any(word in text for word in forbidden), (name,[w for w in forbidden if w in text])
