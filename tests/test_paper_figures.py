from pathlib import Path
import csv
import generate_paper_figures as mainfig
import generate_supplement_figures as suppfig


def test_figure_generators_create_complete_inventories(tmp_path,monkeypatch):
    main=tmp_path/"main";supp=tmp_path/"supplement";monkeypatch.setattr(mainfig,"OUT",main);monkeypatch.setattr(suppfig,"OUT",supp)
    mainfig.main();suppfig.main()
    for root,n in ((main,14),(supp,7)):
        with (root/"figure_inventory.csv").open(newline="") as f:rows=list(csv.DictReader(f))
        assert len(rows)==n
        assert all((root/r["filename_pdf"]).is_file() and (root/r["filename_png"]).is_file() for r in rows)
        assert all(r["short_title"] and r["source_table_or_script"] and r["purpose"] and r["manuscript_location"] for r in rows)
