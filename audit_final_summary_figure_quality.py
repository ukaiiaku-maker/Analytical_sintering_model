#!/usr/bin/env python3
"""File, source, and rendered-text audit for the final ZrO2 figure package."""
from pathlib import Path
import re
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

OUT = Path("results/zro2_forward_final_summary_figures")
REQUIRED = {
    "final_fig1_heating_rate_response": (["time(h)", "temperature(°c)", "relativedensity", "grainsize", "densificationrate", "effectivestress"], ["0.2c/min", "1c/min", "100c/min"]),
    "final_fig2_twostep_vs_isothermal_response": (["physicaltime(h)", "temperature(°c)", "relativedensity", "grainsize", "densityrate", "first-stepgrainsize", "second-steptemperature"], ["lowt2", "centralsuccess", "hight2", "1500°ccomparator"]),
    "final_fig3_standalone_chen_map": (["first-stepgrainsize", "second-steptemperature"], ["success", "growthfailure", "selectedpaths"]),
}

def compact(text):
    return re.sub(r"\s+", "", text).lower().replace("−", "-")

def main():
    inventory = pd.read_csv(OUT / "final_figure_inventory.csv")
    rows = []
    for item in inventory.itertuples(index=False):
        pdf, png, source = OUT/item.pdf_file, OUT/item.png_file, OUT/item.source_table
        pdf_ok, png_ok, source_ok = pdf.is_file(), png.is_file(), source.is_file()
        size_ok = pdf_ok and png_ok and source_ok and pdf.stat().st_size > 5_000 and png.stat().st_size > 20_000 and source.stat().st_size > 100
        variance_ok = False
        if png_ok:
            with Image.open(png) as im:
                sample = np.asarray(im.convert("L").resize((256, 256)), dtype=float)
            variance_ok = bool(np.var(sample) > 10.0 and np.ptp(sample) > 40.0)
        text = ""
        if pdf_ok:
            text = "\n".join(page.extract_text() or "" for page in PdfReader(pdf).pages)
        flat = compact(text)
        expected_labels, expected_legend = REQUIRED.get(item.figure_id, ([], []))
        labels_ok = all(compact(token) in flat for token in expected_labels)
        legend_ok = all(compact(token) in flat for token in expected_legend)
        if item.panel_count == 6:
            labels_ok = labels_ok and all(re.search(rf"(^|\s){letter}(\s|$)", text) for letter in "ABCDEF")
        placeholder_ok = not any(word in text.lower() for word in ("todo", "placeholder", "missing", "tbd"))
        passed = all((pdf_ok, png_ok, source_ok, size_ok, variance_ok, labels_ok, legend_ok, placeholder_ok))
        notes = "all file, variance, rendered-label, legend, panel-label, and placeholder checks passed" if passed else "inspect failed boolean fields"
        rows.append(dict(figure_file=item.pdf_file, source_table=item.source_table, pdf_exists=pdf_ok, png_exists=png_ok, source_exists=source_ok, size_ok=size_ok, pixel_variance_ok=variance_ok, labels_ok=labels_ok, legend_ok=legend_ok, placeholder_text_ok=placeholder_ok, pass_qc=passed, notes=notes))
    report = pd.DataFrame(rows)
    report.to_csv(OUT / "final_figure_qc_report.csv", index=False)
    required_rows = report[report.figure_file.str.contains("final_fig[123]_")]
    if len(required_rows) != 3 or not required_rows.pass_qc.all():
        raise SystemExit("required final figure QC failed\n" + required_rows.to_string(index=False))
    print(report.to_string(index=False))

if __name__ == "__main__":
    main()
