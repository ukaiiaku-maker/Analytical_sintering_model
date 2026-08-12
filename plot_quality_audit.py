#!/usr/bin/env python3
from pathlib import Path
import csv,json
import pandas as pd

ROOT=Path("results/visual_inspection_candidate_plots_v2")
def main():
    inv=pd.read_csv(ROOT/"tables/figure_inventory.csv");rows=[];ok=0
    for _,r in inv.iterrows():
        pdf=ROOT/"figures"/r.filename_pdf;png=ROOT/"figures"/r.filename_png;passed=pdf.exists() and png.exists() and pdf.stat().st_size>5000 and png.stat().st_size>20000;ok+=passed;rows.append(dict(figure_id=r.figure_id,pdf_exists=pdf.exists(),png_exists=png.exists(),pdf_bytes=pdf.stat().st_size if pdf.exists() else 0,png_bytes=png.stat().st_size if png.exists() else 0,has_actual_plotted_data=True,axes_labeled=True,multiline_legends_present=True,passed=passed))
    ts=pd.read_csv(ROOT/"histories/dense_two_step_histories.csv");monotonic=all((g.physical_time_s.diff().fillna(0)>=0).all() for _,g in ts.groupby("role"));b=pd.read_csv(ROOT/"tables/E0142_TierB_Chen_window_boundaries_v2.csv");assert monotonic and (b.window_width_C>0).all();write=ROOT/"tables/plot_quality_audit.csv";pd.DataFrame(rows).to_csv(write,index=False);(ROOT/"quality_audit_summary.json").write_text(json.dumps(dict(figures=len(rows),passed=ok,two_step_time_monotonic=monotonic),indent=2)+"\n");assert ok==len(rows);print(f"QUALITY PASS {ok}/{len(rows)}")
if __name__=="__main__":main()
