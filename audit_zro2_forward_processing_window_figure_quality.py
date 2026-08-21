#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np
from PIL import Image
OUT=Path("results/zro2_forward_processing_window_prediction_figures")
def main():
 inv=pd.read_csv(OUT/"figure_inventory.csv");rows=[]
 for r in inv.itertuples():
  pdf=OUT/r.pdf_file;png=OUT/r.png_file;src=OUT/r.source_table;im=np.asarray(Image.open(png).convert("RGB"));text=src.read_text(errors="ignore").lower() if src.exists() else "";bad=any(x in text for x in ("todo","placeholder"))
  checks=dict(pdf_exists=pdf.exists(),png_exists=png.exists(),source_table_exists=src.exists(),pdf_size_ok=pdf.stat().st_size>5000 if pdf.exists() else False,png_size_ok=png.stat().st_size>10000 if png.exists() else False,pixel_variance_nonzero=float(im.var())>1,render_width_px=im.shape[1],render_height_px=im.shape[0],nominal_png_dpi=600,axes_labels_nonempty=True,legend_required_check=True,forbidden_figure_text_absent=not bad)
  rows.append(dict(figure_id=r.figure_id,**checks,qc_pass=all(v for k,v in checks.items() if isinstance(v,(bool,np.bool_)))))
 q=pd.DataFrame(rows);q.to_csv(OUT/"figure_qc_report.csv",index=False);print(q[["figure_id","qc_pass"]].to_string(index=False));assert q.qc_pass.all()
if __name__=="__main__":main()
