#!/usr/bin/env python3
from pathlib import Path
import json
from pypdf import PdfReader

ROOT=Path("data/targets/mazaheri_8ysz_2008")
def main():
    pdf=ROOT/"1-s2.0-S092150930800302X-main.pdf"
    if not pdf.is_file():
        print(json.dumps({"status":"blocked_missing_pdf","pdf":str(pdf),"curve_points_created":0},indent=2)); return 2
    reader=PdfReader(pdf); text="\n".join(page.extract_text() or "" for page in reader.pages)
    checks={"title":"Processing of nanocrystalline 8 mol% yttria-stabilized zirconia" in text,
            "table1":"41.2 24.5 15–33 10.20 Cubic" in text,
            "table2":all(x in text for x in ["CS 2.14 97.5","LMS 2.35 98","HMS 0.9 98","TSS 0.29 97.6"]),
            "CS_path":"1500" in text and "5" in text,
            "pages":len(reader.pages)}
    required=[ROOT/x for x in ["powder_properties.csv","final_targets_table2.csv","density_vs_temperature_digitized.csv",
              "density_vs_time_digitized.csv","grain_size_vs_temperature_digitized.csv","grain_size_vs_density_digitized.csv"]]
    checks["outputs_nonempty"]=all(p.is_file() and p.stat().st_size>0 for p in required)
    checks["status"]="extracted_and_verified" if all(v for k,v in checks.items() if k!="pages") else "verification_failed"
    print(json.dumps(checks,indent=2)); return 0 if checks["status"]=="extracted_and_verified" else 1
if __name__ == "__main__": raise SystemExit(main())
