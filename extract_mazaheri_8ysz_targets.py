#!/usr/bin/env python3
from pathlib import Path
import json

ROOT=Path("data/targets/mazaheri_8ysz_2008")
def main():
    pdf=ROOT/"1-s2.0-S092150930800302X-main.pdf"
    if not pdf.is_file():
        print(json.dumps({"status":"blocked_missing_pdf","pdf":str(pdf),"curve_points_created":0},indent=2)); return 2
    print(json.dumps({"status":"pdf_present_manual_digitization_required","pdf":str(pdf)},indent=2)); return 0
if __name__ == "__main__": raise SystemExit(main())
