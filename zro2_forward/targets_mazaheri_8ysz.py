from pathlib import Path
import pandas as pd

ROOT = Path("data/targets/mazaheri_8ysz_2008")

def final_targets(): return pd.read_csv(ROOT / "final_targets_table2.csv")

def required_inputs(barrier=Path("data/zro2/bicrystal_creep_barrier_export.json"), pdf=ROOT/"1-s2.0-S092150930800302X-main.pdf"):
    return {"barrier_json": Path(barrier).is_file(), "target_pdf": Path(pdf).is_file()}

def require_full_inputs():
    status=required_inputs()
    if not all(status.values()): raise FileNotFoundError(f"full comparison disabled; missing inputs: {[k for k,v in status.items() if not v]}")
    return status
