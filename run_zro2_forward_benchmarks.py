#!/usr/bin/env python3
from pathlib import Path
from zro2_forward.barrier_json import BarrierModel, BarrierInputError
def main():
    try: BarrierModel.load(Path("data/zro2/bicrystal_creep_barrier_export.json"))
    except BarrierInputError as e: print(f"BLOCKED: {e}"); return 2
    return 0
if __name__ == "__main__": raise SystemExit(main())
