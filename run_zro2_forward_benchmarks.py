#!/usr/bin/env python3
from pathlib import Path
from zro2_forward.barrier_json import BarrierInputError
from zro2_forward.campaign import run_benchmarks
def main():
    try: result=run_benchmarks()
    except BarrierInputError as e: print(f"BLOCKED: {e}"); return 2
    print(result.to_string(index=False))
    return 0
if __name__ == "__main__": raise SystemExit(main())
