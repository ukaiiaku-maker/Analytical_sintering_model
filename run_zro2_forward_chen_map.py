#!/usr/bin/env python3
from zro2_forward.targets_mazaheri_8ysz import require_full_inputs
from zro2_forward.campaign import run_map, write_state
def main():
    try: require_full_inputs()
    except FileNotFoundError as e: print(f"BLOCKED: {e}"); return 2
    frame,best=run_map(); state=write_state()
    print(best.to_string(index=False)); print(state); return 0
if __name__ == "__main__": raise SystemExit(main())
