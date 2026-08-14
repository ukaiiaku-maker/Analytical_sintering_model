#!/usr/bin/env python3
from zro2_forward.targets_mazaheri_8ysz import require_full_inputs
def main():
    try: require_full_inputs()
    except FileNotFoundError as e: print(f"BLOCKED: {e}"); return 2
    return 0
if __name__ == "__main__": raise SystemExit(main())
