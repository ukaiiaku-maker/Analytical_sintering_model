#!/usr/bin/env python3
from pathlib import Path
def main():
    source=Path("results/zro2_forward_natural_pore_evolution_target_8ysz/dense_histories.csv")
    if not source.is_file(): print(f"BLOCKED: no real history file at {source}"); return 2
    return 0
if __name__ == "__main__": raise SystemExit(main())
