#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from zro2_forward.targets_mazaheri_8ysz import require_full_inputs, final_targets
from zro2_forward.campaign import ROOT, run_benchmarks
def main():
    try: require_full_inputs()
    except FileNotFoundError as e: print(f"BLOCKED: {e}"); return 2
    if not (ROOT/"benchmark_summary.csv").is_file(): run_benchmarks()
    model=pd.read_csv(ROOT/"benchmark_summary.csv"); targets=final_targets()
    lookup={"CS":"CS_thermal","HMS":"rate50_thermal"}; rows=[]
    for method,case in lookup.items():
        x=model[model.case.eq(case)].iloc[0]; y=targets[targets.method.eq(method)].iloc[0]
        rows.append({"method":method,"model_case":case,"model_final_rho":x.final_rho,"target_final_rho":y.final_density,
                     "model_final_G_um":x.final_G_um,"target_final_G_um":y.final_G_um,
                     "density_error":x.final_rho-y.final_density,"G_error_um":x.final_G_um-y.final_G_um,
                     "comparison_scope":"thermal calibration" if method=="CS" else "thermal-only comparison to microwave experiment"})
    rows.append({"method":"LMS","model_case":"CS_thermal","model_final_rho":model.loc[model.case.eq('CS_thermal'),'final_rho'].iloc[0],
                 "target_final_rho":.98,"model_final_G_um":model.loc[model.case.eq('CS_thermal'),'final_G_um'].iloc[0],"target_final_G_um":2.35,
                 "density_error":model.loc[model.case.eq('CS_thermal'),'final_rho'].iloc[0]-.98,"G_error_um":model.loc[model.case.eq('CS_thermal'),'final_G_um'].iloc[0]-2.35,
                 "comparison_scope":"same nominal ramp; microwave physics absent"})
    out=pd.DataFrame(rows); out.to_csv(ROOT/"target_comparison_summary.csv",index=False); print(out.to_string(index=False)); return 0
if __name__ == "__main__": raise SystemExit(main())
