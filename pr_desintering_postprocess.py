#!/usr/bin/env python3
"""Keep point-level evidence locally and create review-sized production tables."""
from pathlib import Path
import shutil
import pandas as pd


def compact(path, groups, metrics):
    raw = path.with_name("raw_" + path.name)
    if not raw.exists(): shutil.copyfile(path, raw)
    df = pd.read_csv(raw)
    aggregations = {name: (name, op) for name, op in metrics}
    out = df.groupby(groups, dropna=False).agg(n_cases=("candidate_id", "size"), **aggregations).reset_index()
    out.to_csv(path, index=False, lineterminator="\n")


def main():
    root = Path("results/pr_desintering_fast_firing_memory")
    reduced = pd.read_csv(root/"reduced_joint_screen.csv")
    failed = reduced[~reduced["joint_positive"]]
    failed.to_csv(root/"failed_joint_candidates.csv", index=False, lineterminator="\n")
    failed.to_csv(root/"rejected_parameter_sets.csv", index=False, lineterminator="\n")
    response_groups = ["candidate_id","base_mechanism","variant","rho_target","initial_topology","heating_rate_C_min","peak_T_C","response_class"]
    response_metrics = [("HR_pct","median"),("PR_exposure_difference","median"),("connected_fine_difference","median"),("mean_radius_difference_nm","median")]
    compact(root/"fast_firing_response_map.csv", response_groups, response_metrics)
    compact(root/"fast_firing_response_map_full.csv", response_groups, response_metrics)
    exposure_groups = ["candidate_id","rho_target","initial_topology","heating_rate_C_min"]
    compact(root/"PR_desintering_exposure.csv", exposure_groups, [("HR_pct","median"),("cumulative_PR_desintering_work","median"),("reference_PR_work","median"),("PR_exposure_difference","median")])
    compact(root/"PR_desintering_exposure_full.csv", exposure_groups, [("HR_pct","median"),("cumulative_PR_desintering_work","median"),("reference_PR_work","median"),("PR_exposure_difference","median")])
    memory_metrics=[("connected_fine_pore_fraction","median"),("connected_fine_difference","median"),("pore_mean_radius_nm","median"),("large_pore_fraction","median")]
    compact(root/"connected_fine_pore_memory.csv", exposure_groups, memory_metrics)
    compact(root/"connected_fine_pore_memory_full.csv", exposure_groups, memory_metrics)


if __name__ == "__main__": main()
