#!/usr/bin/env python3
"""Actual size-axis refinement around the calculated two-step onset."""
from pathlib import Path

import expanded_phase_space_exploration as study
from density_window_processing_map import write_csv


def main()->None:
    out=Path("results/expanded_phase_space")
    groups=[("size_onset_refinement",style,g0,.70,1300.,.825,study.T2_COARSE)
            for style in study.STYLES for g0 in (325.,350.,375.,400.,425.,475.,500.,550.)]
    coarse=study.run_groups(groups,4,"size-onset")
    refinement_groups=study.refinement_groups(study.boundary_rows(coarse))
    refined=study.run_groups(refinement_groups,4,"size-onset-refine")
    rows=coarse+refined
    write_csv(out/"size_onset_refinement_trajectories.csv",rows)
    write_csv(out/"size_onset_refinement_boundaries.csv",study.boundary_rows(rows))
    print(f"size onset done coarse={len(coarse)} refined={len(refined)}")


if __name__=="__main__":main()
