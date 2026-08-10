#!/usr/bin/env python3
"""Calculated upper-size continuation for the expanded phase-space audit."""
from pathlib import Path

import expanded_phase_space_exploration as study
from density_window_processing_map import write_csv


def main()->None:
    out=Path("results/expanded_phase_space")
    groups=[]
    for style in study.STYLES:
        for g0 in (3000.,5000.,10000.):
            for T1 in (1200.,1300.,1400.):
                for switch in (.75,.825,.875,.925):
                    groups.append(("upper_size_extension",style,g0,.70,T1,switch,study.T2_COARSE))
    coarse=study.run_groups(groups,4,"upper-size")
    refinement_groups=study.refinement_groups(study.boundary_rows(coarse))
    refined=study.run_groups(refinement_groups,4,"upper-refine")
    rows=coarse+refined
    write_csv(out/"upper_size_extension_trajectories.csv",rows)
    write_csv(out/"upper_size_extension_boundaries.csv",study.boundary_rows(rows))
    print(f"upper done coarse={len(coarse)} refined={len(refined)}")


if __name__=="__main__":main()
