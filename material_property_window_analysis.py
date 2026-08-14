#!/usr/bin/env python3
"""Summarize relative material-property orderings without ratio penalties."""
from pathlib import Path
import pandas as pd

OUT=Path("results/reframe_tierB_experimental_plausibility")


def main():
    d=pd.read_csv(OUT/"relative_material_property_window_reframed.csv")
    groups=[
        ("causal_fast_firing",d.fast_firing_pass),
        ("joint_qualitative",d.joint_qualitative_behavior_pass),
        ("joint_full_evidence",d.joint_full_evidence_pass),
    ]
    rows=[]
    for name,mask in groups:
        q=d[mask]
        for metric in ("Q_nuc_minus_Q_GB_kJ_mol","Q_surface_minus_Q_GB_kJ_mol","Q_transport_minus_Q_exchange_kJ_mol","Q_nuc_over_Q_GB","Q_surface_over_Q_GB","Q_transport_over_Q_exchange"):
            rows.append(dict(group=name,metric=metric,n_materials=len(q),minimum=q[metric].min() if len(q) else pd.NA,
                             maximum=q[metric].max() if len(q) else pd.NA,median=q[metric].median() if len(q) else pd.NA,
                             interpretation=("observed exact-material range; not an identified universal threshold" if len(q)>1
                                             else "single-material observation; no interval identified" if len(q)==1
                                             else "not identified because no fully coupled material passed")))
    pd.DataFrame(rows).to_csv(OUT/"dimensionless_group_thresholds_reframed.csv",index=False)
    print("wrote",OUT/"dimensionless_group_thresholds_reframed.csv")


if __name__=="__main__":main()
