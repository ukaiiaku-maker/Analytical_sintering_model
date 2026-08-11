#!/usr/bin/env python3
"""Compact production PR tables while retaining raw local evidence."""
from pathlib import Path
import shutil
import pandas as pd


def main():
    root=Path("results/production_pr_desintering_assessment");path=root/"failed_or_censored_cases.csv";raw=root/"raw_failed_or_censored_cases.csv"
    if not raw.exists():shutil.copyfile(path,raw)
    d=pd.read_csv(raw,low_memory=False);d.groupby(["candidate_id","map_type","prep_growth_tolerance","second_step_growth_tolerance","boundary_status","G0_nm"],dropna=False).size().reset_index(name="n_cases").to_csv(path,index=False,lineterminator="\n")


if __name__=="__main__":main()
