#!/usr/bin/env python3
"""Compact analysis entry point for completed optimizer outputs."""
from pathlib import Path
import pandas as pd
def main():
 p=Path('results/optimizer_latent_topology_memory_search');d=pd.read_csv(p/'candidate_scorecard.csv');print(d[['candidate_id','tier','score','first_step_divergence','second_step_persistence','window_width_C']].to_string(index=False))
if __name__=='__main__':main()
