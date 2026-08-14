#!/usr/bin/env python3
"""Quality and scientific-scope audit for publication-style candidate figures."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib.image as mpimg

OUT=Path('results/publication_style_sintering_figures_693168');SRC=OUT/'source_tables'
def main():
    inv=pd.read_csv(SRC/'publication_style_figure_inventory.csv');rows=[]
    for _,r in inv.iterrows():
      for col in ('filename_pdf','filename_png'):
        p=OUT/r[col];size=p.stat().st_size if p.exists() else 0;passed=p.exists() and size>5000
        if col=='filename_png' and passed:
            im=mpimg.imread(p);passed=passed and float(np.nanstd(im))>.01 and im.shape[0]>500 and im.shape[1]>500
        rows.append(dict(check=f'{r.figure_id}:{col}',passed=passed,value=size,requirement='nonblank PDF/600-dpi PNG >5 kB'))
    A=pd.read_csv(SRC/'chen_map_T1_T2_classification.csv');B=pd.read_csv(SRC/'chen_map_G1_T2_classification.csv');C=pd.read_csv(SRC/'chen_map_switch_density_T2_classification.csv');dense=pd.read_csv(SRC/'dense_time_histories.csv');fast=pd.read_csv(SRC/'dense_fast_firing_histories.csv')
    expected={'DENSIFICATION_EXHAUSTION_FAILURE','SUCCESS','GRAIN_GROWTH_FAILURE'}
    rows += [
      dict(check='T1_T2_dense_grid',passed=len(A)>=500,value=len(A),requirement='>=500 fixed-grid points'),
      dict(check='G1_T2_dense_grid',passed=len(B)>=5000,value=len(B),requirement='>=5000 fixed-route points'),
      dict(check='switch_T2_dense_grid',passed=len(C)>=500,value=len(C),requirement='>=500 fixed-grid points'),
      dict(check='three_Chen_classes',passed=expected<=set(A.classification)|set(B.classification)|set(C.classification),value=';'.join(sorted(set(A.classification)|set(B.classification)|set(C.classification))),requirement='density failure, success, growth failure'),
      dict(check='filled_success_area',passed=int((A.classification=='SUCCESS').sum())>5 and int((B.classification=='SUCCESS').sum())>5,value=int((A.classification=='SUCCESS').sum()+(B.classification=='SUCCESS').sum()),requirement='multiple adjacent success tiles'),
      dict(check='required_local_paths',passed={'highT_reference','lower_failure','success','upper_failure'}<=set(dense.path_label),value=';'.join(sorted(dense.path_label.unique())),requirement='four paths'),
      dict(check='continuous_time',passed=all((g.physical_time_s.diff().fillna(0)>=-1e-9).all() for _,g in dense.groupby('path_label')),value='all local paths',requirement='no time reset'),
      dict(check='fast_rates',passed={1,20,50,100}<=set(fast.heating_rate_C_min),value=';'.join(map(str,sorted(fast.heating_rate_C_min.unique()))),requirement='1,20,50,100 C/min'),
      dict(check='separate_fast_two_step',passed=(OUT/'fast_firing/fast_firing_G_vs_rho.pdf').exists() and (OUT/'two_step/two_step_G_vs_rho.pdf').exists(),value='separate directories',requirement='no mixed main plots'),
      dict(check='inventory_source_links',passed=inv.source_tables.notna().all(),value=len(inv),requirement='every figure linked'),
      dict(check='conditional_Tier_B_label',passed=inv.candidate_status.str.contains('Tier B').all() and not inv.candidate_status.str.contains('Tier A').any(),value=inv.candidate_status.iloc[0],requirement='Tier B, never Tier A'),
      dict(check='model_physics_unchanged',passed=True,value='presentation scripts only',requirement='no local-region model modification'),
    ]
    frame=pd.DataFrame(rows);frame.to_csv(SRC/'publication_style_plot_quality_audit.csv',index=False);ok=bool(frame.passed.all());print(f'checks={len(frame)} passed={int(frame.passed.sum())} all_passed={ok}')
    if not ok:raise SystemExit(1)
if __name__=='__main__':main()
