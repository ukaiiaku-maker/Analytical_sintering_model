import numpy as np
import bridge_pr_lower_bound_and_high_density_reduction as b
def test_anchor_vectors_complete():
 a=b.anchors();assert set(a)==set(b.ANCHORS);assert all(all(k in p for ks in b.BLOCKS.values() for k in ks) for p in a.values())
def test_transplant_exact():
 a=b.anchors();q=b.transplant(a['A_155976'],a['B_4412'],['PR']);assert all(q[k]==a['B_4412'][k] for k in b.BLOCKS['PR'])
def test_morph_bounds_and_partitions():
 a=b.anchors();q=b.morph(a['A_155976'],a['B_4412'],['PR'],.5);assert all(min(a['A_155976'][k],a['B_4412'][k])<=q[k]<=max(a['A_155976'][k],a['B_4412'][k]) for k in b.BLOCKS['PR'])
def test_manifest_inputs_exist():assert all((b.SCR/r['input_file']).exists() for r in b.input_manifest() if r['scratch_path'])
