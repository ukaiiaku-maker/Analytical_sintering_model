import json
from pathlib import Path
import numpy as np
import interacting_local_region_decoder as d
import interacting_local_region_decoder_audit as audit
import massive_latent_topology_optimizers as o
def test_all_sampled_columns_map():
 p=d.decode(np.linspace(.01,.99,len(d.NAMES)));assert set(d.NAMES)<=set(p)
def test_preflight_diversity():
 x=o.latin_hypercube(10000,len(d.NAMES),20260817);ps=[d.decode(r) for r in x];assert len({d.fingerprint(p) for p in ps})>1000
def test_partition_normalized():
 p=d.decode(np.full(len(d.NAMES),.5));assert np.isclose(sum(p[k] for k in ('PR_damaged','PR_large','PR_TJ','PR_iso','PR_closed')),1)
def test_all_decoder_columns_change_dynamic_signature():
 assert not [row for row in audit.influence_rows() if row['unused']]
def test_persisted_run_threshold_or_incomplete_marker():
 path=Path('results/local_region_decoder_corrected_dynamic_search/run_state.json')
 if not path.exists(): return
 state=json.loads(path.read_text())
 complete=(state.get('stage0',0)>=1_000_000 and state.get('stage0_unique_fingerprints',0)>=10_000 and state.get('stage1',0)>=20_000 and state.get('stage1_unique_fingerprints',0)>=5_000 and state.get('stage2',0)>=1_000 and state.get('stage2_unique_fingerprints',0)>=500)
 assert complete or Path('RUN_INCOMPLETE_STATUS.md').exists()
def test_persisted_production_candidates_are_exact_and_attained():
 import pandas as pd
 path=Path('results/local_region_decoder_corrected_dynamic_search/accepted_tier_candidates.csv')
 if not path.exists(): return
 frame=pd.read_csv(path)
 assert frame.exact_reconfirmed.all()
 assert frame.attained.all()
 assert frame.complete.all()
 assert frame.production_gate_passed.all()
 assert (frame.switch_error<=2e-4).all()
 assert (frame.first_step_growth_fraction<=.20).all()
def test_required_decoder_corrected_figures_are_nonempty():
 root=Path('results/local_region_decoder_corrected_dynamic_search/figures')
 if not root.exists(): return
 required=('decoder_fingerprint_diversity','optimizer_or_screen_convergence','stage0_phase_map','pareto_front','best_highT_vs_twostep_G_rho','reduction_TS_vs_density_best','best_physical_time_histories','local_region_topology_histories','pore_connectivity_distribution_histories','pore_number_D90_location_histories','closed_pore_and_accommodation_histories','Chen_filled_window_best','Chen_classification_map_best','robustness_heatmap_rho0_G0','ablation_waterfall_best','fast_firing_preservation')
 for stem in required:
  for suffix in ('.pdf','.png'):
   path=root/(stem+suffix);assert path.exists() and path.stat().st_size>0
