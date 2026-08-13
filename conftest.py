"""Pytest policy for intentionally archived historical result fixtures."""
from pathlib import Path
import pytest

CURRENT={
 'test_bridge_pr_reduction.py','test_pr_lower_bound_coalescence.py',
 'test_grain_growth_pore_coalescence.py','test_massive_latent_topology_search.py',
 'test_coupled_pr_sweep_state.py',
}
ARCHIVED_MODULES={'test_paper_figures.py'}

def pytest_configure(config):
    config.addinivalue_line('markers','requires_archived_results: requires historical result artifacts that may be stored in results/Archive.zip')

def pytest_collection_modifyitems(config,items):
    """Mark only tests whose modules depend on missing archived result trees."""
    root=Path(__file__).parent
    for item in items:
        path=Path(str(item.fspath))
        if path.name in CURRENT:continue
        try:src=path.read_text()
        except OSError:continue
        if path.name in ARCHIVED_MODULES or ('results/' in src or "'results'" in src or '"results"' in src):
            item.add_marker(pytest.mark.requires_archived_results)
