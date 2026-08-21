from __future__ import annotations

import hashlib
import inspect
from pathlib import Path

import numpy as np
import pandas as pd

import derive_zro2_energy_ledger_and_closed_pore_laws as d

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/zro2_forward_energy_ledger_closed_pore_derivation"


def test_barrier_json_hash_unchanged():
    assert hashlib.sha256(d.BARRIER_PATH.read_bytes()).hexdigest()=="fa7c9a2fb30e55596d9cdf47d2b03d530ef073cbe15f6210282285b2fcee7f37"


def test_diffusivity_laws_unchanged():
    T=1473.15; R=8.31446261815324
    assert np.isclose(d.MAT.D_GB(T),.056*np.exp(-380000/(R*T)))
    assert np.isclose(d.MAT.D_s(T),.10*np.exp(-380000/(R*T)))


def test_failed_global_mobility_fit_not_activated():
    assert d.MAT.mobility_prefactor_status=="calibrated once to conventional-sintering final grain size"
    assert np.isclose(d.MAT.M0_m4_J_s,5.8e-3)


def test_no_physical_qclosed_input():
    sig=inspect.signature(d.closed_state_rate)
    assert not any("q_closed" in n.lower() for n in sig.parameters)
    reg=pd.read_csv(OUT/"closed_pore_law_registry.csv")
    assert not reg.loc[reg.law_id!="empirical_reduced_closure","uses_empirical_Q"].any()


def test_energy_channels_named_and_nonnegative_expenditures():
    reg=pd.read_csv(OUT/"energy_ledger_channel_registry.csv")
    required={"P_open_dens","P_closed_dens","P_surface_smooth","P_pore_coarsen","P_GB_growth","P_drag","P_gas","P_other"}
    assert required.issubset(set(reg.channel))
    h=pd.read_csv(OUT/"energy_ledger_diagnostic_histories.csv")
    for c in required: assert (h[c+"_W_m3"]>=0).all()


def test_conservative_transfers_preserve_pore_volume():
    p=np.array([.2,.1,.05]); q=d.conservative_transfer(p,.002,12,0,2)
    assert np.isclose(p.sum(),q.sum()) and (q>=0).all()


def test_surface_accommodation_does_not_densify():
    z=d.closed_state_rate("surface_diffusion_accommodation_only",1473.15,25e-9,.1,.5,.1)
    assert z["shape_rate_sinv"]>0 and z["rho_dot_closed_sinv"]==0


def test_closed_shrinkage_requires_inventory_and_stress():
    for mode in ("renewal_limited_closed_shrinkage","GB_diffusion_closed_shrinkage","gas_limited_closed_shrinkage"):
        assert d.closed_state_rate(mode,1473.15,25e-9,0,1,1)["rho_dot_closed_sinv"]==0
        assert d.closed_state_rate(mode,1473.15,25e-9,.1,1,1,1.1)["rho_dot_closed_sinv"]==0


def test_gas_pressure_reduces_driving_stress():
    a=d.closed_state_rate("gas_limited_closed_shrinkage",1473.15,25e-9,.1,1,1,0)
    b=d.closed_state_rate("gas_limited_closed_shrinkage",1473.15,25e-9,.1,1,1,.9)
    assert b["sigma_Pa"]<a["sigma_Pa"] and b["rho_dot_closed_sinv"]<a["rho_dot_closed_sinv"]


def test_density_identity():
    assert np.isclose(d.density_identity(np.array([.1,.2]),np.array([.03]),np.array([.07])),.6)


def test_migration_modifiers_do_not_change_density_directly():
    g=pd.read_csv(OUT/"growth_pinning_law_registry.csv")
    modifiers=g[g.term!="intrinsic"]
    assert modifiers.coupling.str.contains("changes migration only").all()
    assert not g.coupling.str.contains("density",case=False).any()


def test_local_law_has_no_processing_labels():
    source=inspect.getsource(d.closed_state_rate).lower()
    forbidden=("cs","lms","hms","tss","fast","slow","protocol","schedule","ramp_rate","target")
    # Token-aware check prevents false positives inside ordinary words.
    import re
    tokens=set(re.findall(r"[a-z_]+",source))
    assert not tokens.intersection(forbidden)


def test_empirical_closure_is_labeled_diagnostic():
    z=d.closed_state_rate("empirical_reduced_closure",1473.15,25e-9,.1,1,1)
    assert z["physical_status"]=="empirical_diagnostic_only"


def test_any_promoted_law_has_registry_and_unit_audit():
    dec=pd.read_csv(OUT/"law_acceptance_decision.csv")
    reg=set(pd.read_csv(OUT/"closed_pore_law_registry.csv").law_id)
    unit=set(pd.read_csv(OUT/"closed_pore_unit_audit.csv").law_id)
    for law in dec.loc[dec.promoted_to_bounded_map,"law_id"]:
        assert law in reg and law in unit


def test_figures_have_sources_and_no_placeholders():
    inv=pd.read_csv(OUT/"figure_inventory.csv")
    assert len(inv)==5 and inv.pdf_nonempty.all() and inv.png_nonempty.all()
    assert not inv.placeholder.any()
    for source in inv.source_table: assert (OUT/source).is_file()


def test_reports_explicitly_disclaim_validation():
    names=["ZRO2_FORWARD_ENERGY_LEDGER_MODEL_REVISION.md","ZRO2_FORWARD_CLOSED_PORE_PHYSICAL_LAW_DERIVATION.md",
           "ZRO2_FORWARD_REDUCED_TO_PHYSICAL_PROPERTY_MAPPING.md","ZRO2_FORWARD_PHYSICAL_LAW_TEST_RESULTS.md",
           "ZRO2_FORWARD_NEXT_DECISION_AFTER_PHYSICAL_LAW_AUDIT.md"]
    for name in names:
        text=(ROOT/"docs"/name).read_text().lower()
        assert "no validation claim" in text or "not validation" in text
