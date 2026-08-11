import inspect
import numpy as np

import persistent_defect_topology_stress_model as persistent
import residual_stress_memory_model as residual
import joint_heterogeneity_residual_stress_search as search
import production_mechanism_assessment as prior


def test_disabled_persistent_mode_exactly_recovers_parent_model():
    p=search.base_params("mech_019_q0",150,.65,"GBseg_rich");rp=residual.ResidualStressParams(mode="mixed_evolving",sigma_res_scale=1,stress_sign="tensile");protocol=prior.FastSchedule(20,1400,2)
    a=persistent.run(p,protocol,rp,persistent.PersistentParams(mode="disabled"),1);b=residual.run(p,protocol,rp,1)
    assert np.array_equal(a["rho"],b["rho"]) and np.array_equal(a["G"],b["G"])


def test_persistent_state_is_bounded_and_density_is_pore_conservative():
    p=search.base_params("mech_019_q0",150,.65,"GBseg_rich");rp=residual.ResidualStressParams(mode="mixed_evolving",sigma_res_scale=1,stress_sign="tensile")
    h=persistent.run(p,prior.FastSchedule(20,1350,0),rp,persistent.PersistentParams(mode="persistent_defect_memory"),1)
    assert np.all((h["f_defect_large_pore"]>=0)&(h["f_defect_large_pore"]<=1));assert np.all((h["defect_connectedness"]>=0)&(h["defect_connectedness"]<=1))
    assert all(np.min(h[k])>=0 for k in ("phi_GBseg","phi_TJ","phi_iso"));assert np.allclose(h["rho"],1-np.sum(h["phi_GBseg"]+h["phi_TJ"]+h["phi_iso"],axis=1))


def test_persistent_local_closures_have_no_schedule_leakage():
    forbidden=("slow","fast","ramp_rate","rho_target","schedule_class")
    for fn in persistent.LOCAL_FUNCTIONS:
        src=inspect.getsource(fn).lower();assert not any(x in src for x in forbidden)


def test_defect_state_has_no_direct_density_term():
    fields=set(persistent.PersistentState.__dataclass_fields__);assert "rho" not in fields and "pore_volume" not in fields

