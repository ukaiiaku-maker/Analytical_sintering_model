import math
from zro2_forward.material_zro2 import MaterialParameters, R
from zro2_forward.barrier_json import BarrierModel, BarrierInputError

def test_diffusivities_and_mobility():
    p=MaterialParameters(); T=1500.
    assert math.isclose(p.D_GB(T),.056*math.exp(-380000/(R*T)))
    assert math.isclose(p.D_s(T),.10*math.exp(-380000/(R*T)))
    assert math.isclose(p.Q_M_J_mol,4.2*96485.33212)

def test_barrier_fails_closed_when_absent(tmp_path):
    try: BarrierModel.load(tmp_path/"absent.json")
    except BarrierInputError: pass
    else: raise AssertionError("missing fitted data must fail closed")

def test_project_barrier_loads_exact_schema_and_slices():
    b=BarrierModel.load("data/zro2/bicrystal_creep_barrier_export.json")
    assert b.schema=="bicrystal_surface_triple_line_arrhenius_EXP_floor_v1"
    assert list(b.temperatures_K)==[1830.15,1953.15,2201.15,2325.15]
    assert b.Gstar(1e8,1953.15)>0
    assert not b.temperature_in_fit_range(1773.15)

def test_triple_line_units_and_formula():
    p=MaterialParameters(); G=1e-7; x=p.triple_line_geometry(G)
    assert math.isclose(x["rho_TL_area_minv"],(p.C_TJ/p.C_GB)/G)
    assert math.isclose(x["eps_event"],p.b_m*x["rho_TL_area_minv"])
