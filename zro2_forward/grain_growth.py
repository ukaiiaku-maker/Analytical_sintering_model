from __future__ import annotations
import numpy as np


def growth_state(G_m: float, radii_m, phi_open, T_K: float, material,
                 C_R: float = 1., C_pd: float = 1e-28, width_fraction: float = .08) -> dict[str, float]:
    moment = float(np.sum(phi_open/np.maximum(radii_m, 1e-30)))
    Rz = C_R/max(moment, 1e-300)
    width = max(width_fraction*G_m, 1e-30)
    Sz = float(1/(1+np.exp(-np.clip((2*Rz-G_m)/width, -60, 60))))
    clean = material.M_GB(T_K)*material.gamma_GB_J_m2/max(G_m, 1e-30)
    Kpd = C_pd*material.D_s(T_K)/max(Rz, 1e-30)**3
    Kclean = material.M_GB(T_K)*material.gamma_GB_J_m2
    mobile = Kpd/(Kclean+Kpd+1e-300)
    factor = mobile+(1-mobile)*Sz
    return {"G_dot_m_s": clean*factor, "G_dot_clean_m_s": clean,
            "R_Z_eff_m": Rz, "S_Z": Sz, "Gamma_growth": factor, "P_Z_Pa": material.gamma_GB_J_m2*moment}
