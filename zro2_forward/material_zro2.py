from __future__ import annotations
from dataclasses import dataclass
import math

R = 8.31446261815324
EV_MOL = 96485.33212


@dataclass(frozen=True)
class MaterialParameters:
    kB: float = 1.380649e-23
    eV_J: float = 1.602176634e-19
    b_m: float = 3.6e-10
    nu0_sinv: float = 1.0e12
    D_GB0_m2_s: float = 0.056
    Q_GB_J_mol: float = 380000.0
    D_s0_m2_s: float = 0.10
    Q_s_J_mol: float = 380000.0
    D_GB0_unc_m2_s: float = 0.05
    Q_GB_unc_J_mol: float = 41000.0
    D_s0_unc_m2_s: float = 0.27
    Q_s_unc_J_mol: float = 28000.0
    Q_M_J_mol: float = 4.2 * EV_MOL
    M0_m4_J_s: float = 5.8e-3
    mobility_prefactor_status: str = "calibrated once to conventional-sintering final grain size"
    gamma_s_J_m2: float = 1.0
    gamma_GB_J_m2: float = 0.5
    Omega_m3: float = 3.35e-29
    C_TJ: float = 10.0
    C_GB: float = 2.0
    site_density_factor: float = 4.0

    def D_GB(self, T_K: float) -> float:
        return self.D_GB0_m2_s * math.exp(-self.Q_GB_J_mol / (R * T_K))

    def D_s(self, T_K: float) -> float:
        return self.D_s0_m2_s * math.exp(-self.Q_s_J_mol / (R * T_K))

    def M_GB(self, T_K: float) -> float:
        return self.M0_m4_J_s * math.exp(-self.Q_M_J_mol / (R * T_K))

    def triple_line_geometry(self, G_m: float) -> dict[str, float]:
        length_per_volume = self.C_TJ / G_m**2
        area_per_volume = self.C_GB / G_m
        density_per_area = length_per_volume / area_per_volume
        return {
            "L_TJ_over_V_m2": length_per_volume,
            "A_GB_over_V_minv": area_per_volume,
            "rho_TL_area_minv": density_per_area,
            "f_TL": length_per_volume * self.b_m**2,
            "eps_event": self.b_m * density_per_area,
            "width_eff_m": self.site_density_factor / density_per_area,
        }
