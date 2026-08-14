from __future__ import annotations


class RampNoHold:
    def __init__(self, rate_C_min, peak_T_C, start_C=25.):
        self.rate = rate_C_min/60; self.peak = peak_T_C; self.start = start_C
        self.t_end_s = (peak_T_C-start_C)/self.rate
    def temperature_K(self, t_s, rho): return min(self.peak, self.start+self.rate*t_s)+273.15


class RampHold(RampNoHold):
    def __init__(self, rate_C_min, peak_T_C, hold_min, start_C=25.):
        super().__init__(rate_C_min, peak_T_C, start_C); self.t_end_s += hold_min*60


class Iso:
    def __init__(self, T_C, hold_h): self.T_C=T_C; self.t_end_s=hold_h*3600
    def temperature_K(self, t_s, rho): return self.T_C+273.15


class TwoStep:
    def __init__(self, T1_C, T2_C, switch_rho, hold_h):
        self.T1_C=T1_C; self.T2_C=T2_C; self.switch_rho=switch_rho; self.t_end_s=hold_h*3600
    def temperature_K(self, t_s, rho): return (self.T2_C if rho >= self.switch_rho else self.T1_C)+273.15


ChenMapProtocol = TwoStep
