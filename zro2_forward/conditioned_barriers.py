from __future__ import annotations
import types
import numpy as np

def barrier_with_mode(barrier,mode):
    allowed={"pchip_extrapolate","fixed_lowT_slope","generic_anchor_barrier"}
    if mode not in allowed: raise ValueError(mode)
    original=barrier.Gstar
    def gstar(self,sigma,T):
        if T>=self.temperatures_K[0]: return original(sigma,T)
        if mode=="generic_anchor_barrier": return max(.25*1.602176634e-19,3.75*1.602176634e-19-4e-29*max(sigma,0))
        # Explicit diagnostic laws operate on the final barrier values at the two lowest slices.
        y0=original(sigma,self.temperatures_K[0]); y1=original(sigma,self.temperatures_K[1]); slope=(y1-y0)/(self.temperatures_K[1]-self.temperatures_K[0])
        if mode=="fixed_lowT_slope": return max(.05*1.602176634e-19,y0+slope*(T-self.temperatures_K[0]))
        # Cubic diagnostic is represented by the local Hermite continuation; flagged in outputs.
        u=(T-self.temperatures_K[0])/(self.temperatures_K[1]-self.temperatures_K[0]); return max(.05*1.602176634e-19,y0+(y1-y0)*(u+u*(u-1)*.5))
    barrier.Gstar=types.MethodType(gstar,barrier); barrier.mode=mode; return barrier
