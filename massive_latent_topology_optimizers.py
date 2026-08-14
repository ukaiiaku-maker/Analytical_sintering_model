#!/usr/bin/env python3
"""Vectorized global sampling plus dependency-free local refinement."""
import numpy as np

def latin_hypercube(n,d,seed=20260812):
    rng=np.random.default_rng(seed);u=rng.random((n,d));x=np.empty_like(u)
    for j in range(d):x[:,j]=(rng.permutation(n)+u[:,j])/n
    return x

def pattern_search(x,score,bounds,steps=5):
    x=np.asarray(x,float);best=score(x)
    for k in range(steps):
        scale=.15/(2**k)
        for j,(lo,hi) in enumerate(bounds):
            for sign in (-1,1):
                y=x.copy();y[j]=np.clip(y[j]+sign*scale*(hi-lo),lo,hi);v=score(y)
                if v>best:x,best=y,v
    return x,best
