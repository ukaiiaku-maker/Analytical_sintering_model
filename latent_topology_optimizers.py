#!/usr/bin/env python3
import numpy as np
def pso(evaluate,bounds,n_particles=64,iterations=25,seed=17,callback=None):
 rng=np.random.default_rng(seed);lo=np.array([x[0] for x in bounds]);hi=np.array([x[1] for x in bounds]);x=lo+(hi-lo)*rng.random((n_particles,len(bounds)));v=np.zeros_like(x);pb=x.copy();ps=np.array([evaluate(z) for z in x]);gb=pb[np.argmax(ps)].copy();trace=[]
 for it in range(iterations):
  v=.65*v+1.4*rng.random(x.shape)*(pb-x)+1.4*rng.random(x.shape)*(gb-x);x=np.clip(x+v,lo,hi);s=np.array([evaluate(z) for z in x]);better=s>ps;pb[better]=x[better];ps[better]=s[better];gb=pb[np.argmax(ps)].copy();trace.append(dict(iteration=it+1,best_score=float(ps.max()),median_score=float(np.median(s))))
  if callback:callback(it+1,x,s,gb,trace)
 return gb,trace
def pattern_search(x,evaluate,bounds,steps=20):
 x=x.copy();best=evaluate(x);delta=np.array([(b-a)*.1 for a,b in bounds])
 for _ in range(steps):
  for j in range(len(x)):
   for sign in (-1,1):
    y=x.copy();y[j]=np.clip(y[j]+sign*delta[j],bounds[j][0],bounds[j][1]);s=evaluate(y)
    if s>best:x,best=y,s
  delta*=.7
 return x,best
