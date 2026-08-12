#!/usr/bin/env python3
"""Strict matched-density scoring for high-density two-step trajectories."""
import numpy as np

def longest_span(rho,mask):
    best=cur=0.;start=None
    for x,ok in zip(rho,mask):
        if ok:
            if start is None:start=x
            cur=x-start;best=max(best,cur)
        else:start=None
    return float(best)

def trajectory_score(high,two,lo=.95,hi=.98):
    top=min(float(np.max(high['rho'])),float(np.max(two['rho'])),hi)
    if top<lo:return dict(attained=False,min_reduction=np.nan,median_reduction=np.nan,max_reduction=np.nan,span20=0.,span30=0.,A_TS=0.,tier='reject',rejection_reason='high_density_unattained')
    rho=np.arange(lo,top+5e-5,1e-4)
    gh=np.interp(rho,high['rho'],high['G_nm']);gt=np.interp(rho,two['rho'],two['G_nm'])
    red=1-gt/gh;span20=longest_span(rho,red>=.2);span30=longest_span(rho,red>=.3)
    full=top>=hi-1e-8;tier='Tier_A' if full and span20>=.03 and np.min(red)>=.2 and np.median(red)>=.25 else ('Tier_B' if span20>=.02 else ('Tier_C' if np.max(red)>0 else 'reject'))
    return dict(attained=full,min_reduction=float(np.min(red)),median_reduction=float(np.median(red)),max_reduction=float(np.max(red)),span20=span20,span30=span30,A_TS=float(np.trapezoid(np.maximum(0,np.log(gh/gt)),rho)),tier=tier,rejection_reason='' if tier!='reject' else 'reduction_below_threshold')

def chen_window(points):
    good=sorted(x['T2_C'] for x in points if x['classification']=='success')
    low=any(x['classification']=='density_exhaustion' for x in points)
    upper=any(x['classification']=='grain_growth' for x in points)
    width=max(good)-min(good) if len(good)>=2 else 0.
    return dict(window_width_C=width,lower_bracketed=low,upper_bracketed=upper,complete=width>=25 and low and upper)
