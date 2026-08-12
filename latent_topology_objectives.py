#!/usr/bin/env python3
import numpy as np
def components(divergence,persistence,width,G1,prep,second,fast=True,complete=True):
    parsimony=.02*3
    tier='Tier_A' if fast and complete and prep<=.05 and second<=.05 and width>=25 and G1<=300 else ('Tier_B' if fast and complete and prep<=.10 and second<=.10 and width>=25 and G1<=450 else ('Tier_C' if complete else 'reject'))
    score=2*divergence+2*persistence+width/50-prep-second-parsimony+(2 if fast else -10)+(2 if complete else -5)
    return dict(score=score,first_step_divergence=divergence,second_step_persistence=persistence,window_width_C=width,G1_nm=G1,prep_growth=prep,second_growth=second,fast_preserved=fast,complete_window=complete,parsimony_penalty=parsimony,tier=tier)
def pareto(rows):
    out=[]
    for r in rows:
      dominated=any(q['first_step_divergence']>=r['first_step_divergence'] and q['second_step_persistence']>=r['second_step_persistence'] and q['window_width_C']>=r['window_width_C'] and q['score']>r['score'] for q in rows)
      if not dominated:out.append(r)
    return out
