#!/usr/bin/env python3
from massive_latent_topology_objectives import trajectory_score,chen_window
def assign_tier(score,window,exact=False):
 if not exact:return 'unscored'
 if score['attained'] and score['span20']>=.03 and score['median_reduction']>=.25 and window['complete']:return 'Tier_A'
 if score['span20']>=.02 and window['complete']:return 'Tier_B'
 return 'Tier_C' if score['max_reduction']>0 or window['complete'] else 'reject'
