#!/usr/bin/env python3
"""Complete mapping from sampled columns to local-region dynamic parameters."""
import hashlib,json,numpy as np
NAMES=('N_regions','weight_sigma','rho_sigma','G_sigma','degree','connected_init','closed_init','cluster','exchange_rate','k_PR','Q_PR','activity_mid','activity_width','PR_damaged','PR_large','PR_TJ','PR_iso','PR_closed','PR_tau','k_sweep_damaged','k_sweep_connected','sweep_exp','coalescence_exp','number_loss','detachment','recapture','closed_transition','k_closed','Q_closed','closed_capacity','capacity_tau','gas_ratio','closed_radius_exp','attached_drag','junction_drag','XJ_prod','XJ_tau','lambda_TJ','K_TJ','q_TJ','pore_relax','pore_drag_fraction','stress_PR','stress_shear','stress_migration','stress_nucleation','stress_tau')
def lg(x,a,b):return 10**(np.log10(a)+x*(np.log10(b)-np.log10(a)))
def decode(x):
 x=np.asarray(x);p=dict(zip(NAMES,map(float,x)));p.update(N_regions=(8,16,32)[min(int(x[0]*3),2)],weight_sigma=.05+1.45*x[1],rho_sigma=.001+.099*x[2],G_sigma=.001+.799*x[3],degree=1+int(x[4]*6),connected_init=.05+.9*x[5],closed_init=.4*x[6],cluster=.6*x[7],exchange_rate=lg(x[8],1e-10,1e-3),k_PR=lg(x[9],1e-9,1e-2),Q_PR=100e3+350e3*x[10],activity_mid=.02+.58*x[11],activity_width=.02+.28*x[12],PR_tau=lg(x[18],1e3,1e8),k_sweep_damaged=lg(x[19],1e-10,1e-2),k_sweep_connected=lg(x[20],1e-10,1e-2),sweep_exp=.5+3.5*x[21],coalescence_exp=.5+5.5*x[22],number_loss=.01+.89*x[23],detachment=lg(x[24],1e-10,1e-2),recapture=lg(x[25],1e-10,1e-2),closed_transition=lg(x[26],1e-10,1e-2),k_closed=lg(x[27],1e-12,1e-3),Q_closed=250e3+500e3*x[28],closed_capacity=.01+.99*x[29],capacity_tau=lg(x[30],1e3,1e8),gas_ratio=x[31],closed_radius_exp=5*x[32],attached_drag=500*x[33],junction_drag=500*x[34],XJ_prod=lg(x[35],1e-4,50),XJ_tau=lg(x[36],1e2,1e9),lambda_TJ=lg(x[37],1e-3,1e3),K_TJ=1+49*x[38],q_TJ=min(int(x[39]*3),2),pore_relax=x[40],pore_drag_fraction=x[41],stress_PR=100*x[42],stress_shear=100*x[43],stress_migration=100*x[44],stress_nucleation=10*x[45],stress_tau=lg(x[46],1e2,1e9))
 parts=np.array([x[i] for i in range(13,18)]);parts/=max(parts.sum(),1e-30)
 for k,v in zip(('PR_damaged','PR_large','PR_TJ','PR_iso','PR_closed'),parts):p[k]=float(v)
 return p
# The dynamic-equivalence fingerprint contains every decoded physical parameter,
# not merely the minimum reporting subset.  This prevents two candidates that
# differ in (for example) PR partition or stress relaxation from being counted
# as the same local dynamics.
FINGERPRINT=NAMES
def fingerprint(p):return hashlib.sha1(json.dumps([round(float(p[k]),12) for k in FINGERPRINT]).encode()).hexdigest()[:16]
