#!/usr/bin/env python3
"""Generate static equation-audit tables without running any simulation."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "equation_functional_form_audit"


def E(i, name, latex, variables, units, source, function, anchor, use, role,
      density=False, migration=False, conservative=False, note="paper-form equivalent of implemented code"):
    return dict(equation_id=i, equation_name=name, equation_latex=latex,
                variable_definitions=variables, units_or_dimensional_status=units,
                source_file=source, source_function=function, source_anchor=anchor,
                branches_results_using_it=use, evidence_role=role,
                changes_density=density, changes_migration=migration,
                conservatively_redistributes_pore_volume=conservative,
                implementation_note=note)


EQUATIONS = [
 E("FF-01","Capillary/pore activation stress",r"\sigma_{\rm loc}={4\gamma_s\over G}[1+c_\sigma{\sum_{r_i>2r_0}\phi_i\over\sum_i\phi_i}]","gamma_s surface energy; G grain size; c_sigma concentration; phi_i pore volume","Pa","separated_fast_chen_model.py","material_rates","T=T_C+273.15;stress=4*p.gamma_s","final fast-firing E0021/E0142","final_evidence"),
 E("FF-02","Stress-assisted nucleation time",r"\tau_{\rm nuc}=\nu_0^{-1}\exp[Q_{\rm nuc}/RT-v^*\sigma_{\rm loc}/(k_BT)]","nu_0 attempt frequency; Q_nuc barrier; v* activation volume","s","separated_fast_chen_model.py","material_rates","tau_nuc=math.exp","final fast-firing E0021/E0142","final_evidence"),
 E("FF-03","Exchange time",r"\tau_{\rm ex}=\tau_{{\rm ex},0}\exp(Q_{\rm ex}/RT)","tau_ex,0 prefactor; Q_ex barrier","s","separated_fast_chen_model.py","material_rates","tau_exchange=p.tau_exchange_prefactor","final fast-firing E0021/E0142","final_evidence"),
 E("FF-04","Transport time",r"\tau_{\rm tr}=\tau_{{\rm tr},0}G^2\exp(Q_{\rm tr}/RT)","tau_tr,0 prefactor; Q_tr barrier; G grain size","s","separated_fast_chen_model.py","material_rates","tau_transport=p.tau_transport_prefactor","final fast-firing E0021/E0142","final_evidence"),
 E("FF-05","Serial cycle time",r"\tau_{\rm cyc}=\tau_{\rm nuc}+\tau_{\rm ex}+\tau_{\rm tr}","three serial waiting times","s","separated_fast_chen_model.py","material_rates","cycle=tau_transport if","full model and transport-only ablation","final_evidence"),
 E("FF-06","Renewal activity",r"a={\tau_{\rm ex}+\tau_{\rm tr}\over\tau_{\rm cyc}}=[1+\tau_{\rm nuc}/(\tau_{\rm ex}+\tau_{\rm tr})]^{-1}","a renewal activity","dimensionless [0,1]","separated_fast_chen_model.py","material_rates","activity=(tau_exchange+tau_transport)/cycle","final fast-firing E0021/E0142","final_evidence"),
 E("FF-07","Connected-fine eligibility",r"\eta_{\rm geo}=\max[1-(\rho-0.82)/0.18,0]\,{\sum_{r_i\le2r_0}\phi_i\over\sum_i\phi_i}","rho density; r_i radius; r_0 reference radius","dimensionless","separated_fast_chen_model.py","material_rates","connected=max(1-(rho-.82)/.18,0)","final fast-firing E0021/E0142","implemented_proxy",note="implemented proxy, not derived physical law"),
 E("FF-08","Serial-cycle densification",r"\dot\rho=\eta_{\rm geo}\epsilon_{\rm event}\zeta_\eta/\tau_{\rm cyc}","epsilon_event density gain; zeta_eta scale ratio","s^-1","separated_fast_chen_model.py","material_rates","rho_dot=geo*p.event_strain","final fast-firing E0021/E0142","final_evidence",density=True),
 E("FF-09","Low-activity surface redistribution",r"J_{\rm PR}=k_{\rm PR}D_{s,0}e^{-Q_s/RT}r_0^{-2}(1-a)^2f_{\rm fine}","surface diffusion and PR parameters","s^-1","separated_fast_chen_model.py","material_rates","pr=p.PR_prefactor*surface*low_activity*fine","PR-off fast-firing ablation","causal_ablation",conservative=True),
 E("FF-10","Intrinsic grain growth",r"\dot G_0=D_{{\rm GB},0}e^{-Q_{\rm GB}/RT}\gamma_{\rm GB}/G","GB mobility/diffusion parameters","m s^-1","separated_fast_chen_model.py","material_rates","growth_base=gb*p.gamma_GB","fast growth and upper-bound channel","final_evidence",migration=True),
 E("FF-11","Migration-only topology factor",r"\Gamma_{\rm top}=[1+D_{\rm TJ}X_{\rm str}/(1+\Lambda/K)+D_p\eta_pC_p]^{-1}","X_str structural state; C_p coverage; Lambda/K completion proxy","dimensionless","separated_fast_chen_model.py","topology_growth_factor","return 1/(1+max(drag,0))","separated topology closures","diagnostic_only",migration=True),
 E("FF-12","Conservative adjacent-bin PR update",r"\phi_i^{n+1}=\phi_i^n-J_i\Delta t+J_{i-1}\Delta t","J_i adjacent-bin volume flux","pore-volume fraction","separated_fast_chen_model.py","run","s[\"phi\"][:-1]-=move","fast PR ablation","causal_ablation",conservative=True),

 E("PR-01","Fine-bin logistic weight",r"w_i^{\rm fine}=[1+\exp((r_i/r_0-r_m)/w_r)]^{-1}","r_m midpoint; w_r width","dimensionless","pr_desintering_memory_model.py","local_competition","fine_weight = 1.0 /","production PR memory","negative_control"),
 E("PR-02","Low-renewal gate",r"g_{\rm low}=\operatorname{sigmoid}[(a_m-a)/w_a](1-a)^{p_a}","a_m,w_a,p_a gate parameters","dimensionless","pr_desintering_memory_model.py","local_competition","low_activity = _sig","production PR memory","negative_control"),
 E("PR-03","Relative PR thermal factor",r"\theta_{\rm PR}=\exp[-Q_{\rm PR}R^{-1}(T^{-1}-T_{\rm ref}^{-1})]","Q_PR barrier; T_ref reference temperature","dimensionless","pr_desintering_memory_model.py","local_competition","thermal = math.exp","production PR memory","negative_control"),
 E("PR-04","Connected fine topology gate",r"g_{\rm top}=[\sum_i(\phi_i^{GB}+\phi_i^{TJ})w_i^{fine}/\sum_i(\phi_i^{GB}+\phi_i^{TJ})]^{p_{\rm top}}","location pore stores; p_top exponent","dimensionless","pr_desintering_memory_model.py","local_competition","topology_gate = (fine_connected / connected_total)","production PR memory","negative_control"),
 E("PR-05","PR source by bin",r"J_i=k_{\rm PR}\theta_{\rm PR}g_{\rm low}g_{\rm top}\alpha_{\rm attr}(\phi_i^{GB}+\phi_i^{TJ})w_i^{fine}(r_0/r_i)^{q_{\rm PR}}","alpha_attr attrition; q_PR size exponent","pore fraction s^-1","pr_desintering_memory_model.py","local_competition","source = J_total * propensity_bins","production PR memory","negative_control",conservative=True),
 E("PR-06","Conservative PR partitions",r"J^{GB}_{i\to i+1}=s_{sm}J_i^{GB},\ J_i^{GB\to TJ}=s_{GT}J_i^{GB},\ J_i^{TJ\to iso}=s_{TI}J_i^{TJ}","three nonnegative shares sum to one","pore fraction s^-1","pr_desintering_memory_model.py","local_competition","pr_smooth[:-1] -= move","production PR memory","negative_control",conservative=True),
 E("PR-07","Densification/PR competition",r"w_{\rm dens}=H_{\rm dens}/(H_{\rm dens}+H_{\rm PR}),\quad w_{\rm PR}=1-w_{\rm dens}","H_dens density hazard; H_PR relocation hazard","dimensionless","pr_desintering_memory_model.py","local_competition","w_dens = H_dens /","production PR memory","negative_control"),
 E("PR-08","Competition-scaled removal",r"\dot\phi_{GB,rem}=w_{\rm dens}\dot\phi^0_{GB,rem},\quad\dot\phi_{TJ,rem}=w_{\rm dens}\dot\phi^0_{TJ,rem}","superscript 0 precompetition flux","s^-1","pr_desintering_memory_model.py","local_competition","d[\"GBseg_remove\"] = old_gb * w_dens","production PR memory","negative_control",density=True),
 E("PR-09","PR/densifying work proxies",r"P_{\rm PR}=\gamma_sJ_{\rm PR}/r_0,\quad P_{\rm dens}=\sigma_{\rm act}\dot\rho","non-densifying and densifying work diagnostics","power proxy","pr_desintering_memory_model.py","local_competition","capillary_work = gamma * pr_flux","production PR diagnostics","diagnostic_only"),

 E("PL-01","Three-store pore conservation",r"\rho=1-\sum_i(\phi_i^{GB}+\phi_i^{TJ}+\phi_i^{iso})","location pore stores","dimensionless","pore_location_topology_model.py","run","s.rho=1-float(np.sum(s.phi_total))","pore-location topology","negative_control",density=True),
 E("PL-02","GB coverage",r"C_{GB}=1-\exp[-\chi_{GB}A_p^{GB}/A_{GB}],\quad A_{GB}=2/G","projected pore area and GB area","dimensionless","pore_location_topology_model.py","topology_diagnostics","C_GB=1-math.exp","pore-location topology","negative_control"),
 E("PL-03","TJ occupancy",r"C_{TJ}=1-\exp[-\chi_{TJ}L_p^{TJ}/L_{TJ}],\quad L_{TJ}=6/G^2","pore-line and TJ-line densities","dimensionless","pore_location_topology_model.py","topology_diagnostics","C_TJ=1-math.exp","pore-location/TJ models","negative_control"),
 E("PL-04","Location densification eligibility",r"\eta_{\rm loc}={C_{GB}\sum_i\phi_i^{GB}(r_0/r_i)^{q_{GB}}+\eta_{TJ}C_{TJ}\sum_i\phi_i^{TJ}(r_0/r_i)^{q_{TJ}}\over\sum_i\phi_i^{tot}}","location stores and removal exponents","dimensionless","pore_location_topology_model.py","instantaneous","eligibility=(gb_elig+tj_elig)","pore-location topology","negative_control"),
 E("PL-05","Location density rate",r"\dot\rho=\min[(1-\rho)\epsilon_{\rm event}\eta_{\rm loc}/\tau_{\rm event},0.05]","tau_event serial event time","s^-1","pore_location_topology_model.py","instantaneous","rho_dot=min((1-s.rho)","pore-location topology","negative_control",density=True),
 E("PL-06","Class-A pore pinning",r"\Gamma_A=(1+R_{GB}+R_{TJ})^{-1}","R_GB,R_TJ nonnegative resistances","dimensionless","pore_location_topology_model.py","instantaneous","mobility=1/(1+Rgb+Rtj)","pore-location topology","negative_control",migration=True),
 E("AL-01","Nonnegative action propensity",r"p_j=\max[w_jA_jg_j(T,a,\mathcal T)/R_j,0]","weight, availability, gate, resistance","relative score","pore_location_agentic_model.py","score_actions","scores={","local action allocation","negative_control",note="paper-form summary of explicit action score dictionary"),
 E("AL-02","Normalized action allocation",r"\omega_j=p_j/\sum_kp_k","nonnegative action scores","dimensionless","pore_location_agentic_model.py","score_actions","a.propensity/total_score","local action allocation","negative_control"),
 E("AL-03","Constrained pair allocation",r"(\omega_a,\omega_b)=(p_a,p_b)/(p_a+p_b)","mutually competing actions","dimensionless","pore_location_agentic_model.py","allocated_fluxes","return (x/z,y/z) if z>0","three action-pair groups","negative_control"),

 E("TJ-01","Persistent-junction production",r"\dot X_J^{prod}=[c_DC_{TJ}\dot\rho+c_R(J_{cap}+J_{reloc})+c_S(\dot G/G)C_{TJ}f_{clean}](1-X_J/X_{cap})","production coefficients and bounded capacity","s^-1","agentic_mechanism_model.py","local_mechanism","production=(p.XJ_prod_TJ","persistent-junction Chen families","negative_control",migration=True),
 E("TJ-02","Junction relaxation",r"\tau_J=\tau_{J,ref}\exp[Q_JR^{-1}(T^{-1}-T_{ref}^{-1})],\quad\dot X_J=\dot X_J^{prod}-X_J/\tau_J","relaxation time and state","s; s^-1","agentic_mechanism_model.py","local_mechanism","tau=p.tau_J_ref_s","persistent-junction Chen families","negative_control",migration=True),
 E("TJ-03","Persistent drag",r"R_J=A_JX_J(G_{ref}/G)^{q_J}","A_J drag; q_J size exponent","dimensionless","agentic_mechanism_model.py","local_mechanism","Rpersistent=p.A_J","persistent-junction Chen families","negative_control",migration=True),
 E("TJ-04","Required TJ packet count",r"K_{TJ}=K_0(G/G_{ref})^{q_{TJ}},\quad q_{TJ}\in\{0,1\}","K_0 count; q_TJ visible variant","event count","agentic_mechanism_model.py","local_mechanism","K=p.K_TJ0","q0/q1 multihit families","negative_control",migration=True),
 E("TJ-05","Pore/structural constraint split",r"C_{\rm constraint}=C_{\rm struct}+\max(1-\eta_{\rm relax}-\eta_{\rm pin},0)C_{\rm pore}","pore, structural, relaxed, pinned TJ fractions","dimensionless","agentic_mechanism_model.py","local_mechanism","C_constraint=C_structural+max","TJ constraint ablations","negative_control",migration=True,note="mixed-mode equivalent; code also implements all-TJ, relaxed, and pinned modes"),
 E("TJ-06","Expected TJ events",r"\Lambda_{TJ}=\lambda_{ref}(C_{\rm constraint}+0.05)\max(f_{clean},0.02)\theta_{TJ}(T)(G_{ref}/G)","availability, thermal and size factors","event count","agentic_mechanism_model.py","local_mechanism","Lambda=p.lambda_TJ_ref","TJ multihit families","negative_control",migration=True),
 E("TJ-07","Poisson completion",r"P_{\rm comp}=1-e^{-\Lambda}\sum_{n=0}^{\lceil K\rceil-1}\Lambda^n/n!,\quad N\sim Poisson(\Lambda)","Lambda expected; K required events","probability","agentic_mechanism_model.py","poisson_completion","return float(np.clip(1-cdf,0,1))","TJ multihit families","negative_control",migration=True),
 E("TJ-08","Coupled TJ mobility",r"\Gamma_{TJ}=P_{\rm comp}/(1+R_J+R_{TJ,pore})","completion and drag resistances","dimensionless","agentic_mechanism_model.py","local_mechanism","mobility_multiplier=gamma/(1+Rpersistent+R_pore)","persistent+multihit families","negative_control",migration=True),

 E("LR-01","Four-store density identity",r"\rho_j=1-(\phi_j^{GB}+\phi_j^{TJ}+\phi_j^{iso}+\phi_j^{closed}),\quad\rho=1-\sum_jw_j\phi_j^{tot}","local stores and region weights","dimensionless","interacting_local_region_model.py","advance","s.rho = 1.0 - (s.phi_GBseg","decoder-corrected local-region; 693168","final_evidence",density=True),
 E("LR-02","Local activity proxy",r"a_j=clip\{sigmoid[(T_C-1180)/70]\exp[-c_\sigma\sigma_j],0,1\}","sigma_j normalized residual stress","dimensionless","interacting_local_region_model.py","local_fluxes","base_activity = sigmoid","candidate 693168","implemented_proxy",note="implemented proxy, distinct from the serial fast-firing law"),
 E("LR-03","Open-pore shrinkage",r"\dot\rho_{open,j}=k_{open}\theta_\rho(T)\phi_j^{GB}C_{rem,j}\eta_{dens,j}","connected/removable eligibility","s^-1","interacting_local_region_model.py","local_fluxes","open_flux = (","candidate 693168","final_evidence",density=True),
 E("LR-04","Closed-radius proxy",r"\xi_{r,j}=clip(\phi_j^{closed}/N_j^{closed},10^{-6},10^6),\quad f_{r,j}=\xi_{r,j}^{-q_r/3}","N_closed pore-number proxy","code-normalized","interacting_local_region_model.py","local_fluxes","radius_proxy = np.clip","candidate 693168","implemented_proxy",note="implemented proxy, not a calibrated physical pore radius"),
 E("LR-05","Closed-pore shrinkage",r"\dot\rho_{closed,j}=k_{closed}\theta_{closed}(T)\phi_j^{closed}A_j(1-r_g)_+f_{r,j}","A accommodation; r_g gas-ratio proxy","s^-1","interacting_local_region_model.py","local_fluxes","closed_flux = (","candidate 693168 high-density support","final_evidence",density=True,note="implemented finite-capacity proxy"),
 E("LR-06","Local PR damage",r"J_{PR,j}=k_{PR}\theta_{PR}(T)sigmoid[(a_m-a_j)/w_a]\phi_j^{GB}","local PR parameters","s^-1","interacting_local_region_model.py","local_fluxes","pr = p[\"k_PR\"]","candidate 693168 preparation","final_evidence",conservative=True),
 E("LR-07","Conservative PR store partition",r"\phi^{GB}\xrightarrow{J_{PR}}(s_d+s_L)\phi^{GB}+s_{TJ}\phi^{TJ}+s_{iso}\phi^{iso}+s_c\phi^{closed},\quad\sum s=1","five decoded PR shares","pore-volume transfer","interacting_local_region_model.py","advance","s.phi_GBseg += (p[\"PR_damaged\"]","candidate 693168 preparation","final_evidence",conservative=True),
 E("LR-08","Conservative closed transition",r"J_{close}=k_{tr}\theta_{closed}(T)(\phi^{iso}+0.5\phi^{TJ})\rho^3","closed transition rate","s^-1","interacting_local_region_model.py","local_fluxes","close_transition = p.get(\"closed_transition\"","candidate 693168 lower boundary","final_evidence",conservative=True),
 E("LR-09","Local TJ completion proxy",r"P_{TJ}=sigmoid[(\Lambda_{eff}-K_{eff})/\sqrt{\Lambda_{eff}+K_{eff}+\epsilon}]","local event-count proxies","dimensionless","interacting_local_region_model.py","local_fluxes","P_comp = sigmoid","candidate 693168 secondary migration","implemented_proxy",migration=True,note="implemented sigmoid proxy, not exact Poisson completion"),
 E("LR-10","Local migration/growth",r"\Gamma_j=clip[P_{TJ,j}/(1+D_j),0,1],\quad\dot G_j=k_g\theta_g(T)\Gamma_j/G_j","D attached/junction/stress drag","dimensionless; grain unit s^-1","interacting_local_region_model.py","local_fluxes","migration_factor = np.clip","candidate 693168 upper boundary","final_evidence",migration=True),
 E("LR-11","Closed accommodation update",r"A_j^{n+1}=clip\{A_j^n+[1-e^{-\Delta t/\tau_A}](A_{cap}-A_j^n)-\Delta\phi_{closed,loss}/\phi_{tot},0,A_{cap}\}","capacity and recovery time","dimensionless bounded state","interacting_local_region_model.py","advance","s.closed_accommodation += cap_relax","candidate 693168 lower boundary","final_evidence",density=True,note="implemented proxy; no closed-pore Poisson Lambda/K law"),
 E("LR-12","Inverse-size exact growth step",r"G^{n+1}=\sqrt{(G^n)^2+2G^n\dot G^n\Delta t}","step-start growth rate","grain-size unit","interacting_local_region_model.py","advance","s.G = np.sqrt","candidate 693168 integration","final_evidence",migration=True),

 E("MET-01","Matched-density interpolation",r"G_r(\rho)=interp[G_r(t),\rho_r(t)],\quad G_c(\rho)=interp[G_c(t),\rho_c(t)]","jointly attained density grid only","grain-size unit","observable_trajectory_effect_audit.py","matched_curve","gr = np.interp(grid","fast and two-step trajectory audits","final_evidence_metric"),
 E("MET-02","Fast-firing ratio and rule",r"R_{fast}=G_{ref}/G_{fast};\quad \max R_{fast}\ge1.5\ \land\ span_\rho(R_{fast}\ge1.5)\ge0.03","matched-density grain sizes and longest span","dimensionless","relative_material_property_window_search.py","fast_metric","ratio.max()>=FAST_RATIO and span>=FAST_SPAN","1000 exact fast promotions","final_evidence_metric"),
 E("MET-03","Second-step growth fraction",r"g_2=(G_2-G_1)/G_1","switch and final grain size","dimensionless","audit_candidate_693168_closed_accommodation.py","classify","growth = float(frame.G_mean_nm.iloc[-1])","candidate 693168 Chen map","final_evidence_metric"),
 E("MET-04","Chen point classification",r"C=\{SUCCESS:D\land B;\ GROWTH:D\land\neg B;\ EXHAUSTION:\neg D\land B;\ MIXED:\neg D\land\neg B\}","D density attained; B growth bounded","categorical","audit_candidate_693168_closed_accommodation.py","classify","if attained and growth <= GROWTH_TOLERANCE","candidate 693168 Chen map","final_evidence_metric"),
 E("MET-05","Complete Chen window",r"W=T_{last\ success}-T_{first\ success}\ge25^\circ C,\quad B_{lower}\land B_{upper},\quad T_2<T_1","lower exhaustion and upper growth brackets","degC","adaptive_T2_boundary_search.py","status","else:label='COMPLETE_WINDOW'","dynamic/production Chen windows","final_evidence_metric"),
 E("MET-06","Two-step matched-density reduction",r"R_{TS}(\rho)=1-G_{two}(\rho)/G_{highT}(\rho)","two-step and high-T grain size","dimensionless","audit_candidate_693168_closed_accommodation.py","score_histories","median_reduction","candidate 693168 trajectory","final_evidence_metric",note="paper-form equivalent of stored reduction columns"),

 E("PROP-01","Activation-energy perturbation",r"Q_x'=Q_x+10^3\Delta Q_x","Delta Q in kJ mol^-1","J mol^-1","relative_material_property_window_search.py","apply_fast","Q_disconnection_nucleation=material","relative-property exact promotion","final_evidence_metric"),
 E("PROP-02","Prefactor perturbation",r"k_x'=f_xk_x,\quad f_x\in\{0.03,0.1,0.3,1,3,10,30\}","registered kinetic prefactors","dimensionless factor","relative_material_property_window_search.py","design","FACTORS=np.array","relative-property exact promotion","final_evidence_metric"),
 E("PROP-03","Latin-hypercube design",r"u_{ij}=[\pi_j(i)+U(0,1)]/N,\quad x_{ij}=x_j^{min}+u_{ij}(x_j^{max}-x_j^{min})","permuted stratified coordinates","dimensionless","relative_material_property_window_search.py","design","u[:,j]=(rng.permutation(n)+rng.random(n))/n","50000-row screen","screening_only"),
 E("PROP-04","Nucleation dominance",r"\Theta_{nuc}=\tau_{nuc}/(\tau_{ex}+\tau_{tr})","serial times","dimensionless","mechanism_dimensionless_groups.py","fast_groups","theta=tau_nuc/np.maximum","screen and exact-success summary","diagnostic_only"),
 E("PROP-05","Closed/growth selectivity",r"S_{closed/growth}=k_{closed}\theta_{closed}(T_2)/[k_g\theta_g(T_2)]","effective rates","dimensionless","mechanism_dimensionless_groups.py","two_step_groups","selectivity=kc/np.maximum(kg","screen and exact-success summary","diagnostic_only"),
 E("PROP-06","Surrogate fast score",r"\widehat R_{fast}=clip[1+0.796162\Delta I_g^{0.65}\Theta_{rel}^{0.12},0.5,10]","normalized exposure and nucleation ratio","surrogate","relative_material_property_window_search.py","screen","Rfast=np.clip(1+.796162","50655-row screen","screening_only",note="screening surrogate; never final evidence"),
 E("PROP-07","Surrogate two-step score",r"\widehat R_{TS}=sigmoid[logit(0.881277)+0.8\ln S_{rel}+0.15\ln\Pi_{PR}]","selectivity and preparation proxies","surrogate","relative_material_property_window_search.py","screen","logit=math.log(.881277","50655-row screen","screening_only",note="screening surrogate; never final evidence"),
 E("PROP-08","Exact joint classification",r"C_{exact}=\{both:F\land T;\ fast:F\land\neg T;\ two:\neg F\land T;\ neither:\neg F\land\neg T\}","F,T exact pass flags","categorical","relative_material_property_window_search.py","main","promotions[\"classification_exact\"]","1903 unique exact-promoted rows","final_evidence_metric"),

 E("NC-01","Nonnegative dissipation partition",r"w_m=q_m/\sum_nq_n,\quad q_m=\max(P_m,0)\max(c_m,0),\quad w_m\ge0,\sum_mw_m=1","power and compatibility proxies","dimensionless","topology_constrained_sintering.py","solve_dissipation_partition","q={n:max(f.power,0)","aggregate architecture","superseded"),
 E("NC-02","Aggregate smoothing redistribution",r"J_i=k_sW_T(T)(1-a)^{p_a}g_s\phi_i(r_0/r_i)^{q_s},\quad\dot\phi_i=-J_i+J_{i-1}","temperature and topology gates","s^-1 per bin","topology_constrained_sintering.py","surface_smoothing_redistribution","J=p.smoothing_rate_s","aggregate memory negative control","superseded",conservative=True),
 E("NC-03","Density efficiency",r"E_G=\dot\rho/(\dot G/G+\epsilon)","density gain per logarithmic growth","dimensionless ratio","topology_constrained_sintering.py","run","'E_G':f.rho_dot/(f.G_dot","cross-model diagnostic","diagnostic_only"),
]


VARIABLES = [
("T","absolute temperature","K","Arrhenius laws"),("T_C","Celsius temperature","degC","protocols"),
("rho","relative density from pore volume","1","all models"),("G","mean/local grain size","m in fast; nm-like local code unit","context required"),
("phi_i","pore-bin volume fraction","1","aggregate/separated"),("phi_GBseg","GB-segment pore store","1","location/local"),
("phi_TJ","TJ pore store","1","location/local"),("phi_iso","isolated pore store","1","location/local"),
("phi_closed","closed-pore store","1","local-region"),("N_i","pore-number proxy","code-normalized","radius proxies"),
("r_i","pore-bin radius","m","bin models"),("r_0","reference pore radius","m","fine weighting"),
("sigma_loc","fast-law activation stress","Pa","serial fast model"),("sigma_res","residual-stress state","normalized proxy","local model"),
("tau_nuc","nucleation waiting time","s","serial law"),("tau_exchange","exchange time","s","serial law"),
("tau_transport","transport time","s","serial law"),("tau_cycle","serial cycle time","s","serial law"),
("activity","renewal activity","1","serial or local proxy"),("eta_geo","connected-fine eligibility","1","fast law"),
("eta_loc","location eligibility","1","pore-location"),("rho_dot","density rate","s^-1","densification"),
("G_dot","grain-size rate","m s^-1 or local unit s^-1","migration"),("E_G","density gain/log grain growth","1","diagnostic"),
("C_GBseg","GB pore coverage","1","topology"),("C_TJ","TJ pore occupancy","1","topology"),
("f_clean_GB","clean boundary fraction","1","migration"),("X_J","persistent junction state","1","TJ closures"),
("C_TJ_pore","pore-occupied TJ fraction","1","TJ ablation"),("C_TJ_constraint","structural constraint fraction","1","TJ multihit"),
("C_TJ_relaxed","relaxed TJ fraction","1","TJ ablation"),("C_TJ_pinned","pinned TJ fraction","1","TJ drag"),
("Lambda_TJ","expected TJ events","1","Poisson TJ closure"),("K_TJ","required TJ events","1","TJ closure"),
("P_comp_TJ","TJ completion probability","1","migration"),("R_J","persistent drag","1","Class A"),
("Gamma","migration multiplier","1","growth"),("J_PR","PR redistribution rate","s^-1","PR"),
("H_dens","densification hazard","s^-1","competition"),("H_PR","PR hazard","s^-1","competition"),
("w_dens","densification weight","1","competition"),("w_PR","PR weight","1","competition"),
("A_closed","available closed accommodation","1","candidate 693168"),("A_cap","accommodation capacity","1","candidate 693168"),
("tau_A","accommodation recovery time","s","candidate 693168"),("gas_ratio","gas/sintering-pressure proxy","1","closed shrinkage"),
("rho_dot_open","open-pore density rate","s^-1","local model"),("rho_dot_closed","closed-pore density rate","s^-1","local model"),
("R_fast","matched-density grain ratio","1","fast criterion"),("g_2","second-step growth fraction","1","Chen criterion"),
("Delta_rho","continuous attained density span","1","trajectory criterion"),("Theta_nuc","nucleation dominance","1","attribution"),
("S_closed_growth","closed/growth selectivity","1","attribution"),("M_PR_closed","PR-closed screening proxy","1","screen only"),
]

# Paper-symbol aliases and derived quantities used explicitly in the numbered
# Methods equations.  They make symbol/unit coverage machine-checkable without
# changing any model equation.
VARIABLES += [
("sigma_j","local residual-stress value in region j","normalized proxy","local activity"),
("phi_tot","total local pore-volume fraction","1","local conservation"),
("w_j","normalized local-region statistical weight","1","local averaging"),
("f_fine","fine-pore volume fraction","1","fast PR eligibility"),
("eta_dens","local densification eligibility","1","open-pore shrinkage"),
("C_rem","connected removable-pore fraction","1","open-pore shrinkage"),
("g_low","low-renewal PR gate","1","PR competition"),
("g_top","connected fine-pore topology gate","1","PR competition"),
("w_i_fine","fine-radius bin weight","1","PR competition"),
("theta_PR","relative PR Arrhenius factor","1","PR laws"),
("theta_rho","relative open-shrinkage Arrhenius factor","1","local densification"),
("theta_closed","relative closed-shrinkage Arrhenius factor","1","closed densification"),
("theta_g","relative grain-growth Arrhenius factor","1","local migration"),
("Gamma_A","Class-A migration multiplier","1","migration closure"),
("Gamma_j","local-region migration multiplier","1","candidate 693168"),
("D","generic nonnegative migration resistance","1","Class-A closure"),
("D_j","local drag-resistance sum","1","local migration"),
("P_TJ","local sigmoid TJ-completion proxy","1 (proxy)","local migration"),
("X_J_prod","persistent-junction production rate","s^-1","TJ closure"),
("J_cap","TJ-to-GB capture flux","s^-1 proxy","junction production"),
("J_reloc","GB-to-TJ relocation flux","s^-1 proxy","junction production"),
("X_cap","persistent-junction capacity","1","TJ closure"),
("R_J","persistent-junction resistance","1","TJ closure"),
("xi_r","closed-pore radius proxy","code-normalized proxy","closed shrinkage"),
("Delta_phi_closed_loss","closed-pore volume removed in one step","1","accommodation consumption"),
("Delta_t","integration time step","s","state update"),
("G_1","grain size at the first-step switch","grain-size unit","Chen classification"),
("G_2","grain size at the end of the second step","grain-size unit","Chen classification"),
("G_ref","matched-density reference grain size","grain-size unit","fast metric"),
("G_fast","matched-density fast-path grain size","grain-size unit","fast metric"),
("W","finite Chen-window width","degC","Chen classification"),
("T_first_success","lowest successful second-step temperature","degC","Chen classification"),
("T_last_success","highest successful second-step temperature","degC","Chen classification"),
("R","molar gas constant","J mol^-1 K^-1","physical constant"),
("k_B","Boltzmann constant","J K^-1","physical constant"),
("a_j","local-region activity proxy","1 (proxy)","candidate 693168"),
("tau_J","persistent-junction relaxation time","s","TJ closure"),
("n","integer summation/event-count index","1","Poisson completion"),
("J_close","conservative isolated/TJ-to-closed transfer rate","s^-1","closed transition"),
]

PARAMETERS = [
("Q_disconnection_nucleation","nucleation barrier","J mol^-1","MaterialKinetics","fast final"),
("v_star","activation volume","m^3","MaterialKinetics","fast final"),("nu0_nucleation","attempt frequency","s^-1","MaterialKinetics","fast final"),
("Q_exchange","exchange barrier","J mol^-1","MaterialKinetics","fast final"),("tau_exchange_prefactor","exchange prefactor","s","MaterialKinetics","fast final"),
("Q_transport","transport barrier","J mol^-1","MaterialKinetics","fast final"),("tau_transport_prefactor","transport prefactor","s m^-2","MaterialKinetics","fast final"),
("event_strain","density gain/event","1","MaterialKinetics","fast final"),("zeta_eta_ratio","density-rate scale ratio","1","MaterialKinetics","fast final"),
("Q_GB_diffusion","growth barrier","J mol^-1","MaterialKinetics","fast/upper bound"),("D_GB_prefactor","growth prefactor","m^2 s^-1","MaterialKinetics","fast/upper bound"),
("Q_surface_diffusion","surface barrier","J mol^-1","MaterialKinetics","PR ablation"),("D_surface_prefactor","surface prefactor","m^2 s^-1","MaterialKinetics","PR ablation"),
("k_PR_ref_s","PR competition rate","s^-1","PRMemoryParams","negative control"),("Q_PR_J_mol","PR barrier","J mol^-1","PRMemoryParams","negative control"),
("renewal_gate_mid","PR activity midpoint","1","PRMemoryParams","negative control"),("renewal_gate_width","PR activity width","1","PRMemoryParams","negative control"),
("smoothing_share","smoothing partition","1","PRMemoryParams","negative control"),("GB_to_TJ_share","GB-to-TJ partition","1","PRMemoryParams","negative control"),
("TJ_to_iso_share","TJ-to-iso partition","1","PRMemoryParams","negative control"),("A_J","junction drag amplitude","1","DiscoveryParams","prior Chen"),
("tau_J_ref_s","junction relaxation time","s","DiscoveryParams","prior Chen"),("Q_relax_J_mol","junction relaxation barrier","J mol^-1","DiscoveryParams","prior Chen"),
("lambda_TJ_ref","expected-event scale","1","DiscoveryParams","prior Chen"),("K_TJ0","required-event reference","1","DiscoveryParams","prior Chen"),
("q_TJ","packet size exponent","1","DiscoveryParams/local","q0/q1"),("k_open","open shrinkage prefactor","code s^-1","local dictionary","693168"),
("Q_density","open shrinkage barrier","J mol^-1","local dictionary","693168"),("k_closed","closed shrinkage prefactor","code s^-1","local dictionary","693168"),
("Q_closed","closed shrinkage barrier","J mol^-1","local dictionary","693168"),("closed_transition","closed-transition rate","code s^-1","local dictionary","693168"),
("closed_capacity","accommodation capacity","1","local dictionary","693168"),("capacity_tau","accommodation recovery time","s","local dictionary","693168"),
("gas_ratio","gas pressure proxy","1","local dictionary","693168"),("closed_radius_exp","closed-radius exponent","1","local dictionary","693168"),
("k_PR","local PR prefactor","code s^-1","local dictionary","693168"),("Q_PR","local PR barrier","J mol^-1","local dictionary","693168"),
("PR_damaged","damaged-GB partition","1","local decoder","693168"),("PR_large","large-attached partition","1","local decoder","693168"),
("PR_TJ","TJ partition","1","local decoder","693168"),("PR_iso","isolated partition","1","local decoder","693168"),
("PR_closed","closed partition","1","local decoder","693168"),("k_growth","local growth prefactor","grain unit^2 s^-1","local dictionary","693168"),
("Q_growth","local growth barrier","J mol^-1","local dictionary","693168"),("lambda_TJ","local TJ proxy scale","1","local dictionary","693168 secondary"),
("K_TJ","local TJ proxy threshold","1","local dictionary","693168 secondary"),
]

PARAMETERS += [
("gamma_s","surface energy","J m^-2","material/local parameters","stress and PR"),
("gamma_GB","grain-boundary energy","J m^-2","MaterialKinetics","growth"),
("stress_concentration","large-pore stress concentration coefficient","1","MaterialKinetics","fast stress"),
("pore_radius0","reference pore radius","m","MaterialKinetics","fine-pore weighting"),
("alpha_attr","PR attrition multiplier","1","PRMemoryParams","negative control"),
("q_PR","PR pore-size exponent","1","PRMemoryParams","negative control"),
("p_top","PR topology-gate exponent","1","PRMemoryParams","negative control"),
("a_mid","activity-gate midpoint","1","PR/local parameters","PR gate"),
("a_width","activity-gate width","1","PR/local parameters","PR gate"),
("p_a","low-activity power","1","PRMemoryParams","PR gate"),
("c_D","TJ densification production coefficient","registered scale","DiscoveryParams","prior Chen"),
("c_R","TJ relocation production coefficient","registered scale","DiscoveryParams","prior Chen"),
("c_S","TJ sweep production coefficient","registered scale","DiscoveryParams","prior Chen"),
("X_cap","persistent-junction capacity","1","DiscoveryParams","prior Chen"),
("Q_J","junction relaxation activation energy","J mol^-1","DiscoveryParams","prior Chen"),
("T_ref","reference temperature","K","multiple parameter sets","relative Arrhenius factors"),
("q_J","persistent-junction size exponent","1","DiscoveryParams","prior Chen"),
("K_0","reference required TJ event count","1","DiscoveryParams","prior Chen"),
("G_ref_parameter","reference grain size in size-scaled closures","m or model grain unit","multiple parameter sets","migration"),
("k_tr","closed-transition prefactor","code s^-1","local parameter dictionary","candidate 693168"),
("k_g","local grain-growth prefactor","code grain unit^2 s^-1","local parameter dictionary","candidate 693168"),
("r_g","gas/sintering-pressure ratio proxy","1 (proxy)","local parameter dictionary","candidate 693168"),
("q_r","closed-radius exponent","1","local parameter dictionary","candidate 693168"),
("c_sigma","stress-inhibition or concentration coefficient","1","material/local parameters","context-specific"),
("g_tol","allowed second-step growth fraction","1","classification constants","Chen classification"),
("rho_target","target relative density","1","classification constants","attainment"),
("f_x","multiplicative kinetic-prefactor perturbation","1","property design","attribution"),
("Delta_Q_x","additive activation-energy perturbation","kJ mol^-1","property design","attribution"),
("Q_x","generic activation energy under perturbation","J mol^-1","property design","attribution"),
("k_x","generic kinetic prefactor under perturbation","process-dependent","property design","attribution"),
]

BRANCH_RESULTS = [
("Aggregate PR/topology baseline","NC-01;NC-02;NC-03","","bounded historical audit","negative fast control","negative control","topology moved boundaries"),
("Pore-location topology","PL-01..PL-06","","bounded historical audit","no nanoscale window or fast trajectory","negative control","placement insufficient"),
("Local action allocation","AL-01..AL-03","","bounded historical audit","no controlling effect","negative control","local competition insufficient"),
("Persistent junction drag","TJ-01..TJ-03","","registered family","insufficient alone","diagnostic","migration-only"),
("TJ multihit","TJ-04..TJ-07","","q0/q1 maps","earlier Chen contributor","diagnostic","not primary for 693168"),
("Persistent junction + TJ multihit","TJ-01..TJ-08","","earlier exact maps","finite Chen windows; no fast trajectory","negative control","prior Chen family"),
("Production PR/de-sintering memory","PR-01..PR-09","","production campaign","internal memory; weak fast trajectory","negative control","conservative PR"),
("Observable trajectory audit","MET-01;MET-02","","matched-density audit","superseded raw HR scoring","controlling metric","censor-aware rule"),
("Heterogeneity/residual stress","LR-02;LR-10","","bounded audit","1.662 ratio over 0.008 density","negative control","failed span"),
("Persistent defect memory","NC-03","","bounded audit","large ratio in unattained interval","negative control","unattained"),
("Local connected-sink mixture","NC-03","","bounded audit","large ratio only censored","negative control","censored"),
("Late-stage closed-pore trajectory","LR-01;LR-03;LR-05","","bounded audit","limited 0.95; no 0.98 support","negative control","precursor"),
("Separated nucleation-limited fast firing","FF-01..FF-10;MET-02","1000","exact fast promotions","55.8% full; 19.3% facile; 71.7% PR-off","controlling fast firing","base ratio 1.796/span 0.17"),
("Grain-growth-driven pore coalescence","PR-06;PL-06","","bounded audit","incomplete joint mechanism","negative control","precursor"),
("PR lower-bound plus coalescence","PR-05;PR-06;MET-05","","anchor audit","lower-bound candidates","diagnostic","bridge precursor"),
("Coupled PR-sweep state","PR-05;LR-07;LR-08","","bounded audit","interaction studied","diagnostic","decoder precursor"),
("Decoder-corrected local-region model","LR-01;LR-03;LR-05..LR-12","6","exact Tier-B base candidates","six conditional Tier-B candidates","controlling two-step","fingerprint-corrected"),
("Candidate 693168 audit","LR-01;LR-05..LR-12;MET-03..MET-06","101","fine T2 points; 12 ablations; 30 IC rows","conditional Tier B; finite window","controlling two-step","closed state uncalibrated"),
("Relative material-property window","PROP-01..PROP-08","50655","screen; 1000+1000 exact promotions","1903 exact: 485 fast; 119 two; 73 both; 1226 neither","final attribution","screen 19880 both is not evidence"),
]

REQUESTED = [
"separated_fast_chen_model.py","separated_fast_chen_search.py","relative_material_property_window_search.py",
"mechanism_dimensionless_groups.py","pr_desintering_memory_model.py","joint_pr_desintering_search.py",
"production_pr_desintering_assessment.py","pore_location_topology_model.py","pore_location_agentic_model.py",
"pore_location_agentic_sensitivity.py","mechanism_registry.py","agentic_mechanism_model.py",
"agentic_mechanism_search.py","adaptive_T2_boundary_search.py","preparation_window_search.py",
"production_mechanism_assessment.py","tj_constraint_ablation.py","interacting_local_region_model.py",
"interacting_local_region_decoder.py","local_region_decoder_corrected_postprocess.py",
"local_region_decoder_corrected_plots.py","audit_candidate_693168_closed_accommodation.py",
"audit_candidate_693168_final_plots.py","final_mechanism_synthesis_and_property_windows.py",
"relative_material_property_window_plots.py","observable_trajectory_effect_audit.py",
"topology_constrained_sintering.py"]


def location(row):
    path = ROOT / row["source_file"]
    if not path.exists():
        return "", ""
    for n, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if row["source_anchor"] in line:
            return n, line.strip()[:500]
    return "", "anchor not found"


def write(name, rows):
    rows = list(rows)
    with (OUT / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    registry = []
    for item in EQUATIONS:
        row = dict(item); row["source_line"], row["code_excerpt"] = location(row); registry.append(row)
    write("equation_registry.csv", registry)
    write("variable_definitions.csv", (dict(variable=v, definition=d, units=u, scope=s) for v,d,u,s in VARIABLES))
    write("parameter_definitions.csv", (dict(parameter=p, definition=d, units=u, source_container=c, evidence_scope=s) for p,d,u,c,s in PARAMETERS))
    keys=("equation_id","equation_name","source_file","source_function","source_line","code_excerpt","branches_results_using_it","evidence_role")
    write("equation_to_source_function.csv", ({k:r[k] for k in keys} for r in registry))
    write("branch_result_to_equation_map.csv", (dict(mechanism_closure=a,key_equation_ids=b,tested_conditions_count=c,count_basis=d,outcome=e,final_role=f,reason=g) for a,b,c,d,e,f,g in BRANCH_RESULTS))
    keys=("equation_id","equation_name","evidence_role","branches_results_using_it","implementation_note")
    write("diagnostic_only_vs_final_equations.csv", ({k:r[k] for k in keys} for r in registry))
    missing=[dict(requested_source=s,status="available_and_audited" if (ROOT/s).exists() else "missing",
                  substitute_or_note="" if (ROOT/s).exists() else "compact reports used for provenance only") for s in REQUESTED]
    missing.append(dict(requested_source="local_region_decoder_corrected_dynamic_search.py",status="missing_as_single_script",
                        substitute_or_note="interacting_local_region_decoder.py + local_region_decoder_corrected_postprocess.py + candidate audit"))
    write("missing_equation_sources.csv", missing)
    final=sum(r["evidence_role"].startswith("final_evidence") for r in registry)
    nonfinal=sum(r["evidence_role"] in {"diagnostic_only","negative_control","superseded","screening_only","implemented_proxy"} for r in registry)
    print(f"registered_equations={len(registry)}")
    print(f"defined_variables={len(VARIABLES)}")
    print(f"defined_parameters={len(PARAMETERS)}")
    print(f"source_files_audited={sum((ROOT/s).exists() for s in REQUESTED)}")
    print(f"final_evidence_equations={final}")
    print(f"diagnostic_negative_screening_equations={nonfinal}")


if __name__ == "__main__":
    main()
