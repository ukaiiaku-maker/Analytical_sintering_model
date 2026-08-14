# Methods text with equations for the paper

## 1. Serial nucleation-limited densification law

The fast-firing calculation used the exact material law in separated_fast_chen_model.material_rates. With \(T=T_C+273.15\),

\[
\sigma_{\rm loc}={4\gamma_s\over G}\left[1+c_\sigma{\sum_{r_i>2r_0}\phi_i\over\sum_i\phi_i}\right], \tag{1}
\]
\[
\tau_{\rm nuc}=\nu_0^{-1}\exp\left({Q_{\rm nuc}\over RT}-{v^*\sigma_{\rm loc}\over k_BT}\right),\quad
\tau_{\rm ex}=\tau_{{\rm ex},0}e^{Q_{\rm ex}/RT},\quad
\tau_{\rm tr}=\tau_{{\rm tr},0}G^2e^{Q_{\rm tr}/RT}. \tag{2}
\]

\(\gamma_s\) is surface energy; \(G\) is grain size; \(c_\sigma\) is the pore-stress concentration; \(r_i,\phi_i\) are pore-bin radius and volume; \(r_0\) is the reference radius; \(\nu_0,Q_{\rm nuc},v^*\) are nucleation parameters; and \(Q_{\rm ex},Q_{\rm tr}\) and their prefactors define exchange and transport. The implemented exponential arguments are clipped to \([-50,50]\).

The cycle is serial:

\[
\tau_{\rm cyc}=\tau_{\rm nuc}+\tau_{\rm ex}+\tau_{\rm tr},\qquad
a={\tau_{\rm ex}+\tau_{\rm tr}\over\tau_{\rm cyc}}. \tag{3}
\]

The implemented geometric eligibility and density rate are

\[
\eta_{\rm geo}=\max\left[1-{\rho-0.82\over0.18},0\right]
{\sum_{r_i\le2r_0}\phi_i\over\sum_i\phi_i},\qquad
\dot\rho={\eta_{\rm geo}\epsilon_{\rm event}\zeta_\eta\over\tau_{\rm cyc}}. \tag{4}
\]

\(\epsilon_{\rm event}\) is density gain per event and \(\zeta_\eta\) is the registered rate-scale ratio. The density gate is an implemented proxy, not a derived percolation law.

## 2. Non-densifying PR/de-sintering competition

The separated model uses conservative adjacent-bin redistribution:

\[
J_{\rm PR}=k_{\rm PR}D_{s,0}e^{-Q_s/RT}r_0^{-2}(1-a)^2f_{\rm fine},\qquad
\dot\phi_i=-J_i+J_{i-1}. \tag{5}
\]

The earlier production PR closure used

\[
g_{\rm low}=\operatorname{sigmoid}\left({a_m-a\over w_a}\right)(1-a)^{p_a},\quad
g_{\rm top}=\left[{\sum_i(\phi_i^{GB}+\phi_i^{TJ})w_i^{fine}\over\sum_i(\phi_i^{GB}+\phi_i^{TJ})}\right]^{p_{\rm top}}, \tag{6}
\]
\[
J_i=k_{\rm PR}\theta_{\rm PR}g_{\rm low}g_{\rm top}\alpha_{\rm attr}
(\phi_i^{GB}+\phi_i^{TJ})w_i^{fine}(r_0/r_i)^{q_{\rm PR}}. \tag{7}
\]

The source is partitioned conservatively among GB smoothing, GB-to-TJ relocation, and TJ-to-isolated relocation. Competition is

\[
w_{\rm dens}={H_{\rm dens}\over H_{\rm dens}+H_{\rm PR}},\qquad w_{\rm PR}=1-w_{\rm dens}. \tag{8}
\]

Only explicit pore-removal fluxes change density. PR is not the controlling fast-firing channel: PR-off frequently preserves that response.

## 3. Pore-state conservation and density identity

The location model satisfies

\[
\rho=1-\sum_i(\phi_i^{GB}+\phi_i^{TJ}+\phi_i^{iso}). \tag{9}
\]

The decoder-corrected model adds a closed store and local-region weights \(w_j\):

\[
\rho_j=1-(\phi_j^{GB}+\phi_j^{TJ}+\phi_j^{iso}+\phi_j^{closed}),\qquad
\rho=1-\sum_jw_j\phi_j^{tot}. \tag{10}
\]

Relocation, PR partition, and closed transition are conservative. Density changes only through open- and closed-pore shrinkage; stores remain nonnegative.

## 4. Grain growth and migration activity

The intrinsic separated-material law is

\[
\dot G_0={D_{{\rm GB},0}e^{-Q_{\rm GB}/RT}\gamma_{\rm GB}\over G}. \tag{11}
\]

Class-A drag uses \(\Gamma_A=(1+D)^{-1}\). Class-B multihit uses the completion probability below. Coupled closures use \(\Gamma=P_{\rm comp}/(1+D)\). These factors modify \(\dot G\), not the material \(\dot\rho\), at a shared state.

The local-region model uses

\[
\Gamma_j=\operatorname{clip}\left[{P_{TJ,j}\over1+D_j},0,1\right],\qquad
\dot G_j=k_g\theta_g(T){\Gamma_j\over G_j}. \tag{12}
\]

The step is integrated as \(G_{n+1}=[G_n^2+2G_n\dot G_n\Delta t]^{1/2}\).

## 5. Persistent junction and TJ multihit mechanisms

\[
\dot X_J^{prod}=\left[c_DC_{TJ}\dot\rho+c_R(J_{cap}+J_{reloc})
+c_S{\dot G\over G}C_{TJ}f_{clean}\right](1-X_J/X_{cap}), \tag{13}
\]
\[
\tau_J=\tau_{J,ref}\exp\left[{Q_J\over R}\left({1\over T}-{1\over T_{ref}}\right)\right],\quad
\dot X_J=\dot X_J^{prod}-{X_J\over\tau_J},\quad
R_J=A_JX_J(G_{ref}/G)^{q_J}. \tag{14}
\]

The multihit closure uses

\[
K_{TJ}=K_0(G/G_{ref})^{q_{TJ}},\quad
P_{\rm comp}=1-e^{-\Lambda_{TJ}}\sum_{n=0}^{\lceil K_{TJ}\rceil-1}{\Lambda_{TJ}^n\over n!}. \tag{15}
\]

\(q_{TJ}=0\) and 1 remain visible. Pore occupancy, structural constraint, relaxed occupancy, and pinned occupancy are separate states. Persistent plus multihit generated prior Chen families but is not primary for candidate 693168.

## 6. Closed-pore shrinkage and finite accommodation

Candidate 693168 uses a local activity proxy, not Eqs. (2)–(3):

\[
a_j=\operatorname{clip}\{\operatorname{sigmoid}[(T_C-1180)/70]e^{-c_\sigma\sigma_j},0,1\}. \tag{16}
\]

The open and closed density rates are

\[
\dot\rho_{open,j}=k_{open}\theta_\rho(T)\phi_j^{GB}C_{rem,j}\eta_{dens,j}, \tag{17}
\]
\[
\xi_{r,j}=\operatorname{clip}(\phi_j^{closed}/N_j^{closed},10^{-6},10^6),\quad
\dot\rho_{closed,j}=k_{closed}\theta_{closed}(T)\phi_j^{closed}A_j(1-r_g)_+\xi_{r,j}^{-q_r/3}. \tag{18}
\]

\(A_j\) is available accommodation and \(r_g\) is a gas-ratio proxy. Accommodation recovers toward a finite capacity and is consumed by closed-pore removal:

\[
A_j^{n+1}=\operatorname{clip}\left\{A_j^n+[1-e^{-\Delta t/\tau_A}](A_{cap}-A_j^n)
-{\Delta\phi_{closed,loss}\over\phi_{tot}},0,A_{cap}\right\}. \tag{19}
\]

Equation (19) is an implemented proxy, not a derived physical law. There is no closed-pore Poisson \(\Lambda/K\) law. The audit therefore leaves closed \(\Lambda\) and \(K\) unavailable.

First-step preparation uses

\[
J_{PR,j}=k_{PR}\theta_{PR}(T)\operatorname{sigmoid}\left({a_m-a_j\over w_a}\right)\phi_j^{GB}, \tag{20}
\]

with conservative partition among GB, TJ, isolated, and closed stores, plus

\[
J_{close,j}=k_{tr}\theta_{closed}(T)(\phi_j^{iso}+0.5\phi_j^{TJ})\rho_j^3. \tag{21}
\]

Removing PR damage, closed transition, closed shrinkage, or finite accommodation destroys the candidate's joint result. Its large closed store is uncalibrated, so the result remains conditional Tier B.

## 7. Chen-window classification

Second-step growth is \(g_2=(G_2-G_1)/G_1\). Let \(D\) denote target attainment and \(B\) denote \(g_2\le g_{tol}\). The mutually exclusive classes are success (\(D\land B\)), grain-growth failure (\(D\land\neg B\)), densification exhaustion (\(\neg D\land B\)), and mixed failure (\(\neg D\land\neg B\)). First-step failure and target-already-reached states are separate ineligible classes.

A complete practical window requires lower exhaustion and upper growth brackets, \(T_2<T_1\), and

\[
W=T_{last\ success}-T_{first\ success}\ge25\,^\circ{\rm C}. \tag{22}
\]

The adaptive map uses a coarse grid, downward/upward boundary extensions, and 10 °C refinement within changing 25 °C intervals.

## 8. Fast-firing matched-density criterion

Interpolation is restricted to jointly attained density:

\[
R_{fast}(\rho)={G_{ref}(\rho)\over G_{fast}(\rho)}. \tag{23}
\]

A pass requires \(R_{fast}\ge1.5\) continuously over \(\Delta\rho\ge0.03\), with both paths attained and uncensored. Nucleation-facile replaces \(\tau_{\rm nuc}\) by \(10^{-6}(\tau_{\rm ex}+\tau_{\rm tr})\); PR-off sets redistribution to zero. Final causal scoring requires the full case to pass and nucleation-facile not to pass.

## 9. Relative material-property attribution

\[
Q_x'=Q_x+10^3\Delta Q_x,\qquad k_x'=f_xk_x,\quad
f_x\in\{0.03,0.1,0.3,1,3,10,30\}. \tag{24}
\]

The bounded design includes OAT rows, eight pair maps, 50,000 Latin-hypercube rows, and diagnostic-only geometry/accommodation rows. Reported groups include

\[
\Theta_{nuc}={\tau_{nuc}\over\tau_{ex}+\tau_{tr}},\qquad
S_{closed/growth}={k_{closed}\theta_{closed}(T_2)\over k_g\theta_g(T_2)}. \tag{25}
\]

These are relative diagnostics, not universal thresholds.

## 10. Numerical guardrails and exact promotion

The 50,655-row algebraic screen was used only for promotion. It predicted 19,880 both-pass rows. Final evidence uses 1,000 exact fast and 1,000 exact two-step promotions: 1,903 unique rows comprising 485 fast-only, 119 two-step-only, 73 both-pass, and 1,226 neither.

Time steps are bounded by temperature, density, grain growth, and store-loss increments. Exponential arguments and bounded states are clipped as recorded in the SI. Unattained or censored intervals are not scored. Topology was frozen during the relative-material perturbation campaign. Candidate 693168 is conditional Tier B, not validation or paper-ready calibration.
