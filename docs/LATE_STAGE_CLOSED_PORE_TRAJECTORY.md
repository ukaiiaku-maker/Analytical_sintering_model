# Late-Stage Closed-Pore Trajectory Audit

## Model and conservation

The extension introduces a distinct `phi_closed_i`/`N_closed_i` store. Pore
isolation, closure, and optional GB/TJ detachment transfer volume
conservatively; they do not change density. Only the named vacancy-transport
closed-pore shrinkage flux removes `phi_closed` and increases density. The
open-pore removal law continues to act only on GB-segment and TJ stores.

The shrinkage drive is capillary pressure minus optional trapped-gas pressure,
plus an ablatable hydrostatic residual-stress term. Gas content is stored
independently and is not consumed by shrinkage. Disabled mode returns the
original open-pore history exactly. The identity

`rho = 1 - sum(phi_GBseg + phi_TJ + phi_iso + phi_closed)`

holds to numerical precision, and all stores remain nonnegative.

## Bounded screen

The 64-case physical design was run uniformly at 96, 240, and 500 h (192
budgeted cases). It includes disabled, vacancy-transport, gas-limited,
detachment, and combined modes. Reference selection is explicit: 0.2 C/min is
used only when attainable; otherwise the recorded 1 C/min fallback is used.

Five fast paths reach density 0.95 at each uniform budget. No fast path reaches
0.98 or 0.99. Longer nominal budgets do not improve these fixed ramp/hold
schedules once their declared hold ends. Thus the module permits limited,
honest high-density sampling but not final-density mapping.

No trajectory is meaningful. Forty-nine cases jointly attain 0.85--0.92, but
none sustains a 1.5 ratio over `Delta rho >= 0.03`. Thirteen jointly attain
0.90--0.95; their maximum ratio is only 1.110. No case jointly attains the
0.95--0.98 or 0.98--0.99 windows. Raw maxima of 2.02 (combined) and 1.58
(gas-limited) occur in short or unattained ranges around density 0.91 and are
rejected.

## Interpretation

1. **Honest sampling at density >=0.95?** Yes, narrowly: five fast paths reach
   0.95. No comparison supports 0.98 or 0.99.
2. **Meaningful fast-firing `G_mean` trajectory?** No.
3. **Early/intermediate or late?** The largest raw signals remain near the
   open/closure transition. The attained 0.90--0.95 window is weak.
4. **Dependence on gas, transport, detachment, stress?** Gas-limited and
   combined modes produce the largest transient ratios; detachment gives the
   highest attained densities (about 0.972 reference, 0.969 fast). Neither
   yields a qualifying trajectory.
5. **Chen preservation?** Disabled mode exactly recovers the established 0.90
   open-pore Chen maps. No late-stage candidate passed the trajectory gate, so
   the adaptive Chen solver was not escalated for active candidates. This
   avoids presenting a nonqualifying parameter set as jointly validated.
6. **Why does the model still fail?** Closed-pore transport improves density
   attainment but does not create enough additional differential grain growth
   across a jointly attained finite interval. Gas limitation can amplify a
   transient response while simultaneously reducing shared attainment.
7. **Calibration measurements?** In-situ or interrupted high-density pore-size
   tomography, closed-pore number/radius distributions, trapped-gas pressure or
   composition, GB detachment statistics, shrinkage rates, and matched-density
   grain distributions at 0.92--0.99 would directly constrain closure and
   vacancy-transport parameters.

Early-stage and late-stage reduced mechanisms tested so far both fail to
produce a sustained experimentally meaningful mean-grain trajectory
separation, although they produce interpretable internal microstructural
memory.

