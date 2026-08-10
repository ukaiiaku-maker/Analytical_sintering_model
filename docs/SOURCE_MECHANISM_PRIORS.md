# Source-grounded mechanism priors

## Sources and modeling boundary

This registry is grounded in `Nucleation_Limited_Sintering_Article9.docx`,
`GrainGrowth_V6.docx`, and the repository's preceding mechanism audits.  The
documents are used as physical priors, not as parameter datasets.  No Y2O3 or
other paper-specific values are fitted here.

## Densification priors

The sintering theory defines densification as interface-normal strain coupled
to pore-volume removal.  Each finite strain quota requires three serial
ingredients: renewal/nucleation of a climb-mediating defect, local
point-defect exchange, and transport between complementary sources and sinks.
Their characteristic cycle is `tau_nuc + tau_exchange + tau_transport`.
Geometrical eligibility (`chi_max`) is distinct from instantaneous renewal
activity.  A boundary must be connected through an admissible transport path
to a pore to remove pore volume.  Pore-disconnected boundaries can migrate,
creep, or relax stress but do not ordinarily densify.

The local activation stress is an evolving internal microstructural state.
Coarsening may continue while densification is arrested and concentrate local
stress until a renewal event occurs.  Exchange resistance and transport
resistance are necessary serial processes; kinetic rate control and energetic
dissipation need not be assigned to the same process.  These principles rule
out a single hidden efficiency factor and require separately reported renewal,
exchange, transport, stress, and power channels.

## Grain-growth priors

The grain-growth theory writes the intrinsic capillarity law as an activity-
modified migration law and organizes closures as:

- **Class A — series-resistance drag:** `Gamma = 1/(1+D)`. Pore, particle,
  solute, and TJ drag dissipate continuously and add as resistances.
- **Class B — hazard-controlled enabling events:** a packet/path step requires
  `K` events from an expected count `Lambda`, giving
  `Gamma = Pr[Poisson(Lambda) >= K]`.  TJ reactions and vacancy accommodation
  are admissible examples.  The source leaves required-count scaling
  unresolved, so both fixed packet (`q_TJ=0`) and accommodation demand
  (`q_TJ=1`) must remain visible ablations.
- **Class C — exchange-limited crossover:** a relaxation timescale produces a
  size-dependent crossover such as `Gamma=1/(1+tau_exchange/tau_available)`.
- **Class D — coupled mechanisms:** necessary serial mechanisms combine by
  reciprocal activities; parallel enabling channels combine by survival
  probabilities.

`Lambda/K` distinguishes smooth (`>>1`), intermittent (`O(1)`), and stagnant
(`<<1`) regimes.  A multihit closure must report `Lambda`, `K`, completion
probability, and its packet/path definition; a constant growth multiplier is
not an acceptable substitute.

## Priors from negative controls

The aggregate, junction-limited, pore-pinning, pore-placement, and constrained
action audits establish the following regression facts:

1. baseline fractional growth scales too strongly at nanoscale;
2. connected pore/TJ drag releases as coverage depletes;
3. explicit location is necessary but does not itself open a nanoscale window;
4. TJ-to-GB capture is observable but 2–3 orders below pore removal;
5. doubled TJ drag reduces growth but still does not overlap density and
   no-growth boundaries;
6. isolated pore volume must not be removed by the open-pore densification
   channel.

## Selected discovery hypotheses

The first bounded discovery iteration selects two migration-only mechanisms:

1. **Persistent junction population (`X_J`, Class A):** local TJ reactions,
   relocation/capture, and boundary sweep create a bounded junction obstacle
   population; thermal relaxation removes it.  Its drag can survive loss of
   connected pore coverage and its unrealized migration power is reported.
2. **TJ multihit reaction (Class B):** migration requires completion of a TJ
   reaction packet.  Fixed-packet and accommodation-demand limits are carried
   separately.  The closure changes migration only and reports
   `Lambda_TJ`, `K_TJ`, `Lambda/K`, completion probability, and regime.

Stress accumulation/release and vacancy-accommodation multihit remain in the
registry but are deliberately deferred.  Adding them simultaneously would
make a success mechanistically non-identifiable.  Closed-pore removal remains
an explicit placeholder and is not used at the 0.90 target.
