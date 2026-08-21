# Reduced-to-physical property mapping

The earlier reduced-model “works” envelope is used only as a diagnostic target. No
candidate-693168 parameter is copied into the ZrO2 model.

| Reduced quantity | Physical interpretation in this audit | Status |
|---|---|---|
| `Q_nuc` | Local effective slope arising from full \(G^*(\sigma,T)\) | diagnostic; not a scalar input |
| `Q_transport` | Fixed \(Q_{GB}=380\) kJ/mol | unchanged physical input |
| `Q_PR` | Fixed \(Q_s=380\) kJ/mol plus radius/topology factors | incomplete geometry mapping |
| `Q_growth` | Current fixed mobility activation, about 405.24 kJ/mol | unchanged; global failed fit inactive |
| `Q_closed` | Post-run apparent slope after inventory/activity effects | not a physical input |
| closed prefactor | inventory × shrinkability × accommodation × dimensional geometry | calibration gap |
| PR prefactor | surface transport × radius moment × topology coefficient | calibration gap |
| lower boundary | finite closed shrinkage at fixed prepared inventory | not preserved as success topology |
| upper boundary | mobility and evolving pore/Zener pinning | growth failure appears, but no success interval |

Across the broad fixed-input rate scan, the GB and gas-limited candidates overlap a
finite-gain diagnostic band for roughly 55% of sampled points, and renewal candidates
for roughly 39–42%. This does not establish a usable processing window: the scan
varies inventory and accommodation independently over broad brackets.

At the actual selected first-step state, \(\phi_c=6.59\times10^{-6}\) and
\(A_c=2.82\times10^{-3}\). Their product caps possible shrinkage far below the
density target regardless of transport rate. The immediate gap is therefore closed
inventory and accommodation/availability. A dimensional GB geometry prefactor,
gas inventory, radius exponent, renewal activity, and transport length remain
secondary unresolved quantities. They must not be fit until the state mapping is
measured or independently constrained.

`Q_closed_app` values in the output are slopes of calculated rates, not material
parameters. The barrier-derived effective slope is also extrapolative because the
fit begins above the processing range.

This comparison is **not validation**. No validation claim is made.
