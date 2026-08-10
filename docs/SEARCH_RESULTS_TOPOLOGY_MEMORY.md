# Topology-memory mechanism result

## Mechanism

`topology_damage` is an accumulated state for non-densifying surface smoothing
during a Gaussian intermediate-temperature window. Its rate is proportional to
low renewal activity, a pre-densification density gate, and the remaining
undamaged fraction. It is bounded between zero and one.

Damage acts only through named topology channels: it lowers removable
pore-boundary coverage and increases isolated-pore fraction. Consequently it
reduces density gain per event and `E_G`; it is not an `eta_total`. The state,
rate, and resulting topology changes are emitted at every integration step and
the mechanism is disabled with `enable_topology_memory=False`.

## Deterministic result at rho = 0.90

| Case | HR_pct | TS_pct | Slow damage | Fast damage | All reached |
|---|---:|---:|---:|---:|---|
| Baseline prototype / disabled | -6.4503 | 11.5345 | 0 | 0 | yes |
| Topology memory | 38.0887 | 11.2949 | 0.8481 | 0.0174 | yes |

The ablation therefore restores the exact prior negative heating-rate result,
while the enabled state flips its sign and preserves the two-step advantage.
No target density, protocol, stress cap, or unrelated kinetic parameter was
changed. Exact generated values are in
`results/initial/topology_memory_ablation.csv`.

## Interpretation and limitation

The result supports the missing-memory diagnosis: a slow ramp spends enough
time in the surface-smoothing window to lose removable pore topology before
densification becomes efficient. The magnitude of the predicted heating-rate
advantage is not calibrated and should not yet be treated quantitatively.
Next work should validate the temperature window and damage rate against pore
size or surface-area observations before any expanded search.
