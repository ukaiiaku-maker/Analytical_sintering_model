# Diagnostic energy-ledger equations

This ledger is an audit of existing trajectories; it does not replace the accepted
forward calculation and is not validation.

Per bulk volume, the stored-energy terms are

\[
E_{surf}=\gamma_s A_{surf},\quad E_{GB}=\gamma_{GB}C_{GB}/G,\quad
E_{pore}=\gamma_s\sum_i3\phi_i/r_i,\quad E_{gas}=\sum_iP_{gas,i}\phi_{c,i}.
\]

The external-area term uses an explicitly provisional compact-scale geometry.
Exported histories do not contain pore-bin arrays, so their pore area is reconstructed
from total pore fraction and D50. Candidate-law tests remain bin-resolved.

Signed rates are \(P_j=-dE_j/dt\), except gas cost \(+dE_{gas}/dt\), and

\[
P_{available}=[P_{surf}]_+ +[P_{GB}]_+ +[P_{pore}]_+-P_{gas,cost}.
\]

Named expenditures are open and closed densification, surface smoothing,
conservative pore redistribution, GB growth, drag, gas work, and an explicit
unresolved channel. The audit records rather than repairs violations of
\(\sum P_{spent}\le P_{available}+\epsilon\). It also reports
\(\Pi_{dens}=P_{dens}/P_{available}\), \(\Pi_{total}=P_{spent}/P_{available}\),
and \(P_{excess}=\max(P_{available}-P_{dens},0)\).

The old equality \(P_{dens}=\gamma_{GB}C_{GB}\dot G/G^2\) is retained only as
`strict_GB_area_loss_balance`, a diagnostic ablation. It is too narrow because it
omits pore/external surface release, gas work, smoothing, redistribution, and drag,
and it feeds the migration law back into activation stress.
