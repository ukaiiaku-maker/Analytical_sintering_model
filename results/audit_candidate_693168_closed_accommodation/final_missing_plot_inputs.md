# Missing or proxy-only plot inputs

The frozen local-region model does not expose the following requested quantities as independent physical channels:

- physical pore `D50_nm` and `D90_nm`; figures use clearly labeled volume/number radius proxies;
- `tau_nuc_s`, `tau_exchange_s`, and `tau_transport_s` inside the local-region closure;
- `sigma_base` and `sigma_act_total`; only the dimensionless residual-stress state is exposed;
- absolute `H_PR`, `H_dens`, `w_PR`, and `w_dens`; plots use named PR redistribution rate and bounded cumulative pore-volume transfer/removal;
- an independent closed-pore `Lambda_closed`, `K_closed`, or Poisson completion probability; the accommodation factor is shown as a proxy and never labeled as a fitted event count;
- `P_persistent_junction_drag` and `P_clean_GB` as separately normalized power channels.

No values were fabricated. Corresponding dense-history columns are `NaN` where requested for schema visibility. Figures requiring unavailable absolute energies or diameters were reformulated around exposed state/flux variables rather than emitted as placeholders.
