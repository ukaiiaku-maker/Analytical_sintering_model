# Codex Task: Initial Architecture Work in `INRL_lambert_onsager`

## Mission

Develop a physically interpretable reduced-order analytical sintering model that can reproduce, in one formulation:

1. fast-heating-rate sintering, where rapid heating gives more density gain per grain-growth increment; and
2. two-step sintering, where a short high-temperature step followed by lower-temperature holding beats a long high-temperature path at matched density.

Do not aim for a one-off parameter set. Build a framework in which the coupling among mechanisms is explicit and searchable.

## Local workspace

The user intends to work in:

```text
INRL_lambert_onsager
```

Use the GitHub repository:

```text
https://github.com/ukaiiaku-maker/Analytical_sintering_model.git
```

Recommended setup:

```bash
mkdir -p INRL_lambert_onsager
cd INRL_lambert_onsager

git clone https://github.com/ukaiiaku-maker/Analytical_sintering_model.git
cd Analytical_sintering_model

python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
python3 -m py_compile *.py
```

Create a working branch:

```bash
git switch -c codex/topology-constrained-mechanisms
```

## First priority

Read:

```text
README.md
docs/CODEX_HANDOFF.md
docs/MECHANISM_OPTIONS_AND_GOALS.md
```

Then inspect the current code:

```text
sinter_reference_model_v3_multibin.py
sinter_reference_model_v6_grainstress_multibin.py
sweep_lambda_window_priority.py
sweep_lambda_window_priority_v4.py
debug_lambda_v6_grainstress.py
```

## Immediate scientific problem

Previous attempts showed that optimizing `Lambda = r_nuc * tau_sink` alone is insufficient. Coarse grains can have high Lambda because `tau_sink` is large, while densification efficiency is poor.

The model must distinguish:

```text
renewal activity
transport/exchange completion rate
density gain per event
grain growth per event
pore topology
energy dissipation
```

Use a combined metric, not raw grain-size difference alone:

```text
HR_pct = 100*(G_slow - G_fast)/G_slow
TS_pct = 100*(G_highT - G_twoStep)/G_highT
E_G    = d rho / d ln G
```

## Architecture target

Create or refactor toward the following components.

### State

```python
@dataclass
class State:
    rho: float
    G: float
    pore_radii: np.ndarray
    pore_phi: np.ndarray
    pore_N: np.ndarray
    topology: TopologyState
    stress: StressState
```

### Topology

```python
@dataclass
class TopologyState:
    f_pore: float
    f_clean: float
    f_PR: float
    f_TL: float
    connectivity: float
    isolated_pore_fraction: float
```

### Mechanism flux

```python
@dataclass
class MechanismFlux:
    rho_dot: float
    G_dot: float
    pore_phi_dot: np.ndarray
    pore_N_dot: np.ndarray
    power: float
    diagnostics: dict
```

### Mechanisms

Implement or stub:

```text
renewal_densification
clean_GB_migration
pore_connected_GB_migration
pore_drag
triple_line_drag
PR_desintering
pore_coarsening
exchange_dissipation
transport_dissipation
```

### Partition

Add:

```python
def solve_dissipation_partition(state, topology, mechanisms, params):
    ...
```

This should enforce nonnegative mechanism weights and conservation constraints.

## Required diagnostics

Every run should report:

```text
rho(t)
G(t)
pore-bin distribution
f_pore, f_clean, f_PR, f_TL
sigma_base, sigma_concentration, sigma_local
r_nuc, tau_exchange, tau_transport, tau_TL
activity
rho_dot
dGdt
E_G = rho_dot/(dGdt/G + eps)
power channels
```

## First deterministic tests

Add or preserve scripts that run:

1. heating-rate comparison:
   ```text
   0.2 C/min vs 20 C/min to the same target density
   ```

2. two-step comparison:
   ```text
   high-T isothermal vs T1 -> T2 to the same target density
   ```

3. second-step activity/efficiency window:
   ```text
   start directly from rho0 = 0.83 and sweep G0 and T
   ```

4. topology sanity:
   ```text
   same pore volume but larger pore sizes gives lower pore-boundary coverage
   ```

5. pore conservation:
   ```text
   rho = 1 - sum(phi_i)
   phi_i >= 0
   N_i >= 0
   ```

## Acceptance criteria for the first Codex run

The first run is successful if it produces:

1. a clean architecture or migration plan;
2. at least one refactored model or new module;
3. deterministic diagnostic plots;
4. a small search using percent metrics;
5. a short report describing what works, what fails, and what to try next;
6. no hidden parameter hacks or unexplained scalar efficiency multipliers;
7. all code compiles.

## Do not do this

Do not simply run a huge brute-force search on the current code without improving the mechanism structure.

Do not optimize raw nanometer grain-size difference.

Do not hide failed densification by changing the target density after the fact.

Do not accept a case where both protocols fail to reach target density.

Do not use one scalar `eta_total` as the whole theory.

## Expected deliverables

Commit with message similar to:

```text
feat: add topology-constrained mechanism architecture prototype
```

Include:

```text
docs/ARCHITECTURE_NOTES.md
docs/SEARCH_RESULTS_INITIAL.md
figures or output summaries in an ignored or small-results folder
updated code
updated README quick-start if needed
```

Before stopping, print:

```bash
git status --short
git log --oneline -5
```
