# Codex Initial Instructions: Unified Analytical Sintering Model

## Working context

Repository:

```text
https://github.com/ukaiiaku-maker/Analytical_sintering_model.git
```

Local working folder requested by the user:

```text
INRL_lambert_onsager
```

This repository contains reduced-order analytical sintering models developed to search for a common physical formulation that can reproduce both:

1. **Fast-heating-rate sintering**: faster heating gives more densification per grain-growth increment than slow heating.
2. **Two-step sintering**: a short high-temperature step followed by a lower-temperature hold gives more efficient density-vs-grain-size evolution than a long high-temperature path.

The current code is exploratory. Do not assume the existing model is correct. Treat it as a compact testbed for mechanism discovery.

## Immediate local setup

From the parent folder that should contain the project:

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

Then run the current baseline diagnostics:

```bash
python3 debug_lambda_v6_grainstress.py
```

Then run a small search:

```bash
python3 sweep_lambda_window_priority_v4.py \
  --model sinter_reference_model_v6_grainstress_multibin \
  --n 500 \
  --rho-target 0.92 \
  --outdir sweep_lambda_window_priority_v6_test
```

## Important instruction

Do **not** begin by blindly expanding the parameter search. The previous searches repeatedly found accidental parameter combinations because the coupling among densification, grain growth, pore topology, pore drag, P-R/de-sintering, and stress generation was not structured enough.

Your first task is to improve the model architecture so the search space is physically interpretable.

## Scientific target

Build a common reduced-order analytical formulation that can satisfy these conditions simultaneously:

```text
Condition A: The model densifies substantially under normal high-temperature sintering schedules.

Condition B: At matched density, fast heating gives smaller grain size than slow heating:
             G_20Cmin(rho*) < G_0.2Cmin(rho*)

Condition C: At matched density, two-step sintering gives smaller grain size than high-T isothermal:
             G_twoStep(rho*) < G_highT(rho*)

Condition D: The successful response comes from interpretable mechanisms:
             renewal-limited densification, topology-dependent eligible boundaries,
             pore-size-distribution evolution, stress generation, and competing dissipation.

Condition E: The model survives held-out schedules not used during search.
```

Use percent metrics rather than raw nanometer differences:

```text
HR_pct = 100 * (G_slow - G_fast) / G_slow
TS_pct = 100 * (G_highT - G_twoStep) / G_highT
```

A useful candidate should eventually achieve at least:

```text
HR_pct >= 5%
TS_pct >= 5%
```

while reaching the target density without unphysical runaway coarsening.

## Known failure mode

Do not optimize `Lambda = r_nuc * tau_sink` alone.

A large Lambda does not necessarily mean efficient densification. Coarse grains can have large `tau_sink`, and therefore large Lambda, while still densifying inefficiently. A good activity-window criterion should include both renewal activity and density gain per grain-growth increment.

Add or use a trajectory-efficiency diagnostic such as:

```text
E_G = d rho / d ln(G)
```

or numerically:

```text
window_efficiency = median(activity * rho_dot / (dGdt/G + eps))
```

A good model should keep the two-step path in a high-efficiency window, not merely in a high-activity window.

## Preferred first milestone

Implement a cleaner mechanism/state architecture before running a very large search.

The first milestone should produce:

1. a refactored or newly added model file with explicit state variables;
2. deterministic diagnostics showing mechanism contributions;
3. unit or smoke tests proving pore volume conservation and finite outputs;
4. a short report summarizing whether the architecture can reproduce the target behaviors;
5. a small search result using percent metrics and held-out checks.

## Required coding style

Use simple, auditable Python. Prefer:

- `dataclasses` for parameters and state;
- pure functions for mechanisms;
- explicit arrays for pore-size bins;
- CSV outputs for search summaries;
- PNG outputs for diagnostic plots;
- deterministic seeds for searches;
- no hidden global state;
- no large binary files committed to git.

## Preserve these principles

1. Densification is an interface-normal strain process mediated by source/sink renewal.
2. Not all boundaries are geometrically eligible for densification.
3. Pore topology determines which boundaries remove pore volume and which only migrate/relax.
4. Pore-size distribution matters: small pores close faster; large pores are harder to remove and may drag/coalesce.
5. Grain growth, clean GB migration, pore drag, P-R/de-sintering, local exchange, and transport may dissipate energy without equal density gain.
6. Stress generation must be separated into baseline capillary stress and microstructure-generated local concentration stress.
7. Searches should be used to compare mechanistic couplings, not to invent unphysical efficiency factors.

## Suggested first branch name

```bash
git switch -c codex/topology-constrained-mechanism-search
```

Commit milestones as small auditable steps.
