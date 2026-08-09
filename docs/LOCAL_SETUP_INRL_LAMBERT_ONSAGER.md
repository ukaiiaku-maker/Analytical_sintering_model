# Local Setup: `INRL_lambert_onsager`

Use this file to initialize the local folder for Codex work.

## Clone

```bash
mkdir -p INRL_lambert_onsager
cd INRL_lambert_onsager

git clone https://github.com/ukaiiaku-maker/Analytical_sintering_model.git
cd Analytical_sintering_model
```

## Environment

```bash
python3 -m venv .venv
source .venv/bin/activate

python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

## Baseline checks

```bash
python3 -m py_compile *.py
python3 debug_lambda_v6_grainstress.py
```

## Start Codex branch

```bash
git switch -c codex/topology-constrained-mechanisms
```

## Give Codex these instructions

Start with:

```text
Read docs/CODEX_INITIAL_INSTRUCTIONS.md, docs/MECHANISM_OPTIONS_AND_GOALS.md, docs/CODEX_TASK_INRL_LAMBERT_ONSAGER.md, and docs/CODEX_HANDOFF.md.

Then inspect the code and propose a staged implementation plan. Do not run a large parameter search until you have made the mechanism coupling and topology state explicit.
```

## Expected first Codex command

```bash
python3 -m py_compile *.py
python3 debug_lambda_v6_grainstress.py
python3 sweep_lambda_window_priority_v4.py \
  --model sinter_reference_model_v6_grainstress_multibin \
  --n 200 \
  --rho-target 0.92 \
  --outdir smoke_lambda_window_priority_v6
```
