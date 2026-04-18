# QAP Scripts

This folder contains the minimal scripts needed to:

- generate train/test instance splits for repeated experiments,
- run the Thomas QAP solver as a standalone black-box evaluator,
- connect the evaluator to a plain irace run.

## Files

- `generate_instance_split.R`: stratified random split by instance family using a seed.
- `qap_runner.R`: standalone QAP runner. Prints one numeric cost.
- `qap_common.R`: shared solver-calling helpers used by the runner and portfolio interface.
- `target-runner`: irace adapter that calls `qap_runner.R`.
- `parameters.txt`: initial irace parameter space.
- `scenario.txt`: minimal single-run irace scenario.
- `qap_portfolio_interface.R`: glue layer that maps QAP instances and `(a,p,l)` into the existing portfolio-irace code.
- `qap_portfolio_run.R`: small CLI entry point for a QAP portfolio-irace run.

## Parameter Space

The current QAP parameter space is:

- `a in {0,1,2,3,4,5,6,7,8,9}`
- `p in {0,1,2,3,4,5,6,7,8,9,10,11}`
- `l in {0,1,2,3,4,5}`

## Split generation

Example:

```bash
Rscript generate_instance_split.R --seed 11
```

This writes:

- `splits/seed_11/train_instances.txt`
- `splits/seed_11/test_instances.txt`
- `splits/seed_11/split_manifest.csv`
- `splits/seed_11/split_summary.csv`

To use the generated split with irace:

```bash
cp splits/seed_11/train_instances.txt train_instances.txt
cp splits/seed_11/test_instances.txt test_instances.txt
```

## Standalone runner

Example:

```bash
Rscript qap_runner.R \
  --instance ../instances/100/GridRandom/GridRandom.974819759.n100.K20.m10.A100.00.B1.00.sp0.00.dat \
  --a 2 --p 3 --l 1 --t 5 --seed 1
```

## Self-test

Example:

```bash
Rscript qap_runner.R --self-test --n 5
```

## Portfolio-irace smoke test

Example:

```bash
Rscript qap_portfolio_run.R --budget 20 --port-size 1 --time-limit 3
```
