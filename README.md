# portfolio-irace-scripts

Scripts for `portfolio-irace` project for synthetic and QAP example portfolio construction and evaluation.

This snapshot keeps only the files needed to:

- run the QAP solver as a black-box evaluator,
- connect the evaluator to the portfolio-irace workflow,
- generate train/test splits,
- evaluate learned portfolios on the QAP test set,
- reproduce the Monte Carlo portfolio baseline calculations used in the 1D plotting code.

## Directory Layout

- `QAP_problem/Prog`
  - C source for the QAP solver.
- `QAP_problem/Scripts`
  - QAP runner, portfolio interface, split generation, and evaluation scripts.
- `QAP_problem/instances`
  - QAP instance files used by the scripts.
- `Multiple_multidim`
  - Core portfolio-irace race logic, contain the synthetic example in run_overall 
  - Also the basic files imported by the QAP interface.
- `Multiple_2combine/performance_eval.R`
  - Monte Carlo helper used by the baseline cost calculations (example scoring and evaluation scripts).

## Included Entry Points

- `QAP_problem/Scripts/generate_instance_split.R`
  - Generate a stratified train/test split.
- `QAP_problem/Scripts/qap_runner.R`
  - Evaluate one QAP configuration on one instance.
- `QAP_problem/Scripts/qap_portfolio_run.R`
  - Run one portfolio-irace experiment.
- `QAP_problem/Scripts/evaluate_qap_config_testset.R`
  - Evaluate configurations or portfolios on the QAP test set.
- `QAP_problem/Scripts/evaluate_qap_result_portfolios_test_instances.R`
  - Evaluate portfolio result CSVs on the QAP test instances.

 
## Setup Notes

This copy does not include compiled binaries. Before running the QAP scripts, build the solver binary so that:

- `QAP_problem/Prog/lsmcqap`

exists and is executable.

The R scripts expect the repository layout in this snapshot to remain unchanged.

## Basic Usage

Generate a split:

```bash
Rscript QAP_problem/Scripts/generate_instance_split.R --seed n (user defined)
```

Run one solver call example:

```bash
Rscript QAP_problem/Scripts/qap_runner.R \
  --instance QAP_problem/instances/100/GridRandom/GridRandom.974819759.n100.K20.m10.A100.00.B1.00.sp0.00.dat \
  --a 2 --p 3 --l 1 --t 5 --seed 1
```

Run one portfolio-irace experiment example:

```bash
Rscript QAP_problem/Scripts/generate_instance_split.R --seed n

Rscript QAP_problem/Scripts/qap_portfolio_run.R \
  --train-split QAP_problem/Scripts/splits/seed_11/train_instances.txt \
  --budget 1000 \
  --port-size 2 \
  --best-k 1 \
  --time-limit 5
```
 