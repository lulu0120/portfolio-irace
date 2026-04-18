## portfolio-irace reproducibility package

This folder contains the code and benchmark data needed for the paper repository.

Included:
- `QAP_problem/NewInstances/`: QAP benchmark instances used in the experiments.
- `QAP_problem/Prog/`: QAP solver source, Makefile, and the compiled `lsmcqap` binary.
- `QAP_problem/Scripts/`: core QAP portfolio-irace runner scripts and the new-instances split generator.
- `QAP-hydra/`: core Hydra runner scripts for QAP.
- `Multiple_multidim/`: core `P-irace` implementation used by the QAP portfolio scripts.
- `Multiple_2combine/performance_eval.R`: shared evaluation utilities sourced by the QAP interface.
- `benchmark/hydra-irace/hydra_wrapper.R`: Hydra wrapper code used by the QAP Hydra interface.

Intentionally excluded:
- result files (`.csv`, `.Rds`, `.RData`)
- plots and derived analysis outputs
- temporary logs, notebooks, and editor history files
- old split outputs under `QAP_problem/Scripts/splits/`

Notes:
- The included QAP scripts use relative repository structure rooted at this folder.
- If fresh train/test splits are needed for `NewInstances`, run `QAP_problem/Scripts/newinstances_generate_split.R`.
- If the solver binary is not executable on the target machine, rebuild it from `QAP_problem/Prog/Makefile`.
