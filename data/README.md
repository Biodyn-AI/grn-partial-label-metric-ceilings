# Data Directory

This directory contains the benchmark evaluation artifacts used as inputs for the partial-label metric ceiling analysis.

## Metric tables

- `score_eval_probe_priors.csv`: Primary evaluation table with F1, AUPR, precision, recall for probe-based methods against regulatory references (TRRUST, DoRothEA, DoRothEA-TRRUST immune union).
- `score_eval_grn_baselines_immune.csv`: Classical GRN inference baselines (GENIE3, GRNBoost2, SCENIC, etc.) evaluated on immune-specific references.
- `score_eval_probe_priors_full_genes.csv`: Extended evaluation with full gene universe (F1 only).
- `score_eval_probe_priors_full_genes_crosswalk.csv`: Extended evaluation with crosswalk gene mapping.
- `score_eval_probe_priors_full_genes_omnipath.csv`: Extended evaluation with OmniPath-expanded references.

## Coverage proxy reports

- `score_eval_probe_priors_missing_report.json`: Gene-mapping coverage report for primary evaluation.
- `score_eval_probe_priors_full_genes_missing_report.json`: Coverage report for full-gene evaluation.
- `score_eval_probe_priors_full_genes_crosswalk_missing_report.json`: Coverage report for crosswalk evaluation.
- `score_eval_probe_priors_full_genes_omnipath_missing_report.json`: Coverage report for OmniPath evaluation.

## Sweep and cross-evaluation artifacts

- `sweep_results_*.json`: Threshold sweep results (top-k and percentile) for various reference configurations.
- `cross_eval_*.json`: Cross-evaluation AUPR and Brier scores between reference families (HPN, BEELINE).

## External validation

- `summary_timeseries_metrics.csv`: Time-series evaluation metrics across methods (HPN, BEELINE AUPR).
- `perturbation_eval_probe.json`: Perturbation validation with per-perturbagen recall statistics.
