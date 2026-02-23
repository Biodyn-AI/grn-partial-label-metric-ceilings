# Partial-Label Metric Ceilings for GRN Evaluation

Companion code and data for the paper: **"Partial-Label Metric Ceilings for GRN Evaluation in Single-Cell Foundation Models"**.

## Overview

Gene regulatory network (GRN) benchmarks are typically interpreted as if curated references were complete. This repository provides a theory-first framework for interpreting observed F1 and AUPR when only a fraction of true edges are labeled, and tests the framework on existing GRN benchmark outputs across 15 methods and 5 references.

Under a missing-at-random (MAR) positive-label model, we derive explicit ceilings for observable F1 and AUPR as functions of label coverage, propagate uncertainty with Beta posteriors, and stress-test assumption violations with adversarial non-MAR simulations.

## Key findings

- Across 39 AUPR-evaluable rows, best ceiling-normalized F1 = 0.137 and AUPR = 0.014; median normalized AUPR = 5.56e-05.
- Only 3/39 rows exceed the observed random baseline after ceiling correction.
- Non-MAR curation shifts AUPR by +0.090 (optimistic) to -0.073 (pessimistic), making curation mechanism as consequential as curation magnitude.
- External perturbation validation recovers signal in only 19/4787 perturbations (0.4%).

## Repository structure

```
.
├── README.md                 # This file
├── LICENSE                   # MIT License
├── requirements.txt          # Python dependencies
├── research_note.md          # Formal derivations and mathematical setup
├── data/                     # Input benchmark data
│   ├── score_eval_*.csv      # Metric tables (F1, AUPR, precision, recall)
│   ├── *_missing_report.json # Coverage proxy reports
│   ├── sweep_results_*.json  # Threshold sweep artifacts
│   ├── cross_eval_*.json     # Cross-evaluation artifacts
│   ├── summary_timeseries_metrics.csv
│   └── perturbation_eval_probe.json
├── scripts/                  # Reproducible analysis pipeline
│   ├── run_analysis.py       # Step 1: Initial ceiling calculations
│   ├── run_paper_pipeline.py # Step 2: Full multi-scenario analysis
│   ├── build_final_paper_markdown.py  # Step 3: Populate manuscript
│   └── run_submission_readiness_analyses.py  # Step 4: Bootstrap + external validation
├── outputs/                  # Pre-computed analysis outputs
│   ├── *.csv                 # Initial ceiling reinterpretation tables
│   ├── figures/              # Initial analysis figures
│   └── paper/
│       ├── tables/           # All supplementary CSV tables (25 files)
│       ├── figures/          # Publication figures (10 panels)
│       └── results_summary.md
```

## Reproduction

### Requirements

Python >= 3.9 with:

```bash
pip install -r requirements.txt
```

### Run the full pipeline

From the repository root:

```bash
# Step 1: Initial ceiling analysis
python scripts/run_analysis.py

# Step 2: Full paper pipeline (scenarios, simulations, figures, tables)
python scripts/run_paper_pipeline.py

# Step 3: Build populated manuscript from pipeline outputs
python scripts/build_final_paper_markdown.py

# Step 4: Bootstrap uncertainty, external validation, submission diagnostics
python scripts/run_submission_readiness_analyses.py
```

All outputs are written to `outputs/`. Pre-computed outputs are included for convenience.

```

## Data description

### Input data (`data/`)

| File pattern | Description |
|---|---|
| `score_eval_probe_priors.csv` | Primary metric table: probe-based methods vs regulatory references |
| `score_eval_grn_baselines_immune.csv` | Classical GRN baselines evaluated on immune references |
| `score_eval_probe_priors_full_genes*.csv` | Extended evaluation tables (F1 only, larger gene universes) |
| `*_missing_report.json` | Gene-mapping coverage reports per evaluation policy |
| `sweep_results_*.json` | Threshold sweep artifacts for external sanity checks |
| `cross_eval_*.json` | Cross-evaluation AUPR/Brier between reference families |
| `summary_timeseries_metrics.csv` | External time-series validation metrics |
| `perturbation_eval_probe.json` | Perturbation validation recall per perturbagen |

### Output tables (`outputs/paper/tables/`)

Key tables include:
- `central_aupr_rows.csv` / `central_f1_rows.csv`: Row-level ceiling reinterpretation
- `reference_summary.csv` / `method_summary.csv`: Aggregate summaries
- `coverage_scenarios_by_reference.csv`: Multi-scenario coverage definitions
- `scenario_method_ranks.csv` / `scenario_rank_stability.csv`: Robustness diagnostics
- `synthetic_mar_summary.csv` / `synthetic_adversarial_summary.csv`: Simulation results
- `bootstrap_global_metrics.csv`: Bootstrap uncertainty for headline statistics
- `family_significance_tests.csv`: Pairwise method-family comparisons

### Figures (`outputs/paper/figures/`)

10 figures covering theoretical curves, empirical reinterpretation, ranking shifts, scenario robustness, synthetic validation, and external validation.

## Mathematical framework

The core derivations are in `research_note.md`. Key results:

- **F1 ceiling:** `F1_max(c) = 2c / (1 + c)`
- **AUPR ceiling:** `AUPR_max(c, t) ≈ c + (1-c) * H_t / t`  (large-t limit: `≈ c`)
- **Coverage uncertainty:** `c ~ Beta(k+1, n-k+1)` with Monte Carlo propagation

where `c` is positive-label coverage, `t` is latent positive count, and `H_t` is the harmonic number.

## License

MIT License. See [LICENSE](LICENSE).
