# GRN Partial-Label Metric Ceilings

Reproducibility artifacts for the paper:

**Partial-Label Metric Ceilings for GRN Evaluation in Single-Cell Foundation Models**

This repository contains:
- manuscript files (`paper/`)
- analysis and figure-generation scripts (`scripts/`)
- generated tables and figures used in the manuscript (`outputs/`)

## Reproducing main results

1. Prepare a Python environment with:
   - `python>=3.10`
   - `numpy`
   - `pandas`
   - `matplotlib`

2. Run the pipeline scripts in order:

```bash
python scripts/run_proposal1_analysis.py
python scripts/run_proposal1_paper_pipeline.py
python scripts/run_submission_readiness_analyses.py
```

3. Build the finalized markdown manuscript:

```bash
python scripts/build_final_paper_markdown.py
```

Main paper file:
- `paper/proposal1_partial_label_metric_ceilings_SUBMISSION_READY.md`

Main artifact directories:
- `outputs/tables/`
- `outputs/figures/`

## Dataset access notes

This project reuses benchmark outputs generated in prior network-inference pipelines.
Please follow dataset and benchmark access instructions from the parent Biodyn-AI resources.
