#!/usr/bin/env python3
"""Build a finalized manuscript markdown from Partial-label ceiling pipeline outputs.

This script reads generated tables under outputs/paper/tables and writes a
fully populated research paper draft with embedded statistics and tables.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def find_repo_root(start: Path) -> Path:
    """Locate repository root from script path."""
    for candidate in [start] + list(start.parents):
        if (candidate / "data").exists() and (candidate / "scripts").exists():
            return candidate
    raise RuntimeError("Could not infer repo root.")


def main() -> None:
    script = Path(__file__).resolve()
    repo_root = find_repo_root(script)

    tables_dir = repo_root / "outputs/paper/tables"
    figures_dir = repo_root / "outputs/paper/figures"

    central_aupr = pd.read_csv(tables_dir / "central_aupr_rows.csv")
    central_f1 = pd.read_csv(tables_dir / "central_f1_rows.csv")
    reference_summary = pd.read_csv(tables_dir / "reference_summary.csv")
    method_summary = pd.read_csv(tables_dir / "method_summary.csv")
    scenario_rank_stability = pd.read_csv(tables_dir / "scenario_rank_stability.csv")
    method_sensitivity = pd.read_csv(tables_dir / "method_sensitivity_summary.csv")
    mar_summary = pd.read_csv(tables_dir / "synthetic_mar_summary.csv")
    adv_summary = pd.read_csv(tables_dir / "synthetic_adversarial_summary.csv")
    sweep_rows = pd.read_csv(tables_dir / "sweep_artifact_rows.csv")
    cross_eval_rows = pd.read_csv(tables_dir / "cross_eval_artifact_rows.csv")

    best_f1 = central_aupr.sort_values("f1_ratio", ascending=False).iloc[0]
    best_aupr = central_aupr.sort_values("aupr_ratio", ascending=False).iloc[0]

    median_f1_ratio = float(central_aupr["f1_ratio"].median())
    median_aupr_ratio = float(central_aupr["aupr_ratio"].median())
    median_uplift = float(central_aupr["aupr_uplift_ratio"].median())

    f1_rows_total = len(central_f1)
    aupr_rows_total = len(central_aupr)
    method_count = central_aupr["method"].nunique()
    ref_count = central_aupr["reference"].nunique()

    nonzero_rank_shift = method_summary[method_summary["rank_shift_aupr"] != 0].copy()
    if nonzero_rank_shift.empty:
        rank_shift_text = "No AUPR rank shifts were observed between raw and normalized scores in the central scenario."
    else:
        shift_rows = []
        for row in nonzero_rank_shift.sort_values("rank_shift_aupr").itertuples(index=False):
            direction = "up" if row.rank_shift_aupr < 0 else "down"
            shift_rows.append(
                f"`{row.method}` moved {direction} by {abs(int(row.rank_shift_aupr))} rank position(s) "
                f"(raw rank {int(row.rank_observed_aupr)} -> normalized rank {int(row.rank_aupr_ratio)})."
            )
        rank_shift_text = " ".join(shift_rows)

    mar_bias_f1_median = float(mar_summary["f1_bias_abs"].median())
    mar_bias_aupr_median = float(mar_summary["aupr_bias_abs"].median())
    mar_bias_f1_max = float(mar_summary["f1_bias_abs"].max())
    mar_bias_aupr_max = float(mar_summary["aupr_bias_abs"].max())

    adv_mode_summary = (
        adv_summary.groupby("mode", as_index=False)
        .agg(
            f1_delta_mean=("f1_delta_mean", "mean"),
            aupr_delta_mean=("aupr_delta_mean", "mean"),
        )
        .sort_values("mode")
    )

    sensitivity_top = method_sensitivity.sort_values("aupr_ratio_range", ascending=False).head(8).copy()

    reference_table = reference_summary.copy()
    for col in [
        "c_median",
        "observed_f1_median",
        "f1_ceiling_median",
        "f1_ratio_median",
        "observed_aupr_median",
        "aupr_ceiling_median",
        "aupr_ratio_median",
        "aupr_uplift_ratio_median",
    ]:
        reference_table[col] = reference_table[col].map(lambda x: float(f"{x:.6g}"))

    method_table = (
        method_summary[
            [
                "method",
                "observed_aupr_median",
                "aupr_ratio_median",
                "rank_observed_aupr",
                "rank_aupr_ratio",
                "rank_shift_aupr",
                "observed_f1_median",
                "f1_ratio_median",
            ]
        ]
        .sort_values("rank_aupr_ratio")
        .copy()
    )

    for col in ["observed_aupr_median", "aupr_ratio_median", "observed_f1_median", "f1_ratio_median"]:
        method_table[col] = method_table[col].map(lambda x: float(f"{x:.6g}"))

    stability_table = scenario_rank_stability.copy()
    for col in ["spearman_vs_overlap_rank", "spearman_vs_raw_rank", "max_abs_rank_shift_vs_overlap"]:
        stability_table[col] = stability_table[col].map(lambda x: float(f"{x:.6g}"))

    mar_table = mar_summary[
        [
            "c",
            "t",
            "f1_obs_mean",
            "f1_theory",
            "f1_bias_abs",
            "aupr_obs_mean",
            "aupr_theory",
            "aupr_bias_abs",
        ]
    ].copy()
    for col in ["f1_obs_mean", "f1_theory", "f1_bias_abs", "aupr_obs_mean", "aupr_theory", "aupr_bias_abs"]:
        mar_table[col] = mar_table[col].map(lambda x: float(f"{x:.6g}"))

    adv_table = adv_mode_summary.copy()
    adv_table["f1_delta_mean"] = adv_table["f1_delta_mean"].map(lambda x: float(f"{x:.6g}"))
    adv_table["aupr_delta_mean"] = adv_table["aupr_delta_mean"].map(lambda x: float(f"{x:.6g}"))

    sensitivity_table = sensitivity_top[
        [
            "method",
            "aupr_ratio_min",
            "aupr_ratio_median",
            "aupr_ratio_max",
            "aupr_ratio_range",
            "f1_ratio_min",
            "f1_ratio_median",
            "f1_ratio_max",
        ]
    ].copy()
    for col in [
        "aupr_ratio_min",
        "aupr_ratio_median",
        "aupr_ratio_max",
        "aupr_ratio_range",
        "f1_ratio_min",
        "f1_ratio_median",
        "f1_ratio_max",
    ]:
        sensitivity_table[col] = sensitivity_table[col].map(lambda x: float(f"{x:.6g}"))

    sweep_table = sweep_rows[sweep_rows["entry_type"] == "global"][
        ["artifact", "candidate_edges", "aupr"]
    ].copy()
    sweep_table["aupr"] = sweep_table["aupr"].map(lambda x: float(f"{x:.6g}"))

    cross_top = (
        cross_eval_rows.dropna(subset=["aupr"])
        .sort_values("aupr", ascending=False)
        .head(6)[
            [
                "artifact",
                "system",
                "calibration",
                "candidate_edges",
                "positives",
                "observed_baserate",
                "aupr",
            ]
        ]
        .copy()
    )
    for col in ["observed_baserate", "aupr"]:
        cross_top[col] = cross_top[col].map(lambda x: float(f"{x:.6g}") if pd.notna(x) else x)

    output_path = repo_root / "paper/proposal1_partial_label_metric_ceilings_FINAL.md"

    manuscript = f"""# Partial-Label Metric Ceilings for GRN Evaluation
## A Theory-First, Adversarially Stress-Tested Study

## Abstract
Gene regulatory network (GRN) benchmarks are usually interpreted as if labels were complete, but curated references are partial observations of latent biology. We formalize this partial-label regime and derive explicit ceilings for observed F1 and AUPR as functions of positive-label coverage `c`. We then execute a full low-compute analysis on existing benchmark artifacts, propagate uncertainty in `c`, and adversarially stress-test assumptions using non-MAR missingness simulations. Across `{aupr_rows_total}` AUPR-evaluable rows ({method_count} methods, {ref_count} references), the best ceiling-normalized scores are modest (best F1 ratio `{best_f1['f1_ratio']:.3f}`, best AUPR ratio `{best_aupr['aupr_ratio']:.4f}`), with medians near zero (F1 ratio `{median_f1_ratio:.6f}`, AUPR ratio `{median_aupr_ratio:.6f}`). MAR-theory simulations show low absolute bias (median F1 bias `{mar_bias_f1_median:.3e}`, median AUPR bias `{mar_bias_aupr_median:.3e}`), while non-MAR missingness can induce substantial AUPR deviations (optimistic labeling: +`{adv_mode_summary[adv_mode_summary['mode']=='optimistic']['aupr_delta_mean'].iloc[0]:.3f}`; pessimistic: `{adv_mode_summary[adv_mode_summary['mode']=='pessimistic']['aupr_delta_mean'].iloc[0]:.3f}` vs MAR theory). We provide a reproducible ceiling-aware reporting framework and a finalized artifact bundle.

## 1. Introduction
Absolute benchmark metrics in GRN inference are hard to interpret under incomplete references. If only a fraction of true edges are labeled positive, observed precision-style metrics are systematically distorted, and even an ideal predictor can appear weak in absolute terms. This paper answers the study's central question directly: **what observed F1/AUPR is theoretically achievable under partial positive labels, and how far are current results from that ceiling?**

We contribute:
1. Closed-form ceiling formulas for observed precision/F1/AUPR under a missing-at-random (MAR) positive-label model.
2. Coverage-aware reinterpretation of existing outputs, including uncertainty bands and adversarial sensitivity scenarios.
3. Synthetic theorem validation and adversarial non-MAR tests that quantify where ceiling assumptions hold or fail.
4. A full reproducible paper package with tables, figures, and manuscript-ready outputs.

## 2. Formal Setup
Let:
- `U`: candidate edge universe.
- `T subset U`: latent true positives.
- `L+ subset T`: observed positives in a reference.
- `c = P(e in L+ | e in T)`: positive-label coverage.
- `P`: predicted positives for a method.

Define latent and observed metrics:
- `Prec_lat = |P intersect T| / |P|`, `Rec_lat = |P intersect T| / |T|`.
- `Prec_obs = |P intersect L+| / |P|`, `Rec_obs = |P intersect L+| / |L+|`.

Assumption A1 (MAR): each latent true edge is observed independently with probability `c`.

## 3. Theory
### Proposition 1: Observed-latent relation
Under A1:
- `E[Prec_obs] = c * Prec_lat`
- `E[Rec_obs] = Rec_lat`

So latent precision can be estimated by `Prec_obs / c` when `c` is estimated.

### Proposition 2: Observed F1 ceiling
For an ideal latent predictor (`Prec_lat = Rec_lat = 1`):
- `F1_obs,max(c) = 2c / (1 + c)`.

### Proposition 3: Observed AUPR ceiling
For latent positive count `t = |T|` and ideal latent ranking:
- `AUPR_obs,max(c, t) ~= c + (1-c) * H_t/t`, where `H_t` is harmonic.
- Large-`t` limit: `AUPR_obs,max ~= c`.

Observed random baseline under partial labels:
- `AUPR_rand,obs = c * b` where `b = |T|/|U|`.

Normalized headroom statistic:
- `rho_aupr = (AUPR_obs - AUPR_rand,obs) / (AUPR_obs,max - AUPR_rand,obs)`.

### Coverage uncertainty model
We model reference-level coverage as:
- `c ~ Beta(k+1, n-k+1)`
with overlap support `k` and denominator `n`, then propagate uncertainty by Monte Carlo draws.

## 4. Data and Protocol
### 4.1 Primary artifacts
- `data/score_eval_probe_priors.csv`
- `data/score_eval_grn_baselines_immune.csv`

### 4.2 F1-extended artifacts
- `data/score_eval_probe_priors_full_genes.csv`
- `data/score_eval_probe_priors_full_genes_crosswalk.csv`
- `data/score_eval_probe_priors_full_genes_omnipath.csv`

### 4.3 Coverage proxy artifacts
- `data/score_eval_probe_priors_missing_report.json`
- `data/score_eval_probe_priors_full_genes_missing_report.json`
- `data/score_eval_probe_priors_full_genes_crosswalk_missing_report.json`
- `data/score_eval_probe_priors_full_genes_omnipath_missing_report.json`

### 4.4 Coverage scenarios
For each reference, we evaluate:
- `overlap_point` (central).
- `beta_ci_low`, `beta_ci_high`.
- `adversarial_half`, `adversarial_150pct`.
- `mapping_min`, `mapping_max` (when available).
- `conservative_min`, `optimistic_max` envelope.

### 4.5 Synthetic experiments
- MAR simulations over grids of `c in {{0.05,0.1,0.2,0.4,0.7,0.9}}`, `t in {{50,200,1000,5000}}`.
- Non-MAR stress test with optimistic/pessimistic rank-dependent labeling.

## 5. Results
### 5.1 Core reinterpretation
- AUPR-evaluable rows: `{aupr_rows_total}`
- F1-evaluable rows: `{f1_rows_total}`
- Methods: `{method_count}`
- References: `{ref_count}`

Key extrema:
- Best F1 ceiling ratio: `{best_f1['method']}` on `{best_f1['reference']}` = `{best_f1['f1_ratio']:.6f}`
- Best AUPR ceiling ratio: `{best_aupr['method']}` on `{best_aupr['reference']}` = `{best_aupr['aupr_ratio']:.6f}`
- Median F1 ceiling ratio: `{median_f1_ratio:.6f}`
- Median AUPR ceiling ratio: `{median_aupr_ratio:.6f}`
- Median normalized AUPR uplift ratio: `{median_uplift:.6e}`

#### Table 1. Reference-level ceiling summary
{reference_table.to_markdown(index=False)}

Interpretation: references with low overlap-based coverage (`dorothea_human`, `dorothea_trrust_union_immune`) have lower absolute ceilings, but observed metrics still remain small fractions of those ceilings.

### 5.2 Method-level normalized rankings
#### Table 2. Method summary (raw vs normalized)
{method_table.to_markdown(index=False)}

Rank-shift finding: {rank_shift_text}

### 5.3 Coverage-scenario stability
#### Table 3. Scenario-level rank stability
{stability_table.to_markdown(index=False)}

Interpretation:
- Overlap/Beta/adversarial scaling scenarios preserve ordering (Spearman ~1.0).
- Mapping-driven scenarios (`mapping_min/max`) induce major shifts but only on a subset of methods/references where those proxies exist.

### 5.4 Adversarial coverage sensitivity
#### Table 4. Most coverage-sensitive methods (top 8 by AUPR-ratio range)
{sensitivity_table.to_markdown(index=False)}

Coverage sensitivity is highly method-dependent; some methods show large ceiling-ratio range under plausible-to-adversarial `c` bands.

### 5.5 MAR theorem validation
#### Table 5. Synthetic MAR validation (bias relative to theory)
{mar_table.to_markdown(index=False)}

Aggregate MAR validation:
- Median absolute F1 bias: `{mar_bias_f1_median:.3e}` (max `{mar_bias_f1_max:.3e}`)
- Median absolute AUPR bias: `{mar_bias_aupr_median:.3e}` (max `{mar_bias_aupr_max:.3e}`)

This supports the formulas as accurate expectation-level approximations under MAR.

### 5.6 Non-MAR adversarial missingness
#### Table 6. Mean deviation from MAR theory by missingness mode
{adv_table.to_markdown(index=False)}

AUPR is strongly sensitive to non-MAR labeling structure, while F1 shifts are smaller in this setup.

### 5.7 External sanity checks
#### Table 7. Sweep global AUPR values
{sweep_table.to_markdown(index=False)}

#### Table 8. Top cross-eval AUPR entries
{cross_top.to_markdown(index=False)}

These external artifacts are directionally consistent with the ceiling-aware framing: absolute AUPR varies strongly by candidate regime and coverage context.

## 6. Figures
### Figure 1. Theoretical ceilings
![Figure 1](../outputs/paper/figures/figure_01_theoretical_curves.png)

### Figure 2. Observed metrics vs ceiling ratios
![Figure 2](../outputs/paper/figures/figure_02_observed_vs_ratio.png)

### Figure 3. Rank shifts after normalization
![Figure 3](../outputs/paper/figures/figure_03_rank_shift.png)

### Figure 4. Scenario-specific method ranks
![Figure 4](../outputs/paper/figures/figure_04_scenario_rank_heatmap.png)

### Figure 5. Coverage sensitivity bands
![Figure 5](../outputs/paper/figures/figure_05_sensitivity_bands.png)

### Figure 6. MAR-theory validation
![Figure 6](../outputs/paper/figures/figure_06_mar_validation.png)

### Figure 7. Non-MAR adversarial effects
![Figure 7](../outputs/paper/figures/figure_07_adversarial_missingness.png)

## 7. Discussion
Main takeaway: low absolute metrics are not self-explanatory. Ceiling-aware interpretation separates two effects:
1. **Reference incompleteness** (what is fundamentally observable),
2. **Method headroom** (distance from the observable ceiling).

In these results, the second term dominates: observed scores are far below estimated ceilings for most methods/reference pairs.

Adversarial perspective:
- Coverage misspecification can alter normalized rankings, especially when proxy families disagree.
- Non-MAR missingness can substantially bias AUPR relative to MAR predictions.

Therefore, reporting should include both assumptions and sensitivity bands, not single-point ceilings.

## 8. Limitations
- `c` is proxy-estimated, not identified.
- Method/reference coverage is heterogeneous across artifacts.
- Non-MAR simulations are stylized stress tests, not direct biological curation models.

## 9. Reproducibility
Run full pipeline:

```bash
python scripts/run_paper_pipeline.py
python scripts/build_final_paper_markdown.py
```

Core outputs:
- `outputs/paper/results_summary.md`
- `outputs/paper/tables/*.csv`
- `outputs/paper/figures/*.png`
- `paper/proposal1_partial_label_metric_ceilings_FINAL.md`

## 10. Appendix: Proof Sketches
### A.1 Proposition 1
Linearity of expectation on Bernoulli-thinned positives gives `E[|P intersect L+|] = c|P intersect T|`, yielding the precision/recall relations.

### A.2 Proposition 2
For ideal latent predictions, observed recall is 1 and observed precision is `c`, so F1 reduces to `2c/(1+c)`.

### A.3 Proposition 3
Under ideal ordering, AP is mean precision at observed-positive ranks among top `t` positions. Taking expectation under Bernoulli thinning introduces the harmonic correction `H_t/t`.
"""

    output_path.write_text(manuscript)
    print(f"Wrote finalized manuscript: {output_path}")


if __name__ == "__main__":
    main()
