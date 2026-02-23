# Research Note: Partial-Label Metric Ceilings for GRN Evaluation

Date: 2026-02-10

## 1. Objective

Execute the partial-label metric ceiling analysis:
- derive partial-label metric ceiling formulas,
- estimate coverage proxies from existing mapping/overlap artifacts,
- reinterpret observed F1/AUPR as fractions of estimated ceilings,
- include uncertainty intervals from uncertainty in coverage `c`.

## 2. Formal Setup and Results

Let:
- `U`: candidate edge universe,
- `T ⊂ U`: latent true edges,
- `L+ ⊂ T`: observed labeled positives,
- `c = P(e ∈ L+ | e ∈ T)`: positive-label coverage,
- `P ⊂ U`: predicted positive edges.

Assumption used in this first-pass execution: labels are missing-at-random within `T` (independent Bernoulli thinning with rate `c`).

### Proposition 1 (Observed vs latent precision/recall)

Define latent and observed precision/recall:

- `Prec_lat = |P ∩ T| / |P|`
- `Rec_lat = |P ∩ T| / |T|`
- `Prec_obs = |P ∩ L+| / |P|`
- `Rec_obs = |P ∩ L+| / |L+|`

Under Bernoulli thinning:

- `E[Prec_obs] = c · Prec_lat`
- `E[Rec_obs] = Rec_lat`

So an estimator for latent precision is `Prec_lat ≈ Prec_obs / c`.

### Proposition 2 (Observed F1 ceiling)

For any method:

- `F1_obs = 2 Prec_obs Rec_obs / (Prec_obs + Rec_obs)`

For an ideal latent predictor (`Prec_lat = Rec_lat = 1`):

- `Prec_obs = c`, `Rec_obs = 1`
- `F1_obs,max(c) = 2c / (1 + c)`

### Proposition 3 (Observed AUPR ceiling)

Let latent positive count `t = |T| = b|U|` where `b` is latent base rate.

Under ideal latent ranking (all true edges above all negatives), expected observed AUPR is approximated by:

- `AUPR_obs,max(c,t) ≈ c + (1-c)·(H_t / t)`

where `H_t` is the harmonic number. Large-`t` limit:

- `AUPR_obs,max ≈ c`

Observed random baseline under partial labels:

- `AUPR_rand,obs = |L+|/|U| = cb`

Useful normalized headroom statistic:

- `rho_aupr = (AUPR_obs - AUPR_rand,obs) / (AUPR_obs,max - AUPR_rand,obs)`

## 3. Coverage Proxy Estimation (Minimal Empirical Component)

### Data used

- Metrics:
  - `data/score_eval_probe_priors.csv`
  - `data/score_eval_grn_baselines_immune.csv`
- Mapping/overlap reports:
  - `data/score_eval_probe_priors_missing_report.json`
  - `data/score_eval_probe_priors_full_genes_missing_report.json`
  - `data/score_eval_probe_priors_full_genes_crosswalk_missing_report.json`
  - `data/score_eval_probe_priors_full_genes_omnipath_missing_report.json`

### Proxy definitions

- Node-overlap proxy (used as default `c` in ceiling calculations):
  - `c_hat = overlap_nodes / ref_nodes`
- Mapping endpoint proxy from missing reports:
  - `c_map = 1 - missing_total / (2 · total_edges_normalized)`

### Estimated proxies

From `outputs/coverage_proxy_summary.csv`:

- Mapping endpoint proxies vary strongly by policy:
  - `hpn_dream`: `0.0496` (`probe_priors`) to `1.0000` (`probe_priors_full_genes`)
  - `beeline_gsd`: `0.0592` (`probe_priors`) to `0.8026` (`probe_priors_full_genes_crosswalk`)
  - `omnipath_interactions`: `0.9948`
- Node-overlap proxies used for reinterpretation:
  - `dorothea_human`: `0.0950`
  - `trrust_human`: `0.1590`
  - `dorothea_trrust_union_immune`: `0.1178`
  - `hpn_dream`: `0.7451`
  - `beeline_gsd`: `0.8824`

Interpretation: `c` is not identifiable from one artifact; mapping and overlap choices create large proxy spread.

## 4. Ceiling Reinterpretation Results

Evaluation rows with non-null observed F1 and AUPR: `39` rows, `15` methods, `5` references.

From `outputs/metric_ceiling_reinterpretation.csv`:

- Best F1 ceiling ratio: `0.137` (genie3 on `dorothea_trrust_union_immune`, 95% CI `[0.133, 0.142]`)
- Best AUPR ceiling ratio: `0.0138` (grnboost2 on `dorothea_trrust_union_immune`, 95% CI `[0.0132, 0.0143]`)
- Global median F1 ceiling ratio: `0.0000`
- Global median AUPR ceiling ratio: `5.56e-05`
- Global median normalized AUPR uplift ratio `rho_aupr`: `-2.90e-06`

Reference-level medians (`outputs/reference_ceiling_summary.csv`):

| reference | c_proxy_median | observed_f1_median | f1_ceiling_median | f1_ratio_median | observed_aupr_median | aupr_ceiling_median | aupr_ratio_median |
|---|---:|---:|---:|---:|---:|---:|---:|
| dorothea_trrust_union_immune | 0.1178 | 0.002179 | 0.210797 | 0.010339 | 0.000343 | 0.117925 | 0.002912 |
| dorothea_human | 0.0950 | 0.000178 | 0.173476 | 0.001027 | 0.000275 | 0.095111 | 0.002889 |
| trrust_human | 0.1590 | 0.000000 | 0.274414 | 0.000000 | 0.000009 | 0.162476 | 0.000056 |
| hpn_dream | 0.7451 | 0.000000 | 0.853933 | 0.000000 | 0.000003 | 0.753650 | 0.000004 |
| beeline_gsd | 0.8824 | 0.000000 | 0.937500 | 0.000000 | 0.000001 | 0.892264 | 0.000001 |

## 5. Uncertainty Intervals

Coverage uncertainty was modeled with Beta posteriors per reference:

- `c ~ Beta(overlap_nodes + 1, ref_nodes - overlap_nodes + 1)`

We propagated this with 12,000 Monte Carlo draws to produce 95% intervals for:
- `c`,
- `F1_obs,max(c)`,
- `AUPR_obs,max(c,t_hat)`,
- ratio metrics (`observed / ceiling`).

These intervals are in `outputs/metric_ceiling_reinterpretation.csv` columns:
- `coverage_proxy_ci_low/high`
- `f1_ceiling_ci_low/high`
- `aupr_ceiling_ci_low/high`
- `f1_ceiling_ratio_ci_low/high`
- `aupr_ceiling_ratio_ci_low/high`

## 6. Figures

- `outputs/figures/observed_metric_vs_ceiling_ratio.png`
  - Observed F1 and AUPR plotted against estimated ceiling ratios with uncertainty bars.
- `outputs/figures/theoretical_ceiling_curves.png`
  - Ceiling functions over coverage `c`, including observed random AUPR baselines (`cb`) for representative latent base rates.

## 7. Caveats and Risk Control

- `c` is structurally non-identifiable from current outputs; proxies are only sensitivity anchors.
- Reported CIs quantify uncertainty conditional on chosen proxy definition, not full epistemic uncertainty.
- The framework should be reported as "ceiling-aware reinterpretation" rather than as recovered latent truth.

## 8. Reproducibility

Command:

```bash
python scripts/run_analysis.py
```

Main script:
- `scripts/run_analysis.py`
