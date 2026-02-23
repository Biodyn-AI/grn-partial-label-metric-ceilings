# Partial-Label Metric Ceilings for GRN Evaluation
## A Theory-First, Adversarially Stress-Tested Study

## Abstract
Gene regulatory network (GRN) benchmarks are usually interpreted as if labels were complete, but curated references are partial observations of latent biology. We formalize this partial-label regime and derive explicit ceilings for observed F1 and AUPR as functions of positive-label coverage `c`. We then execute a full low-compute analysis on existing benchmark artifacts, propagate uncertainty in `c`, and adversarially stress-test assumptions using non-MAR missingness simulations. Across `39` AUPR-evaluable rows (15 methods, 5 references), the best ceiling-normalized scores are modest (best F1 ratio `0.137`, best AUPR ratio `0.0138`), with medians near zero (F1 ratio `0.000000`, AUPR ratio `0.000056`). MAR-theory simulations show low absolute bias (median F1 bias `7.298e-04`, median AUPR bias `1.346e-03`), while non-MAR missingness can induce substantial AUPR deviations (optimistic labeling: +`0.090`; pessimistic: `-0.073` vs MAR theory). We provide a reproducible ceiling-aware reporting framework and a finalized artifact bundle.

## 1. Introduction
Absolute benchmark metrics in GRN inference are hard to interpret under incomplete references. If only a fraction of true edges are labeled positive, observed precision-style metrics are systematically distorted, and even an ideal predictor can appear weak in absolute terms. This paper answers the central question directly: **what observed F1/AUPR is theoretically achievable under partial positive labels, and how far are current results from that ceiling?**

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
- MAR simulations over grids of `c in {0.05,0.1,0.2,0.4,0.7,0.9}`, `t in {50,200,1000,5000}`.
- Non-MAR stress test with optimistic/pessimistic rank-dependent labeling.

## 5. Results
### 5.1 Core reinterpretation
- AUPR-evaluable rows: `39`
- F1-evaluable rows: `117`
- Methods: `15`
- References: `5`

Key extrema:
- Best F1 ceiling ratio: `genie3` on `dorothea_trrust_union_immune` = `0.137393`
- Best AUPR ceiling ratio: `grnboost2` on `dorothea_trrust_union_immune` = `0.013771`
- Median F1 ceiling ratio: `0.000000`
- Median AUPR ceiling ratio: `0.000056`
- Median normalized AUPR uplift ratio: `-2.904159e-06`

#### Table 1. Reference-level ceiling summary
| reference                    |   n_rows |   methods |   c_median |   observed_f1_median |   f1_ceiling_median |   f1_ratio_median |   observed_aupr_median |   aupr_ceiling_median |   aupr_ratio_median |   aupr_uplift_ratio_median |
|:-----------------------------|---------:|----------:|-----------:|---------------------:|--------------------:|------------------:|-----------------------:|----------------------:|--------------------:|---------------------------:|
| dorothea_trrust_union_immune |        9 |         9 |   0.117816 |          0.00217946  |            0.210797 |        0.0103392  |            0.000343395 |             0.117925  |         0.00291199  |               -0.00111666  |
| dorothea_human               |        6 |         6 |   0.094976 |          0.000178214 |            0.173476 |        0.00102731 |            0.000274783 |             0.0951112 |         0.00288907  |               -0.000549558 |
| trrust_human                 |        6 |         6 |   0.159027 |          0           |            0.274414 |        0          |            9.03234e-06 |             0.162476  |         5.55919e-05 |               -2.90038e-05 |
| hpn_dream                    |        9 |         9 |   0.745098 |          0           |            0.853933 |        0          |            3.01472e-06 |             0.75365   |         4.00017e-06 |               -2.90372e-06 |
| beeline_gsd                  |        9 |         9 |   0.882353 |          0           |            0.9375   |        0          |            1.23677e-06 |             0.892264  |         1.38611e-06 |               -8.17867e-07 |

Interpretation: references with low overlap-based coverage (`dorothea_human`, `dorothea_trrust_union_immune`) have lower absolute ceilings, but observed metrics still remain small fractions of those ceilings.

### 5.2 Method-level normalized rankings
#### Table 2. Method summary (raw vs normalized)
| method                     |   observed_aupr_median |   aupr_ratio_median |   rank_observed_aupr |   rank_aupr_ratio |   rank_shift_aupr |   observed_f1_median |   f1_ratio_median |
|:---------------------------|-----------------------:|--------------------:|---------------------:|------------------:|------------------:|---------------------:|------------------:|
| probe_consensus            |            0.000159345 |         0.00165566  |                    1 |                 1 |                 0 |          5.4582e-05  |       5.49029e-05 |
| probe_perturbation         |            0.000144556 |         0.00150017  |                    3 |                 2 |                -1 |          3.92779e-05 |       3.95088e-05 |
| probe_integrated_gradients |            0.000145163 |         0.00148826  |                    2 |                 3 |                 1 |          3.92765e-05 |       3.95074e-05 |
| probe_grad_input           |            0.00014285  |         0.00148224  |                    4 |                 4 |                 0 |          3.92807e-05 |       3.95116e-05 |
| probe_attention            |            0.000122129 |         0.00126438  |                    5 |                 5 |                 0 |          4.67845e-05 |       4.70596e-05 |
| attention_inferred         |            0.000121964 |         0.00126265  |                    6 |                 6 |                 0 |          0           |       0           |
| grnboost2                  |            3.01513e-06 |         4.00071e-06 |                    7 |                 7 |                 0 |          0           |       0           |
| genie3                     |            3.01509e-06 |         4.00065e-06 |                    8 |                 8 |                 0 |          0           |       0           |
| scenic_grnboost2           |            3.01484e-06 |         4.00032e-06 |                    9 |                 9 |                 0 |          0           |       0           |
| scenic_pruned              |            3.01475e-06 |         4.0002e-06  |                   10 |                10 |                 0 |          0           |       0           |
| pidc_full                  |            3.01472e-06 |         4.00017e-06 |                   11 |                11 |                 0 |          0           |       0           |
| pidc_proxy                 |            3.01472e-06 |         4.00016e-06 |                   12 |                12 |                 0 |          0           |       0           |
| spearman                   |            3.01439e-06 |         3.99972e-06 |                   13 |                13 |                 0 |          0           |       0           |
| pearson                    |            3.01439e-06 |         3.99972e-06 |                   14 |                14 |                 0 |          0           |       0           |
| random                     |            3.0138e-06  |         3.99893e-06 |                   15 |                15 |                 0 |          0           |       0           |

Rank-shift finding: `probe_perturbation` moved up by 1 rank position(s) (raw rank 3 -> normalized rank 2). `probe_integrated_gradients` moved down by 1 rank position(s) (raw rank 2 -> normalized rank 3).

### 5.3 Coverage-scenario stability
#### Table 3. Scenario-level rank stability
| scenario           |   n_methods |   spearman_vs_overlap_rank |   spearman_vs_raw_rank |   max_abs_rank_shift_vs_overlap |
|:-------------------|------------:|---------------------------:|-----------------------:|--------------------------------:|
| adversarial_150pct |          15 |                       1    |               0.996429 |                               0 |
| adversarial_half   |          15 |                       1    |               0.996429 |                               0 |
| beta_ci_high       |          15 |                       1    |               0.996429 |                               0 |
| beta_ci_low        |          15 |                       1    |               0.996429 |                               0 |
| conservative_min   |          15 |                       1    |               0.996429 |                               0 |
| mapping_max        |           9 |                       0.5  |               0.5      |                              11 |
| mapping_min        |           9 |                       0.55 |               0.55     |                              10 |
| optimistic_max     |          15 |                       1    |               0.996429 |                               0 |
| overlap_point      |          15 |                       1    |               0.996429 |                               0 |

Interpretation:
- Overlap/Beta/adversarial scaling scenarios preserve ordering (Spearman ~1.0).
- Mapping-driven scenarios (`mapping_min/max`) induce major shifts but only on a subset of methods/references where those proxies exist.

### 5.4 Adversarial coverage sensitivity
#### Table 4. Most coverage-sensitive methods (top 8 by AUPR-ratio range)
| method             |   aupr_ratio_min |   aupr_ratio_median |   aupr_ratio_max |   aupr_ratio_range |   f1_ratio_min |   f1_ratio_median |   f1_ratio_max |
|:-------------------|-----------------:|--------------------:|-----------------:|-------------------:|---------------:|------------------:|---------------:|
| grnboost2          |      1.23789e-06 |         7.0927e-06  |       0.0275386  |         0.0275374  |              0 |       0           |     0.257147   |
| genie3             |      1.2379e-06  |         7.09267e-06 |       0.0221169  |         0.0221157  |              0 |       0           |     0.260305   |
| scenic_grnboost2   |      1.23865e-06 |         7.09456e-06 |       0.00770431 |         0.00770307 |              0 |       0           |     0.0954181  |
| probe_consensus    |      3.71711e-05 |         0.00114075  |       0.00651042 |         0.00647324 |              0 |       0.000600549 |     0.00330375 |
| pidc_proxy         |      1.2388e-06  |         7.09485e-06 |       0.00613842 |         0.00613718 |              0 |       0           |     0.0223148  |
| probe_perturbation |      3.71714e-05 |         0.00103707  |       0.00588855 |         0.00585137 |              0 |       0.00042225  |     0.00232289 |
| spearman           |      1.23786e-06 |         7.09175e-06 |       0.0058233  |         0.00582206 |              0 |       0           |     0.0195886  |
| pearson            |      1.23786e-06 |         7.09175e-06 |       0.00580037 |         0.00579913 |              0 |       0           |     0.0192962  |

Coverage sensitivity is highly method-dependent; some methods show large ceiling-ratio range under plausible-to-adversarial `c` bands.

### 5.5 MAR theorem validation
#### Table 5. Synthetic MAR validation (bias relative to theory)
|    c |    t |   f1_obs_mean |   f1_theory |   f1_bias_abs |   aupr_obs_mean |   aupr_theory |   aupr_bias_abs |
|-----:|-----:|--------------:|------------:|--------------:|----------------:|--------------:|----------------:|
| 0.05 |   50 |     0.100644  |   0.0952381 |   0.00540598  |       0.127244  |     0.135485  |     0.0082409   |
| 0.05 |  200 |     0.0932589 |   0.0952381 |   0.00197922  |       0.0719671 |     0.0779206 |     0.0059535   |
| 0.05 | 1000 |     0.095733  |   0.0952381 |   0.000494944 |       0.0573385 |     0.0571112 |     0.000227278 |
| 0.05 | 5000 |     0.0949053 |   0.0952381 |   0.000332764 |       0.0512554 |     0.051728  |     0.000472586 |
| 0.1  |   50 |     0.180304  |   0.181818  |   0.0015142   |       0.155571  |     0.180986  |     0.0254151   |
| 0.1  |  200 |     0.182858  |   0.181818  |   0.00103984  |       0.123486  |     0.126451  |     0.00296543  |
| 0.1  | 1000 |     0.180888  |   0.181818  |   0.000929746 |       0.10488   |     0.106737  |     0.00185676  |
| 0.1  | 5000 |     0.181712  |   0.181818  |   0.000106567 |       0.101468  |     0.101637  |     0.000168657 |
| 0.2  |   50 |     0.32473   |   0.333333  |   0.00860329  |       0.254503  |     0.271987  |     0.0174844   |
| 0.2  |  200 |     0.332166  |   0.333333  |   0.00116688  |       0.219643  |     0.223512  |     0.00386902  |
| 0.2  | 1000 |     0.33163   |   0.333333  |   0.00170304  |       0.204336  |     0.205988  |     0.00165211  |
| 0.2  | 5000 |     0.333516  |   0.333333  |   0.000183104 |       0.201402  |     0.201455  |     5.2627e-05  |
| 0.4  |   50 |     0.563338  |   0.571429  |   0.00809046  |       0.434901  |     0.45399   |     0.0190897   |
| 0.4  |  200 |     0.569741  |   0.571429  |   0.00168727  |       0.412161  |     0.417634  |     0.005473    |
| 0.4  | 1000 |     0.571482  |   0.571429  |   5.37676e-05 |       0.405464  |     0.404491  |     0.000972989 |
| 0.4  | 5000 |     0.571631  |   0.571429  |   0.000201988 |       0.40085   |     0.401091  |     0.000241273 |
| 0.7  |   50 |     0.82088   |   0.823529  |   0.00264951  |       0.723628  |     0.726995  |     0.00336737  |
| 0.7  |  200 |     0.823298  |   0.823529  |   0.000231415 |       0.707108  |     0.708817  |     0.00170914  |
| 0.7  | 1000 |     0.822737  |   0.823529  |   0.000792673 |       0.701275  |     0.702246  |     0.000970239 |
| 0.7  | 5000 |     0.823673  |   0.823529  |   0.000143324 |       0.701141  |     0.700546  |     0.000595172 |
| 0.9  |   50 |     0.94789   |   0.947368  |   0.000521818 |       0.909245  |     0.908998  |     0.000246889 |
| 0.9  |  200 |     0.948035  |   0.947368  |   0.000667024 |       0.903768  |     0.902939  |     0.000828809 |
| 0.9  | 1000 |     0.94798   |   0.947368  |   0.000612057 |       0.901789  |     0.900749  |     0.00104031  |
| 0.9  | 5000 |     0.947294  |   0.947368  |   7.41732e-05 |       0.900004  |     0.900182  |     0.000178222 |

Aggregate MAR validation:
- Median absolute F1 bias: `7.298e-04` (max `8.603e-03`)
- Median absolute AUPR bias: `1.346e-03` (max `2.542e-02`)

This supports the formulas as accurate expectation-level approximations under MAR.

### 5.6 Non-MAR adversarial missingness
#### Table 6. Mean deviation from MAR theory by missingness mode
| mode        |   f1_delta_mean |   aupr_delta_mean |
|:------------|----------------:|------------------:|
| mar         |    -0.000259052 |       -0.00231265 |
| optimistic  |    -0.000524601 |        0.0897962  |
| pessimistic |    -0.000879459 |       -0.0733929  |

AUPR is strongly sensitive to non-MAR labeling structure, while F1 shifts are smaller in this setup.

### 5.7 External sanity checks
#### Table 7. Sweep global AUPR values
| artifact                                  |   candidate_edges |        aupr |
|:------------------------------------------|------------------:|------------:|
| omnipath                                  |            162835 | 0.000517499 |
| omnipath_dorothea_union_immune_hpn        |            870806 | 0.00127928  |
| omnipath_dorothea_intersection_immune_hpn |             29253 | 0.0251236   |
| regulatory                                |             43198 | 0.00859106  |

#### Table 8. Top cross-eval AUPR entries
| artifact                  | system   | calibration   |   candidate_edges |   positives |   observed_baserate |        aupr |
|:--------------------------|:---------|:--------------|------------------:|------------:|--------------------:|------------:|
| cross_intersection        | hpn      | default       |             29253 |          23 |         0.000786244 | 0.00373156  |
| cross_intersection_scaled | hpn      | default       |             29253 |          23 |         0.000786244 | 0.00350346  |
| cross_union_scaled        | hpn      | default       |            870806 |         123 |         0.000141248 | 0.000677197 |
| cross_union               | hpn      | raw           |            870806 |         nan |       nan           | 0.000674443 |
| cross_union               | hpn      | logistic      |            870806 |         nan |       nan           | 0.000674443 |
| cross_union               | hpn      | isotonic      |            870806 |         nan |       nan           | 0.000645292 |

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
