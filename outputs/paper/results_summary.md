# Full-Pipeline Results Summary

## Core reinterpretation
- AUPR+F1 evaluable rows: `39` across `15` methods and `5` references.
- Best F1 ceiling ratio: `genie3` on `dorothea_trrust_union_immune` = `0.137393`.
- Best AUPR ceiling ratio: `grnboost2` on `dorothea_trrust_union_immune` = `0.013771`.
- Median F1 ceiling ratio: `0.000000`.
- Median AUPR ceiling ratio: `0.000056`.
- Median normalized AUPR uplift ratio: `-2.904159e-06`.

## Reference-level medians
| reference                    |   n_rows |   methods |   c_median |   observed_f1_median |   f1_ceiling_median |   f1_ratio_median |   observed_aupr_median |   aupr_ceiling_median |   aupr_ratio_median |   aupr_uplift_ratio_median |
|:-----------------------------|---------:|----------:|-----------:|---------------------:|--------------------:|------------------:|-----------------------:|----------------------:|--------------------:|---------------------------:|
| dorothea_trrust_union_immune |        9 |         9 |   0.117816 |          0.00217946  |            0.210797 |        0.0103392  |            0.000343395 |             0.117925  |         0.00291199  |               -0.00111666  |
| dorothea_human               |        6 |         6 |   0.094976 |          0.000178214 |            0.173476 |        0.00102731 |            0.000274783 |             0.0951112 |         0.00288907  |               -0.000549558 |
| trrust_human                 |        6 |         6 |   0.159027 |          0           |            0.274414 |        0          |            9.03234e-06 |             0.162476  |         5.55919e-05 |               -2.90038e-05 |
| hpn_dream                    |        9 |         9 |   0.745098 |          0           |            0.853933 |        0          |            3.01472e-06 |             0.75365   |         4.00017e-06 |               -2.90372e-06 |
| beeline_gsd                  |        9 |         9 |   0.882353 |          0           |            0.9375   |        0          |            1.23677e-06 |             0.892264  |         1.38611e-06 |               -8.17867e-07 |

## Method-level ranking shifts
| method                     |   observed_aupr_median |   aupr_ratio_median |   rank_observed_aupr |   rank_aupr_ratio |   rank_shift_aupr |
|:---------------------------|-----------------------:|--------------------:|---------------------:|------------------:|------------------:|
| probe_consensus            |            0.000159345 |         0.00165566  |                    1 |                 1 |                 0 |
| probe_perturbation         |            0.000144556 |         0.00150017  |                    3 |                 2 |                -1 |
| probe_integrated_gradients |            0.000145163 |         0.00148826  |                    2 |                 3 |                 1 |
| probe_grad_input           |            0.00014285  |         0.00148224  |                    4 |                 4 |                 0 |
| probe_attention            |            0.000122129 |         0.00126438  |                    5 |                 5 |                 0 |
| attention_inferred         |            0.000121964 |         0.00126265  |                    6 |                 6 |                 0 |
| grnboost2                  |            3.01513e-06 |         4.00071e-06 |                    7 |                 7 |                 0 |
| genie3                     |            3.01509e-06 |         4.00065e-06 |                    8 |                 8 |                 0 |
| scenic_grnboost2           |            3.01484e-06 |         4.00032e-06 |                    9 |                 9 |                 0 |
| scenic_pruned              |            3.01475e-06 |         4.0002e-06  |                   10 |                10 |                 0 |
| pidc_full                  |            3.01472e-06 |         4.00017e-06 |                   11 |                11 |                 0 |
| pidc_proxy                 |            3.01472e-06 |         4.00016e-06 |                   12 |                12 |                 0 |
| spearman                   |            3.01439e-06 |         3.99972e-06 |                   13 |                13 |                 0 |
| pearson                    |            3.01439e-06 |         3.99972e-06 |                   14 |                14 |                 0 |
| random                     |            3.0138e-06  |         3.99893e-06 |                   15 |                15 |                 0 |

## Scenario rank stability
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

## Most coverage-sensitive method
- `grnboost2` has the widest AUPR-ratio band: `0.027537` (min `0.000001`, max `0.027539`).

## Synthetic theorem validation
- Median MAR absolute F1 bias: `7.298488e-04`.
- Median MAR absolute AUPR bias: `1.346207e-03`.
- Non-MAR missingness shifts AUPR relative to MAR theory:
| mode        |   aupr_delta_mean |   f1_delta_mean |
|:------------|------------------:|----------------:|
| mar         |       -0.00231265 |    -0.000259052 |
| optimistic  |        0.0897962  |    -0.000524601 |
| pessimistic |       -0.0733929  |    -0.000879459 |

## External artifact sanity checks
### Sweep globals
| artifact                                  |   candidate_edges |        aupr |
|:------------------------------------------|------------------:|------------:|
| omnipath                                  |            162835 | 0.000517499 |
| omnipath_dorothea_union_immune_hpn        |            870806 | 0.00127928  |
| omnipath_dorothea_intersection_immune_hpn |             29253 | 0.0251236   |
| regulatory                                |             43198 | 0.00859106  |

### Cross-eval top AUPR entries
| artifact                  | system   | calibration   |   candidate_edges |   positives |   observed_baserate |        aupr |         brier |
|:--------------------------|:---------|:--------------|------------------:|------------:|--------------------:|------------:|--------------:|
| cross_intersection        | hpn      | default       |             29253 |          23 |         0.000786244 | 0.00373156  | nan           |
| cross_intersection_scaled | hpn      | default       |             29253 |          23 |         0.000786244 | 0.00350346  | nan           |
| cross_union_scaled        | hpn      | default       |            870806 |         123 |         0.000141248 | 0.000677197 | nan           |
| cross_union               | hpn      | raw           |            870806 |         nan |       nan           | 0.000674443 |   0.000145614 |
| cross_union               | hpn      | logistic      |            870806 |         nan |       nan           | 0.000674443 |   0.000141782 |
| cross_union               | hpn      | isotonic      |            870806 |         nan |       nan           | 0.000645292 |   0.000142013 |

## Interpretation
- Ceiling-aware reinterpretation remains strongly conservative: observed metrics are tiny fractions of partial-label ceilings.
- Method ranking can shift after coverage normalization because methods are evaluated on different reference mixes and coverage regimes.
- MAR formulas are numerically accurate under their assumptions, but adversarial non-MAR labeling can move AUPR away from MAR-theory ceilings.