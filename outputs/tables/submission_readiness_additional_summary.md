# Submission-Readiness Additional Analyses

## Bootstrap uncertainty (stratified by reference)
| metric                   |   point_estimate |   bootstrap_mean |   bootstrap_std |     ci95_low |    ci95_high |
|:-------------------------|-----------------:|-----------------:|----------------:|-------------:|-------------:|
| median_aupr_ratio        |      5.55915e-05 |      5.56089e-05 |     9.43268e-07 |  5.55915e-05 |  5.5592e-05  |
| median_f1_ratio          |      0           |      0           |     0           |  0           |  0           |
| median_aupr_uplift_ratio |     -2.90416e-06 |     -7.12366e-06 |     9.61001e-06 | -2.90037e-05 | -2.90368e-06 |
| fraction_positive_uplift |      0.0769231   |      0.0782393   |     0.0405328   |  0           |  0.153846    |
| best_aupr_ratio          |      0.0137709   |      0.0120963   |     0.00303087  |  0.00383768  |  0.0137709   |
| best_f1_ratio            |      0.137393    |      0.126963    |     0.0301602   |  0.0117781   |  0.137393    |

## Positive uplift prevalence by reference
| reference                    |   n_rows |   n_positive_uplift |   positive_uplift_fraction |   positive_uplift_ci95_low |   positive_uplift_ci95_high |   median_aupr_ratio |   median_aupr_uplift_ratio |
|:-----------------------------|---------:|--------------------:|---------------------------:|---------------------------:|----------------------------:|--------------------:|---------------------------:|
| dorothea_trrust_union_immune |        9 |                   2 |                   0.222222 |                0.0632251   |                    0.547411 |         0.00291199  |               -0.00111666  |
| trrust_human                 |        6 |                   1 |                   0.166667 |                0.0300534   |                    0.563503 |         5.55919e-05 |               -2.90038e-05 |
| beeline_gsd                  |        9 |                   0 |                   0        |                0           |                    0.299145 |         1.38611e-06 |               -8.17867e-07 |
| dorothea_human               |        6 |                   0 |                   0        |                2.77556e-17 |                    0.390334 |         0.00288907  |               -0.000549558 |
| hpn_dream                    |        9 |                   0 |                   0        |                0           |                    0.299145 |         4.00017e-06 |               -2.90372e-06 |

## Central benchmark family summary
| method_family   |   rows |   methods |   median_aupr |   median_aupr_ratio |   median_f1_ratio |   positive_uplift_fraction |
|:----------------|-------:|----------:|--------------:|--------------------:|------------------:|---------------------------:|
| probe_family    |     12 |         6 |   0.000126162 |         0.00128848  |       6.19291e-05 |                  0.0833333 |
| classical_grn   |     21 |         7 |   3.01472e-06 |         4.00017e-06 |       0           |                  0.0952381 |
| control         |      6 |         2 |   3.01427e-06 |         3.99957e-06 |       0           |                  0         |

## Confound checks (correlation diagnostics)
| target_metric     | candidate_confound   |   n_rows |   pearson_r |   spearman_r |
|:------------------|:---------------------|---------:|------------:|-------------:|
| aupr_ratio        | candidate_edges      |       39 |   0.0397021 |    -0.355427 |
| aupr_ratio        | coverage_scenario_c  |       39 |  -0.532228  |    -0.922158 |
| aupr_ratio        | observed_base_rate   |       39 |   0.738093  |     0.950102 |
| aupr_uplift_ratio | candidate_edges      |       39 |   0.153233  |     0.335681 |
| aupr_uplift_ratio | coverage_scenario_c  |       39 |  -0.063502  |     0.716613 |
| aupr_uplift_ratio | observed_base_rate   |       39 |   0.146182  |    -0.678733 |

## External time-series reinterpretation constants
- HPN random AUPR baseline: `3.013954e-06`
- BEELINE random AUPR baseline: `1.236963e-06`
- HPN coverage proxy (overlap-based): `0.745098`
- BEELINE coverage proxy (overlap-based): `0.882353`

## External time-series family summary
| method_family     |   methods |   median_hpn_aupr |   median_hpn_ratio_to_ceiling_approx |   hpn_positive_uplift_fraction |   median_beeline_aupr |   median_beeline_ratio_to_ceiling_approx |   beeline_positive_uplift_fraction |
|:------------------|----------:|------------------:|-------------------------------------:|-------------------------------:|----------------------:|-----------------------------------------:|-----------------------------------:|
| prior_constrained |        12 |       0.000698699 |                          0.000937727 |                       1        |           2.63232e-05 |                              2.98329e-05 |                           0.916667 |
| classical_grn     |         7 |       3.01472e-06 |                          4.04608e-06 |                       1        |           1.23677e-06 |                              1.40168e-06 |                           0.428571 |
| control           |         3 |       3.01411e-06 |                          4.04526e-06 |                       0.666667 |           1.23724e-06 |                              1.4022e-06  |                           0.666667 |

## Perturbation validation summary
- Overall precision: `3.451251e-03`
- Overall recall: `1.504042e-04`
- Overall F1: `2.882467e-04`
- Perturbations with non-zero recall: `19/4787` (`0.3969%`)

## Cross-assay concordance diagnostic
| comparison                            |   n_methods |     x_range |     y_range |   spearman_r |   pearson_r | informative   | note                                       |
|:--------------------------------------|------------:|------------:|------------:|-------------:|------------:|:--------------|:-------------------------------------------|
| central_hpn_vs_timeseries_hpn         |           9 | 1.77477e-09 | 1.79514e-09 |          nan |         nan | False         | degenerate_dynamic_range_not_interpretable |
| central_beeline_vs_timeseries_beeline |           9 | 1.10483e-09 | 1.11724e-09 |          nan |         nan | False         | degenerate_dynamic_range_not_interpretable |

## Interpretation
- Bootstrap intervals confirm the central claim is robust: ceiling-normalized performance is low even under resampling uncertainty.
- Positive uplift is sparse and reference-specific, indicating biological compatibility matters more than a single global leaderboard.
- Time-series and perturbation validations are directionally consistent with the central benchmark: only constrained/prior-driven settings produce non-trivial uplift, and probe perturbation recall remains sparse.