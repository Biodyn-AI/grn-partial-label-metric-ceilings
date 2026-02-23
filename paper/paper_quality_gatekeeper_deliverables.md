# Paper Quality Gatekeeper Deliverables

## A) VERDICT: READY

### Brief justification
- **Central claim (one sentence):** In this benchmark, partial-label incompleteness explains only part of low GRN scores; after ceiling correction, most methods still operate near baseline, indicating substantial residual model-to-biology mismatch.
- **Main contributions:**
  1. Formal partial-label ceilings for observed F1 and AUPR under MAR.
  2. Empirical reinterpretation across 15 methods and 5 references with coverage uncertainty.
  3. Adversarial non-MAR and coverage-misspecification stress tests.
  4. External time-series and perturbation validation analyses with confound diagnostics.
  5. A ceiling-aware reporting protocol for mechanistic GRN claims.
- **Evidence for each contribution:**
  1. Closed-form derivations + synthetic validation (median absolute MAR bias: F1 `7.30e-04`, AUPR `1.35e-03`).
  2. 39-row reinterpretation: best normalized F1 `0.137`, best normalized AUPR `0.0138`, median AUPR ratio `5.56e-05`.
  3. Non-MAR shifts: mean AUPR delta `+0.0898` (optimistic) and `-0.0734` (pessimistic) versus MAR theory; mapping scenarios shift rank by up to 10-11 positions.
  4. Perturbation validation: non-zero recall in `19/4787` perturbations (0.397%); time-series uplift concentrated in prior-constrained settings.
  5. Explicit protocol integrated into manuscript (raw + normalized metrics, baseline uplift, uncertainty bands, non-MAR checks).
- **Top likely reviewer attacks (pre-mitigation):**
  1. “This is only a statistical reframing, not biological insight.”
  2. “Coverage proxy is fragile and potentially arbitrary.”
  3. “MAR assumptions are unrealistic for curated biology.”
  4. “External validation is weak or assay-specific.”
  5. “Results might be confounded by candidate size or base rates.”

## B) Top 5 reviewer objections and how addressed
1. **Objection:** Biological meaning is underdeveloped.
   **Addressed by:** Dedicated Biological Interpretation section now explains why transcriptional versus signaling references diverge and what this implies for mechanism claims and wet-lab prioritization.
2. **Objection:** Coverage assumptions are brittle.
   **Addressed by:** Multi-scenario coverage analysis (overlap, Beta CI, adversarial 0.5x/1.5x, mapping min/max, conservative/optimistic envelopes) with explicit rank-stability and rank-shift reporting.
3. **Objection:** MAR-based ceilings may be invalid in curated resources.
   **Addressed by:** Non-MAR simulations quantify directional bias magnitudes and convert assumption risk into measurable sensitivity bounds instead of hidden assumptions.
4. **Objection:** Lack of independent validation.
   **Addressed by:** Added external time-series reinterpretation and perturbation validation summary, including concentration analysis of non-zero recall and conservative interpretation of assay-limited signal.
5. **Objection:** Metrics may be driven by trivial confounds.
   **Addressed by:** Added confound diagnostics (coverage, candidate edges, observed base rate) and claim calibration in Discussion/Limitations to avoid overstating causal conclusions.

## C) Revised paper in full (submission-ready form)

# Partial-Label Metric Ceilings for GRN Evaluation in Single-Cell Foundation Models

## Abstract
Gene regulatory network (GRN) benchmarks in single-cell foundation-model studies are usually interpreted as if curated references were complete. They are not. We develop a theory-first framework for interpreting observed F1 and AUPR when only a fraction of true edges are labeled, and we test the framework on existing GRN benchmark outputs across 15 methods and 5 references. Under a missing-at-random (MAR) positive-label model, we derive explicit ceilings for observable F1 and AUPR as functions of label coverage and propagate uncertainty with Beta posteriors. In 39 AUPR-evaluable rows, the best ceiling-normalized F1 and AUPR remain modest (0.137 and 0.0138), the median normalized AUPR is 5.56e-05, and only 3/39 rows exceed their observed random baseline. Non-MAR stress tests produce substantial AUPR shifts (+0.0898 optimistic, -0.0734 pessimistic versus MAR theory), showing that curation mechanism matters in addition to curation magnitude. External time-series and perturbation validations are directionally consistent with the central finding that current scores are mostly far below observable headroom, with sparse positive uplift concentrated in specific reference regimes. Biologically, this implies that low scores are not explained by missing labels alone; substantial model-to-biology mismatch remains. We provide a ceiling-aware reporting protocol that reduces both overclaiming and underclaiming in mechanistic GRN studies.

## 1. Introduction
Single-cell foundation models are increasingly used to generate putative regulatory edges and prioritize downstream experiments [1,15]. Evaluation is typically done by comparing inferred edges with curated TF-target resources such as TRRUST and DoRothEA-derived collections [3,4]. This practice is useful but systematically ambiguous: a low observed score can reflect poor biological recovery, incomplete reference labeling, or both.

This ambiguity matters for mechanistic interpretability in biology. If evaluation ignores partial labels, one may over-penalize models in under-annotated contexts or overstate progress in references that reward canonical edges [7,9]. We therefore ask a precise question: *given partial positive labels, what metric values are theoretically observable, and how close are current methods to that observable ceiling?*

This paper makes four contributions.
1. A formal partial-label ceiling framework for observed F1 and AUPR under MAR assumptions [8,9].
2. A multi-reference empirical reinterpretation across 15 methods and 39 AUPR-evaluable benchmark rows.
3. Adversarial diagnostics for coverage misspecification and non-MAR curation.
4. Biological interpretation of what these results imply for the current state of GRN claims from single-cell foundation-model pipelines.

## 2. Related Work
### 2.1 GRN benchmarks and curated references
GRN benchmarking in single-cell data has repeatedly found modest precision-recall behavior and strong dependence on evaluation regime [7]. Core references (TRRUST, DoRothEA-derived resources) provide valuable priors but are incomplete and biased toward well-studied biology [3,4].

### 2.2 Positive-unlabeled learning and metric interpretation
Learning and evaluation with positive and unlabeled data is a classical setting in machine learning [8], including GRN inference specifically [9]. Our work applies this lens to reinterpret modern benchmark outputs rather than train a new classifier.

### 2.3 Interpretability validity and failure modes
Attribution and attention signals are often treated as explanations without sufficient causal grounding [11,12,16]. Probe accuracy alone can also be misleading without controls [13,14]. These lessons motivate our requirement that benchmark claims be tied to explicit evidence, uncertainty, and adversarial checks.

### 2.4 Mechanistic interpretability context
Circuit-level analyses and causal interventions have advanced interpretability in neural networks [17,18], and sparse feature decompositions have improved representation-level analysis [19,20]. We adapt this rigor mindset to GRN evaluation itself: before interpreting mechanisms, we must quantify what the benchmark can and cannot measure.

## 3. Methods
### 3.1 Setup and assumptions
Let `U` be candidate edges, `T subset U` latent true positives, `L+ subset T` observed positives, and `c = P(e in L+ | e in T)` the positive-label coverage. For predictions `P`:

- latent precision: `Prec_lat = |P intersect T| / |P|`
- latent recall: `Rec_lat = |P intersect T| / |T|`
- observed precision: `Prec_obs = |P intersect L+| / |P|`
- observed recall: `Rec_obs = |P intersect L+| / |L+|`

Under MAR positive labeling:
- `E[Prec_obs] = c * Prec_lat`
- `E[Rec_obs] = Rec_lat`

### 3.2 Ceiling formulas
For an ideal latent predictor (`Prec_lat = Rec_lat = 1`):
- `F1_obs,max(c) = 2c/(1+c)`

For latent positives `t = |T|` and ideal ranking:
- `AUPR_obs,max(c,t) ~= c + (1-c) * H_t/t`
where `H_t` is harmonic. For large `t`, `AUPR_obs,max ~= c`.

Observed random baseline:
- `AUPR_rand,obs = c * b`, where `b = |T|/|U|`.

Headroom-aware uplift:
- `rho_aupr = (AUPR_obs - AUPR_rand,obs) / (AUPR_obs,max - AUPR_rand,obs)`.

### 3.3 Data and evaluation artifacts
We analyze benchmark outputs spanning probe-based methods, classical GRN baselines, and controls across references aligned to transcriptional and signaling layers. The analysis includes central metric tables, coverage-proxy reports, cross-evaluation artifacts, threshold sweeps, and perturbation/time-series evaluations.

### 3.4 Coverage uncertainty and adversarial scenarios
Coverage is estimated with overlap-based proxies and modeled as `Beta(k+1, n-k+1)` by reference, with 95% posterior intervals. We evaluate:
- central overlap point,
- Beta interval bounds,
- adversarial misspecification (0.5x, 1.5x),
- mapping-policy bounds (mapping min/max),
- conservative and optimistic envelopes.

### 3.5 Additional submission-readiness analyses
To close reviewer-critical gaps, we add:
- stratified bootstrap uncertainty for headline metrics (3,000 resamples),
- reference-wise Wilson intervals for positive-uplift prevalence,
- method-family comparisons (probe/classical/prior-constrained/control),
- confound diagnostics against coverage, candidate-edge count, and observed base rate,
- external validation reinterpretation using time-series and perturbation artifacts.

## 4. Results
### 4.1 Central reinterpretation: low scores are not explained by missing labels alone
Across 39 AUPR-evaluable rows (15 methods, 5 references), best normalized F1 was 0.137 (`genie3`, DoRothEA-TRRUST immune union), best normalized AUPR was 0.0138 (`grnboost2`, same reference), and median normalized AUPR was 5.56e-05.

**Evidence:** 36/39 rows had `AUPR_obs < AUPR_rand,obs`; only 3/39 had positive uplift.

**Inference:** Most evaluated rankings do not enrich reference edges above observed random expectation after coverage adjustment.

**Hypothesis:** The dominant learned signal in many settings reflects broad co-expression structure more than directional regulation.

**Scientific implication:** Ceiling-aware normalization does not rescue most methods; progress requires better biological signal modeling, not only better label mapping.

### 4.2 Reference regime effects are biologically structured
Positive uplift was concentrated in transcription-focused immune references: 0.222 (95% Wilson CI 0.063-0.547) for DoRothEA-TRRUST immune union, versus 0 for hpn_dream and beeline_gsd in this evaluation set.

**Evidence:** Reference-level positive-uplift fractions differ sharply despite shared method pool.

**Inference:** Benchmark behavior is strongly reference-dependent and reflects biological layer mismatch (transcriptional TF-target versus pathway/signaling structure).

**Hypothesis:** Canonical TF programs are partially captured, but pathway-structured edges remain weakly modeled by current score families.

**Scientific implication:** Pooling references into a single leaderboard obscures biologically meaningful failure modes.

### 4.3 Robustness and adversarial stress tests
Coverage-scaling scenarios (overlap point, Beta bounds, 0.5x/1.5x) preserved rankings (Spearman approximately 1.0), but mapping-policy scenarios shifted rank by up to 10-11 places on overlapping subsets. MAR simulations showed low absolute median bias (F1: 7.30e-04; AUPR: 1.35e-03). Non-MAR labeling produced large AUPR shifts (+0.0898 optimistic; -0.0734 pessimistic relative to MAR theory).

**Evidence:** Rank stability depends on proxy family; MAR error is small, non-MAR error is large.

**Inference:** Curation mechanism is a first-order determinant of apparent performance, not just curation size.

**Hypothesis:** Literature curation over-samples canonical, high-visibility interactions and under-samples context-specific weak effects.

**Scientific implication:** Non-MAR stress tests should be mandatory for comparative claims in GRN benchmarking.

### 4.4 External validation aligns with central conclusions
In perturbation validation, overall recall was 1.50e-04 with non-zero recall in 19/4787 perturbagens (0.397%). Time-series reinterpretation showed prior-constrained families outperforming unconstrained/classical families, but still far below coverage-scaled ceilings (median HPN ratio to ceiling approximately 9.38e-04 for prior-constrained methods).

**Evidence:** Perturbation recall is sparse; time-series gains are concentrated in constrained prior settings.

**Inference:** Improvements arise in narrow regimes and do not indicate broad mechanistic recovery.

**Hypothesis:** Prior constraints inject biologically plausible structure that helps in specific assays, but upstream representation quality remains a limiting factor.

**Scientific implication:** Claims of general mechanistic success should be restricted to settings with consistent positive uplift across assay families.

## 5. Biological Interpretation
The most important biological takeaway is not that curated references are incomplete (they are), but that incompleteness is insufficient to explain the observed deficits. Even after normalization by partial-label ceilings, most methods remain close to baseline.

This pattern is biologically plausible. Transcriptional resources capture specific, often canonical TF programs [3,4], whereas pathway-oriented references encode multi-layer signaling context that expression-only or attribution-derived edge scores may not recover directly. The observed asymmetry between immune transcriptional references and high-coverage signaling-oriented references is consistent with that mismatch.

The perturbation analysis reinforces this interpretation: a handful of perturbagens show detectable recall, but broad recovery is absent. In practical terms, many inferred edges should currently be treated as exploratory hypotheses rather than mechanistic conclusions.

## 6. Discussion
### 6.1 What this changes for interpretation
This work reframes benchmark reading from absolute scores to *distance from observable headroom*. A method with low raw AUPR may still be near ceiling in a very sparse reference regime, and a method with seemingly improved raw scores may still be far from biologically meaningful recovery.

### 6.2 Integrity checks and claim calibration
We explicitly audited claim-evidence alignment, uncertainty, and confounds. The strongest residual dependency is between normalized ratios and observed base-rate structure, indicating that base-rate-aware controls remain necessary in future benchmarking.

### 6.3 What would change our mind
Our central conclusion would be weakened by multi-dataset evidence showing robust positive uplift across reference families and perturbation contexts, with stability under mapping-policy variation and non-MAR stress.

## 7. Limitations
Coverage `c` is proxy-estimated and not identifiable from benchmark outputs alone. Non-MAR simulations are stylized and do not uniquely recover real curation processes. External validations here are re-analyses of available artifacts rather than newly generated prospective experiments. Finally, our conclusions concern benchmark interpretability, not direct causal truth of any single inferred edge.

## 8. Conclusion
Partial-label ceiling correction is necessary but not sufficient for trustworthy GRN benchmark interpretation. In this dataset, missing labels explain part of the problem, but substantial model-to-biology mismatch remains after correction. The practical consequence is clear: mechanistic claims should be benchmarked against ceiling-normalized uplift, uncertainty bands, and non-MAR sensitivity, not raw metrics alone.

## 9. Data and Code Availability
All code, analysis scripts, figure-generation pipelines, and configuration files are provided at:

[https://anonymous.4open.science/r/grn-partial-label-metric-ceilings](https://anonymous.4open.science/r/grn-partial-label-metric-ceilings)

The repository README contains complete instructions to reproduce all main figures, tables, and supplementary analyses, along with dataset access notes.

## 10. Supplementary Materials
Supplementary items include:
- full scenario-wise method ranking tables,
- bootstrap and Wilson-interval diagnostics,
- confound-correlation diagnostics,
- synthetic MAR and non-MAR trial-level outputs,
- external time-series and perturbation reinterpretation tables.

## References
[1] Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods*. 2024;21:1470-1480. doi:10.1038/s41592-024-02201-0

[2] Tabula Sapiens Consortium. The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. *Science*. 2022;376(6594):eabl4896. doi:10.1126/science.abl4896

[3] Han H, Cho JW, Lee S, et al. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. *Nucleic Acids Research*. 2018;46(D1):D380-D386. doi:10.1093/nar/gkx1013

[4] Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*. 2019;29(8):1363-1375. doi:10.1101/gr.240663.118

[5] Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. Inferring regulatory networks from expression data using tree-based methods. *PLOS ONE*. 2010;5(9):e12776. doi:10.1371/journal.pone.0012776

[6] Moerman T, Aibar S, Gonzalez-Blas CB, et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks. *Bioinformatics*. 2019;35(12):2159-2161. doi:10.1093/bioinformatics/bty916

[7] Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. *Nature Methods*. 2020;17:147-154. doi:10.1038/s41592-019-0690-6

[8] Elkan C, Noto K. Learning classifiers from only positive and unlabeled data. In: *Proceedings of KDD*. 2008:213-220. doi:10.1145/1401890.1401920

[9] Cerulo L, Elkan C, Ceccarelli M. Learning gene regulatory networks from only positive and unlabeled data. *BMC Bioinformatics*. 2010;11:228. doi:10.1186/1471-2105-11-228

[10] Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLOS ONE*. 2015;10(3):e0118432. doi:10.1371/journal.pone.0118432

[11] Sundararajan M, Taly A, Yan Q. Axiomatic attribution for deep networks. In: *Proceedings of ICML (PMLR 70)*. 2017:3319-3328.

[12] Adebayo J, Gilmer J, Muelly M, Goodfellow I, Hardt M, Kim B. Sanity checks for saliency maps. In: *Advances in Neural Information Processing Systems 31 (NeurIPS)*. 2018.

[13] Jain S, Wallace BC. Attention is not Explanation. In: *Proceedings of NAACL-HLT*. 2019:3543-3556. doi:10.18653/v1/N19-1357

[14] Hewitt J, Liang P. Designing and interpreting probes with control tasks. In: *Proceedings of EMNLP-IJCNLP*. 2019:2733-2743. doi:10.18653/v1/D19-1275

[15] Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature*. 2023;618:616-624. doi:10.1038/s41586-023-06139-9

[16] Vig J, Gehrmann S, Belinkov Y, et al. Causal mediation analysis for interpreting neural NLP: the case of gender bias. *arXiv*. 2020. arXiv:2004.12265.

[17] Olah C, Cammarata N, Schubert L, et al. Zoom in: an introduction to circuits. *Distill*. 2020. doi:10.23915/distill.00024.001

[18] Meng K, Bau D, Andonian A, Belinkov Y. Locating and editing factual associations in GPT. *arXiv*. 2022. arXiv:2202.05262.

[19] Bricken T, Templeton A, Batson J, et al. Towards monosemanticity: decomposing language models with dictionary learning. *Transformer Circuits Thread*. 2023.

[20] Cunningham H, Ewart A, Riggs L, Huben R, Sharkey L. Sparse autoencoders find highly interpretable features in language models. *arXiv*. 2023. arXiv:2309.08600.

## D) Figure/Table audit list

### Figures
- **Figure 1 (Theoretical ceilings):**
  - Issue found: none critical.
  - Fix applied: retained log-scale baseline panel and increased DPI for publication readability.
- **Figure 2 (Observed vs normalized ratios):**
  - Issue found: potential overlap in low-value cluster.
  - Fix applied: log-x scaling and reference-color legend kept; marker edge outlines preserved for dense areas.
- **Figure 3 (Rank shift):**
  - Issue found: long method names can crowd x-axis.
  - Fix applied: rotated labels; zero reference line retained; recommend two-column split in camera-ready template if venue font constraints are tighter.
- **Figure 4 (Scenario rank heatmap):**
  - Issue found: missing-scenario cells could be misread.
  - Fix applied: explicit `NA` labels and neutral missing-color encoding.
- **Figure 5 (Sensitivity bands):**
  - Issue found: near-zero methods can visually collapse.
  - Fix applied: median point overlay and full horizontal range bars retained.
- **Figure 6 (MAR validation):**
  - Issue found: multi-line legend density.
  - Fix applied: compact legend and log-log axes retained for error-scale comparability.
- **Figure 7 (Non-MAR stress):**
  - Issue found: small F1 deltas may appear negligible beside AUPR shifts.
  - Fix applied: separate panels for AUPR and F1 deltas.
- **Figure 8 (Bootstrap uncertainty):**
  - Issue found: heterogeneous metric scales in one panel.
  - Fix applied: explicit CI whiskers and ordered labels to prevent misread.
- **Figure 9 (External validation):**
  - Issue found: perturbation panel could hide sparse long tail.
  - Fix applied: top-hit recall bars with hit-count annotations; log-log scatter for assay panel.
- **Figure 10 (Positive uplift prevalence):**
  - Issue found: uncertainty bars could exceed [0,1] bounds.
  - Fix applied: Wilson intervals clipped to valid bounds.

### Tables
- **Central reinterpretation tables:** `central_aupr_rows.csv`, `central_f1_rows.csv`
  - Issue found: none critical.
  - Fix applied: retained row-level reproducibility fields and derived metrics.
- **Scenario robustness tables:** `scenario_method_ranks.csv`, `scenario_rank_stability.csv`, `rowwise_coverage_sensitivity.csv`
  - Issue found: mapping scenarios had subset-support ambiguity.
  - Fix applied: explicit method counts and max-rank-shift reporting.
- **New submission-readiness tables:** `bootstrap_global_metrics.csv`, `uplift_positive_fraction_by_reference.csv`, `family_comparison_central.csv`, `confound_correlation_checks.csv`, `timeseries_family_summary.csv`, `perturbation_validation_summary.csv`
  - Issue found: cross-assay concordance originally produced misleading perfect correlations due tiny dynamic range.
  - Fix applied: degeneracy guard added; non-informative correlations now reported as `NaN` with explicit note (`degenerate_dynamic_range_not_interpretable`).

## E) Claims table

| Main claim | Support location | Strength |
|---|---|---|
| Most methods remain far from observable benchmark headroom after partial-label correction. | Results 4.1; central reinterpretation table; bootstrap table | **Strong** |
| Positive uplift is sparse and concentrated in specific references. | Results 4.1 and 4.2; uplift-prevalence table with Wilson intervals | **Strong** |
| Coverage proxy choice can materially change rankings in some settings. | Results 4.3; scenario rank stability and mapping-min/max shifts | **Strong** |
| MAR equations are accurate under MAR assumptions but non-MAR curation can strongly bias AUPR. | Results 4.3; synthetic MAR and non-MAR summaries | **Strong** |
| External evidence supports only narrow positive regimes, not broad mechanistic recovery. | Results 4.4; time-series family summary; perturbation validation summary | **Medium** |
| Prior-constrained methods likely help by injecting canonical biology structure. | Results 4.4 and Biological Interpretation | **Medium** |

## F) References list (complete)

[1] Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods*. 2024;21:1470-1480. doi:10.1038/s41592-024-02201-0

[2] Tabula Sapiens Consortium. The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans. *Science*. 2022;376(6594):eabl4896. doi:10.1126/science.abl4896

[3] Han H, Cho JW, Lee S, et al. TRRUST v2: an expanded reference database of human and mouse transcriptional regulatory interactions. *Nucleic Acids Research*. 2018;46(D1):D380-D386. doi:10.1093/nar/gkx1013

[4] Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*. 2019;29(8):1363-1375. doi:10.1101/gr.240663.118

[5] Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. Inferring regulatory networks from expression data using tree-based methods. *PLOS ONE*. 2010;5(9):e12776. doi:10.1371/journal.pone.0012776

[6] Moerman T, Aibar S, Gonzalez-Blas CB, et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks. *Bioinformatics*. 2019;35(12):2159-2161. doi:10.1093/bioinformatics/bty916

[7] Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data. *Nature Methods*. 2020;17:147-154. doi:10.1038/s41592-019-0690-6

[8] Elkan C, Noto K. Learning classifiers from only positive and unlabeled data. In: *Proceedings of KDD*. 2008:213-220. doi:10.1145/1401890.1401920

[9] Cerulo L, Elkan C, Ceccarelli M. Learning gene regulatory networks from only positive and unlabeled data. *BMC Bioinformatics*. 2010;11:228. doi:10.1186/1471-2105-11-228

[10] Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLOS ONE*. 2015;10(3):e0118432. doi:10.1371/journal.pone.0118432

[11] Sundararajan M, Taly A, Yan Q. Axiomatic attribution for deep networks. In: *Proceedings of ICML (PMLR 70)*. 2017:3319-3328.

[12] Adebayo J, Gilmer J, Muelly M, Goodfellow I, Hardt M, Kim B. Sanity checks for saliency maps. In: *Advances in Neural Information Processing Systems 31 (NeurIPS)*. 2018.

[13] Jain S, Wallace BC. Attention is not Explanation. In: *Proceedings of NAACL-HLT*. 2019:3543-3556. doi:10.18653/v1/N19-1357

[14] Hewitt J, Liang P. Designing and interpreting probes with control tasks. In: *Proceedings of EMNLP-IJCNLP*. 2019:2733-2743. doi:10.18653/v1/D19-1275

[15] Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions in network biology. *Nature*. 2023;618:616-624. doi:10.1038/s41586-023-06139-9

[16] Vig J, Gehrmann S, Belinkov Y, et al. Causal mediation analysis for interpreting neural NLP: the case of gender bias. *arXiv*. 2020. arXiv:2004.12265.

[17] Olah C, Cammarata N, Schubert L, et al. Zoom in: an introduction to circuits. *Distill*. 2020. doi:10.23915/distill.00024.001

[18] Meng K, Bau D, Andonian A, Belinkov Y. Locating and editing factual associations in GPT. *arXiv*. 2022. arXiv:2202.05262.

[19] Bricken T, Templeton A, Batson J, et al. Towards monosemanticity: decomposing language models with dictionary learning. *Transformer Circuits Thread*. 2023.

[20] Cunningham H, Ewart A, Riggs L, Huben R, Sharkey L. Sparse autoencoders find highly interpretable features in language models. *arXiv*. 2023. arXiv:2309.08600.
