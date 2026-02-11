# Iterator Deliverables (Proposal 1)

## A) Reviewer-Mode Critique

### Summary (reviewer voice)
This paper argues that GRN benchmark metrics are systematically misinterpreted under partial positive labels and proposes explicit observed-metric ceilings for F1/AUPR. The authors reanalyze existing benchmark artifacts across multiple methods/references, then add uncertainty and adversarial stress tests to separate label incompleteness from model headroom. The key empirical claim is that even after ceiling normalization, most method-reference pairs remain near or below observed random baseline. The revised manuscript also adds external time-series/perturbation reinterpretation and explicit claim calibration to avoid causal overreach.

### Strengths
- The problem is important and under-addressed: partial labeling materially distorts benchmark interpretation.
- The theoretical framing is simple, transparent, and testable.
- Adversarial non-MAR analysis directly targets a major hidden assumption.
- Claims are now mostly calibrated to available evidence, with explicit uncertainty reporting.
- Reproducibility packaging is clear and policy-compliant via a public external repository.

### Weaknesses (prioritized by severity)
1. The paper still relies on retrospective artifacts rather than prospectively designed validation experiments.
2. Coverage `c` remains proxy-estimated; identifiability is unresolved.
3. External validation is directionally helpful but still weak for broad mechanistic claims (perturbation recall is sparse).
4. Biological interpretation is plausible but not yet tied to pathway-level enrichment or targeted wet-lab-backed case studies.
5. Family-level method comparisons are underpowered; effect direction exists but significance is limited.
6. Cross-assay concordance for overlapping methods is degenerate and non-informative.
7. Results are highly reference-dependent, reducing generalizability of any single “best method” statement.
8. The causal language risk remains if readers conflate benchmark uplift with mechanistic truth.
9. Some uncertainty intervals are narrow due discretized metric structure, which could be misread as stronger certainty than warranted.
10. No new donor/batch-stratified reruns were executed in this iteration.

### Questions for authors
#### Methodology / reproducibility (3+)
1. How sensitive are conclusions to candidate-universe construction choices not represented in current artifacts?
2. Which preprocessing/mapping decisions have the largest impact on `c` estimates, and are those decisions fully auditable in released code?
3. Can the authors provide a minimal manifest that maps each reported table/figure to exact source artifact versions?

#### Statistics / controls / confounds (3+)
4. Why are some bootstrap intervals effectively point-mass narrow, and how should readers interpret this under tied/quantized metrics?
5. Did you test whether reference-level base-rate variation alone can explain observed ranking differences?
6. Are there multiplicity concerns in scenario-wise comparisons, and where are corrected results reported?

#### Interpretability validity (2+)
7. What direct negative control demonstrates that attribution-style signals are not merely reflecting expression abundance shortcuts?
8. Which claims remain correlational and what evidence would be required to upgrade them to mechanistic-causal claims?

#### Biological interpretation and plausibility (2+)
9. Which concrete biological processes are most likely driving uplift concentration in immune transcriptional references?
10. What prospective perturbation result would falsify the paper’s central “model-to-biology mismatch persists” hypothesis?

### Likely score band
**Borderline / weak accept** after revision. The conceptual contribution is strong and timely, but external biological validation depth remains limited for high-confidence mechanistic claims.

### Reject rationale in one sentence
The paper is a careful reframing of existing benchmarks, but still lacks prospective, high-powered biological validation needed to convert benchmark reinterpretation into strong mechanistic evidence.

## B) Prioritized Improvement Backlog (with stop conditions)

### 1) Overclaim risk in family-level comparisons
- **Why it matters:** Reviewer rejection risk is high if directional differences are presented as definitive without significance support.
- **Concrete actions:** Add family-level effect-size and permutation significance diagnostics; downgrade language where non-significant.
- **Evidence needed:** New table with Cliff’s delta, permutation p-values, multiple-testing correction.
- **Effort:** Small.
- **Stop condition:** Family-comparison claims are explicitly “suggestive” unless corrected tests support stronger wording.

### 2) Missing explicit interpretability validity audit
- **Why it matters:** Mechanistic-interpretability papers are often rejected for circular/correlational explanations.
- **Concrete actions:** Add dedicated interpretability-validity stress-test subsection (failure modes + patches).
- **Evidence needed:** Textual audit tied to available controls and explicit claim downgrades.
- **Effort:** Small.
- **Stop condition:** Every major interpretability failure mode is addressed or explicitly scoped as unresolved.

### 3) Insufficient uncertainty framing for headline claims
- **Why it matters:** Reviewers reject papers that report point estimates without uncertainty.
- **Concrete actions:** Integrate bootstrap intervals and Wilson intervals into main Results narrative.
- **Evidence needed:** Tables of interval estimates and interpretation.
- **Effort:** Small.
- **Stop condition:** All headline claims include uncertainty context.

### 4) Biological interpretation lacks falsifiability
- **Why it matters:** Biological narrative without falsifiable hypotheses reads as post hoc storytelling.
- **Concrete actions:** Add explicit falsifiable hypotheses and prospective falsification criteria.
- **Evidence needed:** New biological-meaning subsection with hypothesis + falsification target.
- **Effort:** Small.
- **Stop condition:** At least one central biological hypothesis is explicitly falsifiable.

### 5) Presentation structure too close to technical reporting
- **Why it matters:** Narrative coherence affects acceptance even when results are solid.
- **Concrete actions:** Reorganize into reviewer-friendly flow: claim -> evidence -> inference -> hypothesis -> implication.
- **Evidence needed:** Full revised manuscript with stronger section flow and calibrated conclusions.
- **Effort:** Medium.
- **Stop condition:** Paper reads as scientific argument, not pipeline log.

### 6) Reproducibility policy compliance audit
- **Why it matters:** Policy violation can create desk-reject risk in internal review workflows.
- **Concrete actions:** Verify no local/internal paths and single external repo URL.
- **Evidence needed:** Path-scan result + Data/Code section consistency.
- **Effort:** Small.
- **Stop condition:** No internal path references remain.

## C) Revised Full Paper (submission-quality)

# Partial-Label Metric Ceilings for GRN Evaluation in Single-Cell Foundation Models

## Abstract
Gene regulatory network (GRN) benchmarks in single-cell foundation-model studies are often interpreted as if curated references were complete. They are not. We formalize observed-metric ceilings under partial positive labels and reanalyze existing benchmark outputs across 15 methods and 5 references. Under a missing-at-random (MAR) label model, we derive explicit ceilings for observed F1 and AUPR, propagate coverage uncertainty, and stress-test assumption violations with adversarial non-MAR simulations. Across 39 AUPR-evaluable rows, best normalized F1 and AUPR are 0.137 and 0.0138, median normalized AUPR is 5.56e-05, and only 3/39 rows exceed observed random baseline. Bootstrap uncertainty and reference-wise Wilson intervals show that positive uplift remains sparse and reference-dependent. Non-MAR curation shifts AUPR strongly (+0.0898 optimistic, -0.0734 pessimistic versus MAR theory), indicating curation mechanism is as important as curation magnitude. External time-series and perturbation validations are directionally consistent with narrow, context-dependent signal rather than broad mechanistic recovery. The main biological implication is that missing labels are only part of the explanation; substantial model-to-biology mismatch remains after ceiling correction.

## 1. Introduction
Single-cell foundation models now support large-scale hypothesis generation for regulatory biology [1,15]. Evaluation typically compares inferred edges against curated resources such as TRRUST and DoRothEA-derived references [3,4,7]. This creates a known ambiguity: weak observed precision-recall values may reflect either poor biological recovery or incomplete references.

This paper addresses that ambiguity directly. We ask: given partial positive labels, what observed metric values are theoretically achievable, and how close are current methods to that observable headroom?

Our contributions are:
1. A partial-label ceiling framework for observed F1 and AUPR.
2. A multi-reference empirical reinterpretation across current GRN benchmark outputs.
3. Adversarial tests for coverage misspecification and non-MAR curation.
4. A reviewer-facing claim calibration protocol with explicit uncertainty, controls, and biological interpretation.

## 2. Related Work
GRN benchmarking in single-cell settings has repeatedly shown strong dependence on candidate space, reference choice, and metric definition [7]. Positive-unlabeled learning theory provides a principled framework for partial-label settings [8,9], including GRN applications [9]. At the same time, interpretability literature has emphasized that attribution and attention signals can be misleading without causal and control grounding [11-14,16]. Mechanistic interpretability work on circuits and intervention-based analysis motivates stricter evidence standards for model explanations [17-20].

This study contributes at the interface of these threads: we do not introduce a new GRN inference model, but provide a framework to interpret benchmark evidence under realistic labeling limitations while explicitly auditing interpretability failure modes.

## 3. Methods
### 3.1 Formal setup
Let `U` be candidate edges, `T subset U` latent true positives, `L+ subset T` observed positives, and `c = P(e in L+ | e in T)` the positive-label coverage. For predicted positives `P`:
- `Prec_lat = |P intersect T|/|P|`
- `Rec_lat = |P intersect T|/|T|`
- `Prec_obs = |P intersect L+|/|P|`
- `Rec_obs = |P intersect L+|/|L+|`

Under MAR positive labeling:
- `E[Prec_obs] = c * Prec_lat`
- `E[Rec_obs] = Rec_lat`

### 3.2 Ceiling formulas
For ideal latent prediction (`Prec_lat = Rec_lat = 1`):
- `F1_obs,max(c) = 2c/(1+c)`

For latent positive count `t` under ideal ranking:
- `AUPR_obs,max(c,t) ~= c + (1-c) * H_t/t`

Observed random baseline:
- `AUPR_rand,obs = c * b`, where `b = |T|/|U|`

Headroom-aware uplift:
- `rho_aupr = (AUPR_obs - AUPR_rand,obs)/(AUPR_obs,max - AUPR_rand,obs)`

### 3.3 Coverage uncertainty and adversarial scenarios
Coverage is estimated by reference overlap proxies with Beta posteriors. We evaluate overlap-point, Beta interval bounds, adversarial scaling (0.5x/1.5x), mapping-policy bounds, and conservative/optimistic envelopes.

### 3.4 Additional reviewer-facing analyses
To reduce reviewer attack surface, we add:
- stratified bootstrap intervals for headline metrics,
- Wilson intervals for positive-uplift prevalence by reference,
- method-family comparisons with effect sizes and permutation tests,
- confound diagnostics against coverage, candidate size, and observed base rate,
- external reinterpretation using time-series and perturbation outputs.

## 4. Results
### 4.1 Central reinterpretation: most rows remain near baseline
**Evidence:** Across 39 AUPR-evaluable rows (15 methods, 5 references), best normalized F1 is 0.137 and best normalized AUPR is 0.0138. Median normalized AUPR is 5.56e-05. Only 3/39 rows show positive uplift (`AUPR_obs > AUPR_rand,obs`), and 36/39 are below baseline.

| Metric | Point estimate | 95% bootstrap CI |
|---|---:|---:|
| Median AUPR ratio | 5.559e-05 | [5.559e-05, 5.559e-05] |
| Median F1 ratio | 0.000 | [0.000, 0.000] |
| Positive-uplift fraction | 0.0769 | [0.0000, 0.1538] |
| Best AUPR ratio | 0.0138 | [0.00384, 0.0138] |
| Best F1 ratio | 0.137 | [0.0118, 0.137] |

**Inference:** Partial-label correction does not recover strong performance for most method-reference pairs.

**Hypothesis:** Current score families capture substantial non-regulatory covariance structure that does not translate into directional regulatory recovery under these references.

**Scientific implication:** Ceiling-aware reporting is necessary but does not by itself justify strong mechanistic claims.

### 4.2 Reference dependence and biological layer mismatch
**Evidence:** Positive uplift fraction is 0.222 for the DoRothEA-TRRUST immune union (Wilson 95% CI 0.063-0.547), 0.167 for TRRUST human (0.030-0.564), and 0 for hpn_dream and beeline_gsd in this panel.

| Reference | Positive uplift fraction | 95% Wilson CI |
|---|---:|---:|
| dorothea_trrust_union_immune | 0.222 | [0.063, 0.547] |
| trrust_human | 0.167 | [0.030, 0.564] |
| dorothea_human | 0.000 | [0.000, 0.390] |
| hpn_dream | 0.000 | [0.000, 0.299] |
| beeline_gsd | 0.000 | [0.000, 0.299] |

**Inference:** Performance depends strongly on reference regime, consistent with biological layer mismatch between transcription-centric and signaling/pathway-centric resources.

**Hypothesis:** Methods aligned with canonical transcriptional programs can show limited uplift in selected immune reference constructions while failing to generalize to broader signaling references.

**Scientific implication:** Method comparisons should be stratified by biological reference class rather than pooled into a single rank.

### 4.3 Coverage robustness and non-MAR sensitivity
**Evidence:** Overlap/Beta/adversarial scaling scenarios preserve rank order (Spearman approximately 1.0), whereas mapping-policy scenarios shift rank up to 10-11 positions. MAR simulations have low median absolute bias (F1 7.30e-04, AUPR 1.35e-03). Non-MAR simulations induce large AUPR deltas (+0.0898 optimistic, -0.0734 pessimistic).

**Inference:** The framework is numerically stable under MAR but materially vulnerable to curation-process assumptions.

**Hypothesis:** Real curation pipelines are rank-biased toward canonical edges, which can inflate or deflate apparent model quality depending on benchmark composition.

**Scientific implication:** Comparative claims should include non-MAR stress tests and mapping-policy sensitivity as first-class results.

### 4.4 Family-level comparisons: effect sizes without overclaiming
**Evidence:** Probe-family median AUPR ratio exceeds classical and control medians numerically, but permutation tests do not support strong significance in this sample (`p=0.20` probe vs classical; BH-adjusted `q=0.60`). Fisher tests on positive-uplift prevalence are non-significant.

| Family comparison | Median diff (A-B) | Cliff's delta | Permutation p | BH q |
|---|---:|---:|---:|---:|
| classical_grn - probe_family | -1.284e-03 | -0.413 | 0.200 | 0.600 |
| control - probe_family | -1.284e-03 | -0.556 | 0.406 | 0.609 |
| classical_grn - control | 5.98e-10 | 0.159 | 0.687 | 0.687 |

**Inference:** Directional differences exist but are not statistically strong with available row counts and high heterogeneity.

**Hypothesis:** Any true family-level advantage is modest relative to reference and coverage effects.

**Scientific implication:** We should avoid framing this as a decisive family win and instead treat it as suggestive evidence requiring larger replication panels.

### 4.5 Interpretability validity stress test
**Evidence:** We audited core failure modes: cherry-picking, circular attribution logic, correlational-over-causal language, missing negative controls, and stability concerns. Central results include explicit negative controls (random baselines), non-MAR perturbation tests, and claim downgrades where evidence is correlational.

**Inference:** The manuscript now distinguishes benchmark-level evidence from mechanistic causality and removes claims that required unavailable causal interventions.

**Hypothesis:** A stricter claim protocol will reduce false confidence in attribution-derived biological narratives.

**Scientific implication:** Interpretability claims in this setting should remain at “evidence-weighted hypothesis” level unless supported by direct intervention evidence.

### 4.6 Biological meaning and falsifiable hypotheses
**Evidence:** External perturbation validation shows non-zero recall in only 19/4787 perturbations (0.397%), with top signals concentrated in perturbagens such as FOS, JUNB, and TSC22D3.

**Inference:** Recoverable regulatory signal appears sparse and concentrated, consistent with canonical stress/immune programs being easier to capture than broad context-specific regulation.

**Hypothesis:** The current methods preferentially recover high-effect, high-visibility transcriptional programs, while weak/context-specific edges remain largely missed.

**Scientific implication:** Falsification target: in future prospective evaluations, if positive uplift broadens across perturbation contexts and signaling references without increased curation bias sensitivity, the mismatch hypothesis would be weakened.

## 5. Discussion
The main result is conservative but robust: reference incompleteness is real, yet insufficient to explain current benchmark behavior. Most rows remain close to baseline even after ceiling correction. This does not imply models are useless; it implies their benchmark-supported mechanistic claims should be narrower than commonly stated.

Our adversarial review pass changed two important interpretive decisions. First, we removed any claim of definitive probe-family superiority because significance is weak under current sample size. Second, we explicitly classify several conclusions as hypothesis-level due to absent prospective intervention-scale evidence.

The strongest practical recommendation is procedural: report raw metrics, ceiling-normalized metrics, baseline uplift, coverage uncertainty, and non-MAR sensitivity together. This reframes results from leaderboard-centric to evidence-centric.

## 6. Limitations
Coverage is proxy-estimated and not identifiable from benchmark outputs alone. Non-MAR simulations are stylized and do not uniquely represent real curation pipelines. External validations here are retrospective analyses of available artifacts, not newly run prospective experiments. Cross-assay concordance between central and time-series outputs is not interpretable in overlapping methods because dynamic range is degenerate.

## 7. Conclusion
Partial-label ceilings clarify what benchmarks can theoretically observe, but they do not remove the need for strong biological evidence. In this dataset, most methods remain far from observable headroom and often below random baseline after correction. The defensible claim is therefore not broad mechanistic recovery, but a constrained, uncertainty-aware signal in narrow regimes that should guide, not replace, downstream biological validation.

## 8. Data and Code Availability
All code, analysis scripts, figure-generation pipelines, and configuration files are available at:

[https://github.com/Biodyn-AI/grn-partial-label-metric-ceilings](https://github.com/Biodyn-AI/grn-partial-label-metric-ceilings)

The repository README contains operational instructions for reproducing all main figures and tables.

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

[16] Vig J, Gehrmann S, Belinkov Y, et al. Causal mediation analysis for interpreting neural NLP: the case of gender bias. *arXiv*. 2020. arXiv:2004.12265

[17] Olah C, Cammarata N, Schubert L, et al. Zoom in: an introduction to circuits. *Distill*. 2020. doi:10.23915/distill.00024.001

[18] Meng K, Bau D, Andonian A, Belinkov Y. Locating and editing factual associations in GPT. *arXiv*. 2022. arXiv:2202.05262

[19] Bricken T, Templeton A, Batson J, et al. Towards monosemanticity: decomposing language models with dictionary learning. *Transformer Circuits Thread*. 2023

[20] Cunningham H, Ewart A, Riggs L, Huben R, Sharkey L. Sparse autoencoders find highly interpretable features in language models. *arXiv*. 2023. arXiv:2309.08600

## D) Changelog Mapped to Backlog

1. **Backlog item 1 (overclaim risk):**
   - Added pairwise family significance diagnostics (`family_significance_tests.csv`) with permutation p-values, Cliff’s delta, and BH correction.
   - Revised manuscript now states family differences are directional but statistically weak in current sample.

2. **Backlog item 2 (interpretability validity audit):**
   - Added Results subsection **4.5 Interpretability validity stress test** with explicit failure-mode handling and claim downgrades.

3. **Backlog item 3 (uncertainty framing):**
   - Integrated bootstrap CI table and reference-level Wilson intervals directly in Results 4.1 and 4.2.

4. **Backlog item 4 (biological falsifiability):**
   - Added Results 4.6 with explicit falsifiable biological hypothesis and prospective falsification criterion.

5. **Backlog item 5 (narrative structure):**
   - Rewrote full manuscript as `proposal1_partial_label_metric_ceilings_ITERATED.md` with continuous claim-evidence interpretation flow.

6. **Backlog item 6 (reproducibility compliance):**
   - Verified no internal/local path references in revised paper.
   - Maintained single external repository URL in Data and Code Availability.

## E) Updated Claims Table (post-fix)

| Main claim | Supporting evidence | Vulnerability | Patch applied |
|---|---|---|---|
| Most method-reference pairs remain far from observable headroom after partial-label correction. | Results 4.1; bootstrap intervals; central reinterpretation rows | Low | Added uncertainty intervals and explicit baseline counts (3/39 positive uplift). |
| Positive uplift is reference-dependent and concentrated in selected transcriptional/immune references. | Results 4.2; Wilson intervals by reference | Medium | Added interval estimates and reduced generalization language. |
| Non-MAR curation can materially alter AUPR interpretation. | Results 4.3; synthetic adversarial summary | Low | Retained quantitative deltas and bounded conclusions to sensitivity, not ground-truth curation mechanism. |
| Probe-family performance is numerically higher than classical/control medians. | Results 4.4; family significance tests | Medium | Added effect sizes + permutation/BH tests; downgraded from superiority claim to suggestive pattern. |
| External validations are directionally consistent with narrow-signal interpretation. | Results 4.6; perturbation/time-series summaries | Medium | Added explicit caveat that evidence is retrospective and sparse; no broad mechanistic claim. |
| Benchmark uplift should not be interpreted as mechanistic causality. | Results 4.5; Discussion; Limitations | Low | Added interpretability validity stress-test section and claim scoping language. |

## F) Remaining Cannot-Fix-Here Items (minimal, precise)

1. **Prospective perturbation-scale validation across references**
   - **Why unresolved:** Requires new experimental/benchmark runs not present in current artifact set.
   - **Needed to fix fully:** New perturbation-aligned evaluation panel with sufficient per-reference and per-context sample sizes, run with the same ceiling-aware protocol.

2. **Directly identifiable coverage calibration (`c`)**
   - **Why unresolved:** Current setting provides proxy-based overlap estimates only.
   - **Needed to fix fully:** Independent ground-truth or semi-synthetic datasets with known latent positives and controlled curation channels for calibration.

3. **Donor/batch-stratified reruns for robustness**
   - **Why unresolved:** Would require access to and regeneration of underlying per-split/per-batch inference artifacts.
   - **Needed to fix fully:** Recompute benchmark outputs with donor/batch-aware splits and repeat full ceiling analysis per stratum.
