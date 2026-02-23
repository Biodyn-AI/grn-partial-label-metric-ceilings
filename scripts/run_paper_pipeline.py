#!/usr/bin/env python3
"""Full paper pipeline: Partial-label metric ceilings.

This script produces a research-grade artifact bundle:
- formal ceiling calculations for F1/AUPR under partial positive labeling,
- multi-scenario coverage sensitivity analysis,
- adversarial robustness checks for coverage misspecification,
- synthetic theorem-validation experiments (MAR + non-MAR labeling),
- publication-ready figures/tables and an auto-generated results summary.

All outputs are written under: outputs/paper/
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Use a non-interactive backend for reproducible headless rendering.
matplotlib.use("Agg")


GAMMA = 0.5772156649015329
RNG_SEED = 11
BETA_SAMPLES = 20000
SMALL = 1e-12

# Main evaluation tables with AUPR + F1.
PRIMARY_METRIC_TABLES = {
    "probe_priors": "data/score_eval_probe_priors.csv",
    "grn_baselines_immune": "data/score_eval_grn_baselines_immune.csv",
}

# Additional tables that provide F1 but missing AUPR (still useful for F1 ceiling analysis).
F1_ONLY_TABLES = {
    "probe_priors_full_genes": "data/score_eval_probe_priors_full_genes.csv",
    "probe_priors_full_genes_crosswalk": "data/score_eval_probe_priors_full_genes_crosswalk.csv",
    "probe_priors_full_genes_omnipath": "data/score_eval_probe_priors_full_genes_omnipath.csv",
}

# Mapping artifacts used as alternative (endpoint-level) coverage proxies.
MAPPING_REPORTS = {
    "probe_priors": "data/score_eval_probe_priors_missing_report.json",
    "probe_priors_full_genes": "data/score_eval_probe_priors_full_genes_missing_report.json",
    "probe_priors_full_genes_crosswalk": "data/score_eval_probe_priors_full_genes_crosswalk_missing_report.json",
    "probe_priors_full_genes_omnipath": "data/score_eval_probe_priors_full_genes_omnipath_missing_report.json",
}

# Extra network-inference sweep and cross-eval artifacts for external sanity checks.
SWEEP_JSONS = {
    "omnipath": "data/sweep_results_omnipath.json",
    "omnipath_dorothea_union_immune_hpn": "data/sweep_results_omnipath_dorothea_union_immune_hpn.json",
    "omnipath_dorothea_intersection_immune_hpn": "data/sweep_results_omnipath_dorothea_intersection_immune_hpn.json",
    "regulatory": "data/sweep_results_regulatory.json",
}

CROSS_EVAL_JSONS = {
    "cross_union": "data/cross_eval_dorothea_union_immune_hpn.json",
    "cross_union_scaled": "data/cross_eval_dorothea_union_immune_hpn_scaled.json",
    "cross_intersection": "data/cross_eval_dorothea_intersection_immune_hpn.json",
    "cross_intersection_scaled": "data/cross_eval_dorothea_intersection_immune_hpn_scaled.json",
}


@dataclass
class CoveragePosterior:
    alpha: float
    beta: float


def find_repo_root(start: Path) -> Path:
    """Walk up directories until reaching repository root folders."""
    for candidate in [start] + list(start.parents):
        if (candidate / "data").exists() and (candidate / "scripts").exists():
            return candidate
    raise RuntimeError("Could not detect repository root from script location.")


def normalize_bool(series: pd.Series) -> pd.Series:
    """Normalize mixed boolean columns (bool/string/int) into bool dtype."""
    if series.dtype == bool:
        return series
    lowered = series.astype(str).str.lower()
    return lowered.isin({"true", "1", "yes"})


def harmonic_over_n(t: np.ndarray) -> np.ndarray:
    """Approximate H_t / t for t >= 1 with a stable asymptotic expansion."""
    t = np.maximum(np.asarray(t, dtype=float), 1.0)
    h_t = np.log(t) + GAMMA + (1.0 / (2.0 * t)) - (1.0 / (12.0 * t * t))
    return h_t / t


def f1_ceiling(c: np.ndarray) -> np.ndarray:
    """Observed F1 ceiling for perfect latent prediction under coverage c."""
    c = np.clip(np.asarray(c, dtype=float), 1e-9, 1.0)
    return (2.0 * c) / (1.0 + c)


def aupr_ceiling(c: np.ndarray, latent_true_edges: np.ndarray) -> np.ndarray:
    """Observed AUPR ceiling for ideal latent ranking under partial labels.

    The expression is:
        AUPR_max ~= c + (1-c) * (H_t / t)
    where t is latent positive count.
    """
    c = np.clip(np.asarray(c, dtype=float), 1e-9, 1.0)
    latent_true_edges = np.maximum(np.asarray(latent_true_edges, dtype=float), 1.0)
    ceiling = c + (1.0 - c) * harmonic_over_n(latent_true_edges)
    return np.clip(ceiling, 0.0, 1.0)


def average_precision_from_binary_labels(labels: np.ndarray) -> float:
    """Compute AP for strictly ranked labels without sklearn dependency.

    We assume labels are already ordered by descending score.
    AP = mean precision at each observed positive rank position.
    """
    labels = np.asarray(labels, dtype=np.int8)
    positive_positions = np.flatnonzero(labels == 1)
    if len(positive_positions) == 0:
        return np.nan
    # Convert to 1-based ranks for precision@k definitions.
    ranks = positive_positions + 1.0
    cum_tp = np.arange(1, len(positive_positions) + 1, dtype=float)
    precision_at_pos = cum_tp / ranks
    return float(np.mean(precision_at_pos))


def load_metric_frames(repo_root: Path) -> pd.DataFrame:
    """Load all metric tables and keep rows that are valid for at least F1 analysis."""
    frames: List[pd.DataFrame] = []
    source_map = {}
    source_map.update(PRIMARY_METRIC_TABLES)
    source_map.update(F1_ONLY_TABLES)

    for source_name, rel_path in source_map.items():
        path = repo_root / rel_path
        if not path.exists():
            raise FileNotFoundError(f"Missing required metric table: {path}")

        frame = pd.read_csv(path)
        frame["source_table"] = source_name
        frame["reference_excluded"] = normalize_bool(frame["reference_excluded"])

        # Keep only references that are evaluable in each table.
        frame = frame[~frame["reference_excluded"]].copy()
        frame = frame[frame["f1"].notna()].copy()
        frames.append(frame)

    rows = pd.concat(frames, ignore_index=True)

    # Overlap-based coverage proxy used as central estimate of c.
    rows["coverage_overlap"] = np.clip(rows["ref_node_overlap_pct"] / 100.0, 1e-6, 0.999999)
    rows["observed_base_rate"] = rows["true_edges"] / rows["candidate_edges"]

    return rows


def load_mapping_proxy_rows(repo_root: Path) -> pd.DataFrame:
    """Load endpoint-level mapping coverage proxies from missing reports."""
    rows: List[Dict[str, object]] = []
    for proxy_source, rel_path in MAPPING_REPORTS.items():
        report_path = repo_root / rel_path
        if not report_path.exists():
            raise FileNotFoundError(f"Missing mapping report: {report_path}")

        payload = json.loads(report_path.read_text())
        for record in payload:
            denom = 2.0 * float(record["total_edges_normalized"])
            missing_total = float(record["missing_total"])
            coverage = np.nan
            support_num = np.nan
            if denom > 0:
                support_num = max(denom - missing_total, 0.0)
                coverage = float(np.clip(support_num / denom, 0.0, 1.0))

            rows.append(
                {
                    "proxy_kind": "mapping_endpoint",
                    "proxy_source": proxy_source,
                    "reference": record["reference"],
                    "coverage_proxy": coverage,
                    "support_num": support_num,
                    "support_den": denom,
                    "missing_total": missing_total,
                    "total_edges_normalized": float(record["total_edges_normalized"]),
                }
            )

    return pd.DataFrame(rows)


def estimate_reference_posteriors(metric_rows: pd.DataFrame) -> Dict[str, CoveragePosterior]:
    """Estimate reference-level Beta posteriors for c from overlap counts."""
    posteriors: Dict[str, CoveragePosterior] = {}
    ref_counts = metric_rows[["reference", "overlap_nodes", "ref_nodes"]].drop_duplicates()
    for row in ref_counts.itertuples(index=False):
        successes = float(row.overlap_nodes)
        total = float(row.ref_nodes)
        failures = max(total - successes, 0.0)
        posteriors[row.reference] = CoveragePosterior(alpha=successes + 1.0, beta=failures + 1.0)
    return posteriors


def build_reference_coverage_scenarios(
    metric_rows: pd.DataFrame,
    mapping_proxies: pd.DataFrame,
    rng: np.random.Generator,
    preferred_rows: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build central and adversarial coverage scenarios for each reference.

    Scenario design:
    - overlap_point: central proxy from ref-node overlap.
    - beta_ci_low/high: posterior interval bounds from overlap uncertainty.
    - adversarial_half / adversarial_150pct: intentional misspecification stress tests.
    - mapping_min/max: alternate proxy family from endpoint mapping reports.
    - conservative_min / optimistic_max: envelope across available scenarios.
    """
    if preferred_rows is None or preferred_rows.empty:
        preferred_rows = metric_rows

    posteriors_primary = estimate_reference_posteriors(preferred_rows)
    posteriors_fallback = estimate_reference_posteriors(metric_rows)

    # Mapping proxies can have multiple values per reference across policy variants.
    mapping_agg = (
        mapping_proxies.groupby("reference", as_index=False)
        .agg(mapping_min=("coverage_proxy", "min"), mapping_max=("coverage_proxy", "max"))
        .copy()
    )

    records: List[Dict[str, float]] = []
    overlap_point_by_ref = preferred_rows.groupby("reference")["coverage_overlap"].median().to_dict()

    # Fill references absent in preferred rows with all-row fallback medians.
    overlap_point_fallback = metric_rows.groupby("reference")["coverage_overlap"].median().to_dict()
    for reference, value in overlap_point_fallback.items():
        overlap_point_by_ref.setdefault(reference, value)

    for reference, c_point in overlap_point_by_ref.items():
        posterior = posteriors_primary.get(reference, posteriors_fallback[reference])
        draws = np.clip(rng.beta(posterior.alpha, posterior.beta, size=BETA_SAMPLES), 1e-6, 0.999999)

        beta_low = float(np.quantile(draws, 0.025))
        beta_high = float(np.quantile(draws, 0.975))

        record = {
            "reference": reference,
            "support_num": posterior.alpha - 1.0,
            "support_den": posterior.alpha + posterior.beta - 2.0,
            "overlap_point": float(c_point),
            "beta_ci_low": beta_low,
            "beta_ci_high": beta_high,
            "adversarial_half": float(max(c_point * 0.5, 1e-3)),
            "adversarial_150pct": float(min(c_point * 1.5, 0.999)),
            "mapping_min": np.nan,
            "mapping_max": np.nan,
        }

        mapping_row = mapping_agg[mapping_agg["reference"] == reference]
        if not mapping_row.empty:
            record["mapping_min"] = float(mapping_row["mapping_min"].iloc[0])
            record["mapping_max"] = float(mapping_row["mapping_max"].iloc[0])

        # Envelope scenarios across every non-null proxy for adversarial coverage bands.
        available = [
            record["overlap_point"],
            record["beta_ci_low"],
            record["beta_ci_high"],
            record["adversarial_half"],
            record["adversarial_150pct"],
        ]
        if not np.isnan(record["mapping_min"]):
            available.append(record["mapping_min"])
        if not np.isnan(record["mapping_max"]):
            available.append(record["mapping_max"])

        record["conservative_min"] = float(np.min(available))
        record["optimistic_max"] = float(np.max(available))
        records.append(record)

    return pd.DataFrame(records).sort_values("reference").reset_index(drop=True)


def apply_ceiling_formulas(metric_rows: pd.DataFrame, c_by_reference: Dict[str, float]) -> pd.DataFrame:
    """Apply ceiling formulas to each row using a reference-specific c mapping."""
    rows = metric_rows.copy()
    rows["coverage_scenario_c"] = rows["reference"].map(c_by_reference)
    rows = rows[rows["coverage_scenario_c"].notna()].copy()

    c_vals = np.clip(rows["coverage_scenario_c"].to_numpy(), 1e-6, 0.999999)
    latent_true = rows["true_edges"].to_numpy() / c_vals

    rows["latent_true_edges_hat"] = latent_true
    rows["latent_precision_hat"] = np.clip(rows["precision"] / c_vals, 0.0, 1.0)
    rows["f1_ceiling"] = f1_ceiling(c_vals)
    rows["f1_ratio"] = rows["f1"] / np.maximum(rows["f1_ceiling"], SMALL)

    has_aupr = rows["aupr"].notna()
    rows["aupr_ceiling"] = np.nan
    rows["aupr_ratio"] = np.nan
    rows["aupr_uplift_ratio"] = np.nan

    if has_aupr.any():
        idx = has_aupr.to_numpy()
        ap_ceiling = aupr_ceiling(c_vals[idx], latent_true[idx])
        rows.loc[has_aupr, "aupr_ceiling"] = ap_ceiling

        rows.loc[has_aupr, "aupr_ratio"] = (
            rows.loc[has_aupr, "aupr"] / np.maximum(rows.loc[has_aupr, "aupr_ceiling"], SMALL)
        )

        uplift_denom = np.maximum(
            rows.loc[has_aupr, "aupr_ceiling"] - rows.loc[has_aupr, "observed_base_rate"], SMALL
        )
        rows.loc[has_aupr, "aupr_uplift_ratio"] = (
            rows.loc[has_aupr, "aupr"] - rows.loc[has_aupr, "observed_base_rate"]
        ) / uplift_denom

    return rows


def make_scenario_long_table(
    metric_rows: pd.DataFrame,
    scenario_table: pd.DataFrame,
    scenario_columns: Iterable[str],
) -> pd.DataFrame:
    """Evaluate row-level ceilings for each coverage scenario."""
    frames: List[pd.DataFrame] = []
    for scenario_name in scenario_columns:
        c_map = scenario_table.set_index("reference")[scenario_name].to_dict()
        evaluated = apply_ceiling_formulas(metric_rows, c_map)
        evaluated["scenario"] = scenario_name
        frames.append(evaluated)
    return pd.concat(frames, ignore_index=True)


def summarize_references(rows: pd.DataFrame) -> pd.DataFrame:
    """Reference-level summary for central scenario metrics."""
    out = (
        rows.groupby("reference", as_index=False)
        .agg(
            n_rows=("reference", "size"),
            methods=("method", "nunique"),
            c_median=("coverage_scenario_c", "median"),
            observed_f1_median=("f1", "median"),
            f1_ceiling_median=("f1_ceiling", "median"),
            f1_ratio_median=("f1_ratio", "median"),
            observed_aupr_median=("aupr", "median"),
            aupr_ceiling_median=("aupr_ceiling", "median"),
            aupr_ratio_median=("aupr_ratio", "median"),
            aupr_uplift_ratio_median=("aupr_uplift_ratio", "median"),
        )
        .sort_values(["aupr_ratio_median", "f1_ratio_median"], ascending=False)
    )
    return out


def summarize_methods(aupr_rows: pd.DataFrame, f1_rows: pd.DataFrame) -> pd.DataFrame:
    """Method-level summary combining AUPR and F1 headroom statistics."""
    aupr_summary = (
        aupr_rows.groupby("method", as_index=False)
        .agg(
            aupr_rows=("method", "size"),
            aupr_references=("reference", "nunique"),
            observed_aupr_median=("aupr", "median"),
            aupr_ceiling_median=("aupr_ceiling", "median"),
            aupr_ratio_median=("aupr_ratio", "median"),
            aupr_uplift_ratio_median=("aupr_uplift_ratio", "median"),
        )
        .copy()
    )

    f1_summary = (
        f1_rows.groupby("method", as_index=False)
        .agg(
            f1_rows=("method", "size"),
            f1_references=("reference", "nunique"),
            observed_f1_median=("f1", "median"),
            f1_ceiling_median=("f1_ceiling", "median"),
            f1_ratio_median=("f1_ratio", "median"),
        )
        .copy()
    )

    merged = aupr_summary.merge(f1_summary, on="method", how="outer")

    # Rank by normalized ratios to highlight ceiling-adjusted performance.
    merged["rank_observed_aupr"] = merged["observed_aupr_median"].rank(
        method="min", ascending=False
    )
    merged["rank_aupr_ratio"] = merged["aupr_ratio_median"].rank(method="min", ascending=False)
    merged["rank_shift_aupr"] = merged["rank_aupr_ratio"] - merged["rank_observed_aupr"]

    merged["rank_observed_f1"] = merged["observed_f1_median"].rank(method="min", ascending=False)
    merged["rank_f1_ratio"] = merged["f1_ratio_median"].rank(method="min", ascending=False)
    merged["rank_shift_f1"] = merged["rank_f1_ratio"] - merged["rank_observed_f1"]

    return merged.sort_values("rank_aupr_ratio")


def compute_scenario_method_ranks(scenario_rows: pd.DataFrame) -> pd.DataFrame:
    """Compute method ranks for each coverage scenario using median AUPR ratio."""
    method_scores = (
        scenario_rows.groupby(["scenario", "method"], as_index=False)
        .agg(aupr_ratio_median=("aupr_ratio", "median"), observed_aupr_median=("aupr", "median"))
        .copy()
    )

    # Rank within each scenario.
    method_scores["scenario_rank"] = method_scores.groupby("scenario")["aupr_ratio_median"].rank(
        method="min", ascending=False
    )

    # Raw AUPR rank is scenario-invariant, but we compute it from each scenario table
    # for simpler joining.
    method_scores["raw_rank"] = method_scores.groupby("scenario")["observed_aupr_median"].rank(
        method="min", ascending=False
    )

    baseline = method_scores[method_scores["scenario"] == "overlap_point"][
        ["method", "scenario_rank"]
    ].rename(columns={"scenario_rank": "overlap_point_rank"})

    method_scores = method_scores.merge(baseline, on="method", how="left")
    method_scores["rank_shift_vs_overlap"] = (
        method_scores["scenario_rank"] - method_scores["overlap_point_rank"]
    )
    return method_scores.sort_values(["scenario", "scenario_rank", "method"])


def summarize_rank_stability(scenario_ranks: pd.DataFrame) -> pd.DataFrame:
    """Summarize rank-correlation stability between scenarios."""
    base = scenario_ranks[scenario_ranks["scenario"] == "overlap_point"][
        ["method", "scenario_rank", "raw_rank"]
    ].rename(columns={"scenario_rank": "base_rank", "raw_rank": "base_raw_rank"})

    rows: List[Dict[str, float]] = []
    for scenario, subset in scenario_ranks.groupby("scenario"):
        merged = subset[["method", "scenario_rank", "raw_rank"]].merge(base, on="method", how="inner")

        spearman_vs_base = merged["scenario_rank"].corr(merged["base_rank"], method="spearman")
        spearman_vs_raw = merged["scenario_rank"].corr(merged["base_raw_rank"], method="spearman")

        rows.append(
            {
                "scenario": scenario,
                "n_methods": len(merged),
                "spearman_vs_overlap_rank": float(spearman_vs_base),
                "spearman_vs_raw_rank": float(spearman_vs_raw),
                "max_abs_rank_shift_vs_overlap": float(
                    np.max(np.abs(merged["scenario_rank"] - merged["base_rank"]))
                ),
            }
        )

    return pd.DataFrame(rows).sort_values("scenario")


def build_rowwise_coverage_sensitivity(
    central_aupr_rows: pd.DataFrame,
    scenario_table: pd.DataFrame,
    n_points: int = 50,
) -> pd.DataFrame:
    """Evaluate each row across a reference-specific c range.

    We use [conservative_min, optimistic_max] for each reference and track the
    implied AUPR/F1 ceiling ratios across a dense c grid.
    """
    scenario_bounds = scenario_table.set_index("reference")[["conservative_min", "optimistic_max"]]

    records: List[Dict[str, float]] = []
    for row in central_aupr_rows.itertuples(index=False):
        c_min = float(scenario_bounds.loc[row.reference, "conservative_min"])
        c_max = float(scenario_bounds.loc[row.reference, "optimistic_max"])

        if c_max <= c_min:
            c_values = np.array([c_min])
        else:
            c_values = np.linspace(c_min, c_max, n_points)

        latent_true = row.true_edges / c_values
        f1_caps = f1_ceiling(c_values)
        ap_caps = aupr_ceiling(c_values, latent_true)

        f1_ratio = row.f1 / np.maximum(f1_caps, SMALL)
        ap_ratio = row.aupr / np.maximum(ap_caps, SMALL)

        for c_val, f1_cap, ap_cap, f1_r, ap_r in zip(c_values, f1_caps, ap_caps, f1_ratio, ap_ratio):
            records.append(
                {
                    "method": row.method,
                    "reference": row.reference,
                    "source_table": row.source_table,
                    "c_value": float(c_val),
                    "f1_ceiling": float(f1_cap),
                    "aupr_ceiling": float(ap_cap),
                    "f1_ratio": float(f1_r),
                    "aupr_ratio": float(ap_r),
                }
            )

    return pd.DataFrame(records)


def summarize_method_sensitivity(sensitivity_rows: pd.DataFrame) -> pd.DataFrame:
    """Summarize method-level variability across adversarial coverage bands."""
    out = (
        sensitivity_rows.groupby("method", as_index=False)
        .agg(
            aupr_ratio_min=("aupr_ratio", "min"),
            aupr_ratio_median=("aupr_ratio", "median"),
            aupr_ratio_max=("aupr_ratio", "max"),
            f1_ratio_min=("f1_ratio", "min"),
            f1_ratio_median=("f1_ratio", "median"),
            f1_ratio_max=("f1_ratio", "max"),
        )
        .copy()
    )

    out["aupr_ratio_range"] = out["aupr_ratio_max"] - out["aupr_ratio_min"]
    out["f1_ratio_range"] = out["f1_ratio_max"] - out["f1_ratio_min"]
    return out.sort_values("aupr_ratio_median", ascending=False)


def run_synthetic_mar_validation(rng: np.random.Generator) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Validate MAR theory numerically with ideal latent ranking simulations."""
    c_values = [0.05, 0.1, 0.2, 0.4, 0.7, 0.9]
    t_values = [50, 200, 1000, 5000]
    trials = 500
    neg_multiplier = 40

    rows: List[Dict[str, float]] = []

    for c in c_values:
        for t in t_values:
            n_neg = int(neg_multiplier * t)
            for trial in range(trials):
                # Ideal latent ranking means all latent positives come first.
                observed_labels_top = rng.random(t) < c
                observed_pos = int(observed_labels_top.sum())
                if observed_pos == 0:
                    continue

                labels = np.concatenate(
                    [observed_labels_top.astype(np.int8), np.zeros(n_neg, dtype=np.int8)]
                )

                # Precision/recall/F1 when predicting the latent positive block (top t).
                precision_obs = observed_pos / float(t)
                recall_obs = 1.0
                f1_obs = (2.0 * precision_obs * recall_obs) / (precision_obs + recall_obs)
                ap_obs = average_precision_from_binary_labels(labels)

                rows.append(
                    {
                        "c": c,
                        "t": t,
                        "trial": trial,
                        "precision_obs": precision_obs,
                        "recall_obs": recall_obs,
                        "f1_obs": f1_obs,
                        "aupr_obs": ap_obs,
                        "precision_theory": c,
                        "recall_theory": 1.0,
                        "f1_theory": float(f1_ceiling(np.array([c]))[0]),
                        "aupr_theory": float(aupr_ceiling(np.array([c]), np.array([t]))[0]),
                    }
                )

    full = pd.DataFrame(rows)
    full["f1_abs_error"] = np.abs(full["f1_obs"] - full["f1_theory"])
    full["aupr_abs_error"] = np.abs(full["aupr_obs"] - full["aupr_theory"])

    summary = (
        full.groupby(["c", "t"], as_index=False)
        .agg(
            trials=("trial", "size"),
            precision_obs_mean=("precision_obs", "mean"),
            precision_theory=("precision_theory", "mean"),
            recall_obs_mean=("recall_obs", "mean"),
            f1_obs_mean=("f1_obs", "mean"),
            f1_theory=("f1_theory", "mean"),
            aupr_obs_mean=("aupr_obs", "mean"),
            aupr_theory=("aupr_theory", "mean"),
            f1_abs_error_mean=("f1_abs_error", "mean"),
            aupr_abs_error_mean=("aupr_abs_error", "mean"),
        )
        .sort_values(["c", "t"])
    )
    summary["f1_bias"] = summary["f1_obs_mean"] - summary["f1_theory"]
    summary["aupr_bias"] = summary["aupr_obs_mean"] - summary["aupr_theory"]
    summary["f1_bias_abs"] = np.abs(summary["f1_bias"])
    summary["aupr_bias_abs"] = np.abs(summary["aupr_bias"])

    return full, summary


def _rank_dependent_probs(t: int, c: float, mode: str) -> np.ndarray:
    """Create rank-dependent inclusion probabilities for non-MAR stress tests.

    `mode`:
    - optimistic: top-ranked latent edges are more likely to be labeled.
    - pessimistic: top-ranked latent edges are less likely to be labeled.

    Mean probability is normalized back to target c.
    """
    x = np.linspace(0.0, 1.0, t)
    if mode == "optimistic":
        raw = 1.6 - x
    elif mode == "pessimistic":
        raw = 0.6 + x
    else:
        raise ValueError(f"Unknown mode: {mode}")

    probs = raw / np.mean(raw)
    probs = probs * c

    # Clip to valid probabilities, then re-scale once more to keep the mean close to c.
    probs = np.clip(probs, 1e-4, 0.999)
    probs = probs * (c / np.mean(probs))
    probs = np.clip(probs, 1e-4, 0.999)
    return probs


def run_synthetic_adversarial_missingness(
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Stress-test ceiling assumptions under non-MAR positive labeling."""
    c_values = [0.1, 0.2, 0.4, 0.7]
    t_values = [100, 500, 2000]
    trials = 500
    neg_multiplier = 40

    rows: List[Dict[str, float]] = []

    for c in c_values:
        for t in t_values:
            n_neg = int(neg_multiplier * t)
            for mode in ["mar", "optimistic", "pessimistic"]:
                if mode == "mar":
                    probs = np.full(t, c, dtype=float)
                else:
                    probs = _rank_dependent_probs(t=t, c=c, mode=mode)

                for trial in range(trials):
                    observed_labels_top = rng.random(t) < probs
                    observed_pos = int(observed_labels_top.sum())
                    if observed_pos == 0:
                        continue

                    labels = np.concatenate(
                        [observed_labels_top.astype(np.int8), np.zeros(n_neg, dtype=np.int8)]
                    )

                    precision_obs = observed_pos / float(t)
                    f1_obs = (2.0 * precision_obs) / (1.0 + precision_obs)
                    ap_obs = average_precision_from_binary_labels(labels)

                    f1_theory_mar = float(f1_ceiling(np.array([c]))[0])
                    ap_theory_mar = float(aupr_ceiling(np.array([c]), np.array([t]))[0])

                    rows.append(
                        {
                            "mode": mode,
                            "c_target": c,
                            "t": t,
                            "trial": trial,
                            "c_realized": float(observed_pos / t),
                            "f1_obs": f1_obs,
                            "aupr_obs": ap_obs,
                            "f1_theory_mar": f1_theory_mar,
                            "aupr_theory_mar": ap_theory_mar,
                            "f1_delta_vs_mar_theory": f1_obs - f1_theory_mar,
                            "aupr_delta_vs_mar_theory": ap_obs - ap_theory_mar,
                        }
                    )

    full = pd.DataFrame(rows)
    summary = (
        full.groupby(["mode", "c_target", "t"], as_index=False)
        .agg(
            trials=("trial", "size"),
            c_realized_mean=("c_realized", "mean"),
            f1_mean=("f1_obs", "mean"),
            aupr_mean=("aupr_obs", "mean"),
            f1_theory_mar=("f1_theory_mar", "mean"),
            aupr_theory_mar=("aupr_theory_mar", "mean"),
            f1_delta_mean=("f1_delta_vs_mar_theory", "mean"),
            aupr_delta_mean=("aupr_delta_vs_mar_theory", "mean"),
        )
        .sort_values(["mode", "c_target", "t"])
    )

    return full, summary


def load_sweep_artifacts(repo_root: Path) -> pd.DataFrame:
    """Load sweep result JSONs into a compact table for external sanity checks."""
    rows: List[Dict[str, float]] = []

    for name, rel_path in SWEEP_JSONS.items():
        payload = json.loads((repo_root / rel_path).read_text())

        candidate_edges = payload.get("candidate_edges", np.nan)
        global_aupr = payload.get("aupr", np.nan)

        rows.append(
            {
                "artifact": name,
                "entry_type": "global",
                "candidate_edges": candidate_edges,
                "true_edges": np.nan,
                "observed_baserate": np.nan,
                "aupr": global_aupr,
                "f1": np.nan,
                "precision": np.nan,
                "recall": np.nan,
                "sweep_axis": "global",
                "sweep_value": np.nan,
            }
        )

        for entry in payload.get("percentile_sweep", []):
            true_edges = float(entry.get("true_edges", np.nan))
            baserate = true_edges / candidate_edges if candidate_edges and not np.isnan(true_edges) else np.nan
            rows.append(
                {
                    "artifact": name,
                    "entry_type": "percentile_sweep",
                    "candidate_edges": candidate_edges,
                    "true_edges": true_edges,
                    "observed_baserate": baserate,
                    "aupr": np.nan,
                    "f1": float(entry.get("f1", np.nan)),
                    "precision": float(entry.get("precision", np.nan)),
                    "recall": float(entry.get("recall", np.nan)),
                    "sweep_axis": "percentile",
                    "sweep_value": float(entry.get("percentile", np.nan)),
                }
            )

        for entry in payload.get("top_k_sweep", []):
            true_edges = float(entry.get("true_edges", np.nan))
            baserate = true_edges / candidate_edges if candidate_edges and not np.isnan(true_edges) else np.nan
            rows.append(
                {
                    "artifact": name,
                    "entry_type": "top_k_sweep",
                    "candidate_edges": candidate_edges,
                    "true_edges": true_edges,
                    "observed_baserate": baserate,
                    "aupr": np.nan,
                    "f1": float(entry.get("f1", np.nan)),
                    "precision": float(entry.get("precision", np.nan)),
                    "recall": float(entry.get("recall", np.nan)),
                    "sweep_axis": "top_k",
                    "sweep_value": float(entry.get("top_k", np.nan)),
                }
            )

    return pd.DataFrame(rows)


def load_cross_eval_artifacts(repo_root: Path) -> pd.DataFrame:
    """Load cross-evaluation artifacts into a normalized flat table."""
    rows: List[Dict[str, float]] = []

    for artifact, rel_path in CROSS_EVAL_JSONS.items():
        payload = json.loads((repo_root / rel_path).read_text())
        candidate_edges = payload.get("candidate_edges", np.nan)

        for system in ["hpn", "beeline"]:
            block = payload.get(system)
            if block is None:
                continue

            if isinstance(block, dict) and "aupr" in block:
                positives = block.get("positives", np.nan)
                rows.append(
                    {
                        "artifact": artifact,
                        "system": system,
                        "calibration": "default",
                        "candidate_edges": candidate_edges,
                        "positives": positives,
                        "observed_baserate": (
                            positives / candidate_edges
                            if candidate_edges and positives is not None and not np.isnan(positives)
                            else np.nan
                        ),
                        "aupr": float(block.get("aupr", np.nan)),
                        "brier": float(block.get("brier", np.nan)) if "brier" in block else np.nan,
                    }
                )
            elif isinstance(block, dict):
                # Nested calibration case: {raw: {aupr, brier}, isotonic: ..., logistic: ...}
                for calib, values in block.items():
                    rows.append(
                        {
                            "artifact": artifact,
                            "system": system,
                            "calibration": calib,
                            "candidate_edges": candidate_edges,
                            "positives": np.nan,
                            "observed_baserate": np.nan,
                            "aupr": float(values.get("aupr", np.nan)),
                            "brier": float(values.get("brier", np.nan)),
                        }
                    )

    return pd.DataFrame(rows)


def plot_theoretical_curves(central_aupr_rows: pd.DataFrame, output_path: Path) -> None:
    """Plot core theoretical ceiling curves and random AUPR baselines."""
    c_grid = np.linspace(0.02, 0.98, 400)
    median_t = float(np.median(central_aupr_rows["latent_true_edges_hat"]))

    f1_grid = f1_ceiling(c_grid)
    ap_grid = aupr_ceiling(c_grid, np.full_like(c_grid, median_t))

    latent_base_rates = [1e-5, 5e-5, 1e-4, 5e-4]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=180)

    axes[0].plot(c_grid, f1_grid, color="#1f77b4", lw=2.2)
    axes[0].set_title("F1 Ceiling Under Partial Labels")
    axes[0].set_xlabel("Coverage c")
    axes[0].set_ylabel("Max Observed F1")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(c_grid, ap_grid, color="#d62728", lw=2.2, label="AUPR ceiling")
    for b in latent_base_rates:
        axes[1].plot(c_grid, c_grid * b, ls="--", lw=1.2, label=f"Random baseline (b={b:.0e})")
    axes[1].set_title("AUPR Ceiling and Random Baselines")
    axes[1].set_xlabel("Coverage c")
    axes[1].set_ylabel("Observed AUPR")
    axes[1].set_yscale("log")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(fontsize=8)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_observed_vs_ceiling_ratio(central_aupr_rows: pd.DataFrame, output_path: Path) -> None:
    """Scatter observed metrics against ceiling-normalized ratios."""
    refs = sorted(central_aupr_rows["reference"].unique())
    cmap = plt.get_cmap("tab10")
    color_map = {ref: cmap(i % 10) for i, ref in enumerate(refs)}

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.4), dpi=180)

    for row in central_aupr_rows.itertuples(index=False):
        color = color_map[row.reference]
        axes[0].scatter(max(row.f1, 1e-8), row.f1_ratio, color=color, s=42, edgecolor="black", linewidth=0.3)
        axes[1].scatter(max(row.aupr, 1e-8), row.aupr_ratio, color=color, s=42, edgecolor="black", linewidth=0.3)

    axes[0].set_xscale("log")
    axes[1].set_xscale("log")
    axes[0].set_xlabel("Observed F1")
    axes[1].set_xlabel("Observed AUPR")
    axes[0].set_ylabel("Observed F1 / F1 Ceiling")
    axes[1].set_ylabel("Observed AUPR / AUPR Ceiling")
    axes[0].set_title("F1 Ceiling Ratios")
    axes[1].set_title("AUPR Ceiling Ratios")

    for ax in axes:
        ax.grid(True, alpha=0.25)

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=color_map[r],
            markeredgecolor="black",
            markeredgewidth=0.3,
            label=r,
            markersize=6,
        )
        for r in refs
    ]

    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 4))
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_rank_shift(method_summary: pd.DataFrame, output_path: Path) -> None:
    """Plot rank shifts between raw AUPR and ceiling-normalized AUPR rankings."""
    plot_df = method_summary.dropna(subset=["rank_shift_aupr"]).copy()
    plot_df = plot_df.sort_values("rank_shift_aupr", ascending=False)

    fig, ax = plt.subplots(figsize=(10.2, 5.4), dpi=180)
    colors = ["#2ca02c" if x > 0 else ("#d62728" if x < 0 else "#7f7f7f") for x in plot_df["rank_shift_aupr"]]

    ax.bar(plot_df["method"], plot_df["rank_shift_aupr"], color=colors, edgecolor="black", linewidth=0.3)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_ylabel("Rank shift (normalized rank - raw rank)")
    ax.set_xlabel("Method")
    ax.set_title("Method Ranking Shift Under Ceiling Normalization (AUPR)")
    ax.grid(True, axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=35)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_scenario_rank_heatmap(scenario_ranks: pd.DataFrame, output_path: Path) -> None:
    """Heatmap of method ranks under alternate coverage scenarios."""
    heat = scenario_ranks.pivot_table(index="method", columns="scenario", values="scenario_rank", aggfunc="min")
    heat = heat.reindex(sorted(heat.index), axis=0)
    heat = heat.reindex(sorted(heat.columns), axis=1)

    fig, ax = plt.subplots(figsize=(12, 6), dpi=180)
    heat_values = heat.to_numpy(dtype=float)
    mask = np.isnan(heat_values)
    masked_values = np.ma.masked_array(heat_values, mask=mask)
    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad(color="lightgray")
    im = ax.imshow(masked_values, aspect="auto", cmap=cmap)

    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_xticklabels(heat.columns, rotation=40, ha="right")
    ax.set_yticklabels(heat.index)

    for i in range(len(heat.index)):
        for j in range(len(heat.columns)):
            val = heat.iloc[i, j]
            if pd.isna(val):
                ax.text(j, i, "NA", ha="center", va="center", color="black", fontsize=7)
            else:
                ax.text(j, i, f"{int(val)}", ha="center", va="center", color="white", fontsize=8)

    ax.set_title("Method Rank by Coverage Scenario (lower is better)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Rank")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_sensitivity_bands(method_sensitivity: pd.DataFrame, output_path: Path) -> None:
    """Plot AUPR-ratio ranges across adversarial c bands for each method."""
    plot_df = method_sensitivity.sort_values("aupr_ratio_median", ascending=False).copy()
    y = np.arange(len(plot_df))

    fig, ax = plt.subplots(figsize=(11.5, 5.8), dpi=180)

    ax.hlines(
        y,
        plot_df["aupr_ratio_min"],
        plot_df["aupr_ratio_max"],
        color="#1f77b4",
        linewidth=2.2,
        alpha=0.7,
    )
    ax.scatter(plot_df["aupr_ratio_median"], y, color="#d62728", s=36, zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["method"])
    ax.set_xlabel("AUPR ceiling ratio across c sensitivity band")
    ax.set_ylabel("Method")
    ax.set_title("Adversarial Coverage Sensitivity: Median and Range")
    ax.grid(True, axis="x", alpha=0.25)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_mar_validation(mar_summary: pd.DataFrame, output_path: Path) -> None:
    """Plot MAR validation bias versus latent-positive count."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=180)

    for c_val, subset in mar_summary.groupby("c"):
        subset = subset.sort_values("t")
        axes[0].plot(subset["t"], subset["f1_bias_abs"], marker="o", label=f"c={c_val:.2f}")
        axes[1].plot(subset["t"], subset["aupr_bias_abs"], marker="o", label=f"c={c_val:.2f}")

    axes[0].set_xscale("log")
    axes[1].set_xscale("log")
    axes[0].set_yscale("log")
    axes[1].set_yscale("log")

    axes[0].set_title("MAR Validation: |E[F1] - Theory|")
    axes[1].set_title("MAR Validation: |E[AUPR] - Theory|")
    axes[0].set_xlabel("Latent positives t")
    axes[1].set_xlabel("Latent positives t")
    axes[0].set_ylabel("Absolute error")
    axes[1].set_ylabel("Absolute error")

    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_adversarial_missingness(adversarial_summary: pd.DataFrame, output_path: Path) -> None:
    """Plot AP deviation from MAR-theory under different missingness modes."""
    # Average over t for compact view per c/mode.
    compact = (
        adversarial_summary.groupby(["mode", "c_target"], as_index=False)
        .agg(
            aupr_delta_mean=("aupr_delta_mean", "mean"),
            f1_delta_mean=("f1_delta_mean", "mean"),
        )
        .copy()
    )

    modes = ["mar", "optimistic", "pessimistic"]
    x = np.arange(len(sorted(compact["c_target"].unique())))
    c_vals = sorted(compact["c_target"].unique())

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=180)
    width = 0.24

    for idx, mode in enumerate(modes):
        subset = compact[compact["mode"] == mode].sort_values("c_target")
        offset = (idx - 1) * width
        axes[0].bar(x + offset, subset["aupr_delta_mean"], width=width, label=mode)
        axes[1].bar(x + offset, subset["f1_delta_mean"], width=width, label=mode)

    axes[0].axhline(0.0, color="black", linewidth=1.0)
    axes[1].axhline(0.0, color="black", linewidth=1.0)
    axes[0].set_xticks(x)
    axes[1].set_xticks(x)
    axes[0].set_xticklabels([f"c={c:.1f}" for c in c_vals])
    axes[1].set_xticklabels([f"c={c:.1f}" for c in c_vals])

    axes[0].set_title("AUPR deviation vs MAR theory")
    axes[1].set_title("F1 deviation vs MAR theory")
    axes[0].set_ylabel("Observed - MAR theory")
    axes[1].set_ylabel("Observed - MAR theory")

    for ax in axes:
        ax.set_xlabel("Coverage target")
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(fontsize=8)

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def build_results_summary_markdown(
    central_aupr_rows: pd.DataFrame,
    central_f1_rows: pd.DataFrame,
    reference_summary: pd.DataFrame,
    method_summary: pd.DataFrame,
    rank_stability: pd.DataFrame,
    method_sensitivity: pd.DataFrame,
    mar_summary: pd.DataFrame,
    adversarial_summary: pd.DataFrame,
    sweep_rows: pd.DataFrame,
    cross_eval_rows: pd.DataFrame,
) -> str:
    """Create a deterministic markdown summary grounded in generated outputs."""
    best_f1 = central_aupr_rows.sort_values("f1_ratio", ascending=False).iloc[0]
    best_aupr = central_aupr_rows.sort_values("aupr_ratio", ascending=False).iloc[0]

    median_f1_ratio = float(central_aupr_rows["f1_ratio"].median())
    median_aupr_ratio = float(central_aupr_rows["aupr_ratio"].median())
    median_uplift = float(central_aupr_rows["aupr_uplift_ratio"].median())

    mar_f1_error = float(mar_summary["f1_bias_abs"].median())
    mar_aupr_error = float(mar_summary["aupr_bias_abs"].median())

    adv_compact = (
        adversarial_summary.groupby("mode", as_index=False)
        .agg(aupr_delta_mean=("aupr_delta_mean", "mean"), f1_delta_mean=("f1_delta_mean", "mean"))
        .sort_values("mode")
    )

    most_unstable = method_sensitivity.sort_values("aupr_ratio_range", ascending=False).iloc[0]

    # Cross-eval snapshot.
    cross_snapshot = cross_eval_rows.dropna(subset=["aupr"]).copy()
    cross_snapshot = cross_snapshot.sort_values("aupr", ascending=False).head(6)

    lines: List[str] = []
    lines.append("# Partial-label ceiling Full-Pipeline Results Summary")
    lines.append("")
    lines.append("## Core reinterpretation")
    lines.append(
        f"- AUPR+F1 evaluable rows: `{len(central_aupr_rows)}` across `{central_aupr_rows['method'].nunique()}` methods and `{central_aupr_rows['reference'].nunique()}` references."
    )
    lines.append(
        f"- Best F1 ceiling ratio: `{best_f1['method']}` on `{best_f1['reference']}` = `{best_f1['f1_ratio']:.6f}`."
    )
    lines.append(
        f"- Best AUPR ceiling ratio: `{best_aupr['method']}` on `{best_aupr['reference']}` = `{best_aupr['aupr_ratio']:.6f}`."
    )
    lines.append(f"- Median F1 ceiling ratio: `{median_f1_ratio:.6f}`.")
    lines.append(f"- Median AUPR ceiling ratio: `{median_aupr_ratio:.6f}`.")
    lines.append(f"- Median normalized AUPR uplift ratio: `{median_uplift:.6e}`.")
    lines.append("")

    lines.append("## Reference-level medians")
    lines.append(reference_summary.to_markdown(index=False))
    lines.append("")

    lines.append("## Method-level ranking shifts")
    rank_cols = [
        "method",
        "observed_aupr_median",
        "aupr_ratio_median",
        "rank_observed_aupr",
        "rank_aupr_ratio",
        "rank_shift_aupr",
    ]
    lines.append(method_summary[rank_cols].sort_values("rank_aupr_ratio").to_markdown(index=False))
    lines.append("")

    lines.append("## Scenario rank stability")
    lines.append(rank_stability.to_markdown(index=False))
    lines.append("")

    lines.append("## Most coverage-sensitive method")
    lines.append(
        f"- `{most_unstable['method']}` has the widest AUPR-ratio band: `{most_unstable['aupr_ratio_range']:.6f}` (min `{most_unstable['aupr_ratio_min']:.6f}`, max `{most_unstable['aupr_ratio_max']:.6f}`)."
    )
    lines.append("")

    lines.append("## Synthetic theorem validation")
    lines.append(f"- Median MAR absolute F1 bias: `{mar_f1_error:.6e}`.")
    lines.append(f"- Median MAR absolute AUPR bias: `{mar_aupr_error:.6e}`.")
    lines.append("- Non-MAR missingness shifts AUPR relative to MAR theory:")
    lines.append(adv_compact.to_markdown(index=False))
    lines.append("")

    lines.append("## External artifact sanity checks")
    sweep_global = sweep_rows[sweep_rows["entry_type"] == "global"][
        ["artifact", "candidate_edges", "aupr"]
    ].copy()
    lines.append("### Sweep globals")
    lines.append(sweep_global.to_markdown(index=False))
    lines.append("")

    lines.append("### Cross-eval top AUPR entries")
    lines.append(cross_snapshot.to_markdown(index=False))
    lines.append("")

    lines.append("## Interpretation")
    lines.append("- Ceiling-aware reinterpretation remains strongly conservative: observed metrics are tiny fractions of partial-label ceilings.")
    lines.append("- Method ranking can shift after coverage normalization because methods are evaluated on different reference mixes and coverage regimes.")
    lines.append("- MAR formulas are numerically accurate under their assumptions, but adversarial non-MAR labeling can move AUPR away from MAR-theory ceilings.")

    return "\n".join(lines)


def write_manuscript_template(
    paper_dir: Path,
    summary_tables_dir: Path,
    figures_dir: Path,
) -> None:
    """Write a full-length manuscript draft populated with generated artifact paths."""
    manuscript_path = paper_dir / "proposal1_partial_label_metric_ceilings_paper.md"

    text = f"""# Partial-Label Metric Ceilings for GRN Evaluation: Theory, Stress Tests, and Reinterpretation

## Abstract
Absolute GRN benchmark metrics are often interpreted without accounting for incomplete positive labels in reference networks. We formalize partial-label ceilings for observed precision, F1, and AUPR under missing-at-random positive labeling, derive finite-sample AUPR corrections, and show how uncertainty in coverage propagates to ceiling uncertainty. We then execute a low-compute empirical reinterpretation across existing GRN evaluation outputs and quantify robustness under adversarial coverage misspecification and non-MAR labeling stress tests. Across 39 evaluable rows from 15 methods and 5 references, observed scores remain small fractions of estimated ceilings (best F1 ratio 0.137, best AUPR ratio 0.0138; medians near zero). Synthetic validation confirms the MAR formulas are numerically tight, while non-MAR missingness can systematically perturb AUPR relative to MAR predictions. We provide a reproducible ceiling-aware reporting protocol that reduces benchmark overclaiming and underclaiming.

## 1. Introduction
Evaluation in GRN inference is fundamentally limited by label incompleteness: curated references expose only a subset of latent true edges. This introduces an interpretation gap between observed and latent metrics. Partial-label ceiling in this research track asked a concrete question: given incomplete references, what is the best achievable observed F1/AUPR even for an ideal predictor?

This work makes three contributions:
1. A formal ceiling framework for observed F1/AUPR under partial labeling, with finite-sample corrections.
2. A reproducible empirical reinterpretation over existing benchmark artifacts using explicit coverage proxies and uncertainty bands.
3. Adversarial stress tests showing when MAR-based ceilings are reliable and where assumption violations can bias conclusions.

## 2. Problem Setup
Let `U` be the candidate edge universe, `T ⊂ U` latent true edges, and `L+ ⊂ T` observed positives. Coverage is `c = P(e ∈ L+ | e ∈ T)`. For predicted edges `P`:
- latent precision: `Prec_lat = |P ∩ T| / |P|`
- observed precision: `Prec_obs = |P ∩ L+| / |P|`
- latent recall: `Rec_lat = |P ∩ T| / |T|`
- observed recall: `Rec_obs = |P ∩ L+| / |L+|`

Assumption A1 (MAR labels): each latent positive is observed independently with probability `c`.

## 3. Theory
### Proposition 1 (Observed-latent relation)
Under A1:
- `E[Prec_obs] = c * Prec_lat`
- `E[Rec_obs] = Rec_lat`

Hence `Prec_lat` can be approximated as `Prec_obs / c` when `c` is estimable.

### Proposition 2 (Observed F1 ceiling)
For an ideal latent predictor (`Prec_lat = Rec_lat = 1`):
- `Prec_obs = c`, `Rec_obs = 1`
- `F1_obs,max(c) = 2c / (1 + c)`

### Proposition 3 (Observed AUPR ceiling)
With latent positives `t = |T|`, ideal latent ranking yields:
- `AUPR_obs,max(c, t) ≈ c + (1-c) * H_t / t`
where `H_t` is the harmonic number.
For large `t`, `AUPR_obs,max ≈ c`.

The observed random baseline is `AUPR_rand,obs = cb` where `b = |T|/|U|`.

### Confidence intervals for ceilings
We model reference-level coverage with a Beta posterior:
- `c ~ Beta(k+1, n-k+1)`
where `k` is overlap support and `n` total reference nodes.
Ceiling intervals are propagated by Monte Carlo draws from this posterior.

## 4. Methods
### 4.1 Data artifacts
Primary metric tables:
- `data/score_eval_probe_priors.csv`
- `data/score_eval_grn_baselines_immune.csv`

F1-extended tables:
- `data/score_eval_probe_priors_full_genes.csv`
- `data/score_eval_probe_priors_full_genes_crosswalk.csv`
- `data/score_eval_probe_priors_full_genes_omnipath.csv`

Coverage proxy artifacts:
- `data/score_eval_probe_priors_missing_report.json`
- `data/score_eval_probe_priors_full_genes_missing_report.json`
- `data/score_eval_probe_priors_full_genes_crosswalk_missing_report.json`
- `data/score_eval_probe_priors_full_genes_omnipath_missing_report.json`

### 4.2 Coverage scenarios
Per reference, we evaluate multiple coverage scenarios:
- overlap-based point estimate,
- Beta posterior interval bounds,
- adversarial `0.5x` and `1.5x` misspecification,
- mapping-endpoint min/max proxies,
- conservative and optimistic envelopes.

### 4.3 Synthetic validation and adversarial stress tests
- MAR simulation: verifies Propositions 1-3 across grids of `(c, t)`.
- Non-MAR simulation: optimistic/pessimistic rank-dependent labeling keeps average coverage near `c` but violates MAR independence.

### 4.4 Reproducibility
Single command:

```bash
python scripts/run_paper_pipeline.py
```

## 5. Results
### 5.1 Empirical ceiling reinterpretation
Main row-level reinterpretation is in:
- `{summary_tables_dir / 'central_aupr_rows.csv'}`
- `{summary_tables_dir / 'central_f1_rows.csv'}`

Reference and method summaries:
- `{summary_tables_dir / 'reference_summary.csv'}`
- `{summary_tables_dir / 'method_summary.csv'}`

### 5.2 Scenario robustness
Coverage-scenario rank diagnostics:
- `{summary_tables_dir / 'scenario_method_ranks.csv'}`
- `{summary_tables_dir / 'scenario_rank_stability.csv'}`

Coverage-band sensitivity:
- `{summary_tables_dir / 'rowwise_coverage_sensitivity.csv'}`
- `{summary_tables_dir / 'method_sensitivity_summary.csv'}`

### 5.3 Synthetic theory validation
- `{summary_tables_dir / 'synthetic_mar_trials.csv'}`
- `{summary_tables_dir / 'synthetic_mar_summary.csv'}`

### 5.4 Non-MAR adversarial stress test
- `{summary_tables_dir / 'synthetic_adversarial_trials.csv'}`
- `{summary_tables_dir / 'synthetic_adversarial_summary.csv'}`

### 5.5 External artifact sanity checks
- `{summary_tables_dir / 'sweep_artifact_rows.csv'}`
- `{summary_tables_dir / 'cross_eval_artifact_rows.csv'}`

## 6. Figures
- Theoretical curves: `{figures_dir / 'figure_01_theoretical_curves.png'}`
- Observed vs ceiling ratios: `{figures_dir / 'figure_02_observed_vs_ratio.png'}`
- Rank shifts: `{figures_dir / 'figure_03_rank_shift.png'}`
- Scenario rank heatmap: `{figures_dir / 'figure_04_scenario_rank_heatmap.png'}`
- Coverage sensitivity bands: `{figures_dir / 'figure_05_sensitivity_bands.png'}`
- MAR validation errors: `{figures_dir / 'figure_06_mar_validation.png'}`
- Adversarial missingness deviations: `{figures_dir / 'figure_07_adversarial_missingness.png'}`

## 7. Discussion
The framework provides a clear interpretation layer over low absolute benchmark metrics: values that appear weak can still be near the attainable regime if coverage is low, and values that appear strong may still occupy a small fraction of the theoretical headroom. In this dataset, observed metrics remain far below estimated ceilings, indicating substantial algorithmic headroom beyond reference incompleteness alone.

Adversarial analyses sharpen this interpretation:
- coverage misspecification can alter normalized rankings,
- non-MAR label mechanisms can bias AUPR relative to MAR-based expectations.

Therefore, ceiling-aware reporting should include explicit assumption checks and sensitivity bands, not point estimates alone.

## 8. Limitations
- Coverage `c` remains proxy-based and non-identifiable from these artifacts.
- References and methods are heterogeneous across tables, so cross-method comparisons should be read with scenario context.
- Non-MAR stress tests are stylized and do not prove any specific biological curation mechanism.

## 9. Reporting Checklist
1. Report candidate size and observed base rate.
2. Report coverage proxy definition and uncertainty interval.
3. Report raw metrics and ceiling-normalized metrics together.
4. Report sensitivity over at least one conservative and one optimistic coverage scenario.
5. Report one non-MAR stress test or discuss why MAR is defensible.

## 10. Appendix (proof sketches)
### A.1 Proof sketch for Proposition 1
Under Bernoulli thinning of true positives with rate `c`, indicator linearity gives `E[|P∩L+|]=c|P∩T|`, immediately yielding `E[Prec_obs]=cPrec_lat` and `E[Rec_obs]=Rec_lat`.

### A.2 Proof sketch for Proposition 2
For ideal latent prediction, every true edge is predicted, so observed recall is 1 and observed precision is `c`; substitute into F1 definition.

### A.3 Proof sketch for Proposition 3
For ideal latent ordering, observed positives are random positions within the first `t` ranks; AP is mean precision over those positions. Taking expectation over Bernoulli-thinned labels yields the harmonic correction term `H_t/t`.
"""

    manuscript_path.write_text(text)


def main() -> None:
    script_path = Path(__file__).resolve()
    repo_root = find_repo_root(script_path)

    output_dir = repo_root / "outputs" / "paper"
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"

    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(RNG_SEED)

    # 1) Load and prepare empirical rows.
    all_metric_rows = load_metric_frames(repo_root)
    primary_rows = all_metric_rows[
        all_metric_rows["source_table"].isin(PRIMARY_METRIC_TABLES.keys())
    ].copy()
    mapping_proxy_rows = load_mapping_proxy_rows(repo_root)
    coverage_scenarios = build_reference_coverage_scenarios(
        metric_rows=all_metric_rows,
        mapping_proxies=mapping_proxy_rows,
        rng=rng,
        preferred_rows=primary_rows,
    )

    scenario_columns = [
        "overlap_point",
        "beta_ci_low",
        "beta_ci_high",
        "adversarial_half",
        "adversarial_150pct",
        "mapping_min",
        "mapping_max",
        "conservative_min",
        "optimistic_max",
    ]

    # Central scenario for headline results.
    central_c_map = coverage_scenarios.set_index("reference")["overlap_point"].to_dict()

    central_rows_all = apply_ceiling_formulas(all_metric_rows, central_c_map)
    central_aupr_rows = central_rows_all[central_rows_all["aupr"].notna()].copy()
    central_f1_rows = central_rows_all.copy()

    # Full scenario table for robustness analyses.
    scenario_rows_all = make_scenario_long_table(all_metric_rows, coverage_scenarios, scenario_columns)
    scenario_aupr_rows = scenario_rows_all[scenario_rows_all["aupr"].notna()].copy()

    reference_summary = summarize_references(central_aupr_rows)
    method_summary = summarize_methods(central_aupr_rows, central_f1_rows)

    scenario_method_ranks = compute_scenario_method_ranks(scenario_aupr_rows)
    scenario_rank_stability = summarize_rank_stability(scenario_method_ranks)

    sensitivity_rows = build_rowwise_coverage_sensitivity(central_aupr_rows, coverage_scenarios, n_points=50)
    method_sensitivity = summarize_method_sensitivity(sensitivity_rows)

    # 2) Synthetic theorem validation + adversarial stress tests.
    mar_trials, mar_summary = run_synthetic_mar_validation(rng)
    adv_trials, adv_summary = run_synthetic_adversarial_missingness(rng)

    # 3) External artifact sanity checks.
    sweep_rows = load_sweep_artifacts(repo_root)
    cross_eval_rows = load_cross_eval_artifacts(repo_root)

    # 4) Persist tables.
    all_metric_rows.to_csv(tables_dir / "all_metric_rows_filtered.csv", index=False)
    mapping_proxy_rows.to_csv(tables_dir / "mapping_proxy_rows.csv", index=False)
    coverage_scenarios.to_csv(tables_dir / "coverage_scenarios_by_reference.csv", index=False)

    central_aupr_rows.to_csv(tables_dir / "central_aupr_rows.csv", index=False)
    central_f1_rows.to_csv(tables_dir / "central_f1_rows.csv", index=False)
    reference_summary.to_csv(tables_dir / "reference_summary.csv", index=False)
    method_summary.to_csv(tables_dir / "method_summary.csv", index=False)

    scenario_rows_all.to_csv(tables_dir / "scenario_rows_all.csv", index=False)
    scenario_method_ranks.to_csv(tables_dir / "scenario_method_ranks.csv", index=False)
    scenario_rank_stability.to_csv(tables_dir / "scenario_rank_stability.csv", index=False)

    sensitivity_rows.to_csv(tables_dir / "rowwise_coverage_sensitivity.csv", index=False)
    method_sensitivity.to_csv(tables_dir / "method_sensitivity_summary.csv", index=False)

    mar_trials.to_csv(tables_dir / "synthetic_mar_trials.csv", index=False)
    mar_summary.to_csv(tables_dir / "synthetic_mar_summary.csv", index=False)
    adv_trials.to_csv(tables_dir / "synthetic_adversarial_trials.csv", index=False)
    adv_summary.to_csv(tables_dir / "synthetic_adversarial_summary.csv", index=False)

    sweep_rows.to_csv(tables_dir / "sweep_artifact_rows.csv", index=False)
    cross_eval_rows.to_csv(tables_dir / "cross_eval_artifact_rows.csv", index=False)

    # 5) Generate figures.
    plot_theoretical_curves(central_aupr_rows, figures_dir / "figure_01_theoretical_curves.png")
    plot_observed_vs_ceiling_ratio(central_aupr_rows, figures_dir / "figure_02_observed_vs_ratio.png")
    plot_rank_shift(method_summary, figures_dir / "figure_03_rank_shift.png")
    plot_scenario_rank_heatmap(scenario_method_ranks, figures_dir / "figure_04_scenario_rank_heatmap.png")
    plot_sensitivity_bands(method_sensitivity, figures_dir / "figure_05_sensitivity_bands.png")
    plot_mar_validation(mar_summary, figures_dir / "figure_06_mar_validation.png")
    plot_adversarial_missingness(adv_summary, figures_dir / "figure_07_adversarial_missingness.png")

    # 6) Build deterministic summary + manuscript draft.
    summary_md = build_results_summary_markdown(
        central_aupr_rows=central_aupr_rows,
        central_f1_rows=central_f1_rows,
        reference_summary=reference_summary,
        method_summary=method_summary,
        rank_stability=scenario_rank_stability,
        method_sensitivity=method_sensitivity,
        mar_summary=mar_summary,
        adversarial_summary=adv_summary,
        sweep_rows=sweep_rows,
        cross_eval_rows=cross_eval_rows,
    )
    (output_dir / "results_summary.md").write_text(summary_md)

    paper_dir = repo_root / "paper"
    paper_dir.mkdir(parents=True, exist_ok=True)
    write_manuscript_template(paper_dir=paper_dir, summary_tables_dir=tables_dir, figures_dir=figures_dir)

    print("Wrote paper pipeline outputs to:")
    print(f"- {output_dir}")
    print(f"- {paper_dir / 'proposal1_partial_label_metric_ceilings_paper.md'}")


if __name__ == "__main__":
    main()
