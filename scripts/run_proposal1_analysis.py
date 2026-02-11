#!/usr/bin/env python3
"""Proposal 1: Partial-label metric ceiling analysis.

This script executes a low-compute analysis for Proposal 1 in
`theory_first_low_compute_proposals.md` by:
1. Deriving metric ceiling quantities from partial-label assumptions.
2. Estimating coverage proxies from existing mapping and overlap reports.
3. Re-interpreting observed F1/AUPR values as fractions of estimated ceilings.
4. Propagating uncertainty in the coverage proxy c via Beta sampling.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


GAMMA = 0.5772156649015329
RNG_SEED = 7
BETA_SAMPLES = 12000


METRIC_TABLES = {
    "probe_priors": "network_inference/outputs/score_eval_probe_priors.csv",
    "grn_baselines_immune": "network_inference/outputs/score_eval_grn_baselines_immune.csv",
}

MAPPING_REPORTS = {
    "probe_priors": "network_inference/outputs/score_eval_probe_priors_missing_report.json",
    "probe_priors_full_genes": "network_inference/outputs/score_eval_probe_priors_full_genes_missing_report.json",
    "probe_priors_full_genes_crosswalk": "network_inference/outputs/score_eval_probe_priors_full_genes_crosswalk_missing_report.json",
    "probe_priors_full_genes_omnipath": "network_inference/outputs/score_eval_probe_priors_full_genes_omnipath_missing_report.json",
}


@dataclass
class CoveragePosterior:
    alpha: float
    beta: float


def find_repo_root(start: Path) -> Path:
    """Find repository root by walking up to directory with known top-level folders."""
    for candidate in [start] + list(start.parents):
        if (candidate / "market_research").exists() and (candidate / "network_inference").exists():
            return candidate
    raise RuntimeError("Could not infer repository root from script location.")


def harmonic_over_n(t: np.ndarray) -> np.ndarray:
    """Return H_t / t using an asymptotic approximation suitable for large t.

    t may be non-integer here because we use latent positive-count estimates.
    For t >= 1, this approximation is accurate enough for our reporting scale.
    """
    t = np.maximum(np.asarray(t, dtype=float), 1.0)
    h_t = np.log(t) + GAMMA + (1.0 / (2.0 * t)) - (1.0 / (12.0 * t * t))
    return h_t / t


def f1_ceiling(c: np.ndarray) -> np.ndarray:
    """Observed F1 ceiling for a perfect latent predictor under coverage c."""
    c = np.clip(np.asarray(c, dtype=float), 1e-9, 1.0)
    return (2.0 * c) / (1.0 + c)


def aupr_ceiling(c: np.ndarray, latent_true_edges: np.ndarray) -> np.ndarray:
    """Observed AUPR ceiling under ideal latent ranking and partial labels.

    Large-sample limit is AUPR_max ~= c.
    The second term is a finite-latent-positive correction through H_t / t.
    """
    c = np.clip(np.asarray(c, dtype=float), 1e-9, 1.0)
    latent_true_edges = np.maximum(np.asarray(latent_true_edges, dtype=float), 1.0)
    ceiling = c + (1.0 - c) * harmonic_over_n(latent_true_edges)
    return np.clip(ceiling, 0.0, 1.0)


def normalize_bool(series: pd.Series) -> pd.Series:
    """Normalize mixed bool/string columns to boolean."""
    if series.dtype == bool:
        return series
    lowered = series.astype(str).str.lower()
    return lowered.isin({"true", "1", "yes"})


def load_metric_rows(repo_root: Path) -> pd.DataFrame:
    """Load and filter rows with observed metrics from existing output tables."""
    frames: List[pd.DataFrame] = []
    for source_name, rel_path in METRIC_TABLES.items():
        path = repo_root / rel_path
        if not path.exists():
            raise FileNotFoundError(f"Required metric table not found: {path}")
        frame = pd.read_csv(path)
        frame["source_table"] = source_name
        frame["reference_excluded"] = normalize_bool(frame["reference_excluded"])
        frame = frame[~frame["reference_excluded"]].copy()
        frame = frame[frame["f1"].notna() & frame["aupr"].notna()].copy()
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    combined["coverage_proxy_c"] = np.clip(combined["ref_node_overlap_pct"] / 100.0, 1e-6, 0.999999)
    combined["observed_base_rate"] = combined["true_edges"] / combined["candidate_edges"]
    combined["latent_true_edges_hat"] = combined["true_edges"] / combined["coverage_proxy_c"]
    combined["latent_precision_hat"] = np.clip(
        combined["precision"] / combined["coverage_proxy_c"], 0.0, 1.0
    )
    return combined


def load_mapping_proxy_rows(repo_root: Path) -> pd.DataFrame:
    """Estimate endpoint-level mapping coverage proxies from missing reports."""
    rows: List[Dict[str, object]] = []
    for policy_name, rel_path in MAPPING_REPORTS.items():
        report_path = repo_root / rel_path
        if not report_path.exists():
            raise FileNotFoundError(f"Required mapping report not found: {report_path}")
        payload = json.loads(report_path.read_text())
        for record in payload:
            denom = 2.0 * float(record["total_edges_normalized"])
            missing_total = float(record["missing_total"])
            endpoint_coverage = np.nan
            if denom > 0:
                endpoint_coverage = float(np.clip(1.0 - (missing_total / denom), 0.0, 1.0))
            rows.append(
                {
                    "proxy_kind": "mapping_endpoint",
                    "proxy_source": policy_name,
                    "reference": record["reference"],
                    "total_edges_normalized": record["total_edges_normalized"],
                    "missing_total": record["missing_total"],
                    "support_num": max(denom - missing_total, 0.0),
                    "support_den": denom,
                    "coverage_proxy_c": endpoint_coverage,
                }
            )
    return pd.DataFrame(rows)


def build_overlap_proxy_rows(metric_rows: pd.DataFrame) -> pd.DataFrame:
    """Extract node-overlap coverage proxies (used as c estimates for ceilings)."""
    dedup = (
        metric_rows[
            [
                "source_table",
                "reference",
                "ref_nodes",
                "overlap_nodes",
                "ref_node_overlap_pct",
                "coverage_proxy_c",
            ]
        ]
        .drop_duplicates()
        .copy()
    )
    dedup = dedup.rename(
        columns={
            "source_table": "proxy_source",
            "ref_nodes": "support_den",
            "overlap_nodes": "support_num",
        }
    )
    dedup["proxy_kind"] = "overlap_node"
    dedup["total_edges_normalized"] = np.nan
    dedup["missing_total"] = np.nan
    return dedup[
        [
            "proxy_kind",
            "proxy_source",
            "reference",
            "total_edges_normalized",
            "missing_total",
            "support_num",
            "support_den",
            "coverage_proxy_c",
        ]
    ]


def estimate_coverage_posterior(metric_rows: pd.DataFrame) -> Dict[str, CoveragePosterior]:
    """Build Beta posteriors for c by reference from overlap counts."""
    posteriors: Dict[str, CoveragePosterior] = {}
    ref_rows = metric_rows[["reference", "overlap_nodes", "ref_nodes"]].drop_duplicates()
    for record in ref_rows.itertuples(index=False):
        overlap_nodes = float(record.overlap_nodes)
        ref_nodes = float(record.ref_nodes)
        failures = max(ref_nodes - overlap_nodes, 0.0)
        posteriors[record.reference] = CoveragePosterior(
            alpha=overlap_nodes + 1.0,
            beta=failures + 1.0,
        )
    return posteriors


def attach_ceiling_estimates(metric_rows: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Compute point estimates and uncertainty intervals for metric ceilings."""
    rows = metric_rows.copy()

    rows["f1_ceiling"] = f1_ceiling(rows["coverage_proxy_c"].to_numpy())
    rows["aupr_ceiling"] = aupr_ceiling(
        rows["coverage_proxy_c"].to_numpy(), rows["latent_true_edges_hat"].to_numpy()
    )

    rows["f1_ceiling_ratio"] = rows["f1"] / rows["f1_ceiling"]
    rows["aupr_ceiling_ratio"] = rows["aupr"] / rows["aupr_ceiling"]

    # AUPR uplift ratio normalizes for observed positive base rate baseline.
    uplift_denom = np.maximum(rows["aupr_ceiling"] - rows["observed_base_rate"], 1e-12)
    rows["aupr_uplift_ratio"] = (rows["aupr"] - rows["observed_base_rate"]) / uplift_denom

    # Monte Carlo uncertainty for each row, using reference-level c posterior.
    posteriors = estimate_coverage_posterior(rows)
    sampled_coverages = {
        reference: np.clip(
            rng.beta(posterior.alpha, posterior.beta, size=BETA_SAMPLES),
            1e-6,
            0.999999,
        )
        for reference, posterior in posteriors.items()
    }
    ci_low = []
    ci_high = []
    ci_low_ap = []
    ci_high_ap = []
    c_ci_low = []
    c_ci_high = []
    ratio_low_f1 = []
    ratio_high_f1 = []
    ratio_low_ap = []
    ratio_high_ap = []

    for record in rows.itertuples(index=False):
        c_samples = sampled_coverages[record.reference]

        latent_true_samples = record.true_edges / c_samples
        f1_samples = f1_ceiling(c_samples)
        ap_samples = aupr_ceiling(c_samples, latent_true_samples)

        ci_low.append(float(np.quantile(f1_samples, 0.025)))
        ci_high.append(float(np.quantile(f1_samples, 0.975)))
        ci_low_ap.append(float(np.quantile(ap_samples, 0.025)))
        ci_high_ap.append(float(np.quantile(ap_samples, 0.975)))
        c_ci_low.append(float(np.quantile(c_samples, 0.025)))
        c_ci_high.append(float(np.quantile(c_samples, 0.975)))

        f1_ratio_samples = record.f1 / f1_samples
        ap_ratio_samples = record.aupr / ap_samples
        ratio_low_f1.append(float(np.quantile(f1_ratio_samples, 0.025)))
        ratio_high_f1.append(float(np.quantile(f1_ratio_samples, 0.975)))
        ratio_low_ap.append(float(np.quantile(ap_ratio_samples, 0.025)))
        ratio_high_ap.append(float(np.quantile(ap_ratio_samples, 0.975)))

    rows["f1_ceiling_ci_low"] = ci_low
    rows["f1_ceiling_ci_high"] = ci_high
    rows["aupr_ceiling_ci_low"] = ci_low_ap
    rows["aupr_ceiling_ci_high"] = ci_high_ap
    rows["coverage_proxy_ci_low"] = c_ci_low
    rows["coverage_proxy_ci_high"] = c_ci_high
    rows["f1_ceiling_ratio_ci_low"] = ratio_low_f1
    rows["f1_ceiling_ratio_ci_high"] = ratio_high_f1
    rows["aupr_ceiling_ratio_ci_low"] = ratio_low_ap
    rows["aupr_ceiling_ratio_ci_high"] = ratio_high_ap

    return rows


def build_reference_summary(rows: pd.DataFrame) -> pd.DataFrame:
    """Aggregate reinterpretation metrics by reference."""
    grouped = (
        rows.groupby("reference", as_index=False)
        .agg(
            rows=("reference", "size"),
            methods=("method", "nunique"),
            c_proxy_median=("coverage_proxy_c", "median"),
            c_proxy_min=("coverage_proxy_c", "min"),
            c_proxy_max=("coverage_proxy_c", "max"),
            observed_f1_median=("f1", "median"),
            observed_aupr_median=("aupr", "median"),
            f1_ceiling_median=("f1_ceiling", "median"),
            aupr_ceiling_median=("aupr_ceiling", "median"),
            f1_ratio_median=("f1_ceiling_ratio", "median"),
            aupr_ratio_median=("aupr_ceiling_ratio", "median"),
            aupr_uplift_ratio_median=("aupr_uplift_ratio", "median"),
            observed_baserate_median=("observed_base_rate", "median"),
        )
        .sort_values(["aupr_ratio_median", "f1_ratio_median"], ascending=False)
    )
    return grouped


def _color_map_for_references(references: Iterable[str]) -> Dict[str, tuple]:
    unique_refs = sorted(set(references))
    cmap = plt.get_cmap("tab10")
    return {ref: cmap(i % 10) for i, ref in enumerate(unique_refs)}


def plot_observed_vs_ceiling_ratio(rows: pd.DataFrame, output_path: Path) -> None:
    """Scatter observed metrics against estimated ceiling ratios with CI bars."""
    colors = _color_map_for_references(rows["reference"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=160)

    # Panel 1: F1 vs F1/ceiling.
    x_f1 = np.maximum(rows["f1"].to_numpy(), 1e-8)
    y_f1 = rows["f1_ceiling_ratio"].to_numpy()
    yerr_f1 = np.vstack(
        [
            y_f1 - rows["f1_ceiling_ratio_ci_low"].to_numpy(),
            rows["f1_ceiling_ratio_ci_high"].to_numpy() - y_f1,
        ]
    )
    axes[0].errorbar(x_f1, y_f1, yerr=yerr_f1, fmt="none", ecolor="0.75", alpha=0.5, zorder=1)
    for record in rows.itertuples(index=False):
        axes[0].scatter(
            max(record.f1, 1e-8),
            record.f1_ceiling_ratio,
            color=colors[record.reference],
            s=36,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.3,
            zorder=2,
        )
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Observed F1")
    axes[0].set_ylabel("Observed F1 / Estimated F1 Ceiling")
    axes[0].set_title("F1 Ceiling Ratio")
    axes[0].grid(True, alpha=0.25)

    # Panel 2: AUPR vs AUPR/ceiling.
    x_ap = np.maximum(rows["aupr"].to_numpy(), 1e-8)
    y_ap = rows["aupr_ceiling_ratio"].to_numpy()
    yerr_ap = np.vstack(
        [
            y_ap - rows["aupr_ceiling_ratio_ci_low"].to_numpy(),
            rows["aupr_ceiling_ratio_ci_high"].to_numpy() - y_ap,
        ]
    )
    axes[1].errorbar(x_ap, y_ap, yerr=yerr_ap, fmt="none", ecolor="0.75", alpha=0.5, zorder=1)
    for record in rows.itertuples(index=False):
        axes[1].scatter(
            max(record.aupr, 1e-8),
            record.aupr_ceiling_ratio,
            color=colors[record.reference],
            s=36,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.3,
            zorder=2,
        )
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Observed AUPR")
    axes[1].set_ylabel("Observed AUPR / Estimated AUPR Ceiling")
    axes[1].set_title("AUPR Ceiling Ratio")
    axes[1].grid(True, alpha=0.25)

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=color,
            markeredgecolor="black",
            markeredgewidth=0.3,
            markersize=6,
            label=ref,
        )
        for ref, color in colors.items()
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=min(len(legend_handles), 4))
    fig.suptitle("Proposal 1: Observed Metrics Relative to Partial-Label Ceilings", y=1.02)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_theoretical_curves(rows: pd.DataFrame, output_path: Path) -> None:
    """Plot theoretical ceiling curves across c and representative base-rate lines."""
    c_grid = np.linspace(0.02, 0.98, 400)
    f1_grid = f1_ceiling(c_grid)

    # Use a representative latent-positive count from the observed tables.
    median_latent_true = float(np.median(rows["latent_true_edges_hat"]))
    ap_grid = aupr_ceiling(c_grid, np.full_like(c_grid, median_latent_true))

    # Random AUPR baselines under observed labels: b_obs = c * b_latent.
    latent_base_rates = [1e-5, 5e-5, 1e-4, 5e-4]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=160)
    axes[0].plot(c_grid, f1_grid, color="#1f77b4", lw=2.0)
    axes[0].set_title("F1 Ceiling vs Coverage c")
    axes[0].set_xlabel("Coverage c")
    axes[0].set_ylabel("Max Observed F1")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(c_grid, ap_grid, color="#d62728", lw=2.0, label="AUPR ceiling")
    for base_rate in latent_base_rates:
        axes[1].plot(
            c_grid,
            c_grid * base_rate,
            linestyle="--",
            lw=1.2,
            label=f"Random baseline (latent b={base_rate:.0e})",
        )
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


def main() -> None:
    script_path = Path(__file__).resolve()
    repo_root = find_repo_root(script_path)
    output_dir = script_path.parent.parent / "outputs"
    figures_dir = output_dir / "figures"

    rng = np.random.default_rng(RNG_SEED)

    metric_rows = load_metric_rows(repo_root)
    metric_rows = attach_ceiling_estimates(metric_rows, rng)

    mapping_proxy_rows = load_mapping_proxy_rows(repo_root)
    overlap_proxy_rows = build_overlap_proxy_rows(metric_rows)
    coverage_proxy_summary = pd.concat(
        [mapping_proxy_rows, overlap_proxy_rows], ignore_index=True
    )

    reference_summary = build_reference_summary(metric_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    metric_rows.to_csv(output_dir / "metric_ceiling_reinterpretation.csv", index=False)
    reference_summary.to_csv(output_dir / "reference_ceiling_summary.csv", index=False)
    coverage_proxy_summary.to_csv(output_dir / "coverage_proxy_summary.csv", index=False)

    plot_observed_vs_ceiling_ratio(
        metric_rows, figures_dir / "observed_metric_vs_ceiling_ratio.png"
    )
    plot_theoretical_curves(metric_rows, figures_dir / "theoretical_ceiling_curves.png")

    print("Wrote:")
    print(f"- {output_dir / 'metric_ceiling_reinterpretation.csv'}")
    print(f"- {output_dir / 'reference_ceiling_summary.csv'}")
    print(f"- {output_dir / 'coverage_proxy_summary.csv'}")
    print(f"- {figures_dir / 'observed_metric_vs_ceiling_ratio.png'}")
    print(f"- {figures_dir / 'theoretical_ceiling_curves.png'}")


if __name__ == "__main__":
    main()
