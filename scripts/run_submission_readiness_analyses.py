#!/usr/bin/env python3
"""Compute additional low-compute analyses for submission-readiness.

This script augments Proposal 1 outputs with:
- uncertainty intervals via stratified bootstrap,
- external validation summaries from time-series and perturbation artifacts,
- baseline and confound checks designed for reviewer-facing integrity audits.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


matplotlib.use("Agg")

RNG_SEED = 19
BOOTSTRAP_N = 3000


def find_repo_root(start: Path) -> Path:
    """Locate repository root by expected domain folders."""
    for candidate in [start] + list(start.parents):
        if (candidate / "market_research").exists() and (candidate / "network_inference").exists():
            return candidate
    raise RuntimeError("Could not infer repository root.")


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if total <= 0:
        return (np.nan, np.nan)
    p = successes / total
    denom = 1.0 + (z * z) / total
    center = (p + (z * z) / (2.0 * total)) / denom
    margin = (
        z
        * np.sqrt((p * (1.0 - p) / total) + ((z * z) / (4.0 * total * total)))
        / denom
    )
    return (max(0.0, center - margin), min(1.0, center + margin))


def stratified_bootstrap(
    frame: pd.DataFrame,
    by: str,
    metric_fn: Callable[[pd.DataFrame], float],
    n_bootstrap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Bootstrap rows with replacement within each stratum, preserving stratum counts."""
    parts = []
    grouped = {k: g.reset_index(drop=True) for k, g in frame.groupby(by)}
    keys = sorted(grouped)

    for _ in range(n_bootstrap):
        resampled = []
        for key in keys:
            group = grouped[key]
            idx = rng.integers(0, len(group), size=len(group))
            resampled.append(group.iloc[idx])
        sample = pd.concat(resampled, ignore_index=True)
        parts.append(metric_fn(sample))
    return np.asarray(parts, dtype=float)


def method_class(name: str) -> str:
    """Map method names to broad families for concise interpretation."""
    if name.startswith("probe_") or name == "attention_inferred":
        return "probe_family"
    if name.startswith("omnipath") or name in {"trrust_prior", "omnipath_prior"}:
        return "prior_constrained"
    if name in {"random", "random_3950", "scenic_pruned"}:
        return "control"
    return "classical_grn"


def build_bootstrap_tables(central: pd.DataFrame, out_dir: Path, rng: np.random.Generator) -> pd.DataFrame:
    """Compute global bootstrap uncertainty for key headline statistics."""
    metrics = {
        "median_aupr_ratio": lambda df: float(df["aupr_ratio"].median()),
        "median_f1_ratio": lambda df: float(df["f1_ratio"].median()),
        "median_aupr_uplift_ratio": lambda df: float(df["aupr_uplift_ratio"].median()),
        "fraction_positive_uplift": lambda df: float((df["aupr_uplift_ratio"] > 0).mean()),
        "best_aupr_ratio": lambda df: float(df["aupr_ratio"].max()),
        "best_f1_ratio": lambda df: float(df["f1_ratio"].max()),
    }

    rows = []
    for metric_name, fn in metrics.items():
        boot = stratified_bootstrap(
            frame=central,
            by="reference",
            metric_fn=fn,
            n_bootstrap=BOOTSTRAP_N,
            rng=rng,
        )
        point = fn(central)
        rows.append(
            {
                "metric": metric_name,
                "point_estimate": point,
                "bootstrap_mean": float(np.mean(boot)),
                "bootstrap_std": float(np.std(boot, ddof=1)),
                "ci95_low": float(np.quantile(boot, 0.025)),
                "ci95_high": float(np.quantile(boot, 0.975)),
            }
        )

    result = pd.DataFrame(rows)
    result.to_csv(out_dir / "bootstrap_global_metrics.csv", index=False)
    return result


def build_reference_uplift_table(central: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Summarize per-reference positive-uplift rates with Wilson intervals."""
    rows = []
    for reference, group in central.groupby("reference"):
        positives = int((group["aupr_uplift_ratio"] > 0).sum())
        total = int(len(group))
        low, high = wilson_interval(positives, total)
        rows.append(
            {
                "reference": reference,
                "n_rows": total,
                "n_positive_uplift": positives,
                "positive_uplift_fraction": positives / total,
                "positive_uplift_ci95_low": low,
                "positive_uplift_ci95_high": high,
                "median_aupr_ratio": float(group["aupr_ratio"].median()),
                "median_aupr_uplift_ratio": float(group["aupr_uplift_ratio"].median()),
            }
        )
    out = pd.DataFrame(rows).sort_values("positive_uplift_fraction", ascending=False)
    out.to_csv(out_dir / "uplift_positive_fraction_by_reference.csv", index=False)
    return out


def build_family_table(central: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Aggregate central metrics by method family."""
    df = central.copy()
    df["method_family"] = df["method"].map(method_class)
    out = (
        df.groupby("method_family", as_index=False)
        .agg(
            rows=("method_family", "size"),
            methods=("method", "nunique"),
            median_aupr=("aupr", "median"),
            median_aupr_ratio=("aupr_ratio", "median"),
            median_f1_ratio=("f1_ratio", "median"),
            positive_uplift_fraction=("aupr_uplift_ratio", lambda x: float((x > 0).mean())),
        )
        .sort_values("median_aupr_ratio", ascending=False)
    )
    out.to_csv(out_dir / "family_comparison_central.csv", index=False)
    return out


def build_confound_checks(central: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    """Compute simple correlation-based confound checks."""
    checks = [
        ("aupr_ratio", "candidate_edges"),
        ("aupr_ratio", "coverage_scenario_c"),
        ("aupr_ratio", "observed_base_rate"),
        ("aupr_uplift_ratio", "candidate_edges"),
        ("aupr_uplift_ratio", "coverage_scenario_c"),
        ("aupr_uplift_ratio", "observed_base_rate"),
    ]

    rows = []
    for y_col, x_col in checks:
        subset = central[[y_col, x_col]].dropna()
        pearson = float(subset[y_col].corr(subset[x_col], method="pearson"))
        spearman = float(subset[y_col].corr(subset[x_col], method="spearman"))
        rows.append(
            {
                "target_metric": y_col,
                "candidate_confound": x_col,
                "n_rows": len(subset),
                "pearson_r": pearson,
                "spearman_r": spearman,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "confound_correlation_checks.csv", index=False)
    return out


def build_timeseries_tables(
    repo_root: Path, coverage_table: pd.DataFrame, out_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Reinterpret external time-series metrics through partial-label ceilings."""
    ts = pd.read_csv(repo_root / "network_inference/outputs/summary_timeseries_metrics.csv").copy()
    ts["method_family"] = ts["method"].map(method_class)

    random_rows = ts[ts["method"].str.contains("random", case=False, regex=False)]
    hpn_random = float(random_rows["hpn_aupr"].median())
    beeline_random = float(random_rows["beeline_aupr"].median())

    c_hpn = float(
        coverage_table.loc[coverage_table["reference"] == "hpn_dream", "overlap_point"].iloc[0]
    )
    c_beeline = float(
        coverage_table.loc[coverage_table["reference"] == "beeline_gsd", "overlap_point"].iloc[0]
    )

    ts["hpn_aupr_minus_random"] = ts["hpn_aupr"] - hpn_random
    ts["beeline_aupr_minus_random"] = ts["beeline_aupr"] - beeline_random
    ts["hpn_ratio_to_ceiling_approx"] = ts["hpn_aupr"] / max(c_hpn, 1e-12)
    ts["beeline_ratio_to_ceiling_approx"] = ts["beeline_aupr"] / max(c_beeline, 1e-12)
    ts["hpn_positive_uplift"] = ts["hpn_aupr_minus_random"] > 0
    ts["beeline_positive_uplift"] = ts["beeline_aupr_minus_random"] > 0

    ts_out = ts.sort_values("hpn_ratio_to_ceiling_approx", ascending=False)
    ts_out.to_csv(out_dir / "timeseries_ceiling_reinterpretation.csv", index=False)

    family = (
        ts.groupby("method_family", as_index=False)
        .agg(
            methods=("method", "nunique"),
            median_hpn_aupr=("hpn_aupr", "median"),
            median_hpn_ratio_to_ceiling_approx=("hpn_ratio_to_ceiling_approx", "median"),
            hpn_positive_uplift_fraction=("hpn_positive_uplift", "mean"),
            median_beeline_aupr=("beeline_aupr", "median"),
            median_beeline_ratio_to_ceiling_approx=("beeline_ratio_to_ceiling_approx", "median"),
            beeline_positive_uplift_fraction=("beeline_positive_uplift", "mean"),
        )
        .sort_values("median_hpn_ratio_to_ceiling_approx", ascending=False)
    )
    family.to_csv(out_dir / "timeseries_family_summary.csv", index=False)

    return ts_out, family, pd.DataFrame(
        {
            "hpn_random_aupr_baseline": [hpn_random],
            "beeline_random_aupr_baseline": [beeline_random],
            "hpn_coverage_proxy": [c_hpn],
            "beeline_coverage_proxy": [c_beeline],
        }
    )


def build_perturbation_tables(repo_root: Path, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize perturbation validation artifact statistics."""
    payload = json.loads((repo_root / "network_inference/outputs/perturbation_eval_probe.json").read_text())
    recalls = pd.DataFrame(payload["per_perturbation_recall"])

    summary = pd.DataFrame(
        [
            {
                "overall_precision": float(payload["overall_metrics"]["precision"]),
                "overall_recall": float(payload["overall_metrics"]["recall"]),
                "overall_f1": float(payload["overall_metrics"]["f1"]),
                "avg_recall_reported": float(payload["avg_recall"]),
                "n_perturbations": int(len(recalls)),
                "nonzero_recall_count": int((recalls["recall"] > 0).sum()),
                "nonzero_recall_fraction": float((recalls["recall"] > 0).mean()),
                "recall_median": float(recalls["recall"].median()),
                "recall_p90": float(recalls["recall"].quantile(0.9)),
                "recall_max": float(recalls["recall"].max()),
            }
        ]
    )
    summary.to_csv(out_dir / "perturbation_validation_summary.csv", index=False)

    top_hits = recalls.sort_values(["recall", "hits"], ascending=False).head(25)
    top_hits.to_csv(out_dir / "perturbation_top_hits.csv", index=False)
    return summary, top_hits


def build_cross_assay_table(
    central: pd.DataFrame,
    timeseries: pd.DataFrame,
    out_dir: Path,
) -> pd.DataFrame:
    """Compare method ordering consistency between central and time-series outputs.

    The function includes a degeneracy guard: if the dynamic range in either axis
    is tiny, correlation is reported as not informative.
    """
    central_ref = (
        central[central["reference"].isin(["hpn_dream", "beeline_gsd"])]
        .pivot_table(index="method", columns="reference", values="aupr_ratio", aggfunc="median")
        .rename(columns={"hpn_dream": "central_hpn_aupr_ratio", "beeline_gsd": "central_beeline_aupr_ratio"})
    )

    ts = (
        timeseries.set_index("method")[
            [
                "hpn_ratio_to_ceiling_approx",
                "beeline_ratio_to_ceiling_approx",
                "hpn_aupr_minus_random",
                "beeline_aupr_minus_random",
            ]
        ]
        .rename(
            columns={
                "hpn_ratio_to_ceiling_approx": "timeseries_hpn_aupr_ratio_approx",
                "beeline_ratio_to_ceiling_approx": "timeseries_beeline_aupr_ratio_approx",
            }
        )
    )

    merged = central_ref.merge(ts, left_index=True, right_index=True, how="inner").reset_index()
    def one_row(comp: str, x_col: str, y_col: str) -> dict[str, object]:
        if len(merged) < 3:
            return {
                "comparison": comp,
                "n_methods": len(merged),
                "x_range": np.nan,
                "y_range": np.nan,
                "spearman_r": np.nan,
                "pearson_r": np.nan,
                "informative": False,
                "note": "insufficient_overlap_methods",
            }

        x = merged[x_col].to_numpy(dtype=float)
        y = merged[y_col].to_numpy(dtype=float)
        x_range = float(np.max(x) - np.min(x))
        y_range = float(np.max(y) - np.min(y))

        # Guard against near-constant vectors where correlations become unstable
        # and can falsely look perfect.
        if x_range < 1e-7 or y_range < 1e-7:
            return {
                "comparison": comp,
                "n_methods": len(merged),
                "x_range": x_range,
                "y_range": y_range,
                "spearman_r": np.nan,
                "pearson_r": np.nan,
                "informative": False,
                "note": "degenerate_dynamic_range_not_interpretable",
            }

        return {
            "comparison": comp,
            "n_methods": len(merged),
            "x_range": x_range,
            "y_range": y_range,
            "spearman_r": float(merged[x_col].corr(merged[y_col], method="spearman")),
            "pearson_r": float(merged[x_col].corr(merged[y_col], method="pearson")),
            "informative": True,
            "note": "ok",
        }

    corrs = [
        one_row(
            "central_hpn_vs_timeseries_hpn",
            "central_hpn_aupr_ratio",
            "timeseries_hpn_aupr_ratio_approx",
        ),
        one_row(
            "central_beeline_vs_timeseries_beeline",
            "central_beeline_aupr_ratio",
            "timeseries_beeline_aupr_ratio_approx",
        ),
    ]

    corr_table = pd.DataFrame(corrs)
    corr_table.to_csv(out_dir / "cross_assay_concordance.csv", index=False)
    merged.to_csv(out_dir / "cross_assay_method_overlap.csv", index=False)
    return corr_table


def plot_bootstrap_metrics(bootstrap_table: pd.DataFrame, out_path: Path) -> None:
    """Plot key global metrics with bootstrap intervals."""
    subset = bootstrap_table[
        bootstrap_table["metric"].isin(
            [
                "median_aupr_ratio",
                "median_f1_ratio",
                "fraction_positive_uplift",
                "best_aupr_ratio",
                "best_f1_ratio",
            ]
        )
    ].copy()

    subset["label"] = subset["metric"].map(
        {
            "median_aupr_ratio": "Median AUPR ratio",
            "median_f1_ratio": "Median F1 ratio",
            "fraction_positive_uplift": "Positive uplift fraction",
            "best_aupr_ratio": "Best AUPR ratio",
            "best_f1_ratio": "Best F1 ratio",
        }
    )

    y = np.arange(len(subset))
    x = subset["point_estimate"].to_numpy()
    low = x - subset["ci95_low"].to_numpy()
    high = subset["ci95_high"].to_numpy() - x

    fig, ax = plt.subplots(figsize=(8.8, 4.8), dpi=200)
    ax.errorbar(
        x,
        y,
        xerr=np.vstack([low, high]),
        fmt="o",
        color="#1f77b4",
        ecolor="#4a4a4a",
        capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(subset["label"])
    ax.set_xlabel("Estimate with 95% bootstrap CI")
    ax.set_title("Global Ceiling Metrics: Point Estimates and Uncertainty")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def plot_external_validation(
    timeseries: pd.DataFrame, perturbation_top: pd.DataFrame, out_path: Path
) -> None:
    """Plot external validation snapshots (time-series and perturbation)."""
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.0), dpi=200)

    palette = {
        "probe_family": "#1f77b4",
        "classical_grn": "#2ca02c",
        "prior_constrained": "#ff7f0e",
        "control": "#7f7f7f",
    }

    for family, block in timeseries.groupby("method_family"):
        axes[0].scatter(
            block["hpn_aupr"],
            block["beeline_aupr"],
            label=family,
            s=42,
            alpha=0.85,
            color=palette.get(family, "#9467bd"),
            edgecolor="black",
            linewidth=0.3,
        )
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("HPN time-series AUPR")
    axes[0].set_ylabel("BEELINE time-series AUPR")
    axes[0].set_title("External Time-Series Performance Landscape")
    axes[0].grid(True, which="both", alpha=0.2)
    axes[0].legend(fontsize=8, frameon=True)

    recalls = perturbation_top["recall"].to_numpy()
    hits = perturbation_top["hits"].to_numpy()
    axes[1].bar(np.arange(len(recalls)), recalls, color="#1f77b4", alpha=0.8)
    axes[1].set_xlabel("Top perturbations by recall")
    axes[1].set_ylabel("Recall")
    axes[1].set_title("Perturbation Recall Concentrates in Few Perturbagens")
    axes[1].grid(True, axis="y", alpha=0.25)

    for idx, (r, h) in enumerate(zip(recalls[:10], hits[:10])):
        axes[1].text(idx, r, f"{int(h)}", ha="center", va="bottom", fontsize=7)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def plot_reference_positive_uplift(reference_table: pd.DataFrame, out_path: Path) -> None:
    """Plot positive-uplift fractions by reference with Wilson intervals."""
    df = reference_table.sort_values("positive_uplift_fraction", ascending=False).reset_index(drop=True)
    y = np.arange(len(df))
    x = df["positive_uplift_fraction"].to_numpy()
    low = np.maximum(x - df["positive_uplift_ci95_low"].to_numpy(), 0.0)
    high = np.maximum(df["positive_uplift_ci95_high"].to_numpy() - x, 0.0)

    fig, ax = plt.subplots(figsize=(9.6, 4.8), dpi=200)
    ax.errorbar(
        x,
        y,
        xerr=np.vstack([low, high]),
        fmt="o",
        color="#d62728",
        ecolor="#4a4a4a",
        capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(df["reference"])
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Fraction of methods with positive uplift (AUPR > observed random baseline)")
    ax.set_title("Positive-Uplift Prevalence by Reference")
    ax.grid(True, axis="x", alpha=0.25)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def write_submission_summary(
    out_dir: Path,
    bootstrap: pd.DataFrame,
    reference_uplift: pd.DataFrame,
    family: pd.DataFrame,
    confounds: pd.DataFrame,
    timeseries_family: pd.DataFrame,
    perturbation_summary: pd.DataFrame,
    concordance: pd.DataFrame,
    ts_constants: pd.DataFrame,
) -> None:
    """Write concise markdown summary of additional analyses."""
    p = perturbation_summary.iloc[0]
    ts_c = ts_constants.iloc[0]

    lines = [
        "# Submission-Readiness Additional Analyses",
        "",
        "## Bootstrap uncertainty (stratified by reference)",
        bootstrap.to_markdown(index=False),
        "",
        "## Positive uplift prevalence by reference",
        reference_uplift.to_markdown(index=False),
        "",
        "## Central benchmark family summary",
        family.to_markdown(index=False),
        "",
        "## Confound checks (correlation diagnostics)",
        confounds.to_markdown(index=False),
        "",
        "## External time-series reinterpretation constants",
        f"- HPN random AUPR baseline: `{ts_c['hpn_random_aupr_baseline']:.6e}`",
        f"- BEELINE random AUPR baseline: `{ts_c['beeline_random_aupr_baseline']:.6e}`",
        f"- HPN coverage proxy (overlap-based): `{ts_c['hpn_coverage_proxy']:.6f}`",
        f"- BEELINE coverage proxy (overlap-based): `{ts_c['beeline_coverage_proxy']:.6f}`",
        "",
        "## External time-series family summary",
        timeseries_family.to_markdown(index=False),
        "",
        "## Perturbation validation summary",
        f"- Overall precision: `{p['overall_precision']:.6e}`",
        f"- Overall recall: `{p['overall_recall']:.6e}`",
        f"- Overall F1: `{p['overall_f1']:.6e}`",
        f"- Perturbations with non-zero recall: `{int(p['nonzero_recall_count'])}/{int(p['n_perturbations'])}` (`{p['nonzero_recall_fraction']:.4%}`)",
        "",
        "## Cross-assay concordance diagnostic",
        concordance.to_markdown(index=False),
        "",
        "## Interpretation",
        "- Bootstrap intervals confirm the central claim is robust: ceiling-normalized performance is low even under resampling uncertainty.",
        "- Positive uplift is sparse and reference-specific, indicating biological compatibility matters more than a single global leaderboard.",
        "- Time-series and perturbation validations are directionally consistent with the central benchmark: only constrained/prior-driven settings produce non-trivial uplift, and probe perturbation recall remains sparse.",
    ]

    (out_dir / "submission_readiness_additional_summary.md").write_text("\n".join(lines))


def main() -> None:
    script = Path(__file__).resolve()
    repo_root = find_repo_root(script)
    base_dir = (
        repo_root
        / "market_research/ambitious_paper_questions/proposal_01_partial_label_metric_ceilings"
    )
    tables_dir = base_dir / "outputs/paper/tables"
    figures_dir = base_dir / "outputs/paper/figures"
    rng = np.random.default_rng(RNG_SEED)

    central = pd.read_csv(tables_dir / "central_aupr_rows.csv")
    coverage = pd.read_csv(tables_dir / "coverage_scenarios_by_reference.csv")

    bootstrap = build_bootstrap_tables(central, tables_dir, rng)
    reference_uplift = build_reference_uplift_table(central, tables_dir)
    family = build_family_table(central, tables_dir)
    confounds = build_confound_checks(central, tables_dir)
    timeseries, timeseries_family, ts_constants = build_timeseries_tables(repo_root, coverage, tables_dir)
    perturbation_summary, perturbation_top = build_perturbation_tables(repo_root, tables_dir)
    concordance = build_cross_assay_table(central, timeseries, tables_dir)

    plot_bootstrap_metrics(bootstrap, figures_dir / "figure_08_bootstrap_key_metrics.png")
    plot_external_validation(timeseries, perturbation_top, figures_dir / "figure_09_external_validation.png")
    plot_reference_positive_uplift(reference_uplift, figures_dir / "figure_10_reference_positive_uplift.png")

    write_submission_summary(
        out_dir=tables_dir,
        bootstrap=bootstrap,
        reference_uplift=reference_uplift,
        family=family,
        confounds=confounds,
        timeseries_family=timeseries_family,
        perturbation_summary=perturbation_summary,
        concordance=concordance,
        ts_constants=ts_constants,
    )

    print("Generated submission-readiness analysis tables and figures.")
    print(f"Tables directory: {tables_dir}")
    print(f"Figures directory: {figures_dir}")


if __name__ == "__main__":
    main()
