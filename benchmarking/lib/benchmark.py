"""Benchmarking helpers — per-gene PCC distributions across methods.

Evaluation utilities for the spatial-transcriptomics method benchmark:
per-gene aggregation, model-by-median ordering for box plots, and a paired
Mann-Whitney U test of Path2Space against each competitor on the matched
per-gene PCC vectors.

This component does **not** train any of the competing methods. The
prediction outputs were generated upstream by the original authors of the
respective methods (see the README). What lives here is the evaluation
layer that consumes those predictions' per-(slide, gene) Pearson correlations
and produces the per-cohort summary plots and Mann-Whitney comparisons.
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

__all__ = ["model_order_by_median",
           "mannwhitney_vs_reference",
           "median_pcc_table"]


def model_order_by_median(df, value_col="cor_pearson", model_col="model",
                          ascending=False):
    """Order models by their median value of `value_col`, descending by default."""
    return (df.groupby(model_col, observed=True)[value_col]
              .median()
              .sort_values(ascending=ascending)
              .index.tolist())


def mannwhitney_vs_reference(df, reference="Path2Space",
                             value_col="cor_pearson", model_col="model",
                             group_cols=("cohort",), alternative="two-sided"):
    """Mann-Whitney U test of each model's per-gene PCC vs the reference model's.

    For each unique combination of `group_cols`, the test compares the
    reference model's per-gene PCC vector against each other model's vector.
    Returns a long DataFrame: one row per (group, model) with U statistic,
    p-value, the reference and competitor medians, and the per-gene n.

    Defaults to a two-sided test (matching the source notebook) — set
    `alternative="greater"` to test whether the reference outperforms the
    competitor.
    """
    group_cols = list(group_cols)
    out = []
    for keys, sub in df.groupby(group_cols, observed=True):
        ref = sub.loc[sub[model_col] == reference, value_col].dropna().to_numpy()
        if len(ref) == 0:
            continue
        for m in sub[model_col].unique():
            if m == reference:
                continue
            x = sub.loc[sub[model_col] == m, value_col].dropna().to_numpy()
            if len(x) == 0:
                continue
            u, p = mannwhitneyu(ref, x, alternative=alternative)
            row = {**dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,))),
                   "model": m, "n_genes_ref": int(len(ref)), "n_genes_model": int(len(x)),
                   "median_ref": float(np.median(ref)),
                   "median_model": float(np.median(x)),
                   "U": float(u), "pvalue": float(p)}
            out.append(row)
    return pd.DataFrame(out)


def median_pcc_table(df, value_col="cor_pearson", model_col="model",
                     group_cols=("cohort",)):
    """Per-(group, model) median PCC — the headline summary table for a panel."""
    return (df.groupby([*group_cols, model_col], observed=True)[value_col]
              .median()
              .unstack(model_col))
