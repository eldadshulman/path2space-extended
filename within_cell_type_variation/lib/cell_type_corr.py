"""Within-cell-type expression-variation analysis for Path2Space (Figures 4E, 4F, S5G, S5H).

Path2Space predicts spatial gene expression from H&E histology. Beyond predicting
cell-type composition, this analysis asks whether the predictions resolve
expression variation *within* individual cell types — by correlating predicted
and measured expression gene-by-gene inside regions that share a cell-type
identity.

Two stages:

  Stage A — `gene_pearson` (the per-region primitive)
    Within one region (a pathologist-annotated area, a transcriptomic
    neighborhood, ...), for every gene, the Pearson correlation between
    predicted and measured expression across the region's spots.

  Stage B — `summarize_by_cell_type`, `concordance_matrix` (this component's
    notebooks). Aggregate the per-gene correlations to per-cell-type medians
    and gene counts, and measure how concordant the gene-wise accuracy is
    across cell types.

Stage A needs the full per-slide ST dataset (predicted and measured expression
matrices, region annotations); its output — the per-gene / per-cell-type
correlation tables — is what the notebooks consume and is bundled in `data/`.
"""

import numpy as np
import pandas as pd

__all__ = ["gene_pearson", "summarize_by_cell_type", "concordance_matrix"]


def gene_pearson(pred_df, obs_df, genes):
    """Per-gene Pearson correlation between predicted and measured expression.

    pred_df, obs_df : spots x genes DataFrames of predicted / measured
        expression for one region. Correlations use the spots common to both.
    genes : sequence of gene names to score.

    Returns a 1-D array of length len(genes); a gene with zero variance in
    either input yields NaN. Vectorized — one pass over the spots x genes
    matrices rather than a Python loop over genes.
    """
    common = pred_df.index.intersection(obs_df.index)
    if len(common) == 0:
        return np.full(len(genes), np.nan)

    pred = pred_df.loc[common, genes].to_numpy(dtype=float)
    obs = obs_df.loc[common, genes].to_numpy(dtype=float)

    pred -= pred.mean(axis=0, keepdims=True)
    obs -= obs.mean(axis=0, keepdims=True)

    numerator = np.sum(pred * obs, axis=0)
    denom = np.sqrt(np.sum(pred ** 2, axis=0) * np.sum(obs ** 2, axis=0))
    return np.where(denom != 0, numerator / denom, np.nan)


def summarize_by_cell_type(df, value="correlation", group="cell_type", threshold=0.4):
    """Per-cell-type summary of a gene-level correlation table.

    df : long DataFrame with one row per (cell type, gene).
    Returns, per cell type, the median correlation, the number of genes whose
    correlation exceeds `threshold`, and the gene count — sorted by descending
    median.
    """
    summary = df.groupby(group)[value].agg(
        median_pcc="median",
        n_genes_over_threshold=lambda x: int((x > threshold).sum()),
        n_genes="size",
    )
    return summary.sort_values("median_pcc", ascending=False)


def concordance_matrix(df, value="correlation", group="cell_type"):
    """Pairwise concordance of gene-wise accuracy across cell types.

    Pivots `df` to a genes x cell-type matrix of correlations, then returns the
    cell-type-by-cell-type Pearson correlation matrix — high values mean a gene
    that is well predicted in one cell type tends to be well predicted in
    another.
    """
    wide = df.pivot_table(index="gene", columns=group, values=value, aggfunc="mean")
    return wide.corr()
