"""Helpers for the SpatioType assignment step.

Used by clustering/notebooks/03_spatiotype_assignment.ipynb. Keeps the
parameter choices (linkage, distance, k, plot ordering, colour palette)
out of the notebook so the notebook reads as a narrative.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import davies_bouldin_score
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# Cluster-validity metrics — matches the paper's main_clustering R script.

def _wss_total(X: np.ndarray, labels: np.ndarray) -> float:
    """Within-cluster sum of squared deviations, summed across features.

    Mirrors `sum((scaled_features - ave(scaled_features, cluster_labels))^2)`
    used in the paper's R notebook.
    """
    total = 0.0
    for k in np.unique(labels):
        pts = X[labels == k]
        total += float(((pts - pts.mean(axis=0)) ** 2).sum())
    return total


def compute_validity_metrics(X: np.ndarray, ks: Sequence[int]) -> pd.DataFrame:
    """Return a (k × {WSS, Davies-Bouldin}) table for hierarchical clusterings."""
    Z = linkage(X, method="ward")
    rows = []
    for k in ks:
        if k == 1:
            labs = np.ones(X.shape[0], dtype=int)
            db = np.nan
        else:
            labs = fcluster(Z, t=k, criterion="maxclust")
            db = float(davies_bouldin_score(X, labs))
        rows.append({"k": int(k), "WSS": _wss_total(X, labs), "DB": db})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# SpatioType assignment

SPATIOTYPE_LEVELS = [
    "Proliferation-Enriched",
    "Immune-Modulated",
    "Immune-Inactive",
    "Lipid-Metabolic",
    "Hedgehog-Active",
]


def assign_spatiotypes(X: np.ndarray) -> pd.DataFrame:
    """Hierarchically cluster `X` at k=5 and label the resulting groups.

    Returns a data frame with columns `Metacluster` and `SpatioType`.
    Mapping is by cluster size, then alphabetical, in the same order as the
    paper:

      biggest → Immune-Modulated
      2nd     → Proliferation-Enriched
      3rd     → Immune-Inactive
      4th     → Hedgehog-Active
      5th     → Lipid-Metabolic
    """
    Z = linkage(X, method="ward")
    labs = fcluster(Z, t=5, criterion="maxclust")

    size_order = pd.Series(labs).value_counts().index.tolist()
    name_order = [
        "Immune-Modulated",
        "Proliferation-Enriched",
        "Immune-Inactive",
        "Hedgehog-Active",
        "Lipid-Metabolic",
    ]
    mapping = dict(zip(size_order, name_order))

    out = pd.DataFrame({
        "Metacluster": labs,
        "SpatioType": pd.Categorical(
            [mapping[m] for m in labs],
            categories=SPATIOTYPE_LEVELS,
            ordered=False,
        ),
    })
    out.attrs["mapping"] = mapping
    out.attrs["linkage"] = Z
    return out


# ---------------------------------------------------------------------------
# Paper colour palette and cluster row ordering

# JAMA palette ("default") used by ggsci::pal_jama in the paper's R notebook.
JAMA_COLORS = [
    "#374E55",  # dark navy   → Proliferation-Enriched
    "#DF8F44",  # orange      → Immune-Modulated
    "#00A1D5",  # light blue  → Immune-Inactive
    "#B24745",  # dark red    → Lipid-Metabolic
    "#79AF97",  # mint green  → Hedgehog-Active
]

SPATIOTYPE_PALETTE = {
    "Proliferation-Enriched": JAMA_COLORS[0],
    "Immune-Modulated":       JAMA_COLORS[1],
    "Immune-Inactive":        JAMA_COLORS[2],
    "Lipid-Metabolic":        JAMA_COLORS[3],
    "Hedgehog-Active":        JAMA_COLORS[4],
}

# Cluster row order in the paper heatmap, frozen by hand.
HEATMAP_ROW_ORDER = [
    "Cluster 9", "Cluster 11", "Cluster 5",
    "Cluster 10", "Cluster 1", "Cluster 3",
    "Cluster 7", "Cluster 6", "Cluster 8",
    "Cluster 2", "Cluster 4",
]


# ---------------------------------------------------------------------------

def scale_features(X) -> np.ndarray:
    """z-score scaling (column-wise)."""
    return StandardScaler().fit_transform(np.asarray(X, dtype=float))
