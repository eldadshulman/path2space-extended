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
    "Metacluster 1",
    "Metacluster 2",
]


def assign_spatiotypes(X: np.ndarray) -> pd.DataFrame:
    """Hierarchically cluster `X` at k=5 and label the resulting groups.

    Returns a data frame with columns `Metacluster` and `SpatioType`.
    The three largest groups get the biologically interpreted names; the two
    smallest are left generic (`Metacluster 1`, `Metacluster 2`) because their
    sample sizes are too small to analyse.

      biggest → Immune-Modulated
      2nd     → Proliferation-Enriched
      3rd     → Immune-Inactive
      4th     → Metacluster 1
      5th     → Metacluster 2
    """
    Z = linkage(X, method="ward")
    labs = fcluster(Z, t=5, criterion="maxclust")

    size_order = pd.Series(labs).value_counts().index.tolist()
    name_order = [
        "Immune-Modulated",
        "Proliferation-Enriched",
        "Immune-Inactive",
        "Metacluster 1",
        "Metacluster 2",
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
# Paper convention (main_clustering_R.ipynb): use JAMA colors 1, 3, 4 for the
# three named SpatioTypes, and grey for the two too-small-to-analyse groups.
JAMA_COLORS = [
    "#374E55",  # 1 dark navy
    "#DF8F44",  # 2 orange
    "#00A1D5",  # 3 light blue
    "#B24745",  # 4 dark red
    "#79AF97",  # 5 mint green
    "#6A6599",  # 6 muted purple
    "#80796B",  # 7 olive
]
_SMALL_GREY = "#908D8B"

SPATIOTYPE_PALETTE = {
    "Proliferation-Enriched": JAMA_COLORS[0],
    "Immune-Modulated":       JAMA_COLORS[2],
    "Immune-Inactive":        JAMA_COLORS[3],
    "Metacluster 1":          _SMALL_GREY,
    "Metacluster 2":          _SMALL_GREY,
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
