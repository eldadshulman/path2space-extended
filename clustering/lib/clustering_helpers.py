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

    Size→name mapping matches the published reference
    (`spatial_clusters/original_meta_tcga.csv`, 562 / 218 / 62 / ... / ...):

      biggest → Immune-Modulated
      2nd     → Immune-Inactive
      3rd     → Proliferation-Enriched
      4th     → Metacluster 1
      5th     → Metacluster 2
    """
    Z = linkage(X, method="ward")
    labs = fcluster(Z, t=5, criterion="maxclust")

    size_order = pd.Series(labs).value_counts().index.tolist()
    name_order = [
        "Immune-Modulated",
        "Immune-Inactive",
        "Proliferation-Enriched",
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


# ---------------------------------------------------------------------------
# SpatioType projection onto an external cohort (Aitchison-distance method).
# Matches METABRIC/scr/spatiotypes_aitchison.r — the script the paper uses
# to assign METABRIC patients to TCGA-derived SpatioTypes.

def _closure(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    return X / X.sum(axis=-1, keepdims=True)


def _clr(X: np.ndarray) -> np.ndarray:
    """Centered log-ratio transform — row-wise."""
    X = np.asarray(X, dtype=float)
    log_X = np.log(X)
    return log_X - log_X.mean(axis=-1, keepdims=True)


def build_spatiotype_centroids(tcga_props: pd.DataFrame,
                                spatiotype_col: str = "SpatioType",
                                cluster_cols: Sequence[str] | None = None,
                                pseudocount: float = 1e-6,
                                keep: Sequence[str] = ("Immune-Modulated",
                                                       "Immune-Inactive",
                                                       "Proliferative")) -> pd.DataFrame:
    """Geometric-mean (CLR) centroids of the reference SpatioTypes.

    Mirrors `prepare_projection_model_aitchison()` in the R script:
    apply a pseudocount, closure-normalize, take the per-SpatioType
    mean in CLR space."""
    if cluster_cols is None:
        cluster_cols = [c for c in tcga_props.columns if c.startswith("Cluster")]
    sub = tcga_props[tcga_props[spatiotype_col].isin(list(keep))].copy()
    X = _clr(_closure(sub[cluster_cols].values + pseudocount))
    sub = sub.assign(**{c: X[:, i] for i, c in enumerate(cluster_cols)})
    cent = sub.groupby(spatiotype_col)[list(cluster_cols)].mean()
    return cent.loc[[s for s in keep if s in cent.index]]


def project_spatiotypes(query_props: pd.DataFrame,
                         centroids: pd.DataFrame,
                         cluster_cols: Sequence[str] | None = None,
                         pseudocount: float = 1e-6) -> pd.DataFrame:
    """Assign each query patient to the nearest TCGA centroid in Aitchison
    (CLR) space. Returns (SpatioType, confidence)."""
    cent_cols = list(centroids.columns)
    if cluster_cols is None:
        # Tolerate "Cluster 1" vs "Cluster_1" by string matching.
        cluster_cols = []
        for cc in cent_cols:
            matches = [c for c in query_props.columns if c.replace("_", " ") == cc.replace("_", " ")]
            if matches:
                cluster_cols.append(matches[0])
        if len(cluster_cols) != len(cent_cols):
            raise ValueError("Could not align query cluster columns to centroid columns.")
    Q = _clr(_closure(query_props[cluster_cols].values + pseudocount))
    C = centroids.values  # already in CLR space
    # Pairwise Euclidean distance Q × C
    dists = np.linalg.norm(Q[:, None, :] - C[None, :, :], axis=-1)
    nearest = dists.argmin(axis=1)
    sorted_d = np.sort(dists, axis=1)
    confidence = (sorted_d[:, 1] - sorted_d[:, 0]) / sorted_d[:, 1]
    return pd.DataFrame({
        "SpatioType":  centroids.index.values[nearest],
        "confidence":  confidence,
    }, index=query_props.index)
