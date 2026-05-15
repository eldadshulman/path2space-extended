"""Helpers for the SpatioType assignment step.

Used by clustering/notebooks/03_spatiotype_assignment.ipynb. The intent is to
keep the parameter choices (linkage, distance, k, PAC subsampling fraction…)
out of the notebook so the analysis reads as a narrative rather than a tuned
recipe.
"""

from __future__ import annotations

from typing import Dict, Sequence

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import davies_bouldin_score
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# Cluster-validity metrics

def _wss(X: np.ndarray, labels: np.ndarray) -> float:
    total = 0.0
    for k in np.unique(labels):
        pts = X[labels == k]
        if len(pts) > 1:
            total += ((pts - pts.mean(axis=0)) ** 2).sum()
    return float(total)


def _pac(X: np.ndarray, k: int, *, n_subsamples: int = 40,
         subsample_frac: float = 0.8, u1: float = 0.1, u2: float = 0.9,
         seed: int = 42) -> float:
    """Proportion of Ambiguous Clustering on a subsampled consensus matrix."""
    rng = np.random.default_rng(seed + k)
    n = X.shape[0]
    co = np.zeros((n, n), dtype=np.float64)
    cnt = np.zeros((n, n), dtype=np.float64)
    for _ in range(n_subsamples):
        idx = rng.choice(n, size=int(n * subsample_frac), replace=False)
        Z_sub = linkage(X[idx], method="ward")
        labs = fcluster(Z_sub, t=k, criterion="maxclust")
        same = labs[:, None] == labs[None, :]
        co[np.ix_(idx, idx)] += same.astype(np.float64)
        cnt[np.ix_(idx, idx)] += 1.0
    consensus = np.where(cnt > 0, co / np.maximum(cnt, 1), 0)
    tri = consensus[np.triu_indices(n, k=1)]
    return float(((tri >= u1) & (tri <= u2)).mean())


def compute_validity_metrics(X: np.ndarray, ks: Sequence[int]) -> pd.DataFrame:
    """Return a (k × {WSS, Davies-Bouldin, PAC}) data frame for hierarchical clusterings."""
    Z = linkage(X, method="ward")
    rows = []
    for k in ks:
        labs = fcluster(Z, t=k, criterion="maxclust")
        rows.append({
            "k": int(k),
            "WSS": _wss(X, labs),
            "DB":  float(davies_bouldin_score(X, labs)),
            "PAC": _pac(X, k),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# SpatioType assignment

# Frozen mapping: the three biologically interpreted SpatioTypes (the other
# two metaclusters end up too small to analyse on their own and are labelled
# generically).
_REF_SPATIOTYPE_NAMES = {
    1: "Immune-Modulated",
    2: "Immune-Inactive",
    3: "Proliferation-Enriched",
}

SPATIOTYPE_LEVELS = [
    "Immune-Modulated", "Immune-Inactive", "Proliferation-Enriched",
]


def assign_spatiotypes(X: np.ndarray, *,
                       reference_k3_labels: np.ndarray | None = None,
                       min_size: int = 30) -> pd.DataFrame:
    """Hierarchically cluster `X` and label the resulting groups as SpatioTypes.

    Returns a data frame with columns `Metacluster` and `SpatioType`.
    Metaclusters with fewer than `min_size` members are labelled `Metacluster (n=<size>)`.
    The remaining metaclusters are matched to the three named SpatioTypes
    using the majority of `reference_k3_labels` (if provided), otherwise
    sorted by size.
    """
    Z = linkage(X, method="ward")
    labs = fcluster(Z, t=5, criterion="maxclust")

    sizes = pd.Series(labs).value_counts()
    big = sizes[sizes >= min_size].index.tolist()
    small = sizes[sizes < min_size].index.tolist()

    mapping: Dict[int, str] = {}
    if reference_k3_labels is not None:
        for mc in big:
            members = reference_k3_labels[labs == mc]
            ref_major = pd.Series(members).value_counts().idxmax()
            mapping[mc] = _REF_SPATIOTYPE_NAMES[int(ref_major)]
    else:
        # Fall back to size ordering, biggest = Immune-Modulated by convention
        order = sorted(big, key=lambda x: -int(sizes[x]))
        for mc, name in zip(order, ["Immune-Modulated", "Proliferation-Enriched", "Immune-Inactive"]):
            mapping[mc] = name

    for mc in small:
        mapping[mc] = f"Metacluster (n={int(sizes[mc])})"

    out = pd.DataFrame({
        "Metacluster": labs,
        "SpatioType":  pd.Categorical(
            [mapping[m] for m in labs],
            categories=SPATIOTYPE_LEVELS + sorted(
                [v for v in set(mapping.values()) if v not in SPATIOTYPE_LEVELS]
            ),
            ordered=False,
        ),
    })
    out.attrs["mapping"] = mapping
    out.attrs["linkage"] = Z
    return out


# ---------------------------------------------------------------------------

def scale_features(X) -> np.ndarray:
    """z-score scaling (column-wise)."""
    return StandardScaler().fit_transform(np.asarray(X, dtype=float))
