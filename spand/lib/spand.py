"""SPAND — Spatial Pattern of Aggregated Neighborhood Diversity.

A per-slide score quantifying how spatially structured a gene-expression-derived
signal is over a spatial-transcriptomics tissue grid. The score is **Global
Moran's I divided by the mean of the signal**:

    SPAND = Moran(signal, lat2W(grid_shape)).I / mean(signal)

The signal can be any per-spot scalar — predicted expression of a single gene
(e.g. ERBB2), a GSEA pathway NES per spot, or a cancer-restricted pathway score
(per-spot NES divided by per-spot cancer-cell fraction, as used in the paper's
HER2 SPAND).

Inputs to `spand_for_slide` are kept minimal — a per-spot signal and the spot
(x, y) coordinates. The spot positions are rasterized into a 2-D grid; the
spatial weights matrix is the queen / rook adjacency on that grid (from
`libpysal.weights.lat2W`); NaN spots are dropped from the Moran calculation.
"""

import numpy as np
import pandas as pd
from esda import Moran
from libpysal import weights

__all__ = ["spand_for_slide", "morans_i_on_grid"]


def _smooth_grid(arr):
    """Fill isolated NaNs by the mean of the four orthogonal neighbours.

    Spots that are NaN AND have at least two NaN neighbours stay NaN, so genuine
    holes in the tissue grid are preserved. This mirrors the smoothing applied
    in the original SPAND pipeline (`SPAND.py`).
    """
    out = arr.astype(float, copy=True)
    R, C = arr.shape
    for i in range(R):
        for j in range(C):
            if np.isnan(arr[i, j]):
                vals = []
                if i > 0:     vals.append(arr[i - 1, j])
                if i < R - 1: vals.append(arr[i + 1, j])
                if j > 0:     vals.append(arr[i, j - 1])
                if j < C - 1: vals.append(arr[i, j + 1])
                vals = np.asarray(vals, dtype=float)
                if np.isnan(vals).sum() > 1:
                    out[i, j] = np.nan
                else:
                    out[i, j] = np.nanmean(vals)
    return out


def morans_i_on_grid(signal, xy):
    """Global Moran's I of a per-spot signal laid out on a 2-D spot grid.

    signal : Series (one value per spot) — values to test for spatial autocorrelation.
    xy     : DataFrame indexed identically to `signal`, with integer-or-castable
             columns ``x`` and ``y`` giving the spot grid coordinates.

    Spots are rasterized into a `(max_x + 1, max_y + 1)` grid; missing positions
    are NaN. The weights matrix is grid adjacency (`libpysal.weights.lat2W`),
    restricted to the non-NaN spots.

    Returns Moran's I as a float, or NaN if the grid has no valid spots.
    """
    df = pd.concat([signal.rename("v"), xy[["x", "y"]]], axis=1).copy()
    df["x"] = df["x"].astype(int)
    df["y"] = df["y"].astype(int)

    grid = df.pivot(index="x", columns="y", values="v").to_numpy(dtype=float)
    grid = _smooth_grid(grid)

    flat = grid.flatten()
    keep = ~np.isnan(flat)
    if keep.sum() == 0:
        return np.nan

    w_full = weights.lat2W(grid.shape[0], grid.shape[1], silence_warnings=True)
    valid_ids = np.where(keep)[0]
    valid_neighs = {i: [j for j in w_full.neighbors[i] if keep[j]] for i in valid_ids}
    w = weights.W(valid_neighs, silence_warnings=True)

    return Moran(flat[keep], w).I


def spand_for_slide(signal, xy):
    """SPAND for one slide: Moran's I divided by the mean of the signal.

    signal : per-spot Series (e.g., predicted ERBB2 expression, or a
             cancer-restricted pathway NES).
    xy     : DataFrame with ``x`` and ``y`` columns, aligned to `signal`'s index.

    Returns SPAND as a float. Negative SPAND values are interpreted as more
    spatially heterogeneous (per the paper's sign convention applied downstream).
    """
    mean_signal = float(signal.mean())
    if mean_signal == 0 or np.isnan(mean_signal):
        return np.nan
    return morans_i_on_grid(signal, xy) / mean_signal
