"""Per-sample QC metrics for spatial-transcriptomics sections.

The six quality-control features of Shulman et al., *Cell* 2026:

Histological image quality (per H&E tile, averaged over a section's tiles)
  1. hematoxylin — mean hematoxylin optical density within tissue
  2. eosin       — mean eosin optical density within tissue
  3. sharpness   — variance of the Laplacian of the grayscale tile

Gene-expression quality (per spot, averaged over a section's spots)
  4. total_counts          — mean total UMI counts per spot
  5. log1p_total_counts_mt — mean log1p mitochondrial counts per spot
  6. log1p_total_counts_hb — mean log1p haemoglobin counts per spot

Image metrics use scikit-image; expression metrics use
`scanpy.pp.calculate_qc_metrics`. Image QC should be run on color-normalized
tiles (see `stain_norm.macenko_normalize`).
"""

import numpy as np
import pandas as pd
import scanpy as sc
from skimage import color, filters

# ---------------------------------------------------------------- image QC ---

def compute_image_qc(tile_rgb):
    """Three histological-image QC metrics for one H&E tile.

    tile_rgb : (H, W, 3) uint8 array — a (color-normalized) H&E tile.
    Returns {'hematoxylin', 'eosin', 'sharpness'}.
    """
    imgf = np.asarray(tile_rgb, dtype=np.float32) / 255.0
    gray = color.rgb2gray(imgf)

    # tissue = pixels darker than the Otsu split (background is bright in H&E)
    tissue = gray < filters.threshold_otsu(gray)

    hed = color.rgb2hed(imgf)                       # H / E / DAB optical density
    H, E = hed[..., 0], hed[..., 1]
    return {
        "hematoxylin": float(H[tissue].mean() if tissue.any() else H.mean()),
        "eosin":       float(E[tissue].mean() if tissue.any() else E.mean()),
        "sharpness":   float(filters.laplace(gray).var()),
    }


def image_qc_per_section(tiles):
    """Mean image QC over a section's tiles.

    tiles : iterable of (H, W, 3) uint8 tiles.
    Returns a per-section {'hematoxylin', 'eosin', 'sharpness'} dict.
    """
    per_tile = pd.DataFrame(compute_image_qc(t) for t in tiles)
    return per_tile.mean().to_dict()


# ----------------------------------------------------------- expression QC ---

_HB_GENES = {"HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2",
             "HBE1", "HBM", "HBZ", "HBQ1"}


def _flag_qc_genes(adata):
    """Tag mitochondrial (`MT-*`) and haemoglobin genes in `adata.var`."""
    sym = pd.Series(adata.var_names.astype(str), index=adata.var_names).str.upper()
    adata.var["mt"] = sym.str.startswith("MT-").values
    adata.var["hb"] = (sym.isin(_HB_GENES) | sym.str.match(r"^HB[ABEDGMQZ]\d?$")).values


def compute_expression_qc(adata):
    """Three gene-expression QC metrics for one ST section.

    adata : AnnData with raw counts in `.X`.
    Returns per-section {'total_counts', 'log1p_total_counts_mt',
    'log1p_total_counts_hb'} — means over spots, via
    `scanpy.pp.calculate_qc_metrics`.
    """
    _flag_qc_genes(adata)
    obs, _ = sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "hb"], percent_top=None, inplace=False)
    return {
        "total_counts":          float(obs["total_counts"].mean()),
        "log1p_total_counts_mt": float(obs["log1p_total_counts_mt"].mean()),
        "log1p_total_counts_hb": float(obs["log1p_total_counts_hb"].mean()),
    }


# ------------------------------------------------------------- thresholds ---
# Paper QC thresholds (Shulman et al., Cell 2026, Methods). A section is
# excluded if a metric falls on the failing side of its threshold.
THRESHOLDS = {
    "hematoxylin":           {"max": 0.053546, "min": 0.029768},  # TCGA mean +/- 2 SD
    "eosin":                 {"max": 0.022495, "min": 0.017486},  # TCGA mean +/- 2 SD
    "sharpness":             {"min": 0.000637},                   # TCGA reference
    "total_counts":          {"min": 6193.7},                     # Bassiouni 1st pctile
    "log1p_total_counts_mt": {"max": 8.3603},                     # Bassiouni 99th pctile
    "log1p_total_counts_hb": {"max": 0.7355},                     # Bassiouni 99th pctile
}


def passes_qc(metrics):
    """Per-metric pass/fail given a {metric: value} dict, using `THRESHOLDS`."""
    out = {}
    for name, value in metrics.items():
        t = THRESHOLDS.get(name, {})
        ok = True
        if "min" in t and value < t["min"]:
            ok = False
        if "max" in t and value > t["max"]:
            ok = False
        out[name] = ok
    return out
