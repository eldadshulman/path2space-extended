"""Macenko H&E stain normalization (numpy-only).

Reproduces the color-normalization step applied before histological-image QC in
Shulman et al., *Cell* 2026. Implements the method of Macenko et al. (2009);
concentrations are estimated by least squares (no `spams` dependency). The
reference stain matrix and target concentrations match the normalizer used in
the paper pipeline.
"""

import numpy as np

# Reference H&E stain matrix (rows = R/G/B optical density, cols = H, E) and
# target 99th-percentile concentrations — the paper-pipeline defaults.
HE_REF = np.array([[0.5626, 0.2159],
                   [0.7201, 0.8012],
                   [0.4062, 0.5581]], dtype=np.float64)
MAXC_REF = np.array([1.9705, 1.0308], dtype=np.float64)


def _standardize_brightness(I):
    """Scale so the 90th-percentile pixel maps to 255 (illumination correction)."""
    p = np.percentile(I, 90)
    return np.clip(I * 255.0 / max(p, 1e-6), 0, 255)


def macenko_normalize(tile_rgb, Io=240, alpha=1, beta=0.15):
    """Macenko-normalize one H&E tile.

    tile_rgb : (H, W, 3) uint8 array.
    Returns a (H, W, 3) uint8 array with stain colour matched to the reference.
    """
    I = _standardize_brightness(np.asarray(tile_rgb, dtype=np.float64))
    h, w, _ = I.shape
    flat = I.reshape(-1, 3)

    # optical density
    OD = -np.log((flat + 1.0) / Io)
    OD_hat = OD[~np.any(OD < beta, axis=1)]
    if OD_hat.shape[0] < 10:                       # too little tissue — return as-is
        return np.clip(I, 0, 255).astype(np.uint8)

    # stain vectors = plane of the two leading eigenvectors of the OD covariance
    _, V = np.linalg.eigh(np.cov(OD_hat.T))
    V = V[:, [2, 1]]
    proj = OD_hat @ V
    phi = np.arctan2(proj[:, 1], proj[:, 0])
    v_min = V @ np.array([np.cos(np.percentile(phi, alpha)),
                          np.sin(np.percentile(phi, alpha))])
    v_max = V @ np.array([np.cos(np.percentile(phi, 100 - alpha)),
                          np.sin(np.percentile(phi, 100 - alpha))])
    HE = np.array([v_min, v_max]).T if v_min[0] > v_max[0] else np.array([v_max, v_min]).T

    # per-pixel stain concentrations (least squares), rescaled to the reference
    C = np.linalg.lstsq(HE, OD.T, rcond=None)[0]
    maxC = np.percentile(C, 99, axis=1)
    C = C * (MAXC_REF / np.where(maxC == 0, 1e-6, maxC))[:, None]

    norm = Io * np.exp(-HE_REF @ C)
    return np.clip(norm.T, 0, 255).reshape(h, w, 3).astype(np.uint8)
