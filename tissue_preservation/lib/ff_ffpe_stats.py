"""Tissue-preservation analysis for Path2Space (Figure 2G).

Tests whether Path2Space — trained only on fresh-frozen (FF) slides — generalizes
to archival FFPE tissue, and asks which slide-level features explain the
remaining preservation effect on prediction accuracy.

Two analyses:

  - `ff_vs_ffpe_test` — linear mixed-effects model of per-sample mean
    predicted-vs-measured PCC against preservation method, with cohort as a
    random intercept. Tests whether FF and FFPE differ in mean accuracy after
    accounting for between-cohort variation.

  - `morphology_associations` — for each of the eight CellPose-derived nuclear
    features (nuclear area, equivalent diameter, circularity, solidity,
    eccentricity, aspect ratio, nearest-neighbor distance, local nuclear
    density), fits a per-sample MixedLM `mean_pcc ~ feature + preservation,
    groups=cohort` and applies Benjamini-Hochberg FDR across the eight
    features.
"""

import warnings

import numpy as np
import pandas as pd
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.stats.multitest import multipletests

__all__ = ["MORPHOLOGY_FEATURES", "ff_vs_ffpe_test", "morphology_associations"]


# Eight CellPose-derived nuclear features used for the morphology analysis.
MORPHOLOGY_FEATURES = [
    "nuclear_area",
    "equivalent_diameter",
    "circularity",
    "solidity",
    "eccentricity",
    "aspect_ratio",
    "nearest_neighbor_distance",
    "nuclear_density_local",
]


def _fit_mixedlm(formula, df, group_col):
    """Fit a MixedLM with cohort as random intercept.

    Tries lbfgs first, then powell as fallback — some narrow predictors (e.g.
    eccentricity over a tight per-sample range) settle on a non-positive-
    definite Hessian under powell, which returns NaN p-values; lbfgs converges
    to a usable fit on the same data. Convergence warnings are silenced.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for method in ("lbfgs", "powell"):
            try:
                fit = MixedLM.from_formula(formula, groups=df[group_col], data=df).fit(
                    reml=True, method=method
                )
                fixed = [c for c in fit.pvalues.index if c not in ("Intercept", "Group Var")]
                if not np.isnan(fit.pvalues[fixed]).any():
                    return fit
            except Exception:
                continue
        # Fall through: return the last fit even if some p-values are NaN
        return fit


def ff_vs_ffpe_test(summary, response="mean_pcc",
                    preservation_col="preservation_method",
                    cohort_col="cohort"):
    """Cohort-adjusted FF vs FFPE test on per-sample prediction accuracy.

    summary : DataFrame with one row per sample; must have `response`,
        `preservation_col` (values 'FF' / 'FFPE') and `cohort_col`.

    Returns a dict with the fitted FF-coefficient, its p-value, n samples used,
    and the median `response` for FF and FFPE.
    """
    d = summary[[response, preservation_col, cohort_col]].dropna().copy()
    d = d[d[preservation_col].isin(["FF", "FFPE"])]
    d["ff"] = (d[preservation_col] == "FF").astype(int)

    fit = _fit_mixedlm(f"{response} ~ ff", d, cohort_col)
    return {
        "n": len(d),
        "median_FF":   float(d.loc[d.ff == 1, response].median()),
        "median_FFPE": float(d.loc[d.ff == 0, response].median()),
        "beta_FF": float(fit.params["ff"]),
        "pvalue":  float(fit.pvalues["ff"]),
    }


def morphology_associations(summary, features=MORPHOLOGY_FEATURES,
                            response="mean_pcc",
                            preservation_col="preservation_method",
                            cohort_col="cohort"):
    """Per-feature MixedLM of accuracy on each nuclear morphology feature.

    For each `feature` in `features`, fits
    `response ~ feature + preservation, groups=cohort`
    on the per-sample table and applies Benjamini-Hochberg FDR across the set.

    Returns a DataFrame with columns: feature, beta, pvalue, fdr, n.
    """
    d0 = summary[[response, preservation_col, cohort_col] + list(features)].copy()
    d0 = d0[d0[preservation_col].isin(["FF", "FFPE"])]
    d0["ff"] = (d0[preservation_col] == "FF").astype(int)

    rows = []
    for f in features:
        d = d0[[response, "ff", cohort_col, f]].dropna()
        try:
            fit = _fit_mixedlm(f"{response} ~ {f} + ff", d, cohort_col)
            rows.append({"feature": f,
                         "beta":   float(fit.params[f]),
                         "pvalue": float(fit.pvalues[f]),
                         "n":      len(d)})
        except Exception:
            rows.append({"feature": f, "beta": np.nan, "pvalue": np.nan, "n": len(d)})

    out = pd.DataFrame(rows)
    pvals = out["pvalue"].to_numpy()
    fdr = np.full(len(out), np.nan)
    mask = ~np.isnan(pvals)
    if mask.any():
        fdr[mask] = multipletests(pvals[mask], method="fdr_bh")[1]
    out["fdr"] = fdr
    return out.sort_values("pvalue").reset_index(drop=True)
