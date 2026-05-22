"""Confounder analysis for Path2Space predictions (Table S2).

Tests whether Path2Space predictions reflect genuine biological variation
rather than technical or compositional artifacts, by checking how much of the
agreement between predicted and measured spatial expression survives once a
potential confounder is controlled for. Three confounders are tested:

  - hematoxylin stain intensity   — directly observable in the input H&E image
  - cancer cell fraction          — SpaCET deconvolution of the measured ST
  - total RNA content per spot    — log10(1 + total counts)

The analysis runs in two stages:

  Stage A — `per_slide_gene_correlations`
    Per slide and per gene: Pearson correlations of predicted and measured
    expression with the confounder, and the partial correlation between
    predicted and measured expression controlling for the confounder
    (`partial_correlation`, residual-based regression).

  Stage B — `aggregate_per_gene` then `confounder_summary`
    Aggregate the per-slide statistics to one row per gene: spot-count-weighted
    mean correlations, p-values combined across slides with Stouffer's method
    (weights proportional to sqrt(n_spots)), and Benjamini-Hochberg FDR across
    genes. `confounder_summary` collapses that to the per-confounder table.

Stage A needs the full per-slide ST dataset (predictions, deconvolution, QC);
its output — the per-slide / per-gene correlation table — is what the
accompanying notebook consumes. Stage B is cheap and fully reproducible from
that table.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

__all__ = ["partial_correlation", "per_slide_gene_correlations",
           "combine_pvalues_stouffer", "aggregate_per_gene",
           "confounder_summary", "METRIC_PAIRS", "SUMMARY_CONFOUNDERS"]


# (statistic column, its p-value column) for every metric in the per-slide
# table. Note partial-correlation columns are `partial_corr_given_<c>` while
# their p-values are `partial_corr_<c>_pvalue` — hence the explicit pairing.
METRIC_PAIRS = [
    ("corr_pred_vs_obs",          "corr_pred_vs_obs_pvalue"),
    ("corr_pred_vs_tumor",        "corr_pred_vs_tumor_pvalue"),
    ("corr_obs_vs_tumor",         "corr_obs_vs_tumor_pvalue"),
    ("partial_corr_given_tumor",  "partial_corr_tumor_pvalue"),
    ("corr_pred_vs_H",            "corr_pred_vs_H_pvalue"),
    ("corr_obs_vs_H",             "corr_obs_vs_H_pvalue"),
    ("partial_corr_given_H",      "partial_corr_H_pvalue"),
    ("corr_pred_vs_counts",       "corr_pred_vs_counts_pvalue"),
    ("corr_obs_vs_counts",        "corr_obs_vs_counts_pvalue"),
    ("partial_corr_given_counts", "partial_corr_counts_pvalue"),
]

# (display name, column suffix) for the three confounders, in the row order
# used by `confounder_summary`.
SUMMARY_CONFOUNDERS = [
    ("Hematoxylin intensity", "H"),
    ("Total RNA content",     "counts"),
    ("Cancer cell fraction",  "tumor"),
]


# ----------------------------------------------------------- Stage A (slide) ---

def partial_correlation(x, y, z):
    """Partial correlation between x and y controlling for z, with p-value.

    Standardizes all three variables, regresses x and y on z separately, and
    Pearson-correlates the residuals. The p-value is two-sided, from a
    t-distribution with n - 3 degrees of freedom.

    Returns (r_partial, p_value), or (nan, nan) if fewer than 4 valid
    (non-NaN) samples remain or any variable is constant.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    valid = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    if valid.sum() < 4:
        return np.nan, np.nan
    x, y, z = x[valid], y[valid], z[valid]
    n = len(x)

    if np.std(x) == 0 or np.std(y) == 0 or np.std(z) == 0:
        return np.nan, np.nan

    x = (x - x.mean()) / x.std()
    y = (y - y.mean()) / y.std()
    z = (z - z.mean()) / z.std()

    # residuals of x and y after regressing out z
    res_x = x - (np.sum(x * z) / np.sum(z * z)) * z
    res_y = y - (np.sum(y * z) / np.sum(z * z)) * z

    denom = np.sqrt(np.sum(res_x ** 2) * np.sum(res_y ** 2))
    if denom == 0 or np.isnan(denom):
        return np.nan, np.nan
    r = np.sum(res_x * res_y) / denom

    if abs(r) >= 1.0:
        return r, (0.0 if abs(r) > 1.0 else np.nan)
    t = r * np.sqrt((n - 3) / (1 - r ** 2))
    p = 2 * (1 - stats.t.cdf(abs(t), df=n - 3))
    return r, p


def per_slide_gene_correlations(pred, obs, confounder, genes):
    """Per-gene correlations within one slide, for one confounder.

    pred, obs   : (n_spots, n_genes) arrays of predicted / measured expression.
    confounder  : (n_spots,) array of the confounder value per spot.
    genes       : sequence of gene names, length n_genes.

    For each gene returns Pearson r (and p-value) of predicted vs measured
    expression, of predicted and measured expression vs the confounder, and the
    partial correlation of predicted vs measured controlling for the confounder.
    Genes with fewer than 3 valid spots or a constant variable yield NaN.
    """
    confounder = np.asarray(confounder, dtype=float)
    rows = []
    for i, gene in enumerate(genes):
        p = np.asarray(pred[:, i], dtype=float)
        o = np.asarray(obs[:, i], dtype=float)
        valid = ~(np.isnan(p) | np.isnan(o) | np.isnan(confounder))

        row = {"gene": gene}
        if valid.sum() < 3:
            rows.append(row)
            continue
        pv, ov, cv = p[valid], o[valid], confounder[valid]
        if np.std(pv) == 0 or np.std(ov) == 0 or np.std(cv) == 0:
            rows.append(row)
            continue

        row["corr_pred_vs_obs"], row["corr_pred_vs_obs_pvalue"] = stats.pearsonr(pv, ov)
        row["corr_pred_vs_confounder"], row["corr_pred_vs_confounder_pvalue"] = stats.pearsonr(pv, cv)
        row["corr_obs_vs_confounder"], row["corr_obs_vs_confounder_pvalue"] = stats.pearsonr(ov, cv)
        row["partial_corr"], row["partial_corr_pvalue"] = partial_correlation(pv, ov, cv)
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------- Stage B (aggregate) ---

def combine_pvalues_stouffer(pvalues, n_spots):
    """Combine per-slide p-values with Stouffer's method.

    Weights are proportional to sqrt(n_spots): the standard error of a
    correlation scales with 1/sqrt(n). NaN p-values are dropped; a one-element
    input is returned unchanged; an empty input returns nan.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    n_spots = np.asarray(n_spots, dtype=float)

    keep = ~np.isnan(pvalues)
    pvalues, n_spots = pvalues[keep], n_spots[keep]
    if len(pvalues) == 0:
        return np.nan
    if len(pvalues) == 1:
        return float(pvalues[0])

    pvalues = np.clip(pvalues, 1e-300, 1 - 1e-10)   # avoid log(0) / log(1)
    return stats.combine_pvalues(pvalues, method="stouffer",
                                 weights=np.sqrt(n_spots))[1]


def aggregate_per_gene(per_slide, pairs=METRIC_PAIRS):
    """Aggregate the per-slide / per-gene table to one row per gene.

    per_slide : DataFrame with a `gene` column, an `n_spots` column, and the
        statistic / p-value columns listed in `pairs`.

    Correlation columns become spot-count-weighted means across slides;
    p-value columns are combined with `combine_pvalues_stouffer`; and a
    Benjamini-Hochberg FDR column (`*_fdr`) is added for each p-value column,
    corrected across genes.
    """
    records = []
    for gene, g in per_slide.groupby("gene", observed=True):
        n_spots = g["n_spots"].to_numpy()
        row = {"gene": gene, "n_slides": len(g), "total_spots": int(n_spots.sum())}
        for val_col, pval_col in pairs:
            vals = g[val_col].to_numpy(dtype=float)
            mask = ~np.isnan(vals)
            row[val_col] = (np.average(vals[mask], weights=n_spots[mask])
                            if mask.any() else np.nan)
            row[pval_col] = combine_pvalues_stouffer(g[pval_col].to_numpy(), n_spots)
        records.append(row)
    genes = pd.DataFrame.from_records(records)

    for _, pval_col in pairs:
        pvals = genes[pval_col].to_numpy(dtype=float)
        fdr = np.full(len(genes), np.nan)
        mask = ~np.isnan(pvals)
        if mask.any():
            fdr[mask] = multipletests(pvals[mask], method="fdr_bh")[1]
        genes[pval_col.replace("pvalue", "fdr")] = fdr
    return genes


def confounder_summary(gene_all, gene_cv, gene_ext):
    """Build the Table S2 confounder summary — one row per confounder.

    gene_all / gene_cv / gene_ext : per-gene tables from `aggregate_per_gene`
    run on all slides / cross-validation slides / external-validation slides.

    A gene "retains" a confounder if its partial correlation is positive and
    significant at FDR < 0.05 after controlling for that confounder.
    """
    rows = []
    n_genes = len(gene_all)
    for name, suffix in SUMMARY_CONFOUNDERS:
        partial = f"partial_corr_given_{suffix}"
        fdr = f"partial_corr_{suffix}_fdr"
        sig_positive = int(((gene_all[fdr] < 0.05) & (gene_all[partial] > 0)).sum())
        rows.append({
            "Confounder": name,
            "Correlation with measured expression":  round(gene_all[f"corr_obs_vs_{suffix}"].mean(), 2),
            "Correlation with predicted expression": round(gene_all[f"corr_pred_vs_{suffix}"].mean(), 2),
            "Genes with significant positive partial correlation (n)": sig_positive,
            "Genes with significant positive partial correlation (%)": round(100 * sig_positive / n_genes, 1),
            "Median partial correlation (all cohorts)":         round(gene_all[partial].median(), 2),
            "Median partial correlation (cross-validation)":    round(gene_cv[partial].median(), 2),
            "Median partial correlation (external validation)": round(gene_ext[partial].median(), 2),
        })
    return pd.DataFrame(rows)
