# confounder_analysis — Table S2

Reproduces the confounder analysis (Table S2) of Shulman et al., *Cell* 2026.

Path2Space predicts spatial gene expression from H&E histology. This component
tests whether those predictions reflect **genuine biological variation** rather
than technical or compositional artifacts: for every gene it measures how much
of the agreement between predicted and measured expression survives once a
potential confounder is controlled for.

## The three confounders

| Confounder | Why it could confound | How it is measured |
|---|---|---|
| **Hematoxylin stain intensity** | Directly visible in the input H&E image — a model could key on staining rather than biology | mean hematoxylin optical density per spot (see the `qc` component) |
| **Cancer cell fraction** | Tumor-dense regions have distinct expression; a model could merely track tumor density | SpaCET deconvolution of the measured spatial transcriptomics |
| **Total RNA content per spot** | Spots with more RNA show higher counts for most genes | `log1p` of total UMI counts per spot |

A gene "retains" a confounder if, after controlling for it, the partial
correlation between predicted and measured expression is still positive and
significant (FDR < 0.05).

## Contents

- `lib/confounder_stats.py` — the analysis functions: `partial_correlation`,
  `per_slide_gene_correlations` (Stage A), `combine_pvalues_stouffer`,
  `aggregate_per_gene`, `confounder_summary` (Stage B).
- `notebooks/01_confounder_analysis.ipynb` — worked example: loads the bundled
  per-slide table, aggregates per gene, and writes the two output tables.
- `data/per_slide_correlations.parquet` — **bundled input** (see below).
- `data/confounder_summary.pkl` — output table 1 (one row per confounder).
- `data/gene_level_confounder_stats.pkl` — output table 2 (one row per gene).

## Method

The analysis runs in two stages.

**Stage A — per slide, per gene.** For each slide and gene: Pearson
correlations of predicted and measured expression with the confounder, and the
partial correlation between predicted and measured expression controlling for
the confounder. Partial correlations use residual-based regression — standardize
the three variables, regress predicted and measured expression on the
confounder, and correlate the residuals.

**Stage B — aggregate across slides.** Per gene: spot-count-weighted mean
correlations; per-slide p-values combined with Stouffer's method (weights
proportional to `sqrt(n_spots)`); Benjamini-Hochberg FDR across genes. Slides
are aggregated three ways — all cohorts, cross-validation only (held-out folds
of the TNBC training cohort), and external validation only (HEST, HTAN,
pierre_martinez).

## Bundled data

Stage A needs the full per-slide spatial-transcriptomics dataset — predicted
and measured expression matrices, SpaCET deconvolution, and per-spot QC — which
is too large to distribute. Its **output** is bundled instead:

`data/per_slide_correlations.parquet` — the per-slide / per-gene correlation
table, 40 slides × 14,068 genes (548,340 rows). Correlations are stored as
`float32` and p-values as `float64`; with zstd compression the file is ~63 MB.
The notebook reproduces both output tables from this file in seconds.

The auxiliary ΔR² columns from the original analysis are not included — the
Table S2 confounder summary is built entirely from partial correlations.

## Running

```bash
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/01_confounder_analysis.ipynb
```

Both output tables are written as pickles; the notebook does not produce Excel.
