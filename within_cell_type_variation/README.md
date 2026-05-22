# within_cell_type_variation — Figures 4E, 4F, S5G, S5H

Reproduces the within-cell-type expression-variation analysis of Shulman et al.,
*Cell* 2026.

Path2Space predicts spatial gene expression from H&E histology. Beyond
predicting *which* cell types are present, this component asks whether the
predictions resolve expression variation **within** individual cell types, and
whether that accuracy is driven by canonical cell-type marker genes.

## The four figures

| Figure | Question | Notebook |
|---|---|---|
| **4E** | Within pathologist-annotated cancer / stromal / lymphocyte regions, how well do predicted and measured expression agree, gene by gene? | `01_within_cell_type_correlations.ipynb` |
| **4F** | The same, inside transcriptomically-defined neighborhoods (spatially coherent regions grouped by dominant cell type). | `01_within_cell_type_correlations.ipynb` |
| **S5G** | Is gene-wise predictive accuracy concordant across cell types — are genes well predicted in one cell type also well predicted in another? | `01_within_cell_type_correlations.ipynb` |
| **S5H** | Is accuracy driven by canonical cell-type **marker** genes, or does it extend to non-marker genes? | `02_marker_vs_nonmarker.ipynb` |

## Contents

- `lib/cell_type_corr.py` — `gene_pearson` (Stage A primitive), `summarize_by_cell_type`, `concordance_matrix`.
- `notebooks/01_within_cell_type_correlations.ipynb` — Figures 4E, 4F, S5G.
- `notebooks/02_marker_vs_nonmarker.ipynb` — Figure S5H.
- `data/pathologist_region_correlations.parquet` — bundled input for 4E.
- `data/neighborhood_correlations.parquet` — bundled input for 4F.
- `data/marker_gene_correlations.parquet` — bundled input for S5H.
- `data/pathologist_summary.pkl`, `data/neighborhood_summary.pkl` — outputs written by notebook 1.

## Method

Two stages.

**Stage A — per-gene correlation within a region.** For each region (a
pathologist-annotated area, a transcriptomic neighborhood, or — for the marker
analysis — a whole slide), and for every gene, the Pearson correlation between
predicted and measured expression is computed across the region's spots
(`gene_pearson`). Per-gene correlations are then aggregated to one value per
(cell type, gene) for 4E/4F, or per gene for S5H.

**Stage B — aggregation and figures** (the notebooks). Per-cell-type median
correlation and the number of genes exceeding PCC > 0.4
(`summarize_by_cell_type`); the pairwise concordance of gene-wise accuracy
across cell types (`concordance_matrix`, S5G); and, for S5H, a two-sided
Mann-Whitney U test of marker vs non-marker gene accuracy.

## Bundled data

Stage A needs the full per-slide ST dataset — predicted and measured expression
matrices, pathologist annotations, deconvolution-based neighborhood labels —
which is too large to distribute. Its **output** is bundled instead, as three
small parquet tables in `data/`:

| File | Rows | Columns | Stage-A source |
|---|---|---|---|
| `pathologist_region_correlations.parquet` | per (cell type, gene) | `cell_type`, `gene`, `correlation` | per-gene Pearson within pathologist-annotated regions on the HEST validation slides; regions outlined from a pathologist annotation table |
| `neighborhood_correlations.parquet` | per (cell type, gene) | `cell_type`, `gene`, `correlation` | per-gene Pearson within transcriptomic neighborhoods — spatially coherent regions grouped by their dominant deconvolved cell type |
| `marker_gene_correlations.parquet` | per gene | `gene`, `correlation`, `gene_group` | per-gene Pearson between predicted and measured whole-tissue expression, averaged across cohorts; `gene_group` labels each gene `Marker genes` or `Non-specific` from an scRNA-seq differential-expression table |

`cell_type` for 4E is `All` (whole tissue), `Cancer`, `Stromal`, `Lympho.`; for
4F it is the dominant deconvolved cell type of each neighborhood.

## Reproduction notes

- **Figure 4F / S5G neighborhood input.** The private analysis notebook
  (`plot_cell_types_new.ipynb`) loads an older SpaGCN-based neighborhood table
  (`plot_df.pkl`, four cell types — Cancer, Macrophage, B cell, CAF) as its
  neighborhood input, whereas the correct Figure 4F data is the five-cell-type
  table (`plot_df_neighborhoods.pkl`, which adds the T cell neighborhood). This
  component uses the correct file — `neighborhood_correlations.parquet` is
  derived from it.
- **S5G concordance range.** On the corrected five-cell-type table, pairwise
  gene-accuracy concordance spans 0.67–0.95 for pathologist regions and
  0.23–0.80 for neighborhoods. The neighborhood low end is the T cell column —
  gene-wise accuracy in T-cell neighborhoods is only weakly concordant with that
  in B-cell or macrophage neighborhoods. The older four-cell-type table, which
  lacks the T cell column, spans a tighter 0.52–0.85.

## Running

```bash
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/01_within_cell_type_correlations.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/02_marker_vs_nonmarker.ipynb
```

Summary tables are written as pickles; the notebooks do not produce Excel.

## Citation

Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023
