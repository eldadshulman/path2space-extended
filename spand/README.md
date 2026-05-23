# spand — SPAND score and the HER2-low / HER2-high PBCP analysis

Reproduces the SPAND (Spatial Pattern of Aggregated Neighborhood Diversity)
analysis of Shulman et al., *Cell* 2026.

## What SPAND is

SPAND is a per-slide score quantifying how spatially structured a gene-
expression-derived signal is over a spatial-transcriptomics tissue grid. The
score is **Global Moran's I divided by the mean of the signal**:

> SPAND(signal) = Moran's I(signal, lat2W(grid)) / mean(signal)

where the weights matrix is the queen/rook adjacency on the spot grid (from
`libpysal.weights.lat2W`). The signal can be any per-spot scalar — predicted
expression of a single gene, a GSEA pathway NES per spot, or a cancer-cell-
restricted pathway score. The paper convention negates the value so positive
scores mean **more** spatially heterogeneous.

The published HER2-SPAND used downstream is **not** raw single-gene Moran's I —
it is Moran's I on a cancer-cell-restricted, GSEA-pathway-derived signature:

```
predicted expression → cancer-cell deconvolution → per-spot pathway NES (GSEA)
                                       │
                                       └─→ NES / cancer_fraction   ─→ SPAND
```

## Contents

- `lib/spand.py` — `spand_for_slide(signal, xy)`: rasterizes the signal onto the
  spot grid and returns Global Moran's I / mean. `morans_i_on_grid(signal, xy)`
  is the Moran-only intermediate. Verbatim port of the canonical `SPAND.py`.
- `notebooks/01_spand_pipeline_example.ipynb` — end-to-end SPAND pipeline on
  one PBCP slide: predicted expression → 5-fold-ensemble cancer-cell
  deconvolution (using `cell_type_deconvolution/models/`) → load per-spot
  GSEA NES (precomputed upstream) → cancer-normalize → SPAND. Reproduces the
  bundled per-patient HER2-SPAND for that slide to ~0.4%.
- `notebooks/02_her2low_vs_her2high_PBCP.ipynb` — PBCP clinical analysis using
  precomputed per-patient HER2-SPAND vs measured bulk ERBB2. ROC AUC with
  bootstrap 95% CIs, split by HER2-low / HER2-high.
- `data/example_slide_predicted_expression.parquet` — one slide (PBCP 732799),
  Path2Space predicted expression restricted to the ~1,400 genes used by the
  deconvolution MLP feature list (~8 MB).
- `data/example_slide_gsea_her2_nes.parquet` — per-(spot, pathway) GSEA NES
  for the six ErbB-related pathways on that slide (~90 KB).
- `data/pbcp_her2_scores.parquet` — per-patient HER2-SPAND, measured ERBB2
  (log2 TPM), and treatment response for 18 PBCP patients (~3 KB).
- `data/erbb_pathway_genes.json` — Wikipathways WP673 ErbB Signaling Pathway
  gene list with provenance; bundled for documentation of the gene-set used
  in the upstream GSEA (not consumed by the notebooks).

## Method

Two stages, mirroring the canonical pipeline:

**Stage A — per-slide HER2-SPAND** (notebook 01 demonstrates on one slide):

1. **Predicted expression.** Path2Space outputs per-spot predicted expression
   for ~14,000 genes. This component bundles a single example slide trimmed to
   the genes consumed by the deconvolution model.
2. **Cancer-cell deconvolution.** Apply the five MLP fold models in
   [`cell_type_deconvolution/models/`](../cell_type_deconvolution/models/) to
   the predicted expression and average across folds. PBCP is not in the
   PanopTILs training set, so no fold is "held out"; averaging is the
   principled default for an out-of-distribution slide.
3. **Per-spot GSEA** (upstream — *not* re-executed in the notebook). For each
   spot, rank the ~14,000 predicted genes by expression and run
   `gseapy.prerank` against the Wikipathways ErbB Signaling Pathway (WP673)
   gene set; collect the Normalized Enrichment Score per spot per pathway.
   The bundled `example_slide_gsea_her2_nes.parquet` is the result of this
   step for the example slide; the recipe is below.
4. **Cancer-restrict.** For each spot, divide the WP673 NES by the
   Epithelial (cancer) fraction from step 2.
5. **SPAND.** Pass the cancer-restricted signal and the spot (x, y) grid
   coordinates to `spand_for_slide`. Negate by convention.

**Stage B — clinical analysis on precomputed per-patient scores**
(notebook 02): the per-patient HER2-SPAND across PBCP is bundled as
`pbcp_her2_scores.parquet`. Computing it at scale across PBCP / IMPRESS /
TransNEO requires the full Path2Space prediction pipeline (companion repo)
and is not reproduced inline here.

## Upstream GSEA recipe

For a reader who wants to recompute the per-spot NES from raw predicted
expression:

```python
import gseapy
# Per-spot rank → GSEA-prerank against Wikipathways WP673
for spot in predicted_expression.index:
    ranked = predicted_expression.loc[spot].sort_values(ascending=False)
    res = gseapy.prerank(rnk=ranked, gene_sets={"WP673": wp673_genes},
                        outdir=None, no_plot=True)
    nes[spot] = res.res2d["NES"].iloc[0]
```

The WP673 gene list is in `data/erbb_pathway_genes.json` (Wikipathways,
CC-BY 3.0; source URL in the JSON).

## Reproduction notes

- **Bundled scope.** Only one example slide of predicted expression ships
  (~8 MB after restriction to deconvolution feature genes). No cohort-level
  predicted-expression data is bundled; the per-patient PBCP table is
  precomputed scores only. The full PBCP predicted-expression matrices live
  upstream in the companion repo and are not redistributed here.
- **Example slide.** PBCP `732799` was chosen as the smallest slide with a
  non-null reference `her2_spand` value in the cohort metadata; the notebook
  output matches that reference to ~0.4%.
- **Sign convention.** Internally `spand_for_slide` returns
  `Moran's I / mean` unmodified (positive = more spatially autocorrelated).
  Notebook 01 negates by `-1` to match the paper / bundled per-patient
  scores convention (positive = more heterogeneous).
- **Deconvolution model ensembling.** Notebook 01 averages all five MLP folds
  for any new slide; the existing `cell_type_deconvolution/` notebooks use
  the held-out-fold pattern for PanopTILs (where every ROI has a designated
  test fold). The averaged ensemble matches the bundled reference deconvolution
  for slide 732799 to correlation 1.000.

## Running

```bash
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/01_spand_pipeline_example.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/02_her2low_vs_her2high_PBCP.ipynb
```

## Citation

Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023
