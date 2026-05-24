# spand — SPAND score and the HER2-low / HER2-high trastuzumab analysis

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
  bundled per-patient HER2-SPAND for that slide to ~0.4%, and contrasts the
  three Moran's I values (raw ERBB2, raw WP673 NES, cancer-restricted WP673)
  to show why cancer restriction matters.
- `notebooks/02_her2low_vs_her2high.ipynb` — three-cohort clinical analysis
  using precomputed per-patient HER2-SPAND vs measured HER2 predictors:
  **PBCP** (18 HER2+ patients, vs bulk ERBB2 log₂ TPM), **IMPRESS** (62
  trastuzumab-arm, vs HER2 FISH ratio + IHC), **TransNEO** (61 trastuzumab,
  vs bulk ERBB2 + bulk ErbB-pathway score). Each cohort: pooled and
  HER2-low / -high split AUC with bootstrap 95% CIs.
- `data/example_slide_predicted_expression.parquet` — one slide (PBCP 732799),
  Path2Space predicted expression restricted to the ~1,400 genes used by the
  deconvolution MLP feature list (~9 MB).
- `data/example_slide_gsea_her2_nes.parquet` — per-(spot, pathway) GSEA NES
  for the six ErbB-related pathways on that slide (~90 KB).
- `data/pbcp_her2_scores.parquet` — 18 PBCP HER2+ patients (~4 KB).
- `data/impress_her2_scores.parquet` — 62 IMPRESS trastuzumab-arm patients
  with FISH + IHC (~6 KB).
- `data/transneo_her2_scores.parquet` — 61 TransNEO trastuzumab patients with
  bulk ERBB2 and bulk pathway scores (~6 KB).
- `data/erbb_pathway_genes.json` — Wikipathways WP673 gene list. **Reference
  material only**: included for readers who want to recompute the upstream
  per-spot GSEA. Notebook 01 does not consume this JSON — it reads the
  precomputed `example_slide_gsea_her2_nes.parquet` directly.

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
(notebook 02): per-patient HER2-SPAND for three trastuzumab-treated cohorts
(PBCP, IMPRESS, TransNEO) is bundled as the three `*_her2_scores.parquet`
files. Computing those scores at scale across all three cohorts requires the
full Path2Space prediction pipeline (companion repo) and is not reproduced
inline here.

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
CC-BY 3.0; source URL in the JSON). The JSON is *reference material*; the
notebooks never read it.

## Reproduction notes

- **Bundled scope.** Only one example slide of predicted expression ships
  (~9 MB after restriction to deconvolution feature genes). No cohort-level
  predicted-expression data is bundled; the three per-patient tables are
  precomputed scores only. The full PBCP / IMPRESS / TransNEO predicted-
  expression matrices live upstream in the companion repo and are not
  redistributed here.
- **Example slide.** PBCP `732799` was chosen as the smallest slide with a
  non-null reference `her2_spand` value in the cohort metadata; the notebook
  output matches that reference to ~0.4%.
- **HER2-SPAND sign convention across cohorts.** All three bundled
  `her2_spand` columns are stored with the convention **higher = more
  predictive of response = 1**. PBCP and TransNEO sources already used this
  convention; the IMPRESS source pipeline produced `her2_spand` with the
  opposite sign, and was **sign-flipped at bundle time** so the three
  cohorts share the same direction. The sign flip is documented in the
  bundling script and surfaced in notebook 02's intro.
- **PBCP subset.** The 18 patients in `pbcp_her2_scores.parquet` are
  exactly the HER2+ half of PBCP's Application Cohort (the paper's
  HER2+ n=18 subset). All have `HER2_status == "POS"`. The matching
  HER2− n=19 subset has no `her2_spand` in the source metadata and is
  not included here.
- **Per-stratum sample sizes are small.** The HER2-low / -high split in
  notebook 02 yields n=9 (PBCP), n=31 (IMPRESS), n=30 (TransNEO) per side.
  Bootstrap 95% CIs are wide — the per-stratum analysis is exploratory,
  not confirmatory. The pooled-cohort AUC bars are the more reliable
  summary at these sample sizes.
- **Deconvolution model ensembling.** Notebook 01 averages all five MLP folds
  for any new slide; the existing `cell_type_deconvolution/` notebooks use
  the held-out-fold pattern for PanopTILs (where every ROI has a designated
  test fold). The averaged ensemble matches the bundled reference deconvolution
  for slide 732799 to correlation 1.000.

## Running

```bash
pip install -r requirements.txt
jupyter nbconvert --to notebook --execute --inplace notebooks/01_spand_pipeline_example.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/02_her2low_vs_her2high.ipynb
```

## Citation

Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023
