# cell_type_deconvolution — cell-type abundance from H&E-inferred expression

Reproduces the cell-type deconvolution analysis from Shulman et al., *Cell* 2026
(Figures 4C, 4D, and S5B). Path2Space-inferred spatial expression is used to
estimate cancer / lymphocyte / stromal abundance, two ways:

- a **supervised MLP regressor** trained on PanopTILs pathologist-annotated
  fractions, and
- **SpaCET**, an unsupervised reference-based deconvolution method, applied to
  the inferred expression.

Both are benchmarked against pathologist annotations on archival TCGA
diagnostic slides (PanopTILs) and on two annotated Visium slides (HEST).

## Notebooks

| Notebook | Language | Figure | What it does |
|---|---|---|---|
| `01_panoptils_deconvolution.ipynb` | Python | 4C, S5B | PanopTILs (TCGA): 5-fold CV of the supervised model + SpaCET on inferred ST; per-slide PCC and AUC. |
| `02_hest_deconvolution.ipynb` | Python | 4D | Annotated HEST: AUC of SpaCET (measured), Path2Space (inferred), SpaCET (inferred) and HoVer-Net. |
| `03_spacet_deconvolution.ipynb` | R | — | Upstream step: runs SpaCET on Path2Space inferred / measured expression to produce the per-spot cell-type proportion matrices. |

The Python notebooks **load the five fold-models and predict live** — no retraining.
The R notebook runs SpaCET (≈ 1 min per slide); it needs the `SpaCET` and `qs`
R packages.

## Layout

```
cell_type_deconvolution/
├── lib/
│   └── new_SpaCET.R              # SpaCET.deconvolution_new() — deconvolves a gene-by-spot matrix
├── models/
│   ├── deconv_fold_{0-4}.joblib   # MinMaxScaler -> MLPRegressor pipelines (5-fold CV)
│   └── best_feat_fold_{0-4}.txt   # the 1,000 genes selected in each fold
├── data/
│   ├── tcga_panoptils_expression.csv.gz   # Path2Space-inferred expression, 1,448 model genes
│   ├── tcga_panoptils_labels.csv          # pathologist fractions + slide + CV fold per ROI
│   ├── panoptils_spacet_inferred.csv      # SpaCET on inferred ST (35 fine cell types)
│   ├── hest_expression.csv.gz             # Path2Space-inferred expression for HEST spots
│   ├── hest_method_predictions.csv        # per-spot scores + labels: SpaCET (measured), HoVer-Net
│   ├── hest_spacet_inferred_TENX{13,39}.csv  # SpaCET on inferred ST, per HEST slide
│   └── spacet/                            # qs inputs/outputs for the R notebook
│       ├── {panoptils,hest_TENX13,hest_TENX39}_*_expression.qs   # gene-by-spot inputs
│       └── {panoptils_spacet_proportions,hest_*_spacet_*}.qs     # SpaCET proportion outputs
└── notebooks/
    ├── 01_panoptils_deconvolution.ipynb
    ├── 02_hest_deconvolution.ipynb
    └── 03_spacet_deconvolution.ipynb
```

## The model

`deconv_fold_{0-4}.joblib` are the five cross-validation fold models behind the
published figures (`mlp_results_1000_MinMaxScaler_all_frac_main`). Each is an
sklearn `Pipeline`: `MinMaxScaler → MLPRegressor`, trained on the top 1,000
genes by ANOVA F-statistic (`SelectKBest`, `f_regression`) within the fold, with
three outputs (TILs / Stromal / Epithelial fractions). Folds keep all ROIs from
a patient together to avoid leakage. For TCGA the fold model predicts on its
held-out ROIs (5-fold CV); for HEST the five models are averaged.

## SpaCET deconvolution (`03`)

`03_spacet_deconvolution.ipynb` runs SpaCET on the Path2Space expression and
writes the per-spot cell-type proportion matrices to `data/spacet/`. Inputs and
outputs are stored as `.qs` (compact R serialization): inferred expression in
log10 space (the notebook applies `10^x − 1`), measured expression as raw
counts. SpaCET returns 35 fine cell types, which the Python notebooks collapse
to cancer / lymphocyte / stromal.

## Outputs

**PanopTILs — Figure S5B / 4C** (per-slide mean over the 20 slides with ≥20 ROIs):

| Cell type | SpaCET PCC | SpaCET AUC | Model PCC | Model AUC |
|---|---|---|---|---|
| Cancer | 0.69 | 0.85 | 0.79 | 0.88 |
| Lympho. | 0.57 | 0.79 | 0.79 | 0.86 |
| Stromal | 0.40 | 0.85 | 0.59 | 0.84 |

**HEST — Figure 4D** (AUC, mean over the two slides):

| Method | Cancer | Stromal | TILs |
|---|---|---|---|
| SpaCET (Measured GE) | 0.931 | 0.937 | 0.823 |
| Path2Space (Inferred GE) | 0.893 | 0.849 | 0.919 |
| SpaCET (Inferred GE) | 0.891 | 0.894 | 0.931 |
| HoVer-Net (H&E) | 0.789 | 0.713 | 0.865 |

## Related

- Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023.
- Sibling components: `../clustering_spatiotypes/`, `../clustering_response_prediction/`.
