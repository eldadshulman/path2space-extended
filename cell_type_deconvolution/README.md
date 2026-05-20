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

| Notebook | Figure | What it does |
|---|---|---|
| `01_panoptils_deconvolution.ipynb` | 4C, S5B | PanopTILs (TCGA): 5-fold CV of the supervised model + SpaCET on inferred ST; per-slide PCC and AUC. |
| `02_hest_deconvolution.ipynb` | 4D | Annotated HEST: AUC of SpaCET (measured), Path2Space (inferred), SpaCET (inferred) and HoVer-Net. |

Both notebooks **load the five fold-models and predict live** — no retraining.

## Layout

```
cell_type_deconvolution/
├── models/
│   ├── deconv_fold_{0-4}.joblib   # MinMaxScaler -> MLPRegressor pipelines (5-fold CV)
│   └── best_feat_fold_{0-4}.txt   # the 1,000 genes selected in each fold
├── data/
│   ├── tcga_panoptils_expression.csv.gz   # Path2Space-inferred expression, 1,448 model genes
│   ├── tcga_panoptils_labels.csv          # pathologist fractions + slide + CV fold per ROI
│   ├── panoptils_spacet_inferred.csv      # SpaCET on inferred ST (35 fine cell types)
│   ├── hest_expression.csv.gz             # Path2Space-inferred expression for HEST spots
│   ├── hest_method_predictions.csv        # per-spot scores + labels: SpaCET (measured), HoVer-Net
│   └── hest_spacet_inferred_TENX{13,39}.csv  # SpaCET on inferred ST, per HEST slide
└── notebooks/
    ├── 01_panoptils_deconvolution.ipynb
    └── 02_hest_deconvolution.ipynb
```

## The model

`deconv_fold_{0-4}.joblib` are the five cross-validation fold models behind the
published figures (`mlp_results_1000_MinMaxScaler_all_frac_main`). Each is an
sklearn `Pipeline`: `MinMaxScaler → MLPRegressor`, trained on the top 1,000
genes by ANOVA F-statistic (`SelectKBest`, `f_regression`) within the fold, with
three outputs (TILs / Stromal / Epithelial fractions). Folds keep all ROIs from
a patient together to avoid leakage. For TCGA the fold model predicts on its
held-out ROIs (5-fold CV); for HEST the five models are averaged.

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
- Sibling components: `../clustering/`, `../clustering_response_prediction/`.
