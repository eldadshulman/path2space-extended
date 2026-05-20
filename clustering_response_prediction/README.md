# clustering_response_prediction — pCR prediction from spatial-cluster composition

Two notebooks that train and evaluate the published pCR-prediction models from
Shulman et al., *Cell* 2026:

- **Chemotherapy** model — trained on TransNEO chemo, evaluated on IMPRESS chemo and PBCP chemo.
- **Trastuzumab** model — trained on TransNEO trastuzumab, evaluated on IMPRESS, PBCP and Cedars-Sinai (Ronai_BRCA). The trastuzumab pipeline additionally combines the cluster-based predictor with the HER2 SPAND heterogeneity score.

Per-fold pipeline: `MinMaxScaler → SelectKBest(f_classif, k='all') → LogisticRegression(L1, C=100, class_weight='balanced')`, wrapped in a 5-fold `StratifiedKFold` and averaged into an `EnsembleModel`. Adapted from
`scr/new_tras_clusters/auc_final_chemo_v3.py` and
`scr/new_tras_clusters/auc_final_trastuzumab_new.py`.

## Layout

```
clustering_response_prediction/
├── lib/
│   └── model_classes.py   # EnsembleModel, FeatureAligner, CombinedModel
├── data/
│   ├── response_labels.csv     # slide_name → Response (pCR=1, non-pCR=0) + Cohort
│   └── her2_spand_scores.csv  # per-slide HER2 SPAND heterogeneity score (sign-flipped)
├── models/
│   ├── chemo_ensemble.joblib              # 5-fold logistic ensemble (chemo)
│   ├── trastuzumab_cluster_ensemble.joblib # 5-fold logistic ensemble (trastuzumab, clusters only)
│   └── trastuzumab_combined.joblib         # CombinedModel: cluster + HER2-SPAND
└── notebooks/
    ├── 01_chemo_response_prediction.ipynb
    └── 02_trastuzumab_response_prediction.ipynb
```

The notebooks read the per-cohort per-patient cluster proportion CSVs from
[`../clustering/data/`](../clustering/data/) (`transneo_chemo_cluster_props.csv`,
`impress_trastuzumab_cluster_props.csv`, …), so this folder is a sibling of the
`clustering/` outputs.

## AUCs from the executed notebooks

| Cohort | n | Prevalence | AUC cluster | AUC combined | Paper combined |
|---|---|---|---|---|---|
| **Chemotherapy** | | | | | |
| TransNEO (5-fold CV) | 93 | 0.226 | 0.753 | — | — |
| IMPRESS | 64 | 0.422 | 0.749 | — | — |
| PBCP | 19 | 0.263 | 0.886 | — | — |
| **Trastuzumab** | | | | | |
| TransNEO (5-fold CV) | 61 | 0.311 | 0.867 | **0.900** | 0.894 |
| IMPRESS | 62 | 0.613 | 0.743 | **0.765** | 0.765 |
| PBCP | 18 | 0.611 | 0.896 | **0.896** | 0.879 |
| Cedars-Sinai | 30 | 0.867 | 0.683 | **0.837** | 0.852 |

IMPRESS trastuzumab AUCs (cluster 0.743 / combined 0.765) reproduce the
paper's Table S4 exactly. TransNEO 5-fold CV and PBCP/Cedars-Sinai combined
are within ~0.02 of the paper; the small remaining deltas come from N
differences (paper N = 25/16/26 for Cedars-Sinai / PBCP cluster / PBCP combined
vs our 30/18 — the paper applied an additional patient-level filter we
haven't reproduced).

## Loading the models in Python

```python
import joblib, sys; sys.path.insert(0, "lib")
from model_classes import EnsembleModel, FeatureAligner, CombinedModel

clf = joblib.load("models/chemo_ensemble.joblib")
proba = clf.predict_proba(X)        # X = per-patient × per-cluster proportion matrix

cmb = joblib.load("models/trastuzumab_combined.joblib")
proba = cmb.predict_proba(X_full)   # X_full = clusters + 'her2_spand' column
```

## Related

- Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023.
- Per-cohort cluster compositions: `../clustering/data/`.
- Projection code that produces those compositions: https://github.com/eldadshulman/path2space-extra
