# path2space-extended

Extended Path2Space codebase — analysis and downstream methods beyond the companion repo ([eldadshulman/path2space-companion](https://github.com/eldadshulman/path2space-companion)).

Code accompanies Shulman et al., *Cell* 2026; licensed under Apache 2.0.

## Components

- [`qc/`](./qc/) — quality control pipeline reproducing the paper's per-sample QC: canonical Macenko stain normalization (vendored from the companion repo), the three image metrics (hematoxylin, eosin, sharpness) and three expression metrics (UMI, mitochondrial, hemoglobin), Table S1 metadata, and a worked example.
- [`clustering_spatiotypes/`](./clustering_spatiotypes/) — spatial clustering pipeline (SpaGCN per slide → Seurat cross-patient → SpatioType). TCGA training + METABRIC validation; per-cohort cluster-projection outputs for Cedars-Sinai, PBCP, TransNEO and IMPRESS.
- [`clustering_response_prediction/`](./clustering_response_prediction/) — pCR-prediction models trained on the ST cluster compositions. Two notebooks (chemotherapy and trastuzumab), the fitted joblib models, per-slide pCR labels, HER2-SPAND scores, and the cross-cohort AUC table.
- [`cell_type_deconvolution/`](./cell_type_deconvolution/) — cell-type abundance from H&E-inferred expression (Figures 4C, 4D, S5B). A supervised MLP regressor and SpaCET, benchmarked on PanopTILs (TCGA) and annotated HEST against pathologist annotations.
- [`confounder_analysis/`](./confounder_analysis/) — per-gene partial correlations of predicted vs. measured expression controlling for hematoxylin stain intensity, cancer cell fraction, and total RNA content, across all four ST cohorts. Stage-A per-slide computation documented; Stage-B aggregation (Stouffer's combination, BH FDR, weighted means) reproducible end-to-end from the bundled parquet.
- [`within_cell_type_variation/`](./within_cell_type_variation/) — within-cell-type expression variation (Figures 4E, 4F, S5G, S5H): per-gene predicted-vs-measured correlation inside pathologist-annotated regions and transcriptomic neighborhoods, gene-wise accuracy concordance across cell types, and a marker vs. non-marker gene comparison.

Future components will be added as additional sibling folders.

## Related

- Companion: https://github.com/eldadshulman/path2space-companion
- Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023
