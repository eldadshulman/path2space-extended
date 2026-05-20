# path2space-extended

Extended Path2Space codebase — analysis and downstream methods beyond the companion repo ([eldadshulman/path2space-companion](https://github.com/eldadshulman/path2space-companion)).

## Components

- [`qc/`](./qc/) — post-QC ST sample metadata reproducing Table S1 of Shulman et al., *Cell* 2026, with paper-defined per-sample QC metrics as columns. Covers 4 ST cohorts (TNBC/Bassiouni, HEST, Martinez, HTAN).
- [`clustering/`](./clustering/) — spatial clustering pipeline (SpaGCN per slide → Seurat cross-patient → SpatioType). TCGA training + METABRIC validation; per-cohort cluster-projection outputs for Cedars-Sinai, PBCP, TransNEO and IMPRESS.
- [`clustering_response_prediction/`](./clustering_response_prediction/) — pCR-prediction models trained on the ST cluster compositions. Two notebooks (chemotherapy and trastuzumab), the fitted joblib models, per-slide pCR labels, HER2-SPAND scores, and the cross-cohort AUC table.
- [`cell_type_deconvolution/`](./cell_type_deconvolution/) — cell-type abundance from H&E-inferred expression (Figures 4C, 4D, S5B). A supervised MLP regressor and SpaCET, benchmarked on PanopTILs (TCGA) and annotated HEST against pathologist annotations.

Future components will be added as additional sibling folders (e.g., `spand/`).

## Related

- Companion: https://github.com/eldadshulman/path2space-companion
- Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023
