# qc — ST sample quality control

The six per-sample QC metrics of Shulman et al., *Cell* 2026: the QC functions,
a worked example on one ST section, and the post-QC metadata for the four paper
ST cohorts.

## Contents

- `lib/qc_metrics.py` — the six QC functions: `compute_image_qc` (hematoxylin, eosin, sharpness per H&E tile), `compute_expression_qc` (total UMI, mitochondrial, hemoglobin counts per spot), the paper `THRESHOLDS`, and `passes_qc`.
- `lib/stain_norm.py` — canonical Macenko H&E color normalization (the `macenko_normalizer` class, SPAMS-based), applied before image QC.
- `notebooks/qc_example.ipynb` — worked example computing all six metrics on one ST section (NCBI776, a breast-cancer Visium section from the public [HEST-1k](https://huggingface.co/datasets/MahmoodLab/hest) dataset).
- `data/example/` — that section's inputs: `NCBI776_he_tiles.npz` (200 raw H&E tiles) and `NCBI776_counts.h5ad` (raw counts).
- `data/metadata/st_samples_all_with_qc.csv` — per-section table for all spatial transcriptomics sections across the four paper cohorts (Bassiouni, HEST, Martinez, HTAN), with the six paper-defined QC metrics (hematoxylin/eosin intensity, sharpness, mean total UMI counts per spot, mean mitochondrial log1p, mean hemoglobin log1p), per-metric pass flags, composite `qc_pass_image` / `qc_pass_expression` / `qc_pass`, and the paper-specific `paper_included` / `paper_final` flags.
- `data/metadata/st_samples_retained.csv` — subset where `paper_final == True`, matching the 40 sections retained in the Cell paper Table S1.

## Validation

`compute_image_qc` reproduces the per-sample hematoxylin, eosin, and sharpness values in `data/metadata/st_samples_all_with_qc.csv` to floating-point exactness on the canonical input tiles. This assumes tiles are 224×224 Macenko-normalized RGB; see the docstrings in `lib/qc_metrics.py` for the resolution-sensitivity caveat on the sharpness metric.

## QC thresholds (Shulman et al., Cell 2026, Methods)

| Metric | Threshold | Exclusion direction | Source |
|---|---|---|---|
| Hematoxylin intensity (`H_mean`) | (mean − 2 SD, mean + 2 SD) = (0.029768, 0.053546)* | exclude if outside the interval | TCGA breast-diagnostic reference distribution |
| Eosin intensity (`E_mean`) | (mean − 2 SD, mean + 2 SD) = (0.017486, 0.022495)* | exclude if outside the interval | TCGA breast-diagnostic reference distribution |
| Sharpness (variance of Laplacian) | 0.000637 | exclude if < | TCGA breast-diagnostic threshold |
| Mean total UMI counts per spot | 6193.7 | exclude if < | Bassiouni training cohort, 1st percentile |
| Mean mitochondrial counts (log1p) | 8.3603 | exclude if > | Bassiouni training cohort, 99th percentile |
| Mean hemoglobin counts (log1p) | 0.7355 | exclude if > | Bassiouni training cohort, 99th percentile |

> Hematoxylin and eosin use two-sided exclusion — a section fails if its value falls outside the reference mean ± 2 SD.

Image-QC metrics are computed on color-normalized 224×224 H&E tiles inside the tissue mask (Otsu on grayscale). Expression QC uses `scanpy.pp.calculate_qc_metrics`. All metrics are aggregated to one value per section by averaging across tiles / spots.

## HTAN cohort definition

The HTAN pool used in the paper is the breast-section subset with `Diagnosis ∈ {"Ductal carcinoma NOS", "Ductal carcinoma in situ NOS"}`. Sections with lobular, mixed lobular/ductal, or infiltrating-duct diagnoses are excluded upstream of QC. This filter is encoded as the `paper_included` column. After the six-threshold QC is applied (`qc_pass`), the final paper set is `paper_final = paper_included & qc_pass`.

## log1p, not log10

The paper Methods text says "log10", but the mitochondrial and hemoglobin metrics reported by `scanpy.pp.calculate_qc_metrics` are its `log1p_*` outputs — natural log of (1 + counts). The thresholds in the code and in the table above (8.36, 0.74) are log1p values.
