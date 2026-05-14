# qc — post-QC ST sample metadata

Post-QC ST sample metadata for the four ST cohorts in Shulman et al., *Cell* 2026, with the six paper-defined per-sample QC metrics as columns.

## Contents

- `data/metadata/st_samples_all_with_qc.csv` — per-section table for all spatial transcriptomics sections across the four paper cohorts (Bassiouni, HEST, Martinez, HTAN), with the six paper-defined QC metrics (hematoxylin/eosin intensity, sharpness, mean total UMI counts per spot, mean mitochondrial log1p, mean hemoglobin log1p), per-metric pass flags, composite `qc_pass_image` / `qc_pass_expression` / `qc_pass`, and the paper-specific `paper_included` / `paper_final` flags.
- `data/metadata/st_samples_retained.csv` — subset where `paper_final == True`, matching the 40 sections retained in the Cell paper Table S1.

## QC thresholds (Shulman et al., Cell 2026, Methods)

| Metric | Threshold | Exclusion direction | Source |
|---|---|---|---|
| Hematoxylin intensity (`H_mean`) | (mean − 2 SD, mean + 2 SD) = (0.029768, 0.053546)* | exclude if outside the interval | TCGA breast-diagnostic reference distribution |
| Eosin intensity (`E_mean`) | (mean − 2 SD, mean + 2 SD) = (0.017486, 0.022495)* | exclude if outside the interval | TCGA breast-diagnostic reference distribution |
| Sharpness (variance of Laplacian) | 0.000637 | exclude if < | TCGA breast-diagnostic threshold |
| Mean total UMI counts per spot | 6193.7 | exclude if < | Bassiouni training cohort, 1st percentile |
| Mean mitochondrial counts (log1p) | 8.3603 | exclude if > | Bassiouni training cohort, 99th percentile |
| Mean hemoglobin counts (log1p) | 0.7355 | exclude if > | Bassiouni training cohort, 99th percentile |

> The paper Methods describes two-sided exclusion (outside mean ± 2 SD). The published implementation in `qc_analysis/collect_image_newest.ipynb` applies only the upper bound. The CSVs here match the published implementation; lower bounds are recorded for completeness.

Image-QC metrics are computed on color-normalized 224×224 H&E tiles inside the tissue mask (Otsu on grayscale). Expression QC uses `scanpy.pp.calculate_qc_metrics`. All metrics are aggregated to one value per section by averaging across tiles / spots.

## HTAN cohort definition

The HTAN pool used in the paper is the breast-section subset with `Diagnosis ∈ {"Ductal carcinoma NOS", "Ductal carcinoma in situ NOS"}`. Sections with lobular, mixed lobular/ductal, or infiltrating-duct diagnoses are excluded upstream of QC. This filter is encoded as the `paper_included` column. After the six-threshold QC is applied (`qc_pass`), the final paper set is `paper_final = paper_included & qc_pass`.

## Note on log1p vs log10

The Cell paper Methods refers to log10; the computed metric and thresholds are log1p (natural log of 1+x), as output by `scanpy.pp.calculate_qc_metrics`. Columns are named `qc_mito_log1p` / `qc_hemo_log1p`, and the thresholds (8.36, 0.74) are log1p values.
