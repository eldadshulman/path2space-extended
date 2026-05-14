# path2space-extended

Private extended Path2Space codebase — analysis and downstream methods beyond the public companion repo ([eldadshulman/path2space-companion](https://github.com/eldadshulman/path2space-companion)).

## Contents

- `data/metadata/st_samples_all_with_qc.csv` — per-section table for all spatial transcriptomics sections across the four paper cohorts (Bassiouni, HEST, Martinez, HTAN), with the six paper-defined QC metrics (hematoxylin/eosin intensity, sharpness, mean total UMI counts per spot, mean mitochondrial log1p, mean hemoglobin log1p), per-metric pass flags, composite `qc_pass_image` / `qc_pass_expression` / `qc_pass`, and the paper-specific `paper_included` / `paper_final` flags.
- `data/metadata/st_samples_retained.csv` — subset where `paper_final == True`, matching the 40 sections retained in the Cell paper Table S1.
- `notebooks/` — analysis notebooks (TBD).

## QC thresholds (Shulman et al., Cell 2026, Methods)

| Metric | Threshold | Exclusion direction | Source |
|---|---|---|---|
| Hematoxylin intensity (`H_mean`) | 0.053546 | exclude if > | TCGA breast-diagnostic mean + 2 SD |
| Eosin intensity (`E_mean`) | 0.022495 | exclude if > | TCGA breast-diagnostic mean + 2 SD |
| Sharpness (variance of Laplacian) | 0.000637 | exclude if < | TCGA breast-diagnostic threshold |
| Mean total UMI counts per spot | 6193.7 | exclude if < | Bassiouni training cohort, 1st percentile |
| Mean mitochondrial counts (log1p) | 8.3603 | exclude if > | Bassiouni training cohort, 99th percentile |
| Mean hemoglobin counts (log1p) | 0.7355 | exclude if > | Bassiouni training cohort, 99th percentile |

Image-QC metrics are computed on color-normalized 224×224 H&E tiles inside the tissue mask (Otsu on grayscale). Expression QC uses `scanpy.pp.calculate_qc_metrics`. All metrics are aggregated to one value per section by averaging across tiles / spots.

## HTAN cohort definition

The HTAN pool used in the paper is the breast-section subset with `Diagnosis ∈ {"Ductal carcinoma NOS", "Ductal carcinoma in situ NOS"}`. Sections with lobular, mixed lobular/ductal, or infiltrating-duct diagnoses are excluded upstream of QC. This filter is encoded as the `paper_included` column. After the six-threshold QC is applied (`qc_pass`), the final paper set is `paper_final = paper_included & qc_pass`.

The paper's pre-QC HTAN pool is defined as 44 slides; our diagnosis-based filter yields 38. The 9 retained slides match the paper exactly. The 6-row denominator gap likely reflects an additional upstream filter (e.g., assay subtype or biopsy timepoint) we did not reverse-engineer.

## Note on log1p vs log10

The Cell paper Methods text refers to log10 for mitochondrial and hemoglobin counts. The actual computed metric and thresholds match log1p (natural log of 1+x), as output by `scanpy.pp.calculate_qc_metrics`. We've named the columns honestly (`qc_mito_log1p`, `qc_hemo_log1p`) and use the log1p thresholds (8.36, 0.74).

## Retention (paper Table S1)

| Cohort | Total sections | After paper inclusion | Final (paper_final) | Paper |
|---|---|---|---|---|
| Bassiouni (TNBC, GSE210616) | 43 | 43 | 22 | 22 / 43 |
| HEST | 5 | 5 | 5 | 5 / 5 |
| Martinez (GSE213688) | 5 | 5 | 4 | 4 / 5 |
| HTAN | 48 | 38 | 9 | 9 / 44 |
| **Total** | **101** | **91** | **40** | **40** |

## Related

- Public companion: https://github.com/eldadshulman/path2space-companion
- Paper: Shulman et al., *Cell* 2026.
