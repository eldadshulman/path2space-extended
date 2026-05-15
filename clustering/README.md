# clustering — spatial cluster derivation pipeline

Reproduces the spatial clustering analysis from Shulman et al., *Cell* 2026, on the TCGA cohort. Three notebooks corresponding to the three steps in the Methods "ST clusters" section.

## Pipeline

| Step | Notebook | Language | What it does |
|---|---|---|---|
| 1. Per-slide spatial domains | `notebooks/01_spagcn_per_slide_example.ipynb` | Python | Runs SpaGCN on the pseudo-spots of one TCGA slide to identify spatially coherent expression domains (~6.6 domains/slide on average). |
| 2. Cross-patient ST clusters | `notebooks/02_seurat_cross_patient_TCGA.ipynb` | R (Seurat) | Log-normalizes the domain-averaged expression profiles across all TCGA slides and jointly clusters them with `FindClusters` to yield the 11 ST clusters that are consistent across patients. |
| 3. SpatioType assignment | `notebooks/03_spatiotype_assignment.ipynb` | R | Aggregates per-domain cluster labels into per-patient SpatioType labels for downstream survival and biomarker analyses. |

## Cohorts

Current: TCGA. Future additions: METABRIC, in-house clinical cohorts.

## Data

Input pseudo-spot expression tables, domain-averaged profiles, and cluster outputs are NOT committed to this repo — they're large and live on Biowulf. The notebooks reference internal paths; replace with your own paths to reuse.

## Related

Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023. The clustering analysis is described in Methods → "ST clusters".
