# clustering — spatial cluster derivation and external-cohort projection

Reproduces the spatial clustering analysis from Shulman et al., *Cell* 2026: per-slide
SpaGCN → cross-patient Seurat clustering → SpatioType assignment on TCGA, plus the
projection of both ST clusters and SpatioTypes onto external cohorts (METABRIC,
Cedars-Sinai, PBCP).

## Pipeline

| Step | Notebook | Language | What it does |
|---|---|---|---|
| 1. Per-slide spatial domains | `notebooks/01_spagcn_per_slide_example.ipynb` | Python | Runs SpaGCN on the pseudo-spots of one TCGA slide to identify spatially coherent expression domains (~6.6 domains/slide on average). |
| 2. Cross-patient ST clusters | `notebooks/02_seurat_cross_patient_TCGA.ipynb` | R (Seurat) | Log-normalizes the domain-averaged expression profiles across all TCGA slides and jointly clusters them with `FindClusters` to yield the 11 ST clusters that are consistent across patients. |
| 3. SpatioType assignment | `notebooks/03_spatiotype_assignment.ipynb` | Python | Hierarchically clusters TCGA patients on their 11-cluster composition vectors (k=5), labels three SpatioTypes (Proliferation-Enriched, Immune-Modulated, Immune-Inactive) plus two small Metaclusters, and validates by KM. |
| 4. Project ST clusters → METABRIC | `notebooks/04_project_clusters_to_METABRIC.ipynb` | Python | Figure-S7A heatmap of per-patient cluster proportions after transferring the 11 TCGA ST clusters onto METABRIC via Seurat anchor-based integration. |
| 5. Project SpatioTypes → METABRIC | `notebooks/05_project_spatiotypes_to_METABRIC.ipynb` | Python | Validates the projected SpatioType labels against METABRIC survival: Figure-S7B main KM, Cox HR for Immune-Inactive (stage + Age>65 adjusted), Figure-S7C KM by clinical subtype, Figure-S7D KM within the poor-NPI subgroup, and the LRT for SpatioType beyond NPI. |
| 6. Project ST clusters → trastuzumab cohorts | `notebooks/06_project_clusters_to_trastuzumab_cohorts.ipynb` | Python | Same anchor-based transfer, cohort-specific parameters, applied to Cedars-Sinai (Ronai_BRCA) and the trastuzumab arm of PBCP. Z-score heatmap per cohort. |
| 7. Project ST clusters → chemotherapy cohorts | `notebooks/07_project_clusters_to_chemotherapy_cohorts.ipynb` | Python | Same anchor-based transfer applied to the chemo arm of PBCP. Set up to extend to TransNEO-Chemo / IMPRESS-Chemo when added. |

The Seurat anchor transfer (notebooks 04, 06, 07) uses **cohort-specific tuned
parameters** chosen by separate grid searches. The presets are kept in one place —
`TRANSFER_PRESETS` in [`lib/transfer_TCGA_to_external.R`](lib/transfer_TCGA_to_external.R) —
and selected by name from the notebook (`METABRIC`, `Cedars_Sinai`,
`PBCP_Trastuzumab`, `PBCP_Chemo`).

## Cohorts

| Cohort | Use | Slides | Notebook |
|---|---|---|---|
| TCGA-BRCA | Training (clusters + SpatioTypes) | 858 | 01–03 |
| METABRIC | External validation | 141 | 04, 05 |
| Cedars-Sinai (Ronai_BRCA) | Trastuzumab response | 30 | 06 |
| PBCP — trastuzumab arm | Trastuzumab response | 18 | 06 |
| TransNEO — trastuzumab arm | Trastuzumab response | 61 | 06 |
| IMPRESS — trastuzumab arm | Trastuzumab response | 62 | 06 |
| PBCP — chemotherapy arm | Chemotherapy response | 19 | 07 |
| TransNEO — chemotherapy arm | Chemotherapy response | 93 | 07 |
| IMPRESS — chemotherapy arm | Chemotherapy response | 64 | 07 |

## Data

Small input artifacts live in `data/`:

| File | Used by |
|---|---|
| `cluster_props_TCGA.csv` | nb 03 (input) |
| `cluster_validity_metrics.csv` | nb 03 (output) |
| `spatiotypes_TCGA_k5.csv` | nb 03 (output) — TCGA SpatioType assignments |
| `tcga_spatiotypes_reference.csv` | nb 05 (input) — TCGA centroids reference for Aitchison projection |
| `tsne.csv` | nb 02 (figure data) |
| `metabric_cluster_props.csv` | nb 04 input / nb 05 input |
| `metabric_spatiotypes.csv` | nb 04 / 05 — METABRIC projected SpatioType labels |
| `metabric_survival.csv` | nb 05 — METABRIC survival + clinical |
| `best_clust.csv` | nb 04 — per-domain METABRIC cluster predictions |
| `cedars_sinai_query.qs` | nb 06 — Cedars-Sinai per-domain Seurat object (13 MB) |
| `cedars_sinai_per_domain.csv` / `cedars_sinai_cluster_props.csv` | nb 06 outputs |
| `pbcp_query.qs` | nb 06 / 07 — PBCP per-domain Seurat object (18 MB) |
| `pbcp_trastuzumab_per_domain.csv` / `pbcp_trastuzumab_cluster_props.csv` | nb 06 outputs |
| `pbcp_chemo_per_domain.csv` / `pbcp_chemo_cluster_props.csv` | nb 07 outputs |

The 533 MB TCGA reference Seurat `.qs` is held on Biowulf and not committed; the
notebooks point at its upstream path. Per-spot expression tables and WSI images
also stay upstream.

## Related

Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023.
The clustering analysis is described in Methods → "ST clusters" and
"External cohort transfer"; the validation results are in the paragraph about
METABRIC and Figure S7A–D.
