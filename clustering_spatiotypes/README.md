# clustering_spatiotypes — spatial cluster derivation and SpatioType assignment

Reproduces the spatial clustering analysis from Shulman et al., *Cell* 2026:
per-slide SpaGCN → cross-patient Seurat clustering → SpatioType assignment on
TCGA → projection onto METABRIC for survival validation.

The **external-cohort projection pipeline** (TCGA → Cedars-Sinai / PBCP /
TransNEO / IMPRESS) has moved to a sibling repo:
[`eldadshulman/path2space-extra`](https://github.com/eldadshulman/path2space-extra).
That repo holds the code and the per-cohort Seurat query objects; this repo
keeps the **outputs** of that pipeline (per-domain and per-patient cluster
proportion CSVs) so they're easy to download without pulling the 150 MB of
`.qs` files.

## Pipeline (notebooks in this repo)

| Step | Notebook | Language | What it does |
|---|---|---|---|
| 1. Per-slide spatial domains | `notebooks/01_spagcn_per_slide_example.ipynb` | Python | SpaGCN on one TCGA slide. |
| 2. Cross-patient ST clusters | `notebooks/02_seurat_cross_patient_TCGA.ipynb` | R (Seurat) | Yields the 11 ST clusters. |
| 3. SpatioType assignment | `notebooks/03_spatiotype_assignment.ipynb` | Python | Hierarchical clustering of TCGA patients into 5 metaclusters → 3 named SpatioTypes + 2 small. |
| 4. Project ST clusters → METABRIC | `notebooks/04_project_clusters_to_METABRIC.ipynb` | Python | Figure-S7A heatmap. |
| 5. Project SpatioTypes → METABRIC | `notebooks/05_project_spatiotypes_to_METABRIC.ipynb` | Python | Figure-S7B–D KM + Cox + LRT. |

## Cluster-projection outputs for the treatment-response cohorts

Produced by [`path2space-extra`](https://github.com/eldadshulman/path2space-extra)
nb 06 (trastuzumab) and nb 07 (chemotherapy). Shipped here as plain CSVs:

| Cohort | Per-domain | Per-patient × per-cluster |
|---|---|---|
| Cedars-Sinai (trastuzumab) | `cedars_sinai_per_domain.csv` | `cedars_sinai_cluster_props.csv` |
| PBCP — trastuzumab | `pbcp_trastuzumab_per_domain.csv` | `pbcp_trastuzumab_cluster_props.csv` |
| PBCP — chemo | `pbcp_chemo_per_domain.csv` | `pbcp_chemo_cluster_props.csv` |
| TransNEO — trastuzumab | `transneo_trastuzumab_per_domain.csv` | `transneo_trastuzumab_cluster_props.csv` |
| TransNEO — chemo | `transneo_chemo_per_domain.csv` | `transneo_chemo_cluster_props.csv` |
| IMPRESS — trastuzumab | `impress_trastuzumab_per_domain.csv` | `impress_trastuzumab_cluster_props.csv` |
| IMPRESS — chemo | `impress_chemo_per_domain.csv` | `impress_chemo_cluster_props.csv` |

## Other data files

| File | Used by |
|---|---|
| `cluster_props_TCGA.csv` | nb 03 input |
| `cluster_validity_metrics.csv` | nb 03 output |
| `spatiotypes_TCGA_k5.csv` | nb 03 output — TCGA SpatioType assignments |
| `tcga_spatiotypes_reference.csv` | nb 05 input — TCGA centroids for Aitchison projection |
| `tsne.csv` | nb 02 figure data |
| `metabric_cluster_props.csv` | nb 04 input |
| `metabric_spatiotypes.csv` | nb 04 / 05 — METABRIC projected SpatioType labels |
| `metabric_survival.csv` | nb 05 — METABRIC survival + clinical |
| `best_clust.csv` | nb 04 — per-domain METABRIC cluster predictions |

The 533 MB TCGA reference Seurat `.qs` is held on Biowulf and not committed.

## Related

- Code + per-cohort query inputs: https://github.com/eldadshulman/path2space-extra
- Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023.
