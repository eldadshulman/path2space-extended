# clustering — spatial cluster derivation and external-cohort projection

Reproduces the spatial clustering analysis from Shulman et al., *Cell* 2026: per-slide
SpaGCN → cross-patient Seurat clustering → SpatioType assignment on TCGA, plus the
projection of both ST clusters and SpatioTypes onto the external METABRIC cohort.

## Pipeline

| Step | Notebook | Language | What it does |
|---|---|---|---|
| 1. Per-slide spatial domains | `notebooks/01_spagcn_per_slide_example.ipynb` | Python | Runs SpaGCN on the pseudo-spots of one TCGA slide to identify spatially coherent expression domains (~6.6 domains/slide on average). |
| 2. Cross-patient ST clusters | `notebooks/02_seurat_cross_patient_TCGA.ipynb` | R (Seurat) | Log-normalizes the domain-averaged expression profiles across all TCGA slides and jointly clusters them with `FindClusters` to yield the 11 ST clusters that are consistent across patients. |
| 3. SpatioType assignment | `notebooks/03_spatiotype_assignment.ipynb` | Python | Hierarchically clusters TCGA patients on their 11-cluster composition vectors (k=5), labels three SpatioTypes (Proliferation-Enriched, Immune-Modulated, Immune-Inactive) plus two small Metaclusters, and validates by KM. |
| 4. Project ST clusters → METABRIC | `notebooks/04_project_clusters_to_METABRIC.ipynb` | Python | Renders the Figure-S7A heatmap of per-patient cluster proportions after transferring the 11 TCGA ST clusters onto METABRIC via Seurat anchor-based integration (`FindTransferAnchors` + `TransferData`; run upstream in R). |
| 5. Project SpatioTypes → METABRIC | `notebooks/05_project_spatiotypes_to_METABRIC.ipynb` | Python | Validates the projected SpatioType labels against METABRIC survival: Figure-S7B main KM, Cox HR for Immune-Inactive (stage + Age>65 adjusted), Figure-S7C KM by clinical subtype, Figure-S7D KM within the poor-NPI subgroup, and the LRT for SpatioType beyond NPI. |

## Cohorts

Current: TCGA (training), METABRIC (external validation). Future additions: TransNEO,
IMPRESS (treatment-response cohorts).

## Data

Small input artifacts live in `data/`:

| File | Used by |
|---|---|
| `cluster_props_TCGA.csv` | nb 03 (input) |
| `cluster_validity_metrics.csv` | nb 03 (output) |
| `spatiotypes_TCGA_k5.csv` | nb 03 (output) — TCGA SpatioType assignments |
| `tsne.csv` | nb 02 (figure data) |
| `metabric_cluster_props.csv` | nb 04 (input) — METABRIC per-patient cluster composition |
| `metabric_survival.csv` | nb 05 (input) — METABRIC survival + SpatioType + clinical |

Per-spot expression tables, per-slide AnnData files, Seurat `.qs` objects, and full
WSI image inputs are NOT committed — they live on Biowulf. The notebooks point at
internal paths for those; replace with your own to re-run end-to-end.

## Related

Paper: Shulman et al., *Cell* 2026. doi:10.1016/j.cell.2026.04.023.
The clustering analysis is described in Methods → "ST clusters" and
"External cohort transfer"; the validation results are in the paragraph about
METABRIC and Figure S7A–D.
