# Helpers for the cross-patient clustering step.
#
# Used by clustering/notebooks/02_seurat_cross_patient_TCGA.ipynb. Keeping
# the parameter choices here lets the notebook read as a narrative rather
# than a tuned recipe.
#
# The Python helper for the SpatioType assignment step (notebook 03) lives
# in clustering_helpers.py — keep this file focused on the Seurat wrappers.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})


# Cross-patient clustering on a per-domain Seurat object.
# Produces the 11 ST clusters used as the unit of cross-cohort transfer.
run_cross_patient_clustering <- function(seurat_obj) {
  params <- list(
    normalization.method   = "LogNormalize",
    selection.method       = "mean.var.plot",
    dims                   = 1:30,
    k.param                = 5,
    resolution             = 0.5
  )

  seurat_obj <- NormalizeData(seurat_obj, normalization.method = params$normalization.method)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = params$selection.method)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = params$dims, k.param = params$k.param)
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution   = params$resolution,
    cluster.name = paste0("clusters_", params$k.param, "_", length(params$dims), "_", params$resolution)
  )
  seurat_obj
}


# Frozen mapping from numeric Seurat IDs to stable Cluster_1..Cluster_11 labels.
# Derived once on TCGA and held fixed so the cluster identity matches the paper.
ST_CLUSTER_MAPPING <- c('0' = 'Cluster_1',
                        '1' = 'Cluster_2',
                        '10' = 'Cluster_3',
                        '2' = 'Cluster_4',
                        '3' = 'Cluster_5',
                        '4' = 'Cluster_6',
                        '5' = 'Cluster_7',
                        '6' = 'Cluster_8',
                        '7' = 'Cluster_9',
                        '8' = 'Cluster_10',
                        '9' = 'Cluster_11')


apply_st_cluster_mapping <- function(seurat_obj) {
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(mapped_clusters = recode(as.character(seurat_clusters), !!!ST_CLUSTER_MAPPING))
  Idents(seurat_obj) <- seurat_obj@meta.data$mapped_clusters
  seurat_obj
}
