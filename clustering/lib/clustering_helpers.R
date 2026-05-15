# Helpers for the spatial-clustering pipeline.
#
# The notebooks in clustering/notebooks/ call these wrappers so the
# parameter choices live in one place. Adjust here, not in the notebooks.

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


# Compute per-patient × cluster composition matrix from a clustered Seurat object.
# Returns a data.frame: rows = patients (tissue_submitter_id), cols = clusters.
patient_cluster_composition <- function(seurat_obj,
                                        patient_col = "tissue_submitter_id",
                                        cluster_col = "seurat_clusters") {
  md <- seurat_obj@meta.data
  fraction_df <- md %>%
    dplyr::count(.data[[patient_col]], .data[[cluster_col]], name = "total_spots") %>%
    dplyr::group_by(.data[[patient_col]]) %>%
    dplyr::mutate(fraction_spots = total_spots / sum(total_spots)) %>%
    dplyr::ungroup()

  wide <- fraction_df %>%
    dplyr::select(dplyr::all_of(c(patient_col, cluster_col)), fraction_spots) %>%
    tidyr::pivot_wider(names_from  = dplyr::all_of(cluster_col),
                       values_from = fraction_spots,
                       values_fill = 0)
  as.data.frame(wide)
}


# Hierarchical clustering of patients on their cluster-composition vectors.
# Returns a data.frame with tissue_submitter_id, Metacluster (1..k) and SpatioType label.
assign_spatiotypes <- function(composition_df,
                                patient_col = "tissue_submitter_id",
                                k = 3) {
  mat <- as.matrix(composition_df[, setdiff(colnames(composition_df), patient_col)])
  rownames(mat) <- composition_df[[patient_col]]

  scaled <- scale(mat)
  linkage <- hclust(dist(scaled), method = "ward.D2")
  meta    <- cutree(linkage, k = k)

  out <- data.frame(
    setNames(list(rownames(mat)), patient_col),
    Metacluster = unname(meta),
    stringsAsFactors = FALSE
  )

  spatiotype_levels <- c('Proliferative', 'Immune-Modulated', 'Immune-Inactive')
  out$SpatioType <- dplyr::case_when(
    out$Metacluster == 3 ~ spatiotype_levels[1],
    out$Metacluster == 1 ~ spatiotype_levels[2],
    out$Metacluster == 2 ~ spatiotype_levels[3],
    TRUE                 ~ NA_character_
  )
  out$SpatioType <- factor(out$SpatioType, levels = spatiotype_levels, ordered = TRUE)

  attr(out, "linkage") <- linkage
  attr(out, "scaled_features") <- scaled
  out
}
