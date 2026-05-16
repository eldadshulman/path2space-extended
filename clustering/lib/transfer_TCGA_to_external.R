#!/usr/bin/env Rscript
# Project the 11 TCGA ST-cluster labels onto an external cohort's per-domain
# Seurat object via anchor-based integration. Functions are factored out so
# the tuned parameters live in one place; the notebook only sees the
# high-level pipeline.
#
# Pipeline:
#   reference.qs ----+
#                    +--> transfer_cluster_labels()  --> per-domain predictions
#   query.qs --------+
#                                                     |
#                                                     v
#                              aggregate_to_patient_proportions()
#                                                     |
#                                                     v
#                                  per-patient × per-cluster CSV
#
# Usage:
#   Rscript transfer_TCGA_to_external.R \
#       <reference.qs> <query.qs> <out_per_patient_csv> [<refdata_field>]

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(data.table)
  library(dplyr)
  library(tidyr)
})


# ----- Step 1: Seurat anchor-based label transfer ----------------------------
# All tuning parameters are kept here. They were chosen by the grid search
# in spatial_clusters/3_transfer_check_param_METABRIC.r; the column-name
# encoding in clusters_reviewer_3/analysis_row_1_data.rds pins them down
# (k_anchor_40 / k_score_30 / l2_norm_FALSE / k_filter_50 / nn_method_hnsw /
# k_weight_15 / reduction_pcaproject). dims and npcs are fixed in the
# original loop body, not in the column name.
transfer_cluster_labels <- function(reference, query, refdata_col = "mapped_clusters") {
  anchors <- FindTransferAnchors(
    reference = reference,
    query     = query,
    dims      = 1:20,
    npcs      = 20,
    reduction = "pcaproject",
    k.anchor  = 40,
    k.score   = 30,
    k.filter  = 50,
    l2.norm   = FALSE,
    nn.method = "hnsw"
  )
  preds <- TransferData(
    anchorset        = anchors,
    refdata          = reference@meta.data[[refdata_col]],
    weight.reduction = "pcaproject",
    k.weight         = 15
  )
  data.table(
    domain_id           = rownames(preds),
    slide_name          = query@meta.data$slide_name,
    proportion_of_spots = query@meta.data$proportion_of_spots,
    predicted_cluster   = preds$predicted.id,
    prediction_score    = preds$prediction.score.max
  )
}


# ----- Step 2: aggregate per-domain predictions to per-patient proportions ---
aggregate_to_patient_proportions <- function(per_domain, n_clusters = 11) {
  wide <- per_domain[, .(prop = sum(proportion_of_spots, na.rm = TRUE)),
                     by = .(slide_name, predicted_cluster)] |>
          tidyr::pivot_wider(names_from = predicted_cluster,
                             values_from = prop,
                             values_fill = 0)
  cluster_cols <- paste0("Cluster_", seq_len(n_clusters))
  for (c in cluster_cols) if (!c %in% colnames(wide)) wide[[c]] <- 0
  wide <- wide[, c("slide_name", cluster_cols)]
  wide
}


# ----- Main script -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  reference_qs <- args[1]
  query_qs     <- args[2]
  out_csv      <- args[3]
  refdata_col  <- ifelse(length(args) >= 4, args[4], "mapped_clusters")

  cat(sprintf("Loading reference: %s\n", reference_qs))
  reference <- qread(reference_qs)
  cat(sprintf("Loading query:     %s\n", query_qs))
  query     <- qread(query_qs)

  cat("Running label transfer...\n")
  per_domain <- transfer_cluster_labels(reference, query, refdata_col = refdata_col)

  cat("Aggregating to per-patient × per-cluster proportions...\n")
  per_patient <- aggregate_to_patient_proportions(per_domain)

  fwrite(per_patient, out_csv)
  cat(sprintf("Wrote %d patients to %s\n", nrow(per_patient), out_csv))
}
