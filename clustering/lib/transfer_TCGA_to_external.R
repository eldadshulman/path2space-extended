#!/usr/bin/env Rscript
# Project the 11 TCGA ST-cluster labels onto an external cohort's per-domain
# Seurat object via anchor-based integration.
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
# Per-cohort tuned parameter presets are kept in TRANSFER_PRESETS below.
# They were chosen by the grid searches in
# spatial_clusters/3_transfer_check_param_METABRIC.r (for METABRIC) and
# spatial_clusters/2_Ronai_Seurat_winningc.r (for Cedars-Sinai / Ronai_BRCA).
#
# Usage:
#   Rscript transfer_TCGA_to_external.R <reference.qs> <query.qs> <out_csv> [cohort]
#     cohort: METABRIC (default) or Cedars_Sinai

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(data.table)
  library(dplyr)
  library(tidyr)
})


# ----- Cohort-specific tuned parameter presets -------------------------------
TRANSFER_PRESETS <- list(
  # Recovered from clusters_reviewer_3/analysis_row_1_data.rds; the dims/npcs
  # values are fixed in the source loop, not the column-name encoding.
  METABRIC = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 40,    k.score  = 30,    k.filter  = 50,
    l2.norm   = FALSE, nn.method = "hnsw",
    weight.reduction = "pcaproject", k.weight = 15
  ),
  # Hardcoded params at the top of spatial_clusters/2_Ronai_Seurat_winningc.r.
  Cedars_Sinai = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 195,   k.score  = 200,   k.filter  = 50,
    l2.norm   = TRUE,  nn.method = "rann",
    weight.reduction = "pcaproject", k.weight = 15
  ),
  # Column-name encoded in PBCP/scr/get_aucs.ipynb:
  #   pred_id_k_anchor_75_k_score_30_l2_norm_TRUE_k_filter_200_
  #   nn_method_rann_k_weight_115_reduction_pcaproject
  PBCP_Trastuzumab = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 75,    k.score  = 30,    k.filter  = 200,
    l2.norm   = TRUE,  nn.method = "rann",
    weight.reduction = "pcaproject", k.weight = 115
  ),
  # Column-name encoded in PBCP/scr/get_aucs.ipynb:
  #   pred_id_k_anchor_145_k_score_135_l2_norm_FALSE_k_filter_100_
  #   nn_method_hnsw_k_weight_55_reduction_pcaproject
  PBCP_Chemo = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 145,   k.score  = 135,   k.filter  = 100,
    l2.norm   = FALSE, nn.method = "hnsw",
    weight.reduction = "pcaproject", k.weight = 55
  ),
  # Column-name encoded in scr/new_tras_clusters/auc_final_chemo_v3.py:
  #   pred_id_k_anchor_20_k_score_10_l2_norm_TRUE_k_filter_50_
  #   nn_method_rann_k_weight_25_reduction_pcaproject
  TransNEO_Chemo = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 20,    k.score  = 10,    k.filter  = 50,
    l2.norm   = TRUE,  nn.method = "rann",
    weight.reduction = "pcaproject", k.weight = 25
  ),
  # Column-name encoded in scr/new_tras_clusters/auc_final_chemo_v3.py:
  #   pred_id_k_anchor_3_k_score_2_l2_norm_TRUE_k_filter_200_
  #   nn_method_hnsw_k_weight_5_reduction_pcaproject
  IMPRESS_Chemo = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 3,     k.score  = 2,     k.filter  = 200,
    l2.norm   = TRUE,  nn.method = "hnsw",
    weight.reduction = "pcaproject", k.weight = 5
  ),
  # Column-name encoded in scr/new_tras_clusters/auc_final_trastuzumab_new.py:
  #   pred_id_k_anchor_120_k_score_200_l2_norm_FALSE_k_filter_NA_
  #   nn_method_annoy_k_weight_15_reduction_pcaproject
  # NA in k_filter means anchor filtering disabled (k.filter = NA).
  TransNEO_Trastuzumab = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 120,   k.score  = 200,   k.filter  = NA,
    l2.norm   = FALSE, nn.method = "annoy",
    weight.reduction = "pcaproject", k.weight = 15
  ),
  # Column-name encoded in scr/new_tras_clusters/auc_final_trastuzumab_new.py:
  #   pred_id_k_anchor_50_k_score_10_l2_norm_TRUE_k_filter_100_
  #   nn_method_annoy_k_weight_20_reduction_pcaproject
  IMPRESS_Trastuzumab = list(
    dims      = 1:20,  npcs     = 20,    reduction = "pcaproject",
    k.anchor  = 50,    k.score  = 10,    k.filter  = 100,
    l2.norm   = TRUE,  nn.method = "annoy",
    weight.reduction = "pcaproject", k.weight = 20
  )
)


# ----- Step 1: Seurat anchor-based label transfer ----------------------------
transfer_cluster_labels <- function(reference, query, params,
                                    refdata_col = "mapped_clusters") {
  anchors <- FindTransferAnchors(
    reference = reference,
    query     = query,
    dims      = params$dims,
    npcs      = params$npcs,
    reduction = params$reduction,
    k.anchor  = params$k.anchor,
    k.score   = params$k.score,
    k.filter  = params$k.filter,
    l2.norm   = params$l2.norm,
    nn.method = params$nn.method
  )
  preds <- TransferData(
    anchorset        = anchors,
    refdata          = reference@meta.data[[refdata_col]],
    weight.reduction = params$weight.reduction,
    k.weight         = params$k.weight
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
  cohort       <- ifelse(length(args) >= 4, args[4], "METABRIC")
  refdata_col  <- ifelse(length(args) >= 5, args[5], "mapped_clusters")

  if (!cohort %in% names(TRANSFER_PRESETS))
    stop(sprintf("Unknown cohort '%s'. Known: %s", cohort,
                 paste(names(TRANSFER_PRESETS), collapse = ", ")))
  params <- TRANSFER_PRESETS[[cohort]]
  cat(sprintf("Cohort: %s\n", cohort))

  cat(sprintf("Loading reference: %s\n", reference_qs))
  reference <- qread(reference_qs)
  cat(sprintf("Loading query:     %s\n", query_qs))
  query     <- qread(query_qs)

  cat("Running label transfer...\n")
  per_domain <- transfer_cluster_labels(reference, query, params,
                                        refdata_col = refdata_col)

  cat("Aggregating to per-patient × per-cluster proportions...\n")
  per_patient <- aggregate_to_patient_proportions(per_domain)

  fwrite(per_patient, out_csv)
  cat(sprintf("Wrote %d patients to %s\n", nrow(per_patient), out_csv))
}
