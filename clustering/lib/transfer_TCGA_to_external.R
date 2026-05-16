#!/usr/bin/env Rscript
# Project the 11 TCGA ST-cluster labels onto an external cohort's per-domain
# Seurat object via anchor-based integration (`FindTransferAnchors` +
# `TransferData`). Adapted from
# spatial_clusters/3_transfer_check_param_METABRIC.r — the parameter grid in
# that script was pruned to one tuned set here.
#
# Usage:
#   Rscript transfer_TCGA_to_external.R \
#       <reference.qs> <query.qs> <out_per_domain.csv> [<refdata_field>]
#
#   reference.qs : Seurat object with per-domain mapped TCGA clusters (label
#                  column named in `refdata_field`, default "mapped_clusters").
#   query.qs     : Seurat object for the external cohort, one "cell" per
#                  SpaGCN domain, with `slide_name` and `proportion_of_spots`
#                  in @meta.data.
#   out_per_domain.csv : output CSV — one row per query domain with the
#                  predicted TCGA cluster label.

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript transfer_TCGA_to_external.R <reference.qs> <query.qs> <out.csv> [refdata_field]")
}
reference_qs <- args[1]
query_qs     <- args[2]
out_csv      <- args[3]
refdata_col  <- ifelse(length(args) >= 4, args[4], "mapped_clusters")

cat(sprintf("Loading reference: %s\n", reference_qs))
target <- qread(reference_qs)
cat(sprintf("Loading query:     %s\n", query_qs))
query  <- qread(query_qs)

# Winning parameter combination from the METABRIC grid search in
# spatial_clusters/3_transfer_check_param_METABRIC.r.
# Recovered from clusters_reviewer_3/analysis_row_1_data.rds, which stored
# the top-ranked column name as
#   pred_id_k_anchor_40_k_score_30_l2_norm_FALSE_k_filter_50_
#   nn_method_hnsw_k_weight_15_reduction_pcaproject
params <- list(
  dims              = 1:20,    # fixed in source loop, not in column name
  npcs              = 20,      # fixed in source loop, not in column name
  reduction         = "pcaproject",
  k.anchor          = 40,
  k.score           = 30,
  k.filter          = 50,
  l2.norm           = FALSE,
  nn.method         = "hnsw",
  weight.reduction  = "pcaproject",
  k.weight          = 15
)

cat("Finding transfer anchors...\n")
anchors <- FindTransferAnchors(
  reference = target,
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

cat("Transferring labels...\n")
preds <- TransferData(
  anchorset       = anchors,
  refdata         = target@meta.data[[refdata_col]],
  weight.reduction = params$weight.reduction,
  k.weight         = params$k.weight
)

# Build per-domain output: one row per query domain.
md <- query@meta.data
out <- data.table(
  domain_id           = rownames(md),
  slide_name          = md$slide_name,
  proportion_of_spots = if ("proportion_of_spots" %in% colnames(md)) md$proportion_of_spots else NA_real_,
  predicted_cluster   = preds$predicted.id,
  prediction_score    = preds$prediction.score.max
)

fwrite(out, file = out_csv)
cat(sprintf("Wrote %d rows to %s\n", nrow(out), out_csv))
