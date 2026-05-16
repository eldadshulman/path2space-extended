#!/usr/bin/env Rscript
# Quick summary of a Seurat .qs object: dimensions, metadata columns,
# and head() of the first few cells. Used by notebook 04 to show the reader
# what `target` and `query` look like.

suppressPackageStartupMessages({
  library(Seurat); library(qs)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript inspect_seurat.R <object.qs> [label]")
qs_path <- args[1]
label   <- ifelse(length(args) >= 2, args[2], basename(qs_path))

obj <- qread(qs_path)
cat(sprintf("=== %s ===\n", label))
cat(sprintf("class:      %s\n", paste(class(obj), collapse = ", ")))
cat(sprintf("cells:      %d\n", ncol(obj)))
cat(sprintf("features:   %d\n", nrow(obj)))
cat(sprintf("assays:     %s\n", paste(Assays(obj), collapse = ", ")))
cat(sprintf("meta.data cols (%d):\n  %s\n", ncol(obj@meta.data),
            paste(colnames(obj@meta.data), collapse = ", ")))
cat("\nhead(@meta.data):\n")
print(head(obj@meta.data, 4))

# If a `mapped_clusters` or similar label column exists, show its distribution.
for (col in c("mapped_clusters", "predicted_cluster", "Cluster", "slide_name")) {
  if (col %in% colnames(obj@meta.data)) {
    cat(sprintf("\nDistribution of %s:\n", col))
    print(table(obj@meta.data[[col]]))
    break
  }
}
