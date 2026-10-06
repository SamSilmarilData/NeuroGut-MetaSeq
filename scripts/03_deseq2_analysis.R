#!/usr/bin/env Rscript
# scripts/03_deseq2_analysis.R
# Differential gene expression analysis using Bioconductor DESeq2.
# Usage: Rscript scripts/03_deseq2_analysis.R [--cohort GSE107925]

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
cohorts <- if (length(args) > 0) args else c("GSE107925", "GSE108045", "GSE266602")

for (cohort in cohorts) {
  counts_file <- sprintf("data/processed/%s_counts.csv", cohort)
  meta_file <- sprintf("data/metadata/%s_metadata.csv", cohort)
  out_file <- sprintf("results/de_results/%s_deg.csv", cohort)

  if (!file.exists(counts_file) || !file.exists(meta_file)) {
    cat(sprintf("[SKIP] Files for %s not found. Skipping.\n", cohort))
    next
  }

  cat(sprintf("[*] Running DESeq2 on %s...\n", cohort))
  counts_df <- read_csv(counts_file, show_col_types = FALSE)
  meta_df <- read_csv(meta_file, show_col_types = FALSE)

  # Prepare count matrix
  gene_symbols <- counts_df$gene_symbol
  gene_ids <- counts_df$gene_id
  sample_names <- meta_df$sample_id

  count_mat <- as.matrix(counts_df[, sample_names])
  rownames(count_mat) <- gene_symbols

  # Metadata
  coldata <- as.data.frame(meta_df)
  rownames(coldata) <- coldata$sample_id
  coldata$condition <- factor(coldata$condition, levels = c("reference", "perturbed"))

  # Construct DESeqDataSet
  design_formula <- if ("sex" %in% colnames(coldata) && length(unique(coldata$sex)) > 1) {
    coldata$sex <- factor(coldata$sex)
    ~ sex + condition
  } else {
    ~ condition
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = coldata,
    design = design_formula
  )

  # Filter low counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  # Run DESeq
  dds <- DESeq(dds, quiet = TRUE)

  # Extract results
  res <- results(dds, contrast = c("condition", "perturbed", "reference"))
  res_df <- as.data.frame(res)
  res_df$gene_symbol <- rownames(res_df)
  
  # Map back gene_id
  gene_map <- data.frame(gene_id = gene_ids, gene_symbol = gene_symbols)
  gene_map <- distinct(gene_map, gene_symbol, .keep_all = TRUE)
  res_df <- left_join(res_df, gene_map, by = "gene_symbol")

  res_df <- res_df %>%
    select(gene_id, gene_symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    arrange(padj)

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  write_csv(res_df, out_file)
  cat(sprintf("[SUCCESS] Saved %s DESeq2 results to %s\n", cohort, out_file))
}
