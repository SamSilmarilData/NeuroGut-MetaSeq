#!/usr/bin/env Rscript
# scripts/04_meta_analysis.R
# Multi-cohort RNA-seq meta-analysis using R and metafor.
# Usage: Rscript scripts/04_meta_analysis.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(metafor)
})

deg_files <- list.files("results/de_results", pattern = "_deg.csv$", full.names = TRUE)
if (length(deg_files) == 0) {
  stop("No DEG results found in results/de_results/")
}

cat(sprintf("[*] Found %d DEG files for meta-analysis.\n", length(deg_files)))

# Read and combine data
all_degs <- list()
for (f in deg_files) {
  cname <- gsub("_deg.csv", "", basename(f))
  df <- read_csv(f, show_col_types = FALSE)
  df$cohort <- cname
  all_degs[[cname]] <- df
}

combined_df <- bind_rows(all_degs)
genes <- unique(combined_df$gene_symbol)

cat(sprintf("[*] Meta-analyzing %d unique genes...\n", length(genes)))

meta_rows <- list()
for (g in genes) {
  sub <- combined_df %>% filter(gene_symbol == g, !is.na(log2FoldChange), !is.na(lfcSE))
  if (nrow(sub) < 2) next

  yi <- sub$log2FoldChange
  vi <- sub$lfcSE^2

  # DerSimonian-Laird Random Effects Model via rma()
  res <- tryCatch({
    rma(yi = yi, vi = vi, method = "DL")
  }, error = function(e) NULL)

  if (is.null(res)) next

  # Fisher's combination of p-values
  pvals <- pmax(sub$pvalue, 1e-15)
  fisher_stat <- -2 * sum(log(pvals))
  p_fisher <- pchisq(fisher_stat, df = 2 * length(pvals), lower.tail = FALSE)

  meta_rows[[length(meta_rows) + 1]] <- data.frame(
    gene_symbol = g,
    n_cohorts = nrow(sub),
    meta_log2fc = as.numeric(res$b),
    meta_se = res$se,
    ci_lower = res$ci.lb,
    ci_upper = res$ci.ub,
    cochran_q = res$QE,
    i2_heterogeneity = res$I2,
    p_random_effects = res$pval,
    fisher_stat = fisher_stat,
    p_fisher = p_fisher,
    stringsAsFactors = FALSE
  )
}

res_table <- bind_rows(meta_rows)
res_table$fdr_fisher <- p.adjust(res_table$p_fisher, method = "BH")
res_table$fdr_random_effects <- p.adjust(res_table$p_random_effects, method = "BH")
res_table$significance_flag <- (res_table$fdr_fisher < 0.05) & (abs(res_table$meta_log2fc) >= 0.5)

res_table <- res_table %>% arrange(fdr_fisher)

dir.create("results/meta_results", recursive = TRUE, showWarnings = FALSE)
out_file <- "results/meta_results/microglia_meta_analysis_summary_R.csv"
write_csv(res_table, out_file)
cat(sprintf("[SUCCESS] R meta-analysis completed and saved to %s\n", out_file))
