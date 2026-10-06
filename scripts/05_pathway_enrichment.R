#!/usr/bin/env Rscript
# scripts/05_pathway_enrichment.R
# Biological pathway analysis using clusterProfiler and org.Mm.eg.db.
# Usage: Rscript scripts/05_pathway_enrichment.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

meta_file <- "results/meta_results/microglia_meta_analysis_summary.csv"
if (!file.exists(meta_file)) {
  stop("Meta-analysis results not found. Please run meta-analysis first.")
}

meta_df <- read_csv(meta_file, show_col_types = FALSE)
sig_genes <- meta_df %>% filter(fdr_fisher < 0.05, abs(meta_log2fc) >= 0.5)

cat(sprintf("[*] Running clusterProfiler on %d consensus significant genes...\n", nrow(sig_genes)))

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(sig_genes$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO Enrichment (Biological Process)
go_res <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# KEGG Enrichment
kegg_res <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

dir.create("results/pathways", recursive = TRUE, showWarnings = FALSE)
if (!is.null(go_res)) {
  write_csv(as.data.frame(go_res), "results/pathways/go_bp_enrichment_R.csv")
}
if (!is.null(kegg_res)) {
  write_csv(as.data.frame(kegg_res), "results/pathways/kegg_enrichment_R.csv")
}

cat("[SUCCESS] R pathway enrichment completed and saved to results/pathways/\n")
