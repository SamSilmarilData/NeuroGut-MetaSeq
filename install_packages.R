#!/usr/bin/env Rscript
# NeuroGut-MetaSeq: Automated Bioconductor and CRAN Package Installer

cat("========================================================\n")
cat("NeuroGut-MetaSeq: Installing R Dependencies\n")
cat("========================================================\n")

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("[+] Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran_packages <- c(
  "yaml", "dplyr", "readr", "tidyr", "ggplot2", 
  "pheatmap", "metafor", "metap", "optparse"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[+] Installing CRAN package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("[OK] CRAN package already installed: %s\n", pkg))
  }
}

bioc_packages <- c(
  "DESeq2", "apeglm", "clusterProfiler", "org.Mm.eg.db", "EnhancedVolcano"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[+] Installing Bioconductor package: %s\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat(sprintf("[OK] Bioconductor package already installed: %s\n", pkg))
  }
}

cat("\n[SUCCESS] All R packages verified successfully.\n")
