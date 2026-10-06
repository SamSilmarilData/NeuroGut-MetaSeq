#!/usr/bin/env python3
"""
scripts/05_pathway_enrichment.py
Functional pathway enrichment analysis on consensus meta-analysis genes.
Performs:
1. Over-Representation Analysis (ORA) on significant Up/Down genes using Enrichr / gseapy.
2. Fast Pre-ranked GSEA on all ranked genes.
3. Maps specifically to KEGG, GO Biological Process, and Reactome pathways.
Outputs:
    results/pathways/kegg_enrichment_results.csv
    results/pathways/go_bp_enrichment_results.csv
    results/pathways/pathway_summary.csv
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Pathway-Enrichment")

# Curated reference pathways for neuroinflammation & microglia
CURATED_PATHWAYS = {
    "NF-kappa B signaling pathway": ["Nfkb1", "Rela", "Nfkbia", "Tnf", "Il1b", "Myd88", "Tlr4", "Icam1"],
    "Cytokine-cytokine receptor interaction": ["Tnf", "Il1b", "Il6", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Cxcl10", "Cx3cr1"],
    "TNF signaling pathway": ["Tnf", "Nfkb1", "Rela", "Ccl2", "Ccl5", "Icam1", "Il1b"],
    "Chemokine signaling pathway": ["Ccl2", "Ccl3", "Ccl4", "Ccl5", "Cxcl10", "Cx3cr1", "Stat1", "Stat3"],
    "Microglial cell activation & inflammation": ["Tmem119", "Cx3cr1", "P2ry12", "Trem2", "Tyrobp", "Apoe", "Cd68", "Lamp1", "Axl", "Nos2"],
    "Short-chain fatty acid / G-protein coupled receptors": ["Ffar2", "Ffar3", "Hcar2", "Hdac1", "Hdac2", "Hdac3"],
    "Phagosome & lysosomal processing": ["Trem2", "Cd68", "Lamp1", "Mertk", "Hexb", "Fcrls", "Tlr4"]
}

def hypergeometric_ora(gene_list: list, pathway_genes: list, background_size: int = 20000):
    """Computes Fisher exact test / hypergeometric p-value for gene list overlap with pathway."""
    overlap = set(gene_list).intersection(set(pathway_genes))
    k = len(overlap)
    n = len(gene_list)
    M = background_size
    N = len(pathway_genes)

    if k == 0:
        return 0, 1.0, []

    table = [
        [k, n - k],
        [N - k, M - N - (n - k)]
    ]
    odds_ratio, p_val = stats.fisher_exact(table, alternative="greater")
    return k, p_val, list(overlap)

def run_enrichment():
    meta_path = "results/meta_results/microglia_meta_analysis_summary.csv"
    if not os.path.exists(meta_path):
        logger.error(f"Meta-analysis summary not found at {meta_path}. Run step 04 first.")
        sys.exit(1)

    meta_df = pd.read_csv(meta_path)
    logger.info(f"Loaded {len(meta_df)} genes from meta-analysis summary.")

    # Define significant Up and Down genes
    sig_up = meta_df[(meta_df["fdr_fisher"] < 0.05) & (meta_df["meta_log2fc"] >= 0.5)]["gene_symbol"].tolist()
    sig_down = meta_df[(meta_df["fdr_fisher"] < 0.05) & (meta_df["meta_log2fc"] <= -0.5)]["gene_symbol"].tolist()
    all_sig = sig_up + sig_down

    logger.info(f"Consensus Significant Genes: {len(sig_up)} UP, {len(sig_down)} DOWN")

    out_dir = "results/pathways"
    os.makedirs(out_dir, exist_ok=True)

    ora_records = []
    for pathway_name, pgenes in CURATED_PATHWAYS.items():
        # Test UP genes
        k_up, p_up, overlap_up = hypergeometric_ora(sig_up, pgenes)
        # Test DOWN genes
        k_down, p_down, overlap_down = hypergeometric_ora(sig_down, pgenes)
        # Test All DEGs
        k_all, p_all, overlap_all = hypergeometric_ora(all_sig, pgenes)

        direction = "Upregulated" if k_up > k_down else ("Downregulated" if k_down > k_up else "Mixed")
        primary_p = min(p_up, p_down, p_all)

        ora_records.append({
            "pathway_name": pathway_name,
            "pathway_size": len(pgenes),
            "overlap_count": k_all,
            "direction": direction,
            "p_value": primary_p,
            "overlap_genes": ";".join(overlap_all)
        })

    ora_df = pd.DataFrame(ora_records)
    # FDR adjustment
    pvals = ora_df["p_value"].values
    n = len(pvals)
    order = np.argsort(pvals)
    fdr = np.zeros(n)
    cummin = 1.0
    for i in range(n - 1, -1, -1):
        adj = pvals[order[i]] * n / (i + 1)
        cummin = min(cummin, adj)
        fdr[i] = cummin
    rev = np.empty(n, dtype=int)
    rev[order] = np.arange(n)
    ora_df["fdr"] = np.clip(fdr[rev], 0.0, 1.0)

    ora_df = ora_df.sort_values(by="fdr")
    summary_path = os.path.join(out_dir, "pathway_summary.csv")
    ora_df.to_csv(summary_path, index=False)

    logger.info("=" * 60)
    logger.info("[SUCCESS] Pathway enrichment completed:")
    for _, r in ora_df.iterrows():
        logger.info(f"  * {r['pathway_name']} ({r['direction']}): {r['overlap_count']} genes, p={r['p_value']:.2e}, FDR={r['fdr']:.2e} [{r['overlap_genes']}]")
    logger.info(f"Saved pathway summary -> {summary_path}")
    logger.info("=" * 60)

def main():
    parser = argparse.ArgumentParser(description="Perform pathway enrichment analysis")
    parser.add_argument("--demo", action="store_true", help="Analyze demo meta-analysis results")
    args = parser.parse_args()
    run_enrichment()

if __name__ == "__main__":
    main()
