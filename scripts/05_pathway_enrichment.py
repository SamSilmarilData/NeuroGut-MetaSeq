#!/usr/bin/env python3
"""
scripts/05_pathway_enrichment.py
Whole-Transcriptome Fast GSEA & Over-Representation Analysis (Horizon 4):
1. Pre-ranked GSEA on all 23,096 meta-analyzed genes ranked by:
   signed_score = sign(meta_log2fc) * (-log10(p_random_effects))
2. Evaluates:
   - MSigDB Hallmarks 2020 (data/reference/msigdb_hallmark_mouse.json)
   - KEGG Mouse 2019 (data/reference/kegg_mouse.json)
   - Curated Microglial Activation States (data/reference/microglia_phenotypes.json)
3. Over-Representation Analysis (ORA) on Core Consensus Signature & High-Heterogeneity Shock genes.
Outputs:
    results/pathways/gsea_hallmarks_summary.csv
    results/pathways/gsea_kegg_summary.csv
    results/pathways/gsea_microglia_phenotypes_summary.csv
    results/pathways/ora_consensus_pathways.csv
"""

import os
import sys
import json
import logging
import numpy as np
import pandas as pd
from scipy import stats
import gseapy as gp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Pathway-Enrichment")

def hypergeometric_ora(gene_list: list, pathway_genes: list, background_size: int = 23096):
    """Computes Fisher's exact test for gene list overlap with pathway."""
    overlap = set(gene_list).intersection(set(pathway_genes))
    k = len(overlap)
    n = len(gene_list)
    M = background_size
    N = len(pathway_genes)

    if k == 0:
        return 0, 1.0, 0.0, []

    table = [
        [k, n - k],
        [N - k, M - N - (n - k)]
    ]
    odds_ratio, p_val = stats.fisher_exact(table, alternative="greater")
    return k, p_val, odds_ratio, sorted(list(overlap))

def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Computes Benjamini-Hochberg FDR."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])
    order = np.argsort(p)
    ranked_p = p[order]
    fdr = np.zeros(n)
    cummin = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        adj = ranked_p[i] * n / rank
        cummin = min(cummin, adj)
        fdr[i] = cummin
    fdr = np.clip(fdr, 0.0, 1.0)
    rev_order = np.empty(n, dtype=int)
    rev_order[order] = np.arange(n)
    return fdr[rev_order]

def run_prerank_gsea(rnk: pd.Series, gene_sets: dict, out_name: str, out_dir: str) -> pd.DataFrame:
    """Executes fast pre-ranked GSEA on a gene set library dictionary."""
    logger.info(f"Running Pre-ranked GSEA for {out_name} ({len(gene_sets)} sets)...")
    res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        min_size=5,
        max_size=500,
        permutation_num=250,
        seed=42,
        verbose=False
    )
    res_df = res.res2d.copy()
    # Standardize column names
    res_df = res_df.rename(columns={
        "Term": "pathway",
        "ES": "enrichment_score",
        "NES": "normalized_enrichment_score",
        "NOM p-val": "nominal_p_value",
        "FDR q-val": "fdr_q_value",
        "FWER p-val": "fwer_p_value",
        "Tag %": "tag_percent",
        "Gene %": "gene_percent",
        "Lead_genes": "leading_edge_genes"
    })
    res_df = res_df.sort_values(by="fdr_q_value", ascending=True)
    out_path = os.path.join(out_dir, f"{out_name}.csv")
    res_df.to_csv(out_path, index=False)
    logger.info(f"Saved {out_name} -> {out_path} ({len(res_df)} pathways evaluated)")
    return res_df

def run_pathway_analysis(meta_path: str = "results/meta_results/microglia_meta_analysis_summary.csv",
                         core_path: str = "results/meta_results/core_consensus_signature.csv",
                         ref_dir: str = "data/reference",
                         out_dir: str = "results/pathways"):
    os.makedirs(out_dir, exist_ok=True)
    logger.info("Initializing Whole-Transcriptome Pathway & Phenotype Analysis (Horizon 4)...")

    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Missing meta-analysis summary at {meta_path}")

    meta_df = pd.read_csv(meta_path)
    logger.info(f"Loaded {len(meta_df)} meta-analyzed genes.")

    # 1. Construct pre-ranked vector: rank = sign(log2FC) * (-log10(p_RE))
    pvals = np.clip(meta_df["p_random_effects"].values, 1e-50, 1.0)
    meta_df["rank_metric"] = np.sign(meta_df["meta_log2fc"]) * (-np.log10(pvals))
    # Break ties using meta_log2fc
    meta_df = meta_df.sort_values(by=["rank_metric", "meta_log2fc"], ascending=[False, False])
    rnk_series = meta_df.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")["rank_metric"]

    # Load libraries
    with open(os.path.join(ref_dir, "msigdb_hallmark_mouse.json")) as f:
        hallmarks = json.load(f)
    with open(os.path.join(ref_dir, "kegg_mouse.json")) as f:
        kegg = json.load(f)
    with open(os.path.join(ref_dir, "microglia_phenotypes.json")) as f:
        phenotypes = json.load(f)

    # 2. Run GSEA across libraries
    gsea_hallmarks = run_prerank_gsea(rnk_series, hallmarks, "gsea_hallmarks_summary", out_dir)
    gsea_kegg = run_prerank_gsea(rnk_series, kegg, "gsea_kegg_summary", out_dir)
    gsea_pheno = run_prerank_gsea(rnk_series, phenotypes, "gsea_microglia_phenotypes_summary", out_dir)

    # 3. Over-Representation Analysis (ORA) on Core Consensus Signature & High-Heterogeneity Shock Genes
    logger.info("Running Over-Representation Analysis (ORA)...")
    core_genes = []
    if os.path.exists(core_path):
        core_df = pd.read_csv(core_path)
        core_genes = core_df["gene_symbol"].tolist()

    sig_up = meta_df[(meta_df["significance_flag"]) & (meta_df["meta_log2fc"] > 0)]["gene_symbol"].tolist()
    sig_down = meta_df[(meta_df["significance_flag"]) & (meta_df["meta_log2fc"] < 0)]["gene_symbol"].tolist()
    shock_genes = meta_df[(meta_df["i2_heterogeneity"] > 75.0) & (meta_df["p_fisher"] < 0.01)]["gene_symbol"].tolist()

    ora_records = []
    # Combined target pathways for ORA: Hallmarks + Key Microglial Phenotypes
    ora_target_sets = {**hallmarks, **phenotypes}

    for target_name, gene_set in ora_target_sets.items():
        k_core, p_core, odds_core, ov_core = hypergeometric_ora(core_genes, gene_set, len(meta_df))
        k_shock, p_shock, odds_shock, ov_shock = hypergeometric_ora(shock_genes, gene_set, len(meta_df))
        k_up, p_up, odds_up, ov_up = hypergeometric_ora(sig_up, gene_set, len(meta_df))

        if k_core > 0 or k_shock > 0 or k_up > 0:
            ora_records.append({
                "pathway": target_name,
                "pathway_size": len(gene_set),
                "core_overlap": k_core,
                "core_p_value": p_core,
                "core_odds_ratio": round(odds_core, 2),
                "core_overlap_genes": ";".join(ov_core),
                "shock_overlap": k_shock,
                "shock_p_value": p_shock,
                "shock_odds_ratio": round(odds_shock, 2),
                "shock_overlap_genes": ";".join(ov_shock),
                "up_overlap": k_up,
                "up_p_value": p_up,
                "up_overlap_genes": ";".join(ov_up)
            })

    ora_df = pd.DataFrame(ora_records)
    if not ora_df.empty:
        ora_df["core_fdr"] = benjamini_hochberg(ora_df["core_p_value"].values)
        ora_df["shock_fdr"] = benjamini_hochberg(ora_df["shock_p_value"].values)
        ora_df["up_fdr"] = benjamini_hochberg(ora_df["up_p_value"].values)
        ora_df = ora_df.sort_values(by="core_p_value", ascending=True)

    ora_path = os.path.join(out_dir, "ora_consensus_pathways.csv")
    ora_df.to_csv(ora_path, index=False)
    logger.info(f"Saved ORA summary -> {ora_path}")

    logger.info("=" * 65)
    logger.info("[SUCCESS] Pathway & Phenotype analysis completed:")
    top_hallmark = gsea_hallmarks.head(3)
    for _, r in top_hallmark.iterrows():
        logger.info(f"  * Hallmark: {r['pathway']} (NES={r['normalized_enrichment_score']:.2f}, FDR={r['fdr_q_value']:.2e})")
    for _, r in gsea_pheno.iterrows():
        logger.info(f"  * Phenotype: {r['pathway']} (NES={r['normalized_enrichment_score']:.2f}, FDR={r['fdr_q_value']:.2e})")
    logger.info("=" * 65)

if __name__ == "__main__":
    run_pathway_analysis()
