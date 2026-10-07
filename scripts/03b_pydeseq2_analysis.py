#!/usr/bin/env python3
"""
scripts/03b_pydeseq2_analysis.py
Performs full-transcriptome differential expression analysis using PyDESeq2 Negative Binomial GLMs.
Supports multi-factor design (~ sex + condition) to partition baseline sex variance,
pre-filters low-count genes (>=10 counts across samples), and applies empirical Bayes shrinkage.
Outputs full statistical tables to results/de_results/<cohort>_deg.csv.

Usage:
    python scripts/03b_pydeseq2_analysis.py [--demo] [--cohort COHORT_ID]
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("PyDESeq2-Analysis")

def run_pydeseq2(cohort: str, counts_path: str, meta_path: str, out_path: str, min_counts: int = 10, n_cpus: int = 2):
    logger.info(f"=== Starting PyDESeq2 GLM Analysis: {cohort} ===")
    if not os.path.exists(counts_path) or not os.path.exists(meta_path):
        raise FileNotFoundError(f"Missing input files for {cohort}: {counts_path} or {meta_path}")

    counts_df = pd.read_csv(counts_path)
    meta_df = pd.read_csv(meta_path)

    # Clean string identifiers
    counts_df.columns = counts_df.columns.astype(str)
    meta_df["sample_id"] = meta_df["sample_id"].astype(str)
    sample_cols = meta_df["sample_id"].tolist()

    # Gene mapping
    gene_info = counts_df[["gene_id", "gene_symbol"]].copy()

    # Matrix: samples as rows, genes as columns
    counts_matrix = counts_df.set_index("gene_symbol")[sample_cols].T
    counts_matrix = counts_matrix.round().astype(int)

    # Filter lowly expressed genes
    gene_totals = counts_matrix.sum(axis=0)
    keep_genes = gene_totals >= min_counts
    c_filt = counts_matrix.loc[:, keep_genes].copy()
    logger.info(f"Filtered {counts_matrix.shape[1]:,} genes -> {c_filt.shape[1]:,} expressed genes (>={min_counts} total counts)")

    # Metadata alignment
    meta_df = meta_df.set_index("sample_id").loc[c_filt.index]

    # Check whether sex varies to formulate multi-factor design
    has_sex_var = meta_df["sex"].nunique() > 1
    design = "~ sex + condition" if has_sex_var else "~ condition"
    logger.info(f"Design formula for {cohort}: {design}")
    logger.info(f"Samples ({len(meta_df)}): {meta_df['condition'].value_counts().to_dict()}, Sex: {meta_df['sex'].value_counts().to_dict()}")

    # Initialize PyDESeq2
    dds = DeseqDataSet(
        counts=c_filt,
        metadata=meta_df,
        design=design,
        refit_cooks=True,
        n_cpus=n_cpus
    )

    logger.info(f"Fitting Negative Binomial dispersions and GLMs for {cohort}...")
    dds.deseq2()

    # Statistical test: perturbed vs reference
    logger.info("Executing Wald hypothesis tests (perturbed vs. reference)...")
    stat_res = DeseqStats(
        dds,
        contrast=["condition", "perturbed", "reference"],
        alpha=0.05,
        n_cpus=n_cpus
    )
    stat_res.summary()

    res_df = stat_res.results_df.copy()
    res_df = res_df.reset_index().rename(columns={"index": "gene_symbol"})

    # Merge gene_id
    res_df = pd.merge(res_df, gene_info.drop_duplicates(subset=["gene_symbol"]), on="gene_symbol", how="left")

    # Standardize column schema
    cols_order = ["gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    existing_cols = [c for c in cols_order if c in res_df.columns]
    res_df = res_df[existing_cols]

    # Handle uncalculated p-values / padj (outliers or zero counts)
    res_df["padj"] = res_df["padj"].fillna(1.0)
    res_df["pvalue"] = res_df["pvalue"].fillna(1.0)

    # Sort by significance (padj ascending, then absolute log2FC descending)
    res_df["abs_lfc"] = res_df["log2FoldChange"].abs()
    res_df = res_df.sort_values(by=["padj", "abs_lfc"], ascending=[True, False]).drop(columns=["abs_lfc"])

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    res_df.to_csv(out_path, index=False)

    sig_up = ((res_df["padj"] < 0.05) & (res_df["log2FoldChange"] >= 0.5)).sum()
    sig_down = ((res_df["padj"] < 0.05) & (res_df["log2FoldChange"] <= -0.5)).sum()
    total_sig = sig_up + sig_down

    logger.info(f"[SUCCESS] {cohort}: {total_sig} significant DEGs (padj < 0.05, |log2FC| >= 0.5)")
    logger.info(f"  -> Upregulated:   {sig_up}")
    logger.info(f"  -> Downregulated: {sig_down}")
    logger.info(f"Saved: {out_path} ({len(res_df):,} genes)")

    # Print top 5 landmark DEGs
    top5 = res_df.head(5)[["gene_symbol", "baseMean", "log2FoldChange", "padj"]]
    logger.info(f"Top 5 DEGs in {cohort}:\n{top5.to_string(index=False)}")
    return res_df

def main():
    parser = argparse.ArgumentParser(description="Run PyDESeq2 differential expression")
    parser.add_argument("--demo", action="store_true", help="Analyze demo cohorts")
    parser.add_argument("--cohort", type=str, default=None, help="Specific cohort to analyze")
    parser.add_argument("--cpus", type=int, default=2, help="Number of CPU cores")
    args = parser.parse_args()

    out_dir = "results/de_results"
    os.makedirs(out_dir, exist_ok=True)

    if args.demo:
        cohorts = ["GSE107925", "GSE108045", "GSE266602"]
        data_dir = "data/demo"
        meta_dir = "data/demo"
    else:
        cohorts = ["GSE107925", "GSE108045", "GSE266602", "GSE186210"]
        data_dir = "data/processed"
        meta_dir = "data/metadata"

    if args.cohort:
        if args.cohort not in cohorts:
            logger.error(f"Unknown cohort: {args.cohort}")
            sys.exit(1)
        cohorts = [args.cohort]

    for cid in cohorts:
        counts_path = f"{data_dir}/{cid}_counts.csv"
        meta_path = f"{meta_dir}/{cid}_metadata.csv"
        out_path = os.path.join(out_dir, f"{cid}_deg.csv")
        try:
            run_pydeseq2(cid, counts_path, meta_path, out_path, n_cpus=args.cpus)
        except Exception as e:
            logger.error(f"Failed PyDESeq2 analysis on {cid}: {e}")
            raise

    logger.info("[SUCCESS] Full PyDESeq2 differential expression completed across all cohorts.")

if __name__ == "__main__":
    main()
