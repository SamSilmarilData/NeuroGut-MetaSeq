#!/usr/bin/env python3
"""
scripts/03b_pydeseq2_analysis.py
Performs differential gene expression analysis on RNA-seq counts using PyDESeq2.
Outputs results to results/de_results/<cohort>_deg.csv.
Usage:
    python scripts/03b_pydeseq2_analysis.py [--demo]
"""

import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
import yaml

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("PyDESeq2-Analysis")

def run_pydeseq2(cohort: str, counts_path: str, meta_path: str, out_path: str):
    logger.info(f"Running PyDESeq2 analysis on cohort: {cohort}")
    counts_df = pd.read_csv(counts_path)
    meta_df = pd.read_csv(meta_path)

    # Gene mapping
    gene_info = counts_df[["gene_id", "gene_symbol"]].copy()
    
    # Filter count matrix to sample columns only
    sample_cols = meta_df["sample_id"].tolist()
    counts_matrix = counts_df.set_index("gene_symbol")[sample_cols].T

    # Ensure integer counts
    counts_matrix = counts_matrix.round().astype(int)

    # Prepare metadata
    meta_df = meta_df.set_index("sample_id")
    # Verify alignment
    counts_matrix = counts_matrix.loc[meta_df.index]

    logger.info(f"Design matrix for {cohort}: {counts_matrix.shape[0]} samples x {counts_matrix.shape[1]} genes")
    logger.info(f"Conditions: {meta_df['condition'].value_counts().to_dict()}")

    # Initialize PyDESeq2
    dds = DeseqDataSet(
        counts=counts_matrix,
        metadata=meta_df,
        design_factors="condition",
        refit_cooks=True,
        n_cpus=1
    )

    logger.info(f"Fitting Negative Binomial dispersions and GLMs for {cohort}...")
    dds.deseq2()

    # Statistical test: perturbed vs reference
    stat_res = DeseqStats(
        dds,
        contrast=["condition", "perturbed", "reference"],
        alpha=0.05,
        n_cpus=1
    )
    stat_res.summary()

    res_df = stat_res.results_df.copy()
    res_df = res_df.reset_index().rename(columns={"index": "gene_symbol"})

    # Merge gene_id
    res_df = pd.merge(res_df, gene_info.drop_duplicates(subset=["gene_symbol"]), on="gene_symbol", how="left")

    # Reorder columns
    cols_order = ["gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    existing_cols = [c for c in cols_order if c in res_df.columns]
    res_df = res_df[existing_cols]

    # Fill NaNs in padj with 1.0
    res_df["padj"] = res_df["padj"].fillna(1.0)
    res_df["pvalue"] = res_df["pvalue"].fillna(1.0)

    # Sort by padj
    res_df = res_df.sort_values(by="padj", ascending=True)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    res_df.to_csv(out_path, index=False)
    
    sig_up = ((res_df["padj"] < 0.05) & (res_df["log2FoldChange"] > 0.5)).sum()
    sig_down = ((res_df["padj"] < 0.05) & (res_df["log2FoldChange"] < -0.5)).sum()
    logger.info(f"[SUCCESS] {cohort} DEG completed: {sig_up} upregulated, {sig_down} downregulated (padj < 0.05, |log2FC| > 0.5)")
    logger.info(f"Saved results -> {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Run PyDESeq2 differential expression")
    parser.add_argument("--demo", action="store_true", help="Analyze demo cohorts")
    args = parser.parse_args()

    out_dir = "results/de_results"
    os.makedirs(out_dir, exist_ok=True)

    cohorts = ["GSE107925", "GSE108045", "GSE266602"]
    for cohort in cohorts:
        counts_path = f"data/processed/{cohort}_counts.csv"
        meta_path = f"data/metadata/{cohort}_metadata.csv"
        out_path = os.path.join(out_dir, f"{cohort}_deg.csv")

        if not os.path.exists(counts_path) or not os.path.exists(meta_path):
            logger.warning(f"Skipping {cohort}: Missing {counts_path} or {meta_path}")
            continue

        run_pydeseq2(cohort, counts_path, meta_path, out_path)

if __name__ == "__main__":
    main()
