#!/usr/bin/env python3
"""
scripts/02_curate_metadata.py
Standardizes sample metadata and count tables across cohorts.
Ensures uniform schema: sample_id, cohort, condition, group_label, sex, tissue.
Usage:
    python scripts/02_curate_metadata.py [--demo]
"""

import os
import sys
import argparse
import logging
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Curate-Metadata")

def curate_demo_data():
    logger.info("Curating bundled demo cohorts...")
    demo_dir = "data/demo"
    proc_dir = "data/processed"
    meta_dir = "data/metadata"
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    cohorts = ["GSE107925", "GSE108045", "GSE266602"]
    for cohort in cohorts:
        counts_path = os.path.join(demo_dir, f"{cohort}_counts.csv")
        meta_path = os.path.join(demo_dir, f"{cohort}_metadata.csv")

        if not os.path.exists(counts_path) or not os.path.exists(meta_path):
            raise FileNotFoundError(f"Missing demo files for {cohort} in {demo_dir}")

        counts_df = pd.read_csv(counts_path)
        meta_df = pd.read_csv(meta_path)

        # Validate schema
        required_meta_cols = ["sample_id", "cohort", "condition", "group_label", "sex", "tissue", "sequencing_type"]
        for col in required_meta_cols:
            if col not in meta_df.columns:
                raise ValueError(f"Missing required column '{col}' in {meta_path}")

        # Check that samples in metadata match columns in counts
        sample_ids = meta_df["sample_id"].tolist()
        count_sample_cols = [c for c in counts_df.columns if c not in ["gene_id", "gene_symbol"]]
        if set(sample_ids) != set(count_sample_cols):
            raise ValueError(f"Sample ID mismatch between counts and metadata for {cohort}")

        # Verify integer counts
        for s in sample_ids:
            if not pd.api.types.is_numeric_dtype(counts_df[s]):
                raise TypeError(f"Non-numeric count column {s} in {counts_path}")
            if (counts_df[s] < 0).any():
                raise ValueError(f"Negative counts detected in column {s} in {counts_path}")

        dest_counts = os.path.join(proc_dir, f"{cohort}_counts.csv")
        dest_meta = os.path.join(meta_dir, f"{cohort}_metadata.csv")
        counts_df.to_csv(dest_counts, index=False)
        meta_df.to_csv(dest_meta, index=False)
        logger.info(f"[OK] Curated {cohort}: {counts_df.shape[0]} genes, {len(sample_ids)} samples")

    logger.info("[SUCCESS] All demo cohorts curated successfully.")

def curate_raw_geo():
    logger.info("Processing raw downloaded NCBI GEO data...")
    proc_dir = "data/processed"
    meta_dir = "data/metadata"
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)
    # Checks if raw data is downloaded, otherwise falls back to demo
    if not os.path.exists("data/raw/GSE107925/GSE107925_readCount_geneName.txt.gz"):
        logger.warning("Raw GEO files not downloaded yet. Falling back to demo data.")
        curate_demo_data()
        return

    # Process GSE107925 raw file
    logger.info("Processing GSE107925 raw counts...")
    raw_path = "data/raw/GSE107925/GSE107925_readCount_geneName.txt.gz"
    df = pd.read_csv(raw_path, sep="\t", compression="gzip")
    df = df.rename(columns={df.columns[0]: "gene_id"})
    
    # Filter adult samples (A2M)
    meta_rows = []
    adult_cols = []
    for col in df.columns:
        if col.startswith("A2M_"):
            adult_cols.append(col)
            cond = "perturbed" if "GF" in col else "reference"
            grp = "GF" if "GF" in col else "SPF"
            sex = "Male" if "Male" in col else ("Female" if "Female" in col else "Unspecified")
            meta_rows.append({
                "sample_id": col, "cohort": "GSE107925", "condition": cond,
                "group_label": grp, "sex": sex, "tissue": "Microglia", "sequencing_type": "Bulk RNA-seq"
            })

    counts_out = df[["gene_id", "gene_name"] + adult_cols].rename(columns={"gene_name": "gene_symbol"})
    counts_out.to_csv("data/processed/GSE107925_counts.csv", index=False)
    pd.DataFrame(meta_rows).to_csv("data/metadata/GSE107925_metadata.csv", index=False)
    logger.info(f"[OK] Processed GSE107925: {counts_out.shape[0]} genes, {len(adult_cols)} samples")

def main():
    parser = argparse.ArgumentParser(description="Curate sample metadata and RNA-seq counts")
    parser.add_argument("--demo", action="store_true", help="Curate bundled demo dataset")
    args = parser.parse_args()

    if args.demo:
        curate_demo_data()
    else:
        curate_raw_geo()

if __name__ == "__main__":
    main()
