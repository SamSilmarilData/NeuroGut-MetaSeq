#!/usr/bin/env python3
"""
scripts/02_curate_metadata.py
Standardizes sample metadata and RNA-seq count tables across all cohorts.
Ensures uniform schema:
  - Metadata: sample_id, cohort, condition, group_label, sex, tissue, sequencing_type
  - Counts: gene_id, gene_symbol, <sample_1>, <sample_2>, ...
Supports:
  - Bundled demo cohorts (--demo)
  - Full-scale raw NCBI GEO datasets: GSE107925, GSE108045, GSE266602, GSE186210

Usage:
    python scripts/02_curate_metadata.py [--demo] [--cohort COHORT_ID]
"""

import os
import sys
import gzip
import argparse
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Curate-Metadata")

def load_reference_mapping(ref_path: str = "data/reference/mouse_ensembl_to_symbol.tsv.gz"):
    """Load cached Ensembl ID <-> MGI Symbol mapping."""
    if not os.path.exists(ref_path):
        logger.warning(f"Reference mapping not found at {ref_path}. Ensembl IDs will serve as fallback.")
        return {}, {}
    df_ref = pd.read_csv(ref_path, sep="\t")
    ens2sym = dict(zip(df_ref["Gene stable ID"].dropna(), df_ref["Gene name"].dropna()))
    sym2ens = dict(zip(df_ref["Gene name"].dropna(), df_ref["Gene stable ID"].dropna()))
    return ens2sym, sym2ens

def validate_and_save(counts_df: pd.DataFrame, meta_df: pd.DataFrame, cohort_id: str, proc_dir: str, meta_dir: str):
    """Ensure strict schema adherence and write to CSV."""
    required_meta_cols = ["sample_id", "cohort", "condition", "group_label", "sex", "tissue", "sequencing_type"]
    for col in required_meta_cols:
        if col not in meta_df.columns:
            raise ValueError(f"Missing required column '{col}' in metadata for {cohort_id}")

    sample_ids = meta_df["sample_id"].tolist()
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"Duplicate sample IDs detected in metadata for {cohort_id}")

    if not set(meta_df["condition"]).issubset({"reference", "perturbed"}):
        invalid_conds = set(meta_df["condition"]) - {"reference", "perturbed"}
        raise ValueError(f"Invalid condition values {invalid_conds} in {cohort_id}")

    count_samples = [c for c in counts_df.columns if c not in ["gene_id", "gene_symbol"]]
    if set(sample_ids) != set(count_samples):
        mismatch_counts = set(count_samples) - set(sample_ids)
        mismatch_meta = set(sample_ids) - set(count_samples)
        raise ValueError(f"Sample mismatch in {cohort_id}: extra in counts: {mismatch_counts}, extra in meta: {mismatch_meta}")

    # Reorder count columns to match metadata sample_id order
    counts_df = counts_df[["gene_id", "gene_symbol"] + sample_ids]

    # Validate non-negative integers
    for s in sample_ids:
        col_vals = counts_df[s].values
        if np.isnan(col_vals).any():
            raise ValueError(f"NaN counts detected in sample {s} for {cohort_id}")
        if (col_vals < 0).any():
            raise ValueError(f"Negative counts detected in sample {s} for {cohort_id}")
        counts_df[s] = counts_df[s].round().astype(int)

    dest_counts = os.path.join(proc_dir, f"{cohort_id}_counts.csv")
    dest_meta = os.path.join(meta_dir, f"{cohort_id}_metadata.csv")
    counts_df.to_csv(dest_counts, index=False)
    meta_df.to_csv(dest_meta, index=False)
    logger.info(f"[OK] Curated {cohort_id}: {counts_df.shape[0]:,} genes, {len(sample_ids)} samples -> {dest_counts}")

def curate_demo_data():
    """Curates bundled demo cohorts."""
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
        validate_and_save(counts_df, meta_df, cohort, proc_dir, meta_dir)

    logger.info("[SUCCESS] All demo cohorts curated successfully.")

def curate_gse107925(proc_dir: str, meta_dir: str):
    """Curate GSE107925: Adult SPF vs. GF microglia."""
    raw_path = "data/raw/GSE107925/GSE107925_readCount_geneName.txt.gz"
    if not os.path.exists(raw_path):
        raise FileNotFoundError(f"Missing raw file: {raw_path}")

    logger.info(f"Curating GSE107925 from {raw_path}...")
    df = pd.read_csv(raw_path, sep="\t", compression="gzip")
    df = df.rename(columns={df.columns[0]: "gene_id", "gene_name": "gene_symbol"})

    # Filter adult samples (A2M)
    adult_cols = [c for c in df.columns if c.startswith("A2M_")]
    meta_rows = []
    for col in adult_cols:
        cond = "perturbed" if "GF" in col else "reference"
        grp = "GF" if "GF" in col else "SPF"
        sex = "Male" if "Male" in col else ("Female" if "Female" in col else "Unspecified")
        meta_rows.append({
            "sample_id": col, "cohort": "GSE107925", "condition": cond,
            "group_label": grp, "sex": sex, "tissue": "Microglia", "sequencing_type": "Bulk RNA-seq"
        })

    # Group duplicate gene symbols by taking sum of counts
    count_data = df[["gene_symbol"] + adult_cols].groupby("gene_symbol").sum().reset_index()
    # Map first Ensembl ID for each symbol
    gene_id_map = df.groupby("gene_symbol")["gene_id"].first().to_dict()
    count_data.insert(0, "gene_id", count_data["gene_symbol"].map(gene_id_map))

    meta_df = pd.DataFrame(meta_rows)
    validate_and_save(count_data, meta_df, "GSE107925", proc_dir, meta_dir)

def curate_gse108045(proc_dir: str, meta_dir: str):
    """Curate GSE108045: Adult CTR vs. ABX microglia."""
    raw_path = "data/raw/GSE108045/GSE108045_readCount_geneName_exp2.txt.gz"
    if not os.path.exists(raw_path):
        raise FileNotFoundError(f"Missing raw file: {raw_path}")

    logger.info(f"Curating GSE108045 from {raw_path}...")
    df = pd.read_csv(raw_path, sep="\t", compression="gzip")
    df = df.rename(columns={df.columns[0]: "gene_id", "gene_name": "gene_symbol"})

    sample_cols = [c for c in df.columns if c.startswith("2_CTR") or c.startswith("2_ABX")]
    meta_rows = []
    for col in sample_cols:
        cond = "reference" if "CTR" in col else "perturbed"
        grp = "CTR" if "CTR" in col else "ABX"
        # Standardize sample id string (clean spaces for safety)
        clean_id = col.replace(" ", "_")
        sex = "Female" if " F" in col or "_F" in clean_id else "Male"
        meta_rows.append({
            "sample_id": clean_id, "cohort": "GSE108045", "condition": cond,
            "group_label": grp, "sex": sex, "tissue": "Microglia", "sequencing_type": "Bulk RNA-seq"
        })

    rename_map = {col: col.replace(" ", "_") for col in sample_cols}
    df = df.rename(columns=rename_map)
    clean_sample_cols = [rename_map[c] for c in sample_cols]

    count_data = df[["gene_symbol"] + clean_sample_cols].groupby("gene_symbol").sum().reset_index()
    gene_id_map = df.groupby("gene_symbol")["gene_id"].first().to_dict()
    count_data.insert(0, "gene_id", count_data["gene_symbol"].map(gene_id_map))

    meta_df = pd.DataFrame(meta_rows)
    validate_and_save(count_data, meta_df, "GSE108045", proc_dir, meta_dir)

def curate_gse266602(proc_dir: str, meta_dir: str, ens2sym: dict):
    """Curate GSE266602: Adult sham baseline microglia (SPF vs. GF and ABX)."""
    raw_path = "data/raw/GSE266602/GSE266602_gene_count_matrix.txt.gz"
    if not os.path.exists(raw_path):
        raise FileNotFoundError(f"Missing raw file: {raw_path}")

    logger.info(f"Curating GSE266602 from {raw_path}...")
    df = pd.read_csv(raw_path, sep="\t", compression="gzip")
    df["gene_id_clean"] = df["gene_id"].str.split(".").str[0]
    df["gene_symbol"] = df["gene_id_clean"].map(ens2sym).fillna(df["gene_id"])

    # Extract sham baseline samples (SPF_Sham, GF_Sham, ABX_Sham)
    sham_cols = [c for c in df.columns if c.startswith("SPF_Sham") or c.startswith("GF_Sham") or c.startswith("ABX_Sham")]
    meta_rows = []
    for col in sham_cols:
        cond = "reference" if col.startswith("SPF_Sham") else "perturbed"
        grp = col.split("Sham")[0] + "Sham"
        meta_rows.append({
            "sample_id": col, "cohort": "GSE266602", "condition": cond,
            "group_label": grp, "sex": "Male", "tissue": "Microglia", "sequencing_type": "Bulk RNA-seq"
        })

    count_data = df[["gene_symbol"] + sham_cols].groupby("gene_symbol").sum().reset_index()
    gene_id_map = df.groupby("gene_symbol")["gene_id"].first().to_dict()
    count_data.insert(0, "gene_id", count_data["gene_symbol"].map(gene_id_map))

    meta_df = pd.DataFrame(meta_rows)
    validate_and_save(count_data, meta_df, "GSE266602", proc_dir, meta_dir)

def curate_gse186210(proc_dir: str, meta_dir: str, sym2ens: dict):
    """Curate GSE186210: Adult dietary fiber perturbation (Normal Fiber vs. Zero Fiber, WT)."""
    raw_path = "data/raw/GSE186210/GSE186210_table.tsv.gz"
    series_path = "data/raw/GSE186210/GSE186210_series_matrix.txt.gz"
    if not os.path.exists(raw_path) or not os.path.exists(series_path):
        raise FileNotFoundError(f"Missing raw files for GSE186210: {raw_path} or {series_path}")

    logger.info(f"Curating GSE186210 from {raw_path} and {series_path}...")
    # Parse series matrix to extract metadata
    titles, accessions, chars = [], [], []
    with gzip.open(series_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"\n ') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accessions = [x.strip('"\n ') for x in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip('"\n ') for x in line.split("\t")[1:]])

    meta_dict = {"sample_title": titles, "geo_accession": accessions}
    for row in chars:
        if not row:
            continue
        first = row[0]
        if ":" in first:
            key = first.split(":")[0].strip().lower()
            meta_dict[key] = [r.split(":")[1].strip() if ":" in r else r for r in row]

    s_meta = pd.DataFrame(meta_dict)
    s_meta["sample_id"] = s_meta["sample_title"].str.replace("Microglia, ", "", regex=False).str.strip()

    # Filter to WT mice for core dietary contrast (Normal Fiber vs. Zero Fiber)
    # Diet values: 'Control', 'zero fiber'
    wt_meta = s_meta[(s_meta["genotype"] == "WT") & (s_meta["diet"].isin(["Control", "zero fiber"]))].copy()

    meta_rows = []
    for _, row in wt_meta.iterrows():
        sid = f"S_{row['sample_id']}"
        diet = row["diet"]
        cond = "reference" if diet == "Control" else "perturbed"
        grp = "Fiber_Normal" if diet == "Control" else "Zero_Fiber"
        sex = "Female" if row.get("sex", "").upper() == "F" else "Male"
        meta_rows.append({
            "sample_id": sid, "cohort": "GSE186210", "condition": cond,
            "group_label": grp, "sex": sex, "tissue": "Microglia", "sequencing_type": "Bulk RNA-seq"
        })

    meta_df = pd.DataFrame(meta_rows)
    selected_raw_ids = wt_meta["sample_id"].tolist()

    df_counts = pd.read_csv(raw_path, sep="\t", compression="gzip")
    df_counts = df_counts.rename(columns={"Geneid": "gene_symbol"})
    df_counts.columns = df_counts.columns.astype(str)

    # Ensure selected sample columns exist in count table
    avail_cols = [c for c in selected_raw_ids if c in df_counts.columns]
    if len(avail_cols) != len(selected_raw_ids):
        missing = set(selected_raw_ids) - set(avail_cols)
        raise ValueError(f"Sample columns missing in GSE186210 count table: {missing}")

    rename_counts = {c: f"S_{c}" for c in avail_cols}
    df_counts = df_counts.rename(columns=rename_counts)
    clean_sample_cols = [rename_counts[c] for c in avail_cols]

    count_data = df_counts[["gene_symbol"] + clean_sample_cols].groupby("gene_symbol").sum().reset_index()
    count_data.insert(0, "gene_id", count_data["gene_symbol"].map(sym2ens).fillna(count_data["gene_symbol"]))

    validate_and_save(count_data, meta_df, "GSE186210", proc_dir, meta_dir)

def curate_raw_geo(cohort: str = None):
    """Processes real downloaded NCBI GEO data for all 4 cohorts."""
    logger.info("Processing raw downloaded NCBI GEO datasets across cohorts...")
    proc_dir = "data/processed"
    meta_dir = "data/metadata"
    os.makedirs(proc_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    ens2sym, sym2ens = load_reference_mapping()

    cohort_runners = {
        "GSE107925": lambda: curate_gse107925(proc_dir, meta_dir),
        "GSE108045": lambda: curate_gse108045(proc_dir, meta_dir),
        "GSE266602": lambda: curate_gse266602(proc_dir, meta_dir, ens2sym),
        "GSE186210": lambda: curate_gse186210(proc_dir, meta_dir, sym2ens)
    }

    targets = [cohort] if cohort else list(cohort_runners.keys())
    for cid in targets:
        if cid in cohort_runners:
            try:
                cohort_runners[cid]()
            except Exception as e:
                logger.error(f"Failed to curate {cid}: {e}")
                raise
        else:
            logger.error(f"Unknown cohort ID: {cid}")

    logger.info("[SUCCESS] All targeted raw GEO cohorts curated into data/processed and data/metadata.")

def main():
    parser = argparse.ArgumentParser(description="Curate sample metadata and RNA-seq counts")
    parser.add_argument("--demo", action="store_true", help="Curate bundled demo dataset")
    parser.add_argument("--cohort", type=str, default=None, help="Process specific cohort only")
    args = parser.parse_args()

    if args.demo:
        curate_demo_data()
    else:
        curate_raw_geo(cohort=args.cohort)

if __name__ == "__main__":
    main()
