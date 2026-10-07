#!/usr/bin/env python3
"""
scripts/05b_tf_regulon_analysis.py
Upstream Transcription Factor (TF) Regulon Deconvolution (Horizon 4):
1. Ingests consensus meta-analysis summary statistics.
2. Ingests TRRUST Mouse TF regulon database (data/reference/trrust_mouse.json).
3. Quantifies Regulon Activity Z-scores and two-sample test statistics for all TFs with >= 5 targets.
4. Computes Fisher's exact test for target over-representation.
5. Applies Benjamini-Hochberg FDR correction.
Outputs:
    results/pathways/tf_regulon_activity_summary.csv
"""

import os
import json
import logging
import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("TF-Regulon")

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

def run_tf_regulon_analysis(meta_path: str = "results/meta_results/microglia_meta_analysis_summary.csv",
                            trrust_path: str = "data/reference/trrust_mouse.json",
                            out_dir: str = "results/pathways"):
    os.makedirs(out_dir, exist_ok=True)
    logger.info("Initializing Upstream TF Regulon Deconvolution (Horizon 4)...")

    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Missing meta-analysis summary at {meta_path}")
    if not os.path.exists(trrust_path):
        raise FileNotFoundError(f"Missing TRRUST library at {trrust_path}")

    meta_df = pd.read_csv(meta_path)
    logger.info(f"Loaded {len(meta_df)} meta-analyzed genes.")

    with open(trrust_path) as f:
        trrust_raw = json.load(f)

    # Build mapping from uppercase gene symbols to dataset index
    expressed_upper = {g.upper(): g for g in meta_df["gene_symbol"].values}
    gene_to_lfc = meta_df.set_index("gene_symbol")["meta_log2fc"].to_dict()
    gene_to_p = meta_df.set_index("gene_symbol")["p_random_effects"].to_dict()
    all_lfcs = np.array(list(gene_to_lfc.values()))
    bg_mean = np.mean(all_lfcs)
    bg_var = np.var(all_lfcs)

    # Clean and consolidate TRRUST TF keys
    tf_regulons = {}
    for k, targets in trrust_raw.items():
        # Standardize TF symbol to capitalized name
        tf_name = k.split()[0].capitalize()
        matched_targets = [expressed_upper[t] for t in targets if t.upper() in expressed_upper]
        if len(matched_targets) >= 5:
            if tf_name not in tf_regulons or len(matched_targets) > len(tf_regulons[tf_name]):
                tf_regulons[tf_name] = matched_targets

    logger.info(f"Identified {len(tf_regulons)} TFs with >= 5 measured target genes in transcriptome.")

    # High-confidence hits for Fisher overlap test
    sig_hits = set(meta_df[meta_df["significance_flag"]]["gene_symbol"].values)
    total_genes = len(meta_df)

    tf_records = []
    for tf, targets in tf_regulons.items():
        k_targets = len(targets)
        target_lfcs = np.array([gene_to_lfc[g] for g in targets])
        bg_other_lfcs = np.array([gene_to_lfc[g] for g in gene_to_lfc if g not in set(targets)])

        # Target metrics
        mean_target_lfc = float(np.mean(target_lfcs))
        median_target_lfc = float(np.median(target_lfcs))

        # Standardized Activity Z-score: (mean_target - bg_mean) / SE_bg
        se_target = float(np.sqrt(bg_var / k_targets))
        z_activity = float((mean_target_lfc - bg_mean) / se_target) if se_target > 0 else 0.0

        # Two-sample Welch t-test
        t_stat, p_welch = stats.ttest_ind(target_lfcs, bg_other_lfcs, equal_var=False)

        # Mann-Whitney U test
        u_stat, p_mwu = stats.mannwhitneyu(target_lfcs, bg_other_lfcs, alternative="two-sided")

        # Kolmogorov-Smirnov test (distribution shift)
        ks_stat, p_ks = stats.ks_2samp(target_lfcs, bg_other_lfcs)

        # Fisher's exact test for overlap with significant hits
        overlap_sig = set(targets).intersection(sig_hits)
        n_overlap = len(overlap_sig)
        table = [
            [n_overlap, k_targets - n_overlap],
            [len(sig_hits) - n_overlap, total_genes - len(sig_hits) - (k_targets - n_overlap)]
        ]
        odds_ratio, p_fisher = stats.fisher_exact(table, alternative="greater")

        tf_records.append({
            "tf_symbol": tf,
            "target_count": k_targets,
            "mean_target_log2fc": round(mean_target_lfc, 4),
            "median_target_log2fc": round(median_target_lfc, 4),
            "activity_z_score": round(z_activity, 3),
            "p_welch": float(p_welch),
            "p_mann_whitney": float(p_mwu),
            "p_ks_test": float(p_ks),
            "p_fisher_overlap": float(p_fisher),
            "sig_overlap_count": n_overlap,
            "sig_overlap_genes": ";".join(sorted(list(overlap_sig)))
        })

    tf_df = pd.DataFrame(tf_records)

    # Compute False Discovery Rates
    tf_df["fdr_welch"] = benjamini_hochberg(tf_df["p_welch"].values)
    tf_df["fdr_mann_whitney"] = benjamini_hochberg(tf_df["p_mann_whitney"].values)
    tf_df["fdr_ks"] = benjamini_hochberg(tf_df["p_ks_test"].values)
    tf_df["fdr_fisher"] = benjamini_hochberg(tf_df["p_fisher_overlap"].values)

    # Assign Regulon Status
    tf_df["regulon_status"] = "Unchanged"
    tf_df.loc[(tf_df["activity_z_score"] >= 1.96) & (tf_df["fdr_welch"] < 0.05), "regulon_status"] = "Activated"
    tf_df.loc[(tf_df["activity_z_score"] <= -1.96) & (tf_df["fdr_welch"] < 0.05), "regulon_status"] = "Repressed"

    # Sort primarily by absolute Activity Z-score
    tf_df["abs_z"] = tf_df["activity_z_score"].abs()
    tf_df = tf_df.sort_values(by="abs_z", ascending=False).drop(columns=["abs_z"])

    out_path = os.path.join(out_dir, "tf_regulon_activity_summary.csv")
    tf_df.to_csv(out_path, index=False)
    logger.info(f"Saved TF regulon summary -> {out_path} ({len(tf_df)} TFs evaluated)")

    n_act = (tf_df["regulon_status"] == "Activated").sum()
    n_rep = (tf_df["regulon_status"] == "Repressed").sum()

    logger.info("=" * 65)
    logger.info(f"[SUCCESS] TF Regulon Deconvolution completed across {len(tf_df)} TFs:")
    logger.info(f"  - Significantly Activated Regulons: {n_act}")
    logger.info(f"  - Significantly Repressed Regulons: {n_rep}")

    # Log landmark TFs
    landmarks = ["Fosb", "Fos", "Jun", "Rela", "Nfkb1", "Spi1", "Stat1", "Cebpb"]
    for tf in landmarks:
        sub = tf_df[tf_df["tf_symbol"] == tf]
        if not sub.empty:
            r = sub.iloc[0]
            logger.info(f"  * {r['tf_symbol']}: Targets={r['target_count']}, Mean LFC={r['mean_target_log2fc']:+.3f}, Z={r['activity_z_score']:+.2f}, Status={r['regulon_status']}")
    logger.info("=" * 65)

if __name__ == "__main__":
    run_tf_regulon_analysis()
