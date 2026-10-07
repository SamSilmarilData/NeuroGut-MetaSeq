#!/usr/bin/env python3
"""
scripts/04_meta_analysis.py
Performs statistical meta-analysis across multiple RNA-seq cohorts:
- DerSimonian-Laird Random-Effects inverse-variance pooling (log2FC, SE & 95% CI)
- Heterogeneity quantification: Cochran's Q, between-study variance tau^2, and Higgins I^2
- Complementary non-parametric p-value combination: Fisher's chi-square and sample-size weighted Stouffer's Z
- Multiple hypothesis correction: Benjamini-Hochberg FDR
- Leave-One-Out (LOO) sensitivity analysis for multi-cohort robustness
- Heterogeneity tiering and directional concordance classification
- Two-tiered core consensus signature extraction (Tier 1: Omnipresent k=4, Tier 2: Robust Broad k=3)

Outputs:
    results/meta_results/microglia_meta_analysis_summary.csv
    results/meta_results/microglia_meta_analysis_loo.csv
    results/meta_results/core_consensus_signature.csv
"""

import os
import sys
import argparse
import logging
import glob
import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Meta-Analysis")

def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Computes Benjamini-Hochberg False Discovery Rate (FDR)."""
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
    # Restore original order
    rev_order = np.empty(n, dtype=int)
    rev_order[order] = np.arange(n)
    return fdr[rev_order]

def dersimonian_laird(thetas: np.ndarray, variances: np.ndarray):
    """
    Computes DerSimonian-Laird random-effects pooled effect size,
    between-study variance tau^2, Cochran's Q, and Higgins I^2.
    """
    k = len(thetas)
    if k == 1:
        se = float(np.sqrt(variances[0]))
        return float(thetas[0]), se, float(thetas[0] - 1.96 * se), float(thetas[0] + 1.96 * se), 0.0, 0.0, 0.0, 1.0

    w_fe = 1.0 / variances
    sum_w = np.sum(w_fe)
    theta_fe = np.sum(w_fe * thetas) / sum_w

    # Cochran's Q
    Q = float(np.sum(w_fe * (thetas - theta_fe)**2))
    df = k - 1
    p_Q = 1.0 - stats.chi2.cdf(Q, df) if df > 0 else 1.0

    # Higgins I^2 (%)
    I2 = float(max(0.0, (Q - df) / Q * 100.0)) if Q > 0 else 0.0

    # Between-study variance tau^2
    denom = sum_w - (np.sum(w_fe**2) / sum_w)
    tau2 = float(max(0.0, (Q - df) / denom)) if denom > 0 else 0.0

    # Random-effects weights
    w_re = 1.0 / (variances + tau2)
    sum_w_re = np.sum(w_re)
    theta_re = float(np.sum(w_re * thetas) / sum_w_re)
    se_re = float(np.sqrt(1.0 / sum_w_re))

    ci_lower = theta_re - 1.96 * se_re
    ci_upper = theta_re + 1.96 * se_re

    # Z-statistic and two-tailed p-value
    z = theta_re / se_re if se_re > 0 else 0.0
    p_val = float(2.0 * (1.0 - stats.norm.cdf(abs(z))))

    return theta_re, se_re, ci_lower, ci_upper, tau2, Q, I2, p_val

def combine_pvalues_fisher(pvals: np.ndarray) -> tuple:
    """Fisher's method: -2 * sum(ln(p)) ~ chi2(2k)."""
    p_clipped = np.clip(pvals, 1e-15, 1.0)
    chi2_stat = float(-2.0 * np.sum(np.log(p_clipped)))
    df = 2 * len(pvals)
    p_combined = float(1.0 - stats.chi2.cdf(chi2_stat, df))
    return chi2_stat, p_combined

def combine_pvalues_stouffer(pvals: np.ndarray, thetas: np.ndarray, sample_sizes: np.ndarray) -> tuple:
    """Sample-size weighted Stouffer's Z-transform."""
    p_clipped = np.clip(pvals, 1e-15, 1.0 - 1e-15)
    z_scores = stats.norm.ppf(1.0 - p_clipped / 2.0) * np.sign(thetas)
    weights = np.sqrt(sample_sizes)
    denom = np.sqrt(np.sum(weights**2))
    z_weighted = float(np.sum(weights * z_scores) / denom) if denom > 0 else 0.0
    p_stouffer = float(2.0 * (1.0 - stats.norm.cdf(abs(z_weighted))))
    return z_weighted, p_stouffer

def load_cohort_data(deg_dir: str = "results/de_results", meta_dir: str = "data/metadata") -> dict:
    """Loads and standardizes DEG tables and sample sizes for all cohorts."""
    deg_files = sorted(glob.glob(os.path.join(deg_dir, "*_deg.csv")))
    if not deg_files:
        raise FileNotFoundError(f"No DEG files found in {deg_dir}")

    cohort_data = {}
    for f in deg_files:
        cname = os.path.basename(f).replace("_deg.csv", "")
        df = pd.read_csv(f)
        df = df.dropna(subset=["gene_symbol", "log2FoldChange", "lfcSE", "pvalue"])
        # Keep positive finite SE
        df = df[np.isfinite(df["log2FoldChange"]) & np.isfinite(df["lfcSE"]) & (df["lfcSE"] > 0) & np.isfinite(df["pvalue"])]
        # Deduplicate gene symbols by prioritizing largest baseMean
        if "baseMean" in df.columns:
            df = df.sort_values(by="baseMean", ascending=False).drop_duplicates(subset=["gene_symbol"])
        else:
            df = df.drop_duplicates(subset=["gene_symbol"])

        meta_f = os.path.join(meta_dir, f"{cname}_metadata.csv")
        n_samples = len(pd.read_csv(meta_f)) if os.path.exists(meta_f) else 10
        cohort_data[cname] = {
            "df": df.set_index("gene_symbol"),
            "n": n_samples
        }
        logger.info(f"Loaded {cname}: {len(df)} expressed genes (N={n_samples} biological samples)")

    return cohort_data

def run_meta_analysis(deg_dir: str = "results/de_results",
                      meta_dir: str = "data/metadata",
                      out_dir: str = "results/meta_results"):
    """Executes full random-effects meta-analysis and LOO sensitivity analysis."""
    os.makedirs(out_dir, exist_ok=True)
    logger.info("Initializing Cross-Study Statistical Meta-Analysis (Horizon 3)...")

    cohort_data = load_cohort_data(deg_dir=deg_dir, meta_dir=meta_dir)
    all_cohort_names = sorted(list(cohort_data.keys()))
    logger.info(f"Analyzing {len(all_cohort_names)} cohorts: {all_cohort_names}")

    # Map all gene symbols present in at least 2 cohorts
    gene_to_cohorts = {}
    for cname, cinfo in cohort_data.items():
        for gene in cinfo["df"].index:
            if gene not in gene_to_cohorts:
                gene_to_cohorts[gene] = []
            gene_to_cohorts[gene].append(cname)

    evaluated_genes = sorted([g for g, clist in gene_to_cohorts.items() if len(clist) >= 2])
    logger.info(f"Total genes detected in >= 2 cohorts: {len(evaluated_genes)} loci")

    # 1. Main Meta-Analysis Computation
    records = []
    loo_records = []

    for gene in evaluated_genes:
        thetas = []
        variances = []
        pvals = []
        ns = []
        cnames = []

        for cname in all_cohort_names:
            cinfo = cohort_data[cname]
            if gene in cinfo["df"].index:
                row = cinfo["df"].loc[gene]
                lfc = float(row["log2FoldChange"])
                se = float(row["lfcSE"])
                p = float(row["pvalue"])
                thetas.append(lfc)
                variances.append(se**2)
                pvals.append(p)
                ns.append(cinfo["n"])
                cnames.append(cname)

        k = len(thetas)
        thetas_arr = np.array(thetas)
        vars_arr = np.array(variances)
        pvals_arr = np.array(pvals)
        ns_arr = np.array(ns)

        # DerSimonian-Laird RE
        theta_re, se_re, ci_low, ci_high, tau2, Q, I2, p_re = dersimonian_laird(thetas_arr, vars_arr)
        chi2_f, p_f = combine_pvalues_fisher(pvals_arr)
        z_s, p_s = combine_pvalues_stouffer(pvals_arr, thetas_arr, ns_arr)

        # Directional Concordance
        if (thetas_arr > 0).all():
            concordance = "Concordant Up"
        elif (thetas_arr < 0).all():
            concordance = "Concordant Down"
        else:
            concordance = "Mixed"

        # Heterogeneity Tier
        if I2 < 25.0:
            het_tier = "Low (<25%)"
        elif I2 <= 75.0:
            het_tier = "Moderate (25-75%)"
        else:
            het_tier = "High (>75%)"

        record = {
            "gene_symbol": gene,
            "n_cohorts": k,
            "cohorts_detected": ";".join(cnames),
            "meta_log2fc": round(theta_re, 4),
            "meta_se": round(se_re, 4),
            "ci_lower": round(ci_low, 4),
            "ci_upper": round(ci_high, 4),
            "tau2": round(tau2, 4),
            "cochran_q": round(Q, 3),
            "i2_heterogeneity": round(I2, 1),
            "heterogeneity_tier": het_tier,
            "p_random_effects": p_re,
            "fisher_stat": round(chi2_f, 3),
            "p_fisher": p_f,
            "stouffer_z": round(z_s, 3),
            "p_stouffer": p_s,
            "direction_concordance": concordance
        }
        records.append(record)

        # 2. Leave-One-Out (LOO) Sensitivity calculations for k >= 3
        if k >= 3:
            for i, omitted in enumerate(cnames):
                loo_thetas = np.delete(thetas_arr, i)
                loo_vars = np.delete(vars_arr, i)
                loo_ns = np.delete(ns_arr, i)
                loo_theta_re, loo_se_re, _, _, _, loo_q, loo_i2, loo_p_re = dersimonian_laird(loo_thetas, loo_vars)
                loo_shift = abs(theta_re - loo_theta_re)
                loo_records.append({
                    "gene_symbol": gene,
                    "omitted_cohort": omitted,
                    "k_remaining": k - 1,
                    "loo_log2fc": round(loo_theta_re, 4),
                    "loo_se": round(loo_se_re, 4),
                    "loo_cochran_q": round(loo_q, 3),
                    "loo_i2": round(loo_i2, 1),
                    "loo_p_value": loo_p_re,
                    "shift_from_full": round(loo_shift, 4)
                })

    meta_df = pd.DataFrame(records)

    # Compute False Discovery Rates
    meta_df["fdr_random_effects"] = benjamini_hochberg(meta_df["p_random_effects"].values)
    meta_df["fdr_fisher"] = benjamini_hochberg(meta_df["p_fisher"].values)
    meta_df["fdr_stouffer"] = benjamini_hochberg(meta_df["p_stouffer"].values)

    # Process Leave-One-Out table and calculate per-slice FDR
    loo_df = pd.DataFrame(loo_records)
    if not loo_df.empty:
        # Group by omitted cohort to calculate BH FDR per omitted slice
        loo_df["loo_fdr"] = 1.0
        for omitted_c in loo_df["omitted_cohort"].unique():
            mask = loo_df["omitted_cohort"] == omitted_c
            loo_df.loc[mask, "loo_fdr"] = benjamini_hochberg(loo_df.loc[mask, "loo_p_value"].values)

        # Aggregate LOO robustness per gene
        loo_summary = {}
        for gene, group in loo_df.groupby("gene_symbol"):
            total_folds = len(group)
            sig_folds = (group["loo_fdr"] < 0.05).sum()
            robustness = round(sig_folds / total_folds, 3)
            max_shift = round(group["shift_from_full"].max(), 4)
            # Find vulnerable study causing max shift
            vulnerable = group.sort_values(by="shift_from_full", ascending=False).iloc[0]["omitted_cohort"]
            loo_summary[gene] = {
                "robustness_score": robustness,
                "max_lfc_shift": max_shift,
                "loo_vulnerable_study": vulnerable
            }

        loo_sum_df = pd.DataFrame.from_dict(loo_summary, orient="index")
        meta_df = meta_df.set_index("gene_symbol").join(loo_sum_df).reset_index()
    else:
        meta_df["robustness_score"] = np.nan
        meta_df["max_lfc_shift"] = np.nan
        meta_df["loo_vulnerable_study"] = None

    # Flag primary consensus significance: FDR_RE < 0.05 and |meta_log2fc| >= 0.5
    meta_df["significance_flag"] = (meta_df["fdr_random_effects"] < 0.05) & (meta_df["meta_log2fc"].abs() >= 0.5)

    # Sort primarily by random-effects FDR, then by absolute effect size
    meta_df = meta_df.sort_values(by=["fdr_random_effects", "meta_log2fc"], ascending=[True, False])

    # Save complete meta-analysis table
    summary_path = os.path.join(out_dir, "microglia_meta_analysis_summary.csv")
    meta_df.to_csv(summary_path, index=False)

    # Save detailed LOO table
    loo_path = os.path.join(out_dir, "microglia_meta_analysis_loo.csv")
    loo_df.to_csv(loo_path, index=False)

    # 3. Stratified Core Consensus Signature Extraction
    # Requirements: significance_flag == True, n_cohorts >= 3, i2_heterogeneity < 50.0
    core_mask = (meta_df["significance_flag"]) & (meta_df["n_cohorts"] >= 3) & (meta_df["i2_heterogeneity"] < 50.0)
    core_df = meta_df[core_mask].copy()

    core_df["consensus_tier"] = np.where(
        core_df["n_cohorts"] == 4,
        "Tier 1: Omnipresent Core (k=4)",
        "Tier 2: Robust Broad Core (k=3)"
    )

    core_path = os.path.join(out_dir, "core_consensus_signature.csv")
    core_df.to_csv(core_path, index=False)

    # Summary Statistics
    n_total = len(meta_df)
    n_sig = meta_df["significance_flag"].sum()
    n_up = ((meta_df["significance_flag"]) & (meta_df["meta_log2fc"] > 0)).sum()
    n_down = ((meta_df["significance_flag"]) & (meta_df["meta_log2fc"] < 0)).sum()
    n_tier1 = (core_df["consensus_tier"] == "Tier 1: Omnipresent Core (k=4)").sum()
    n_tier2 = (core_df["consensus_tier"] == "Tier 2: Robust Broad Core (k=3)").sum()

    logger.info("=" * 70)
    logger.info(f"[SUCCESS] Meta-analysis completed across {n_total} genes.")
    logger.info(f"Total Consensus Significant Genes (FDR_RE < 0.05, |LFC| >= 0.5): {n_sig}")
    logger.info(f"  - Consensus Upregulated:   {n_up}")
    logger.info(f"  - Consensus Downregulated: {n_down}")
    logger.info(f"Core Invariant Consensus Signature (I² < 50%, k >= 3): {len(core_df)}")
    logger.info(f"  - Tier 1: Omnipresent Core (k = 4): {n_tier1}")
    logger.info(f"  - Tier 2: Robust Broad Core (k = 3): {n_tier2}")
    logger.info(f"Master Summary -> {summary_path}")
    logger.info(f"LOO Sensitivity -> {loo_path}")
    logger.info(f"Core Consensus Signature -> {core_path}")
    logger.info("=" * 70)

def main():
    parser = argparse.ArgumentParser(description="Cross-study statistical meta-analysis for microglia")
    parser.add_argument("--deg_dir", default="results/de_results", help="Directory with DEG CSV files")
    parser.add_argument("--meta_dir", default="data/metadata", help="Directory with sample metadata")
    parser.add_argument("--out_dir", default="results/meta_results", help="Output directory")
    args = parser.parse_args()

    run_meta_analysis(deg_dir=args.deg_dir, meta_dir=args.meta_dir, out_dir=args.out_dir)

if __name__ == "__main__":
    main()
