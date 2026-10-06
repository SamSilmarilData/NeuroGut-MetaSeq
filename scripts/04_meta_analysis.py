#!/usr/bin/env python3
"""
scripts/04_meta_analysis.py
Performs statistical meta-analysis across multiple RNA-seq cohorts:
- DerSimonian-Laird Random-Effects inverse-variance pooling (log2FC & 95% CI)
- Heterogeneity testing: Cochran's Q and Higgins I²
- Non-parametric p-value combination: Fisher's chi-square and Stouffer's Z-test
- Multiple hypothesis correction: Benjamini-Hochberg FDR
Outputs:
    results/meta_results/microglia_meta_analysis_summary.csv
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
        return thetas[0], np.sqrt(variances[0]), thetas[0] - 1.96*np.sqrt(variances[0]), thetas[0] + 1.96*np.sqrt(variances[0]), 0.0, 0.0, 0.0, 1.0

    w_fe = 1.0 / variances
    sum_w = np.sum(w_fe)
    theta_fe = np.sum(w_fe * thetas) / sum_w

    # Cochran's Q
    Q = np.sum(w_fe * (thetas - theta_fe)**2)
    df = k - 1
    p_Q = 1.0 - stats.chi2.cdf(Q, df) if df > 0 else 1.0

    # Higgins I^2 (%)
    I2 = max(0.0, (Q - df) / Q * 100.0) if Q > 0 else 0.0

    # Between-study variance tau^2
    denom = sum_w - (np.sum(w_fe**2) / sum_w)
    tau2 = max(0.0, (Q - df) / denom) if denom > 0 else 0.0

    # Random-effects weights
    w_re = 1.0 / (variances + tau2)
    sum_w_re = np.sum(w_re)
    theta_re = np.sum(w_re * thetas) / sum_w_re
    se_re = np.sqrt(1.0 / sum_w_re)

    ci_lower = theta_re - 1.96 * se_re
    ci_upper = theta_re + 1.96 * se_re

    # Z-statistic and p-value for pooled effect
    z = theta_re / se_re if se_re > 0 else 0.0
    p_val = 2.0 * (1.0 - stats.norm.cdf(abs(z)))

    return theta_re, se_re, ci_lower, ci_upper, tau2, Q, I2, p_val

def combine_pvalues_fisher(pvals: np.ndarray) -> tuple:
    """Fisher's method: -2 * sum(ln(p)) ~ chi2(2k)."""
    p_clipped = np.clip(pvals, 1e-15, 1.0)
    chi2_stat = -2.0 * np.sum(np.log(p_clipped))
    df = 2 * len(pvals)
    p_combined = 1.0 - stats.chi2.cdf(chi2_stat, df)
    return chi2_stat, p_combined

def combine_pvalues_stouffer(pvals: np.ndarray, thetas: np.ndarray, sample_sizes: np.ndarray) -> tuple:
    """Sample-size weighted Stouffer's Z-transform."""
    p_clipped = np.clip(pvals, 1e-15, 1.0 - 1e-15)
    # One-tailed standard normal quantiles signed by direction of effect
    z_scores = stats.norm.ppf(1.0 - p_clipped / 2.0) * np.sign(thetas)
    weights = np.sqrt(sample_sizes)
    z_weighted = np.sum(weights * z_scores) / np.sqrt(np.sum(weights**2))
    p_stouffer = 2.0 * (1.0 - stats.norm.cdf(abs(z_weighted)))
    return z_weighted, p_stouffer

def run_meta_analysis():
    logger.info("Starting cross-study statistical meta-analysis...")
    deg_files = sorted(glob.glob("results/de_results/*_deg.csv"))
    if not deg_files:
        logger.error("No DEG results found in results/de_results/. Please run DESeq2 first.")
        sys.exit(1)

    logger.info(f"Found {len(deg_files)} DEG cohort files: {[os.path.basename(f) for f in deg_files]}")

    cohort_data = {}
    for f in deg_files:
        cname = os.path.basename(f).replace("_deg.csv", "")
        df = pd.read_csv(f)
        # Filter valid rows
        df = df.dropna(subset=["gene_symbol", "log2FoldChange", "lfcSE", "pvalue"])
        # Estimate sample size from metadata if available
        meta_f = f"data/metadata/{cname}_metadata.csv"
        n_samples = len(pd.read_csv(meta_f)) if os.path.exists(meta_f) else 10
        cohort_data[cname] = {"df": df.set_index("gene_symbol"), "n": n_samples}

    # Identify all genes across studies
    all_genes = sorted(list(set.union(*[set(cd["df"].index) for cd in cohort_data.values()])))
    logger.info(f"Total union of unique genes across cohorts: {len(all_genes)}")

    results = []
    for gene in all_genes:
        thetas = []
        variances = []
        pvals = []
        ns = []
        cnames = []

        for cname, cinfo in cohort_data.items():
            if gene in cinfo["df"].index:
                row = cinfo["df"].loc[gene]
                # Handle possible duplicate gene rows by taking first
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                lfc = float(row["log2FoldChange"])
                se = float(row["lfcSE"])
                p = float(row["pvalue"])
                if np.isfinite(lfc) and np.isfinite(se) and se > 0 and np.isfinite(p):
                    thetas.append(lfc)
                    variances.append(se**2)
                    pvals.append(p)
                    ns.append(cinfo["n"])
                    cnames.append(cname)

        if len(thetas) < 2:
            continue

        thetas = np.array(thetas)
        variances = np.array(variances)
        pvals = np.array(pvals)
        ns = np.array(ns)

        # Meta-analysis calculations
        theta_re, se_re, ci_low, ci_high, tau2, Q, I2, p_re = dersimonian_laird(thetas, variances)
        chi2_f, p_f = combine_pvalues_fisher(pvals)
        z_s, p_s = combine_pvalues_stouffer(pvals, thetas, ns)

        # Direction concordance
        all_up = (thetas > 0).all()
        all_down = (thetas < 0).all()
        concordance = "Concordant Up" if all_up else ("Concordant Down" if all_down else "Mixed")

        results.append({
            "gene_symbol": gene,
            "n_cohorts": len(thetas),
            "cohorts_detected": ";".join(cnames),
            "meta_log2fc": round(theta_re, 4),
            "meta_se": round(se_re, 4),
            "ci_lower": round(ci_low, 4),
            "ci_upper": round(ci_high, 4),
            "cochran_q": round(Q, 3),
            "i2_heterogeneity": round(I2, 1),
            "p_random_effects": p_re,
            "fisher_stat": round(chi2_f, 3),
            "p_fisher": p_f,
            "stouffer_z": round(z_s, 3),
            "p_stouffer": p_s,
            "direction_concordance": concordance
        })

    meta_df = pd.DataFrame(results)

    # Compute FDR adjustments
    meta_df["fdr_random_effects"] = benjamini_hochberg(meta_df["p_random_effects"].values)
    meta_df["fdr_fisher"] = benjamini_hochberg(meta_df["p_fisher"].values)
    meta_df["fdr_stouffer"] = benjamini_hochberg(meta_df["p_stouffer"].values)

    # Flag consensus significant genes: FDR < 0.05 and |meta_log2fc| >= 0.5
    meta_df["significance_flag"] = (meta_df["fdr_fisher"] < 0.05) & (meta_df["meta_log2fc"].abs() >= 0.5)

    # Sort by Fisher FDR
    meta_df = meta_df.sort_values(by="fdr_fisher", ascending=True)

    out_dir = "results/meta_results"
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "microglia_meta_analysis_summary.csv")
    meta_df.to_csv(out_file, index=False)

    n_sig_up = ((meta_df["significance_flag"]) & (meta_df["meta_log2fc"] > 0)).sum()
    n_sig_down = ((meta_df["significance_flag"]) & (meta_df["meta_log2fc"] < 0)).sum()

    logger.info("=" * 60)
    logger.info(f"[SUCCESS] Meta-analysis completed on {len(meta_df)} genes.")
    logger.info(f"Consensus Significant Upregulated:   {n_sig_up} genes (FDR < 0.05, Log2FC >= +0.5)")
    logger.info(f"Consensus Significant Downregulated: {n_sig_down} genes (FDR < 0.05, Log2FC <= -0.5)")
    logger.info(f"Summary written -> {out_file}")
    logger.info("=" * 60)

def main():
    parser = argparse.ArgumentParser(description="Run cross-study statistical meta-analysis")
    parser.add_argument("--demo", action="store_true", help="Analyze demo DEG results")
    args = parser.parse_args()
    run_meta_analysis()

if __name__ == "__main__":
    main()
