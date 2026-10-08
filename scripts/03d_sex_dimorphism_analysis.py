#!/usr/bin/env python3
"""
scripts/03d_sex_dimorphism_analysis.py
Factorial Sex-by-Condition Interaction Meta-Regression Analysis:
Models ~ sex + condition + sex:condition across 51 sex-informative biological samples
(GSE107925, GSE108045, GSE186210).

Computes:
1. Male-specific log2FC and standard errors per cohort.
2. Female-specific log2FC and standard errors per cohort.
3. Interaction effect sizes (theta_int = theta_female - theta_male) per cohort.
4. Random-effects meta-pooling across cohorts with I2_sex heterogeneity.
5. Three-Tier classification: Sex-Shared, Female-Biased, Male-Biased.
6. Publication figures: fig_sex_stratified_forest (PNG + SVG), fig_sex_concordance_scatter (PNG + SVG).

Outputs:
    results/meta_results/sex_dimorphism_meta_analysis.csv
    results/meta_results/figures/fig_sex_stratified_forest.png / .svg
    results/meta_results/figures/fig_sex_concordance_scatter.png / .svg
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Sex-Dimorphism")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "axes.edgecolor": "#333333",
    "axes.linewidth": 1.1,
    "grid.color": "#EAEAEA",
    "grid.linestyle": "--",
    "grid.linewidth": 0.6,
    "figure.titlesize": 13,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5
})

SEX_COHORTS = ["GSE107925", "GSE108045", "GSE186210"]

def compute_cpm(counts: np.ndarray) -> np.ndarray:
    total = np.sum(counts, axis=0)
    total[total == 0] = 1.0
    return (counts / total) * 1e6

def fit_cohort_interaction(cohort_id: str, counts_df: pd.DataFrame, meta_df: pd.DataFrame):
    """Fits factorial linear model per gene: y = b0 + b_sex*Sex + b_cond*Cond + b_int*(Sex*Cond)."""
    logger.info(f"Fitting sex-by-condition interaction for cohort: {cohort_id}...")
    sample_ids = meta_df["sample_id"].tolist()
    counts_sub = counts_df[sample_ids].values.astype(float)
    cpm = compute_cpm(counts_sub)
    log_cpm = np.log2(cpm + 1.0)
    
    # Coding: Sex: Male = 0, Female = 1; Condition: Reference = 0, Perturbed = 1
    sex_vec = np.array([1.0 if str(s).strip().lower() == "female" else 0.0 for s in meta_df["sex"]])
    cond_vec = np.array([1.0 if str(c).strip().lower() == "perturbed" else 0.0 for c in meta_df["condition"]])
    
    n_samples = len(sample_ids)
    X = np.column_stack([np.ones(n_samples), sex_vec, cond_vec, sex_vec * cond_vec])
    
    # Solve via normal equations (X^T X)^{-1} X^T Y
    XtX = X.T @ X
    inv_XtX = np.linalg.pinv(XtX)
    beta = inv_XtX @ (X.T @ log_cpm.T)  # Shape: (4, n_genes)
    
    # Residuals & standard errors
    y_pred = (X @ beta).T  # Shape: (n_genes, n_samples)
    residuals = log_cpm - y_pred
    dof = n_samples - 4
    if dof <= 0:
        dof = 1
    sigma2 = np.sum(residuals**2, axis=1) / dof
    
    # Covariance diagonals
    # beta_0: intercept, beta_1: sex (F vs M in ref), beta_2: cond (Male effect), beta_3: interaction
    se_beta2 = np.sqrt(sigma2 * inv_XtX[2, 2])  # SE Male effect
    se_beta3 = np.sqrt(sigma2 * inv_XtX[3, 3])  # SE Interaction
    
    # Var(beta_2 + beta_3) = Var(beta_2) + Var(beta_3) + 2*Cov(beta_2, beta_3)
    var_female = sigma2 * (inv_XtX[2, 2] + inv_XtX[3, 3] + 2.0 * inv_XtX[2, 3])
    var_female[var_female < 0] = 1e-6
    se_female = np.sqrt(var_female)
    
    theta_male = beta[2, :]
    theta_female = beta[2, :] + beta[3, :]
    theta_int = beta[3, :]
    
    res_df = pd.DataFrame({
        "gene_symbol": counts_df["gene_symbol"].values,
        f"{cohort_id}_theta_male": theta_male,
        f"{cohort_id}_se_male": np.maximum(se_beta2, 0.01),
        f"{cohort_id}_theta_female": theta_female,
        f"{cohort_id}_se_female": np.maximum(se_female, 0.01),
        f"{cohort_id}_theta_int": theta_int,
        f"{cohort_id}_se_int": np.maximum(se_beta3, 0.01)
    })
    return res_df

def dersimonian_laird_pool(thetas: np.ndarray, ses: np.ndarray):
    """Vectorized DL random-effects pooling across K studies."""
    k = thetas.shape[1]
    vars_ = ses**2
    w_fe = 1.0 / np.maximum(vars_, 1e-6)
    w_sum = np.sum(w_fe, axis=1)
    theta_fe = np.sum(w_fe * thetas, axis=1) / w_sum
    
    # Cochran's Q
    q = np.sum(w_fe * (thetas - theta_fe[:, None])**2, axis=1)
    df = k - 1
    denom = w_sum - np.sum(w_fe**2, axis=1) / w_sum
    denom[denom <= 0] = 1e-6
    tau2 = np.maximum(0.0, (q - df) / denom)
    
    # Random-effects weights
    w_re = 1.0 / (vars_ + tau2[:, None])
    w_re_sum = np.sum(w_re, axis=1)
    theta_re = np.sum(w_re * thetas, axis=1) / w_re_sum
    se_re = np.sqrt(1.0 / w_re_sum)
    
    # Heterogeneity I^2
    q_adj = np.maximum(q, 1e-6)
    i2 = np.maximum(0.0, (q - df) / q_adj) * 100.0
    
    # Z and p-value
    z = theta_re / np.maximum(se_re, 1e-6)
    p_val = 2.0 * (1.0 - stats.norm.cdf(np.abs(z)))
    
    return theta_re, se_re, p_val, i2, q

def run_sex_dimorphism_analysis():
    logger.info("Initializing Factorial Sex-by-Condition Interaction Meta-Regression...")
    
    cohort_dfs = []
    for cid in SEX_COHORTS:
        counts_path = f"data/processed/{cid}_counts.csv"
        meta_path = f"data/metadata/{cid}_metadata.csv"
        if not os.path.exists(counts_path) or not os.path.exists(meta_path):
            raise FileNotFoundError(f"Missing data files for {cid}")
            
        cdf = pd.read_csv(counts_path)
        mdf = pd.read_csv(meta_path)
        c_res = fit_cohort_interaction(cid, cdf, mdf)
        cohort_dfs.append(c_res)
        
    # Merge across cohorts on gene_symbol
    merged = cohort_dfs[0]
    for c_res in cohort_dfs[1:]:
        merged = pd.merge(merged, c_res, on="gene_symbol", how="inner")
        
    logger.info(f"Retained {len(merged)} common genes across 3 sex-informative cohorts.")
    
    # Extract matrices
    male_thetas = merged[[f"{c}_theta_male" for c in SEX_COHORTS]].values
    male_ses = merged[[f"{c}_se_male" for c in SEX_COHORTS]].values
    
    female_thetas = merged[[f"{c}_theta_female" for c in SEX_COHORTS]].values
    female_ses = merged[[f"{c}_se_female" for c in SEX_COHORTS]].values
    
    int_thetas = merged[[f"{c}_theta_int" for c in SEX_COHORTS]].values
    int_ses = merged[[f"{c}_se_int" for c in SEX_COHORTS]].values
    
    # Pool across cohorts
    th_m, se_m, p_m, i2_m, _ = dersimonian_laird_pool(male_thetas, male_ses)
    th_f, se_f, p_f, i2_f, _ = dersimonian_laird_pool(female_thetas, female_ses)
    th_int, se_int, p_int, i2_int, q_int = dersimonian_laird_pool(int_thetas, int_ses)
    
    # Benjamini-Hochberg FDR
    def bh_fdr(p_vals):
        n = len(p_vals)
        order = np.argsort(p_vals)
        ranked_p = p_vals[order]
        fdr = np.zeros(n)
        cummin = 1.0
        for i in range(n - 1, -1, -1):
            val = ranked_p[i] * n / (i + 1)
            cummin = min(cummin, val)
            fdr[order[i]] = min(1.0, cummin)
        return fdr
        
    fdr_m = bh_fdr(p_m)
    fdr_f = bh_fdr(p_f)
    fdr_int = bh_fdr(p_int)
    
    # Assign Three-Tier Classification
    # 1. Sex-Shared: |theta_int| < 0.585 and (p_int >= 0.05 or I2_int < 25%)
    # 2. Female-Biased: theta_int >= 0.585 and p_int < 0.05
    # 3. Male-Biased: theta_int <= -0.585 and p_int < 0.05
    tiers = []
    for ti, pi, i2 in zip(th_int, p_int, i2_int):
        if abs(ti) < 0.585 and (pi >= 0.05 or i2 < 25.0):
            tiers.append("Sex-Shared")
        elif ti >= 0.585 and pi < 0.05:
            tiers.append("Female-Biased Vulnerability")
        elif ti <= -0.585 and pi < 0.05:
            tiers.append("Male-Biased Vulnerability")
        else:
            tiers.append("Model-Discordant / Ambiguous")
            
    res_df = pd.DataFrame({
        "gene_symbol": merged["gene_symbol"],
        "pooled_male_log2fc": np.round(th_m, 4),
        "pooled_male_se": np.round(se_m, 4),
        "pooled_male_pval": p_m,
        "pooled_male_fdr": fdr_m,
        "pooled_female_log2fc": np.round(th_f, 4),
        "pooled_female_se": np.round(se_f, 4),
        "pooled_female_pval": p_f,
        "pooled_female_fdr": fdr_f,
        "interaction_log2fc": np.round(th_int, 4),
        "interaction_se": np.round(se_int, 4),
        "interaction_pval": p_int,
        "interaction_fdr": fdr_int,
        "i2_sex_heterogeneity": np.round(i2_int, 2),
        "cochran_q_sex": np.round(q_int, 2),
        "sex_dimorphism_tier": tiers
    })
    
    out_dir = "results/meta_results"
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "sex_dimorphism_meta_analysis.csv")
    res_df.to_csv(out_csv, index=False)
    logger.info(f"[SUCCESS] Saved sex dimorphism meta-analysis table -> {out_csv} ({len(res_df)} genes)")
    
    tier_counts = res_df["sex_dimorphism_tier"].value_counts()
    print("\n=== Microglial Sex Dimorphism Classification Summary ===")
    print(tier_counts.to_string())
    
    # Plotting
    fig_dir = "results/meta_results/figures"
    os.makedirs(fig_dir, exist_ok=True)
    plot_figures(res_df, fig_dir)
    
    return res_df

def plot_figures(res_df: pd.DataFrame, fig_dir: str):
    logger.info("Generating Publication Figures (PNG + SVG) for Sex Dimorphism...")
    
    # -------------------------------------------------------------------
    # Figure 1: Sex Concordance Scatter Plot
    # -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 7.0), dpi=300)
    
    # Correlation across all genes
    r_val, _ = stats.pearsonr(res_df["pooled_male_log2fc"], res_df["pooled_female_log2fc"])
    rho_val, _ = stats.spearmanr(res_df["pooled_male_log2fc"], res_df["pooled_female_log2fc"])
    
    # Background points
    ax.scatter(
        res_df["pooled_male_log2fc"], res_df["pooled_female_log2fc"],
        c="#D0D0D0", alpha=0.35, s=12, edgecolors="none", rasterized=True,
        label=f"All Genes (N={len(res_df):,}; r={r_val:.2f}, rho={rho_val:.2f})"
    )
    
    # Highlight Sex-Biased
    fem_biased = res_df[res_df["sex_dimorphism_tier"] == "Female-Biased Vulnerability"]
    male_biased = res_df[res_df["sex_dimorphism_tier"] == "Male-Biased Vulnerability"]
    
    ax.scatter(
        fem_biased["pooled_male_log2fc"], fem_biased["pooled_female_log2fc"],
        c="#E15759", alpha=0.75, s=24, edgecolors="none",
        label=f"Female-Biased (N={len(fem_biased)})"
    )
    ax.scatter(
        male_biased["pooled_male_log2fc"], male_biased["pooled_female_log2fc"],
        c="#4E79A7", alpha=0.75, s=24, edgecolors="none",
        label=f"Male-Biased (N={len(male_biased)})"
    )
    
    # Key Landmark Genes
    landmarks = ["Irf1", "Stat1", "Oas1a", "Slfn2", "Llgl2", "Clu", "Plin3", "Tnf", "Fosb", "Sap30", "Tsc22d3"]
    landmark_df = res_df[res_df["gene_symbol"].isin(landmarks)]
    
    for _, row in landmark_df.iterrows():
        g = row["gene_symbol"]
        mx, fy = row["pooled_male_log2fc"], row["pooled_female_log2fc"]
        tier = row["sex_dimorphism_tier"]
        color = "#2CA02C" if tier == "Sex-Shared" else "#E15759" if "Female" in tier else "#4E79A7"
        ax.scatter(mx, fy, color=color, s=65, edgecolors="#111111", lw=1.2, zorder=5)
        ax.annotate(
            g, (mx, fy), textcoords="offset points", xytext=(5, 4),
            fontsize=9, fontweight="bold", color="#111111", zorder=6
        )
        
    ax.axhline(0, color="#888888", linestyle=":", lw=0.9)
    ax.axvline(0, color="#888888", linestyle=":", lw=0.9)
    # Diagonal 1:1 line
    lims = [-4.5, 4.5]
    ax.plot(lims, lims, color="#222222", linestyle="--", lw=1.0, alpha=0.7, label="1:1 Concordance")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    ax.set_xlabel("Male Pooled Effect Size (log2FC, Perturbed vs Reference)", fontweight="bold")
    ax.set_ylabel("Female Pooled Effect Size (log2FC, Perturbed vs Reference)", fontweight="bold")
    ax.set_title("Cross-Sex Concordance of Microglial Response to Depletion\n(N=51 Samples, 3 Factorial Cohorts)", fontweight="bold", pad=12)
    ax.legend(loc="upper left", frameon=True, fontsize=8.5)
    ax.grid(True, alpha=0.3)
    
    plt.savefig(os.path.join(fig_dir, "fig_sex_concordance_scatter.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(fig_dir, "fig_sex_concordance_scatter.svg"), format="svg", bbox_inches="tight")
    plt.close()
    logger.info("Saved fig_sex_concordance_scatter (PNG + SVG)")
    
    # -------------------------------------------------------------------
    # Figure 2: Sex-Stratified Forest Plot for Key Targets
    # -------------------------------------------------------------------
    forest_genes = ["Slfn2", "Llgl2", "Clu", "Plin3", "Irf1", "Stat1", "Oas1a", "Tnf", "Fosb", "Sap30"]
    forest_sub = res_df[res_df["gene_symbol"].isin(forest_genes)].copy()
    forest_sub["sort_val"] = forest_sub["pooled_female_log2fc"]
    forest_sub = forest_sub.sort_values(by="sort_val", ascending=True)
    
    fig, ax = plt.subplots(figsize=(8.5, 6.5), dpi=300)
    y_pos = np.arange(len(forest_sub))
    offset = 0.18
    
    # Male estimates
    ax.errorbar(
        forest_sub["pooled_male_log2fc"], y_pos - offset,
        xerr=1.96 * forest_sub["pooled_male_se"], fmt="o", color="#4E79A7",
        ecolor="#4E79A7", elinewidth=1.5, capsize=3.5, label="Male (Pooled)", markersize=6
    )
    # Female estimates
    ax.errorbar(
        forest_sub["pooled_female_log2fc"], y_pos + offset,
        xerr=1.96 * forest_sub["pooled_female_se"], fmt="s", color="#E15759",
        ecolor="#E15759", elinewidth=1.5, capsize=3.5, label="Female (Pooled)", markersize=6
    )
    
    ax.axvline(0, color="#444444", linestyle="--", lw=1.0)
    ax.set_yticks(y_pos)
    labels = [f"{g} [{tier.replace(' Vulnerability', '')}]" for g, tier in zip(forest_sub["gene_symbol"], forest_sub["sex_dimorphism_tier"])]
    ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
    ax.set_xlabel("Pooled Effect Size (log2FC ± 95% CI)", fontweight="bold")
    ax.set_title("Sex-Stratified Random-Effects Meta-Estimates of Landmark Genes", fontweight="bold", pad=12)
    ax.legend(loc="lower right", frameon=True)
    ax.grid(True, axis="x", alpha=0.3)
    
    plt.savefig(os.path.join(fig_dir, "fig_sex_stratified_forest.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(fig_dir, "fig_sex_stratified_forest.svg"), format="svg", bbox_inches="tight")
    plt.close()
    logger.info("Saved fig_sex_stratified_forest (PNG + SVG)")

if __name__ == "__main__":
    run_sex_dimorphism_analysis()
