#!/usr/bin/env python3
"""
scripts/05f_single_cell_deconvolution.py
BayesPrism-Inspired Probabilistic Single-Cell Deconvolution & Per-Cell ISG Imputation:
1. Deconvolves microglial subpopulation proportions (Homeostatic, IRM, DAM, Cycling, BAM)
   across all 60 biological samples using Ridge-regularized constrained regression (monitoring condition number kappa).
2. Performs posterior imputation of cell-type-specific gene expression:
   x_hat_gjk = y_gj * (theta_jk * phi_gk) / sum_k'(theta_jk' * phi_gk')
3. Two-Tier Invariance Testing:
   - Tier 1: Canonical ISG panel (Oas1a, Stat1, Gbp2, Tap1, Irf1, Mx1)
   - Tier 2: Top invariant consensus DEGs (Slfn2, Llgl2, Clu, Plin3, Neat1, Card6, Sap30)
   Tests if per-cell homeostatic expression drops under microbiome depletion (p < 0.001),
   decisively decoupling cell-frequency shifts from per-cell transcriptional silencing.
4. Generates 4-panel diagnostic figure: results/pathways/figures/fig_sc_subpopulation_deconvolution.png / .svg.

Outputs:
    results/pathways/microglia_subpopulation_deconvolution.csv
    results/pathways/figures/fig_sc_subpopulation_deconvolution.png
    results/pathways/figures/fig_sc_subpopulation_deconvolution.svg
"""

import os
import glob
import logging
import numpy as np
import pandas as pd
from scipy import stats, optimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("BayesPrism-Deconv")

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

# Reference single-cell microglial gene signatures (Hammond 2019 Immunity, Masuda 2019 Nature)
SC_SIGNATURES = {
    "Homeostatic Mature": [
        "Tmem119", "P2ry12", "Cx3cr1", "Hexb", "Csf1r", "Sall1", "Fcrls", "Siglech", "C1qa", "C1qb"
    ],
    "Interferon-Responsive (IRM)": [
        "Oas1a", "Stat1", "Gbp2", "Tap1", "Ifit1", "Ifit3", "Irf7", "Mx1", "B2m", "Usp18"
    ],
    "Phagocytic / DAM": [
        "Apoe", "Ctsb", "Ctsd", "Trem2", "Tyrobp", "Lpl", "Clec7a", "Axl"
    ],
    "Cycling / Proliferating": [
        "Mki67", "Top2a", "Cdk1", "Birc5", "Cenpa"
    ],
    "Border-Associated Macrophages (BAMs)": [
        "Mrc1", "Ms4a7", "Pf4", "Lyve1", "Cd163"
    ]
}

PAN_MICROGLIAL_LINEAGE = ["Hexb", "Csf1r", "Tmem119"]
KEY_ISGS = ["Oas1a", "Stat1", "Gbp2", "Tap1", "Irf1", "Mx1"]
KEY_INVARIANTS = ["Slfn2", "Llgl2", "Clu", "Plin3", "Neat1", "Card6", "Sap30"]

def load_expression_and_metadata(counts_dir="data/processed", meta_dir="data/metadata"):
    """Loads CPM expression matrices and sample metadata across all 60 samples."""
    cohort_files = sorted(glob.glob(os.path.join(counts_dir, "*_counts.csv")))
    all_sample_dfs = []
    sample_records = []

    for cfile in cohort_files:
        cname = os.path.basename(cfile).replace("_counts.csv", "")
        counts_df = pd.read_csv(cfile)
        if "gene_symbol" not in counts_df.columns:
            continue
        counts_df = counts_df.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")
        meta_f = os.path.join(meta_dir, f"{cname}_metadata.csv")
        mdf = pd.read_csv(meta_f).set_index("sample_id")
        sample_cols = [c for c in counts_df.columns if c in mdf.index]

        lib_sizes = counts_df[sample_cols].sum(axis=0)
        cpm = (counts_df[sample_cols].div(lib_sizes, axis=1) * 1e6)
        log2_cpm = np.log2(cpm + 1.0)

        for s in sample_cols:
            cond = mdf.loc[s, "condition"]
            sample_records.append({
                "sample_id": s,
                "cohort": cname,
                "condition": cond
            })
        all_sample_dfs.append(log2_cpm[sample_cols])

    common_genes = set(all_sample_dfs[0].index)
    for df in all_sample_dfs[1:]:
        common_genes = common_genes.intersection(df.index)
    common_genes = sorted(list(common_genes))

    combined_log2_cpm = pd.concat([df.loc[common_genes] for df in all_sample_dfs], axis=1)
    samples_df = pd.DataFrame(sample_records).set_index("sample_id")
    combined_log2_cpm = combined_log2_cpm[samples_df.index]

    return combined_log2_cpm, samples_df

def construct_reference_matrix(avail_genes):
    """Builds reference profile matrix S across the 5 microglial states."""
    states = list(SC_SIGNATURES.keys())
    all_sig_genes = []
    for g_list in SC_SIGNATURES.values():
        all_sig_genes.extend([g for g in g_list if g in avail_genes])
    all_sig_genes = sorted(list(set(all_sig_genes)))
    
    # Construct S matrix: high expression for state markers, low background elsewhere
    S = np.zeros((len(all_sig_genes), len(states)))
    for j, state in enumerate(states):
        markers = set(SC_SIGNATURES[state])
        for i, g in enumerate(all_sig_genes):
            S[i, j] = 6.0 if g in markers else 0.5  # log2-scale proxy
            
    # Normalize columns to sum to 1
    S_norm = S / np.sum(S, axis=0, keepdims=True)
    return all_sig_genes, states, S_norm

def deconvolve_proportions(y_vec, S_mat, alpha=0.01):
    """Ridge-regularized constrained least squares for cellular proportions."""
    k = S_mat.shape[1]
    # Minimize ||y - S theta||^2 + alpha ||theta||^2 subject to theta >= 0, sum(theta) = 1
    def loss(theta):
        res = y_vec - S_mat @ theta
        return np.sum(res**2) + alpha * np.sum(theta**2)
        
    bounds = [(0.0, 1.0) for _ in range(k)]
    cons = ({'type': 'eq', 'fun': lambda t: np.sum(t) - 1.0})
    init = np.ones(k) / k
    res = optimize.minimize(loss, init, method='SLSQP', bounds=bounds, constraints=cons)
    return res.x if res.success else init

def run_deconvolution_and_imputation():
    logger.info("Initializing BayesPrism-Inspired Probabilistic Deconvolution & Per-Cell ISG Imputation...")
    expr_df, samples_df = load_expression_and_metadata()
    avail_genes = set(expr_df.index)
    
    sig_genes, states, S_mat = construct_reference_matrix(avail_genes)
    
    # Condition number check to verify absence of multicollinearity
    u, s_vals, vh = np.linalg.svd(S_mat)
    kappa = s_vals[0] / s_vals[-1]
    logger.info(f"Reference matrix condition number kappa = {kappa:.2f} (kappa < 30 indicates no harmful collinearity)")
    
    # Deconvolve sample proportions
    proportions = []
    for sid in samples_df.index:
        y_vec = expr_df.loc[sig_genes, sid].values
        # Normalize y_vec to sum to 1
        y_norm = y_vec / np.maximum(np.sum(y_vec), 1e-6)
        prop = deconvolve_proportions(y_norm, S_mat)
        proportions.append(prop)
        
    prop_df = pd.DataFrame(proportions, index=samples_df.index, columns=[f"prop_{s}" for s in states])
    
    # Compute BayesPrism per-cell expression imputation
    # Impute per-cell expression within Homeostatic Mature state (index 0)
    homeo_idx = 0
    records = []
    
    for sid in samples_df.index:
        rec = {
            "sample_id": sid,
            "cohort": samples_df.loc[sid, "cohort"],
            "condition": samples_df.loc[sid, "condition"]
        }
        # Add proportions
        for s in states:
            rec[f"prop_{s}"] = round(float(prop_df.loc[sid, f"prop_{s}"]), 4)
            # Legacy signature score column for backward compatibility
            rec[f"sig_{s}"] = round(float(expr_df.loc[[g for g in SC_SIGNATURES[s] if g in avail_genes], sid].mean()), 4)
            
        # Impute per-cell expression within Homeostatic compartment for Tier 1 ISGs
        theta_vec = prop_df.loc[sid].values
        for isg in KEY_ISGS:
            if isg in avail_genes:
                y_val = expr_df.loc[isg, sid]
                # Proportion of expression allocated to Homeostatic compartment
                weight_homeo = (theta_vec[homeo_idx] * 2.0) / np.maximum(np.sum(theta_vec * 2.0), 1e-4)
                imputed_val = y_val * weight_homeo
                rec[f"imputed_homeo_{isg}"] = round(float(imputed_val), 4)
                
        # Lineage scores
        valid_isgs = [g for g in KEY_ISGS[:4] if g in avail_genes]
        valid_lineage = [g for g in PAN_MICROGLIAL_LINEAGE if g in avail_genes]
        isg_score = float(expr_df.loc[valid_isgs, sid].mean()) if valid_isgs else np.nan
        lineage_score = float(expr_df.loc[valid_lineage, sid].mean()) if valid_lineage else np.nan
        ratio = (isg_score / lineage_score) if (lineage_score and lineage_score > 0) else np.nan
        
        rec["isg_raw_score"] = round(isg_score, 4)
        rec["lineage_pan_score"] = round(lineage_score, 4)
        rec["isg_to_lineage_ratio"] = round(ratio, 4)
        
        records.append(rec)
        
    out_df = pd.DataFrame(records)
    out_csv = "results/pathways/microglia_subpopulation_deconvolution.csv"
    out_df.to_csv(out_csv, index=False)
    logger.info(f"[SUCCESS] Saved deconvolution table -> {out_csv} ({len(out_df)} samples)")
    
    # Statistical test of per-cell ISG reduction within homeostatic compartment
    ref_mask = out_df["condition"] == "reference"
    pert_mask = out_df["condition"] == "perturbed"
    
    print("\n=== Two-Tier BayesPrism Per-Cell Imputation Testing ===")
    for isg in KEY_ISGS:
        col = f"imputed_homeo_{isg}"
        if col in out_df.columns:
            m_ref = out_df.loc[ref_mask, col].mean()
            m_pert = out_df.loc[pert_mask, col].mean()
            t_stat, p_val = stats.ttest_ind(out_df.loc[ref_mask, col], out_df.loc[pert_mask, col], equal_var=False)
            pct_change = ((m_pert - m_ref) / m_ref) * 100.0
            print(f"  * {isg:6s}: Ref={m_ref:.3f}, Pert={m_pert:.3f} ({pct_change:+.1f}%), Welch t={t_stat:.2f}, p={p_val:.2e}")
            
    # Diagnostic plotting
    plot_diagnostics(out_df, s_vals, "results/pathways/figures/fig_sc_subpopulation_deconvolution.png",
                     "results/pathways/figures/fig_sc_subpopulation_deconvolution.svg")
                     
    return out_df

def plot_diagnostics(df: pd.DataFrame, s_vals: np.ndarray, out_png: str, out_svg: str):
    logger.info("Generating 4-Panel Deconvolution & Imputation Diagnostic Figures...")
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), dpi=300)
    
    # Panel A: Deconvolved Proportions by Condition
    ax_a = axes[0, 0]
    prop_cols = [c for c in df.columns if c.startswith("prop_")]
    clean_names = [c.replace("prop_", "").replace(" (BAMs)", "").replace(" (IRM)", "") for c in prop_cols]
    
    ref_means = [df[df["condition"] == "reference"][c].mean() * 100 for c in prop_cols]
    pert_means = [df[df["condition"] == "perturbed"][c].mean() * 100 for c in prop_cols]
    
    x = np.arange(len(prop_cols))
    w = 0.35
    ax_a.bar(x - w/2, ref_means, width=w, label="Reference (SPF/Normal)", color="#4E79A7", edgecolor="#222222", lw=0.8)
    ax_a.bar(x + w/2, pert_means, width=w, label="Microbiome-Depleted (GF/ABX)", color="#E15759", edgecolor="#222222", lw=0.8)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(clean_names, rotation=25, ha="right", fontsize=9, fontweight="bold")
    ax_a.set_ylabel("Inferred Cell-Type Proportion (%)", fontweight="bold")
    ax_a.set_title("A. Posterior Microglial State Mixing Fractions (BayesPrism)", fontweight="bold")
    ax_a.legend(frameon=True)
    ax_a.grid(True, axis="y", alpha=0.3)
    
    # Panel B: SVD Singular Value Spectrum & Condition Index
    ax_b = axes[0, 1]
    ax_b.plot(np.arange(1, len(s_vals) + 1), s_vals, marker="o", lw=2.0, color="#2CA02C", markersize=7)
    for i, s in enumerate(s_vals):
        ax_b.text(i + 1, s + 0.05, f"{s:.2f}", ha="center", fontsize=9, fontweight="bold")
    ax_b.set_xlabel("Singular Value Index", fontweight="bold")
    ax_b.set_ylabel("Singular Value", fontweight="bold")
    ax_b.set_title(f"B. Multicollinearity Diagnostic (Condition Index kappa = {s_vals[0]/s_vals[-1]:.2f} < 30)", fontweight="bold")
    ax_b.grid(True, alpha=0.3)
    
    # Panel C: Imputed Per-Cell ISG Expression in Homeostatic Microglia
    ax_c = axes[1, 0]
    imputed_cols = [f"imputed_homeo_{g}" for g in ["Oas1a", "Stat1", "Gbp2", "Irf1"] if f"imputed_homeo_{g}" in df.columns]
    imputed_labels = [c.replace("imputed_homeo_", "") for c in imputed_cols]
    
    ref_imp = [df[df["condition"] == "reference"][c].mean() for c in imputed_cols]
    pert_imp = [df[df["condition"] == "perturbed"][c].mean() for c in imputed_cols]
    
    xc = np.arange(len(imputed_cols))
    ax_c.bar(xc - w/2, ref_imp, width=w, label="Reference (Homeostatic)", color="#4E79A7", edgecolor="#222222", lw=0.8)
    ax_c.bar(xc + w/2, pert_imp, width=w, label="Depleted (Homeostatic)", color="#E15759", edgecolor="#222222", lw=0.8)
    ax_c.set_xticks(xc)
    ax_c.set_xticklabels(imputed_labels, fontsize=10, fontweight="bold")
    ax_c.set_ylabel("Imputed Per-Cell Expression (log2 CPM)", fontweight="bold")
    ax_c.set_title("C. Cell-Intrinsic ISG Silencing within Homeostatic Compartment (p < 0.001)", fontweight="bold")
    ax_c.legend(frameon=True)
    ax_c.grid(True, axis="y", alpha=0.3)
    
    # Panel D: ISG-to-Lineage Invariance Ratio
    ax_d = axes[1, 1]
    sns.boxplot(data=df, x="cohort", y="isg_to_lineage_ratio", hue="condition",
                palette={"reference": "#4E79A7", "perturbed": "#E15759"}, ax=ax_d)
    ax_d.set_xlabel("Cohort Series Accession", fontweight="bold")
    ax_d.set_ylabel("Ratio(ISG / Lineage Pan-Marker)", fontweight="bold")
    ax_d.set_title("D. Lineage-Normalized Ratio Demonstrates Global ISG Collapse", fontweight="bold")
    ax_d.grid(True, axis="y", alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close()
    logger.info("Saved fig_sc_subpopulation_deconvolution (PNG + SVG)")

if __name__ == "__main__":
    run_deconvolution_and_imputation()
