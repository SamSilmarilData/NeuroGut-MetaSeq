#!/usr/bin/env python3
"""
scripts/06_generate_figures.py
Generates publication-quality figures (300 DPI PNG and vector SVG):
- Figure 1: PCA plot across all cohorts and conditions
- Figure 2: Consensus Meta-Analysis Volcano Plot
- Figure 3: Forest Plots of landmark microglial and inflammatory genes
- Figure 4: Hierarchically clustered heatmap of top consensus DEGs
- Figure 5: KEGG and GO pathway enrichment dotplot
Outputs to results/figures/ and copies to docs/assets/ for the web paper.
"""

import os
import sys
import glob
import logging
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Generate-Figures")

# Set publication style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "axes.edgecolor": "#2c3e50",
    "axes.linewidth": 1.2,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "figure.titlesize": 14,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 300
})

def save_fig(fig, base_name, out_dir, web_dir):
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(web_dir, exist_ok=True)
    for ext in ["png", "svg"]:
        fpath = os.path.join(out_dir, f"{base_name}.{ext}")
        fig.savefig(fpath, dpi=300, bbox_inches="tight")
        # Copy to web paper assets
        web_path = os.path.join(web_dir, f"{base_name}.{ext}")
        fig.savefig(web_path, dpi=300, bbox_inches="tight")
    logger.info(f"[FIGURE] Saved {base_name}.png and {base_name}.svg")

def plot_fig1_pca(out_dir, web_dir):
    """Figure 1: PCA across samples and cohorts."""
    count_files = glob.glob("data/processed/*_counts.csv")
    meta_files = glob.glob("data/metadata/*_metadata.csv")
    if not count_files or not meta_files:
        return

    # Merge samples
    all_counts = []
    meta_list = []
    for mf in meta_files:
        meta_list.append(pd.read_csv(mf))
    meta_df = pd.concat(meta_list, ignore_index=True).set_index("sample_id")

    for cf in count_files:
        df = pd.read_csv(cf).set_index("gene_symbol")
        sample_cols = [c for c in df.columns if c in meta_df.index]
        all_counts.append(df[sample_cols])

    combined_counts = pd.concat(all_counts, axis=1, join="inner")
    # Log-transform normalized
    log_c = np.log2(combined_counts + 1).T
    # PCA
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(log_c)
    pca_df = pd.DataFrame(pca_coords, index=log_c.index, columns=["PC1", "PC2"])
    pca_df = pca_df.join(meta_df)

    fig, ax = plt.subplots(figsize=(7, 5))
    palette = {"reference": "#2b5c8f", "perturbed": "#d9534f"}
    markers = {"GSE107925": "o", "GSE108045": "s", "GSE266602": "^"}

    for (cohort, cond), grp in pca_df.groupby(["cohort", "condition"]):
        color = palette.get(cond, "#555555")
        marker = markers.get(cohort, "o")
        ax.scatter(grp["PC1"], grp["PC2"], label=f"{cohort} ({cond})", s=80, alpha=0.85, edgecolors="#111111", c=color, marker=marker)

    ax.set_xlabel(f"Principal Component 1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
    ax.set_ylabel(f"Principal Component 2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
    ax.set_title("Cross-Cohort Principal Component Analysis (Microglial Transcriptomes)")
    ax.grid(True)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    save_fig(fig, "fig1_cohort_pca", out_dir, web_dir)
    plt.close(fig)

def plot_fig2_volcano(meta_df, out_dir, web_dir):
    """Figure 2: Consensus Meta-Analysis Volcano Plot."""
    fig, ax = plt.subplots(figsize=(8, 6))

    x = meta_df["meta_log2fc"].values
    y = -np.log10(np.clip(meta_df["fdr_fisher"].values, 1e-30, 1.0))

    # Color classification
    is_up = (meta_df["fdr_fisher"] < 0.05) & (meta_df["meta_log2fc"] >= 0.5)
    is_down = (meta_df["fdr_fisher"] < 0.05) & (meta_df["meta_log2fc"] <= -0.5)
    is_ns = ~is_up & ~is_down

    ax.scatter(x[is_ns], y[is_ns], c="#95a5a6", alpha=0.45, s=25, label=f"Not Significant (n={is_ns.sum()})")
    ax.scatter(x[is_up], y[is_up], c="#e74c3c", alpha=0.9, s=45, label=f"Consensus Upregulated (n={is_up.sum()})")
    ax.scatter(x[is_down], y[is_down], c="#3498db", alpha=0.9, s=45, label=f"Consensus Downregulated (n={is_down.sum()})")

    ax.axvline(0.5, color="#7f8c8d", linestyle="--", alpha=0.7)
    ax.axvline(-0.5, color="#7f8c8d", linestyle="--", alpha=0.7)
    ax.axhline(-np.log10(0.05), color="#7f8c8d", linestyle="--", alpha=0.7)

    # Label top landmark genes
    highlight_genes = ["Tnf", "Il1b", "Nfkb1", "Ccl2", "Cx3cr1", "Tmem119", "P2ry12", "Ffar2", "Trem2", "Apoe", "Stat1", "Nos2"]
    for _, row in meta_df[meta_df["gene_symbol"].isin(highlight_genes)].iterrows():
        gx = row["meta_log2fc"]
        gy = -np.log10(max(row["fdr_fisher"], 1e-30))
        ax.annotate(row["gene_symbol"], xy=(gx, gy), xytext=(gx + 0.08, gy + 0.3),
                    fontsize=9, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="#2c3e50", lw=0.8))

    ax.set_xlabel("Pooled Effect Size (DerSimonian-Laird log2FoldChange)")
    ax.set_ylabel("-log10(Fisher FDR Adjusted p-value)")
    ax.set_title("Cross-Cohort Meta-Analysis Volcano Plot (Microglia Perturbation)")
    ax.legend(frameon=True, loc="upper right")
    ax.grid(True)
    save_fig(fig, "fig2_meta_volcano", out_dir, web_dir)
    plt.close(fig)

def plot_fig3_forest(out_dir, web_dir):
    """Figure 3: Forest plots of landmark inflammatory and microglial genes."""
    deg_files = sorted(glob.glob("results/de_results/*_deg.csv"))
    meta_df = pd.read_csv("results/meta_results/microglia_meta_analysis_summary.csv").set_index("gene_symbol")

    target_genes = ["Tnf", "Il1b", "Nfkb1", "Ccl2", "Cx3cr1", "Tmem119", "Ffar2", "Trem2"]
    cohort_data = {}
    for f in deg_files:
        cname = os.path.basename(f).replace("_deg.csv", "")
        cohort_data[cname] = pd.read_csv(f).set_index("gene_symbol")

    fig, axes = plt.subplots(4, 2, figsize=(11, 12), sharex=True)
    axes = axes.flatten()

    for idx, gene in enumerate(target_genes):
        ax = axes[idx]
        y_pos = []
        labels = []
        thetas = []
        lows = []
        highs = []

        # Cohort points
        for c_idx, (cname, cdf) in enumerate(cohort_data.items()):
            if gene in cdf.index:
                row = cdf.loc[gene]
                lfc = float(row["log2FoldChange"])
                se = float(row["lfcSE"])
                y_pos.append(c_idx)
                labels.append(cname)
                thetas.append(lfc)
                lows.append(lfc - 1.96 * se)
                highs.append(lfc + 1.96 * se)

        # Meta summary
        if gene in meta_df.index:
            mrow = meta_df.loc[gene]
            m_lfc = float(mrow["meta_log2fc"])
            m_low = float(mrow["ci_lower"])
            m_high = float(mrow["ci_upper"])
            y_pos.append(len(labels))
            labels.append("Pooled Meta (RE)")
            thetas.append(m_lfc)
            lows.append(m_low)
            highs.append(m_high)

        # Plot error bars
        for yp, th, lo, hi, lbl in zip(y_pos, thetas, lows, highs, labels):
            is_meta = "Pooled" in lbl
            color = "#e74c3c" if th > 0 else "#3498db"
            if is_meta:
                ax.plot([lo, hi], [yp, yp], color="#2c3e50", lw=3.0)
                ax.scatter(th, yp, marker="D", color="#2c3e50", s=90, zorder=5)
            else:
                ax.plot([lo, hi], [yp, yp], color=color, lw=1.8)
                ax.scatter(th, yp, marker="o", color=color, s=50, zorder=5)

        ax.axvline(0, color="#7f8c8d", linestyle="--", lw=1.0)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.set_title(f"Gene: {gene}", fontweight="bold", loc="left")
        ax.grid(True, axis="x")
        ax.invert_yaxis()

    fig.text(0.5, 0.04, "Effect Size (log2FoldChange ± 95% Confidence Interval)", ha="center", fontsize=12)
    fig.suptitle("Forest Plots: Cross-Study Reproducibility of Landmark Microglial Genes", fontsize=14, y=0.99)
    plt.tight_layout(rect=[0, 0.05, 1, 0.97])
    save_fig(fig, "fig3_forest_plots", out_dir, web_dir)
    plt.close(fig)

def plot_fig4_heatmap(meta_df, out_dir, web_dir):
    """Figure 4: Clustered Heatmap of Top Consensus Genes."""
    top_genes = meta_df.head(28)["gene_symbol"].tolist()
    deg_files = sorted(glob.glob("results/de_results/*_deg.csv"))

    heat_data = pd.DataFrame(index=top_genes)
    for f in deg_files:
        cname = os.path.basename(f).replace("_deg.csv", "")
        df = pd.read_csv(f).set_index("gene_symbol")
        heat_data[cname] = df.loc[top_genes, "log2FoldChange"]

    heat_data["Pooled Meta"] = meta_df.set_index("gene_symbol").loc[top_genes, "meta_log2fc"]

    fig, ax = plt.subplots(figsize=(8, 10))
    sns.heatmap(heat_data, cmap="vlag", center=0, annot=True, fmt=".2f", cbar_kws={"label": "log2FoldChange"},
                linewidths=0.5, linecolor="#eeeeee", ax=ax)
    ax.set_title("Top Consensus Meta-Analysis Microglial Genes Across Cohorts")
    ax.set_ylabel("Gene Symbol")
    save_fig(fig, "fig4_consensus_heatmap", out_dir, web_dir)
    plt.close(fig)

def plot_fig5_pathways(out_dir, web_dir):
    """Figure 5: Pathway Enrichment Dotplot."""
    pathway_file = "results/pathways/pathway_summary.csv"
    if not os.path.exists(pathway_file):
        return

    df = pd.read_csv(pathway_file)
    df = df.sort_values(by="fdr", ascending=False)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    y_pos = np.arange(len(df))
    colors = ["#e74c3c" if d == "Upregulated" else "#3498db" for d in df["direction"]]

    bars = ax.barh(y_pos, -np.log10(df["fdr"]), color=colors, alpha=0.85, edgecolor="#2c3e50", height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["pathway_name"], fontsize=10)
    ax.set_xlabel("-log10(FDR Adjusted Enrichment p-value)")
    ax.set_title("Enriched Inflammatory & Microglial Functional Pathways")
    ax.axvline(-np.log10(0.05), color="#7f8c8d", linestyle="--", label="FDR = 0.05")
    ax.grid(True, axis="x")

    # Add gene count badges
    for idx, row in enumerate(df.iterrows()):
        r = row[1]
        ax.text(-np.log10(r["fdr"]) + 0.1, idx, f"{r['overlap_count']} genes ({r['direction']})", va="center", fontsize=9, fontweight="bold")

    ax.legend(loc="lower right")
    save_fig(fig, "fig5_pathway_enrichment", out_dir, web_dir)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Generate publication figures")
    parser.add_argument("--demo", action="store_true", help="Generate figures from demo data")
    args = parser.parse_args()

    out_dir = "results/figures"
    web_dir = "docs/assets"
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(web_dir, exist_ok=True)

    meta_file = "results/meta_results/microglia_meta_analysis_summary.csv"
    if not os.path.exists(meta_file):
        logger.error(f"Meta-analysis summary not found at {meta_file}. Run step 04 first.")
        sys.exit(1)

    meta_df = pd.read_csv(meta_file)
    logger.info("Generating publication figures...")

    plot_fig1_pca(out_dir, web_dir)
    plot_fig2_volcano(meta_df, out_dir, web_dir)
    plot_fig3_forest(out_dir, web_dir)
    plot_fig4_heatmap(meta_df, out_dir, web_dir)
    plot_fig5_pathways(out_dir, web_dir)

    logger.info(f"[SUCCESS] All figures generated in {out_dir} and copied to {web_dir}")

if __name__ == "__main__":
    main()
