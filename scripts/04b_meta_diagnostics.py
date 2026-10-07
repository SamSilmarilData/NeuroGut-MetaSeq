#!/usr/bin/env python3
"""
scripts/04b_meta_diagnostics.py
Generates publication-grade diagnostic graphics (300 DPI) for Horizon 3:
1. Whole-Transcriptome Meta-Analysis Volcano Plot (fig_meta_volcano.png)
2. Multi-Study Forest Plots for Top Consensus Hits & Tracked Blooms (fig_forest_plots_top.png)
3. Leave-One-Out (LOO) Sensitivity Stability Scatter Plot (fig_loo_stability.png)
4. Heterogeneity Distribution & Directional Concordance Breakdown (fig_heterogeneity_distribution.png)
5. Cross-Study Normalized Expression Heatmap across 60 biological samples (fig_consensus_heatmap.png)
"""

import os
import glob
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Meta-Diagnostics")

# Global style settings
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "axes.edgecolor": "#333333",
    "axes.linewidth": 1.1,
    "grid.color": "#E5E5E5",
    "grid.linestyle": "--",
    "grid.linewidth": 0.6,
    "figure.titlesize": 14,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9.5
})

def plot_meta_volcano(meta_df: pd.DataFrame, out_path: str):
    """Whole-transcriptome meta-analysis volcano plot colored by heterogeneity tier."""
    logger.info("Generating Meta-Analysis Volcano Plot...")
    fig, ax = plt.subplots(figsize=(10, 7.5), dpi=300)

    # Clip p-values to prevent inf in -log10
    pvals = np.clip(meta_df["p_random_effects"].values, 1e-50, 1.0)
    neg_log10_p = -np.log10(pvals)
    lfcs = meta_df["meta_log2fc"].values

    # Determine p-value threshold corresponding to FDR_RE = 0.05
    sig_re = meta_df[meta_df["fdr_random_effects"] < 0.05]
    fdr_thresh_p = sig_re["p_random_effects"].max() if not sig_re.empty else 0.001
    neg_log10_thresh = -np.log10(fdr_thresh_p)

    # Background non-significant points
    non_sig_mask = (meta_df["fdr_random_effects"] >= 0.05) | (meta_df["meta_log2fc"].abs() < 0.5)
    ax.scatter(lfcs[non_sig_mask], neg_log10_p[non_sig_mask],
               c="#D0D4DC", alpha=0.35, s=16, edgecolors="none", rasterized=True, label="Non-Significant")

    # Significant points colored by heterogeneity tier
    sig_low = (meta_df["significance_flag"]) & (meta_df["heterogeneity_tier"] == "Low (<25%)")
    sig_mod = (meta_df["significance_flag"]) & (meta_df["heterogeneity_tier"] == "Moderate (25-75%)")
    sig_high = (meta_df["significance_flag"]) & (meta_df["heterogeneity_tier"] == "High (>75%)")

    ax.scatter(lfcs[sig_low], neg_log10_p[sig_low],
               c="#0072B2", alpha=0.9, s=48, edgecolors="#003B6F", linewidth=0.8,
               label=f"Low I² (<25%, Core Invariant: n={sig_low.sum()})", zorder=4)

    ax.scatter(lfcs[sig_mod], neg_log10_p[sig_mod],
               c="#E69F00", alpha=0.9, s=48, edgecolors="#8C5C00", linewidth=0.8,
               label=f"Moderate I² (25-75%: n={sig_mod.sum()})", zorder=4)

    ax.scatter(lfcs[sig_high], neg_log10_p[sig_high],
               c="#CC79A7", alpha=0.9, s=48, edgecolors="#662244", linewidth=0.8,
               label=f"High I² (>75%, Perturbation-Specific: n={sig_high.sum()})", zorder=4)

    # Threshold guidelines
    ax.axhline(neg_log10_thresh, color="#D55E00", linestyle="--", linewidth=1.1, alpha=0.8,
               label=f"FDR = 0.05 (p = {fdr_thresh_p:.1e})")
    ax.axvline(0.5, color="#555555", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.axvline(-0.5, color="#555555", linestyle=":", linewidth=1.0, alpha=0.7)

    # Landmark and top consensus gene annotations
    labels_to_show = ["Fosb", "Slfn2", "Clu", "Neat1", "Llgl2", "Sap30", "Card6",
                      "Tsc22d3", "Ddit4", "Plin3", "Tnf", "Sall1", "Ffar2", "Cx3cr1", "Tmem119"]

    annot_df = meta_df[meta_df["gene_symbol"].isin(labels_to_show)]
    for _, row in annot_df.iterrows():
        g = row["gene_symbol"]
        x = row["meta_log2fc"]
        p_val = max(row["p_random_effects"], 1e-50)
        y = -np.log10(p_val)

        # Highlight point
        ax.scatter(x, y, s=75, facecolors="none", edgecolors="#111111", linewidths=1.5, zorder=6)

        # Adjust text offset
        dx = 0.25 if x >= 0 else -0.3
        dy = 0.8
        if g in ["Tsc22d3", "Ddit4"]:
            dx = -0.5
            dy = -1.2
        elif g in ["Fosb"]:
            dx = 0.3
            dy = 0.2
        elif g in ["Slfn2"]:
            dx = -0.5
            dy = 0.5
        elif g in ["Tnf"]:
            dx = 0.3
            dy = 0.5

        ax.annotate(g, xy=(x, y), xytext=(x + dx, y + dy),
                    fontsize=9.5, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8, shrinkA=3, shrinkB=3),
                    zorder=7)

    ax.set_xlabel("Pooled Effect Size (REML / Hartung-Knapp log₂ Fold Change)")
    ax.set_ylabel("-log₁₀ (Random Effects p-value)")
    ax.set_title("Cross-Study Microglial Meta-Analysis Landscape (N=60, 4 Cohorts, 23,096 Genes)",
                 pad=12, fontweight="bold")
    ax.legend(frameon=True, facecolor="white", edgecolor="#DDDDDD", loc="upper left")
    ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    svg_path = out_path.replace(".png", ".svg")
    plt.savefig(svg_path, format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and {svg_path}")

def plot_forest_plots(meta_df: pd.DataFrame, cohort_data: dict, out_path: str):
    """Multi-study forest plots for top consensus hits and tracked blooms."""
    logger.info("Generating Multi-Study Forest Plots...")

    # Selected genes (12 genes in 4x3 grid)
    selected_genes = [
        # Consensus Hits
        "Fosb", "Clu", "Llgl2", "Neat1",
        "Slfn2", "Sap30", "Card6", "Tnf",
        # Tracked Biological Blooms & Landmarks
        "Tsc22d3", "Ddit4", "Plin3", "Sall1"
    ]

    fig, axes = plt.subplots(4, 3, figsize=(15, 14), dpi=300)
    axes = axes.flatten()

    cohort_order = ["GSE107925", "GSE108045", "GSE186210", "GSE266602"]
    cohort_labels = ["GSE107925 (GF)", "GSE108045 (ABX)", "GSE186210 (Fiber)", "GSE266602 (Sham)"]

    for idx, gene in enumerate(selected_genes):
        ax = axes[idx]
        if gene not in meta_df["gene_symbol"].values:
            ax.set_visible(False)
            continue

        gene_meta = meta_df[meta_df["gene_symbol"] == gene].iloc[0]

        studies = []
        thetas = []
        ci_lows = []
        ci_highs = []

        # Extract study stats
        for cname, clabel in zip(cohort_order, cohort_labels):
            if gene in cohort_data[cname]["df"].index:
                row = cohort_data[cname]["df"].loc[gene]
                lfc = float(row["log2FoldChange"])
                se = float(row["lfcSE"])
                thetas.append(lfc)
                ci_lows.append(lfc - 1.96 * se)
                ci_highs.append(lfc + 1.96 * se)
                studies.append(clabel)

        # Append pooled summary
        studies.append("RE Pooled")
        thetas.append(gene_meta["meta_log2fc"])
        ci_lows.append(gene_meta["ci_lower"])
        ci_highs.append(gene_meta["ci_upper"])

        y_pos = np.arange(len(studies))[::-1]

        # Reference null line
        ax.axvline(0, color="#888888", linestyle="--", linewidth=1.0, alpha=0.8)

        # Plot study effect sizes
        for i in range(len(studies) - 1):
            ax.plot([ci_lows[i], ci_highs[i]], [y_pos[i], y_pos[i]], color="#2B5B84", lw=1.6)
            ax.scatter(thetas[i], y_pos[i], s=55, marker="s", color="#1F77B4", edgecolors="#0D3B66", zorder=3)

        # Plot pooled diamond
        i_pool = len(studies) - 1
        diamond_y = y_pos[i_pool]
        diamond_x = [ci_lows[i_pool], thetas[i_pool], ci_highs[i_pool], thetas[i_pool]]
        diamond_ys = [diamond_y, diamond_y + 0.22, diamond_y, diamond_y - 0.22]
        ax.fill(diamond_x, diamond_ys, color="#D95F02", alpha=0.9, edgecolor="#8B3E00", zorder=4)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(studies, fontsize=8.5)

        # Annotations (I² and p-value)
        i2_val = gene_meta["i2_heterogeneity"]
        p_val = gene_meta["p_random_effects"]
        fdr_val = gene_meta["fdr_random_effects"]
        ax.set_title(f"{gene} (I²={i2_val:.1f}%, FDR={fdr_val:.2e})", fontsize=10.5, fontweight="bold")
        ax.grid(True, axis="x", alpha=0.3)
        ax.set_xlabel("log₂ Fold Change", fontsize=8.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def plot_loo_stability(loo_df: pd.DataFrame, meta_df: pd.DataFrame, out_path: str):
    """Leave-One-Out (LOO) stability scatter across omitted studies."""
    logger.info("Generating LOO Stability Scatter Plot...")

    cohorts = sorted(loo_df["omitted_cohort"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(12, 11), dpi=300)
    axes = axes.flatten()

    meta_dict = meta_df.set_index("gene_symbol")["meta_log2fc"].to_dict()

    for idx, omitted in enumerate(cohorts):
        ax = axes[idx]
        sub = loo_df[loo_df["omitted_cohort"] == omitted].copy()
        sub["full_log2fc"] = sub["gene_symbol"].map(meta_dict)
        sub = sub.dropna(subset=["full_log2fc", "loo_log2fc"])

        x = sub["full_log2fc"].values
        y = sub["loo_log2fc"].values

        # Compute correlations
        r_val = float(np.corrcoef(x, y)[0, 1])
        rho_val = float(pd.Series(x).corr(pd.Series(y), method="spearman"))

        # Scatter
        ax.scatter(x, y, c="#4E79A7", alpha=0.3, s=14, edgecolors="none", rasterized=True)

        # Diagonal unity line
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0], -3),
                max(ax.get_xlim()[1], ax.get_ylim()[1], 3)]
        ax.plot(lims, lims, color="#E15759", linestyle="--", lw=1.2, alpha=0.8, label="Unity (y = x)")

        # Highlight key genes
        key_genes = ["Fosb", "Slfn2", "Clu", "Llgl2", "Tsc22d3", "Ddit4", "Plin3", "Tnf"]
        for g in key_genes:
            if g in sub["gene_symbol"].values:
                grow = sub[sub["gene_symbol"] == g].iloc[0]
                gx, gy = grow["full_log2fc"], grow["loo_log2fc"]
                ax.scatter(gx, gy, color="#E15759", s=45, edgecolors="#661111", zorder=5)
                ax.annotate(g, (gx, gy), textcoords="offset points", xytext=(5, 5),
                            fontsize=8.5, fontweight="bold", zorder=6)

        ax.set_title(f"Omit {omitted} (r = {r_val:.3f}, ρ = {rho_val:.3f})", fontweight="bold", fontsize=11)
        ax.set_xlabel("Full Meta-Analysis log₂FC")
        ax.set_ylabel(f"LOO Pooled log₂FC (without {omitted})")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def plot_heterogeneity_distribution(meta_df: pd.DataFrame, out_path: str):
    """Distribution of Higgins I² and Directional Concordance proportions."""
    logger.info("Generating Heterogeneity & Concordance Distribution Plot...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)

    # Panel A: Higgins I² distribution
    i2_vals = meta_df["i2_heterogeneity"].values
    ax1.hist(i2_vals, bins=40, color="#59A14F", edgecolor="#2B5B24", alpha=0.75, density=True)
    sns.kdeplot(i2_vals, ax=ax1, color="#1B4D14", linewidth=2.0)

    ax1.axvline(25.0, color="#E69F00", linestyle="--", linewidth=1.5, label="Low Threshold (25%)")
    ax1.axvline(75.0, color="#CC79A7", linestyle="--", linewidth=1.5, label="High Threshold (75%)")

    # Percentage breakdown
    p_low = (i2_vals < 25.0).mean() * 100
    p_mod = ((i2_vals >= 25.0) & (i2_vals <= 75.0)).mean() * 100
    p_high = (i2_vals > 75.0).mean() * 100

    ax1.text(0.05, 0.90, f"Low (<25%): {p_low:.1f}%\nMod (25-75%): {p_mod:.1f}%\nHigh (>75%): {p_high:.1f}%",
             transform=ax1.transAxes, fontsize=10, bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8))

    ax1.set_xlabel("Higgins I² Heterogeneity (%)")
    ax1.set_ylabel("Density")
    ax1.set_title("A. Higgins I² Distribution Across Common Transcriptome", fontweight="bold")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)

    # Panel B: Directional Concordance Breakdown
    conc_counts = meta_df["direction_concordance"].value_counts()
    colors = ["#4E79A7", "#F28E2B", "#E15759"]
    bars = ax2.bar(conc_counts.index, conc_counts.values, color=colors, edgecolor="#333333", width=0.55)

    for bar in bars:
        h = bar.get_height()
        pct = (h / len(meta_df)) * 100
        ax2.text(bar.get_x() + bar.get_width()/2., h + 200, f"{h:,}\n({pct:.1f}%)",
                 ha="center", va="bottom", fontsize=9.5, fontweight="bold")

    ax2.set_ylabel("Number of Genes")
    ax2.set_title("B. Directional Concordance Across Cohorts", fontweight="bold")
    ax2.grid(True, axis="y", alpha=0.3)
    ax2.set_ylim(0, max(conc_counts.values) * 1.15)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def plot_consensus_heatmap(meta_df: pd.DataFrame, counts_dir: str, meta_dir: str, out_path: str):
    """Normalized cross-study z-score heatmap of top consensus & landmark genes across 60 samples."""
    logger.info("Generating Cross-Study Consensus Heatmap...")

    # Key landmark & top consensus genes to include
    target_genes = [
        "Fosb", "Clu", "Llgl2", "Neat1", "Ppif", "Slfn2", "Sap30", "Card6",
        "Hes1", "Slc25a19", "Osbpl11", "Zbtb48", "Tnf", "Tsc22d3", "Ddit4",
        "Plin3", "Sall1", "Ffar2", "Cx3cr1", "Tmem119", "Fos"
    ]

    # Load normalized expression (CPM) per cohort
    cohort_files = sorted(glob.glob(os.path.join(counts_dir, "*_counts.csv")))
    all_sample_dfs = []
    sample_annotations = []

    for cfile in cohort_files:
        cname = os.path.basename(cfile).replace("_counts.csv", "")
        counts_df = pd.read_csv(cfile)
        if "gene_symbol" not in counts_df.columns:
            continue

        counts_df = counts_df.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")
        sample_cols = [c for c in counts_df.columns if c != "gene_id"]

        # Compute log2(CPM + 1)
        lib_sizes = counts_df[sample_cols].sum(axis=0)
        cpm = np.log2((counts_df[sample_cols].div(lib_sizes, axis=1) * 1e6) + 1.0)

        # Standardize Z-score per gene within cohort (removes baseline inter-cohort shifts)
        cpm_z = cpm.sub(cpm.mean(axis=1), axis=0).div(cpm.std(axis=1).replace(0, 1), axis=0)

        # Load metadata
        meta_f = os.path.join(meta_dir, f"{cname}_metadata.csv")
        mdf = pd.read_csv(meta_f).set_index("sample_id")

        for s in sample_cols:
            cond = mdf.loc[s, "condition"] if s in mdf.index else "unknown"
            sample_annotations.append({
                "sample_id": s,
                "cohort": cname,
                "condition": cond
            })

        all_sample_dfs.append(cpm_z)

    # Intersect available target genes
    avail_genes = [g for g in target_genes if all(g in df.index for df in all_sample_dfs)]
    logger.info(f"Heatmap target genes available across all cohorts: {len(avail_genes)} / {len(target_genes)}")

    # Combine matrices
    combined_expr = pd.concat([df.loc[avail_genes] for df in all_sample_dfs], axis=1)
    annot_df = pd.DataFrame(sample_annotations).set_index("sample_id")

    # Sort columns by Cohort, then Condition
    annot_df = annot_df.sort_values(by=["cohort", "condition"])
    combined_expr = combined_expr[annot_df.index]

    # Set up colors for metadata
    cohort_colors = {
        "GSE107925": "#4E79A7",
        "GSE108045": "#F28E2B",
        "GSE186210": "#E15759",
        "GSE266602": "#76B7B2"
    }
    condition_colors = {
        "reference": "#59A14F",
        "perturbed": "#E15759"
    }

    col_cohort = annot_df["cohort"].map(cohort_colors)
    col_condition = annot_df["condition"].map(condition_colors)
    col_colors = pd.DataFrame({"Cohort": col_cohort, "Condition": col_condition})

    g = sns.clustermap(combined_expr,
                       col_cluster=False,
                       row_cluster=True,
                       col_colors=col_colors,
                       cmap="vlag",
                       center=0.0,
                       vmin=-2.5,
                       vmax=2.5,
                       figsize=(14, 10),
                       dendrogram_ratio=(0.15, 0.05),
                       cbar_pos=(0.02, 0.8, 0.02, 0.15),
                       cbar_kws={"label": "Within-Cohort Relative Expression (Z-score)"})

    g.ax_heatmap.set_xlabel("Biological Samples (N=60 across 4 Cohorts)", fontweight="bold", fontsize=11)
    g.ax_heatmap.set_ylabel("Consensus & Landmark Microglial Genes", fontweight="bold", fontsize=11)
    g.ax_heatmap.set_xticklabels([])

    # Add custom legend for Cohort and Condition
    legend_patches = [
        mpatches.Patch(color=c, label=f"Cohort: {k}") for k, c in cohort_colors.items()
    ] + [
        mpatches.Patch(color=c, label=f"Condition: {k.capitalize()}") for k, c in condition_colors.items()
    ]
    g.ax_heatmap.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1.0), loc="upper left", frameon=True)

    g.savefig(out_path, dpi=300)
    g.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def main():
    meta_path = "results/meta_results/microglia_meta_analysis_summary.csv"
    loo_path = "results/meta_results/microglia_meta_analysis_loo.csv"
    fig_dir = "results/meta_results/figures"
    os.makedirs(fig_dir, exist_ok=True)

    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"{meta_path} does not exist. Run scripts/04_meta_analysis.py first.")

    meta_df = pd.read_csv(meta_path)
    loo_df = pd.read_csv(loo_path) if os.path.exists(loo_path) else pd.DataFrame()

    import sys
    import importlib.util
    script_path = os.path.join(os.path.dirname(__file__), "04_meta_analysis.py")
    spec = importlib.util.spec_from_file_location("meta_module", script_path)
    meta_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(meta_mod)
    cohort_data = meta_mod.load_cohort_data()

    # 1. Volcano Plot
    plot_meta_volcano(meta_df, os.path.join(fig_dir, "fig_meta_volcano.png"))

    # 2. Forest Plots
    plot_forest_plots(meta_df, cohort_data, os.path.join(fig_dir, "fig_forest_plots_top.png"))

    # 3. LOO Stability Plot
    if not loo_df.empty:
        plot_loo_stability(loo_df, meta_df, os.path.join(fig_dir, "fig_loo_stability.png"))

    # 4. Heterogeneity Distribution Plot
    plot_heterogeneity_distribution(meta_df, os.path.join(fig_dir, "fig_heterogeneity_distribution.png"))

    # 5. Consensus Heatmap
    plot_consensus_heatmap(meta_df, "data/processed", "data/metadata", os.path.join(fig_dir, "fig_consensus_heatmap.png"))

    logger.info("All Horizon 3 diagnostic figures successfully generated!")

if __name__ == "__main__":
    main()
