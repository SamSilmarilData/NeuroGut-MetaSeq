#!/usr/bin/env python3
"""
scripts/03c_cohort_phenotyping.py
Cross-cohort phenotypic comparison, effect size correlation, and publication graphics.
Outputs:
  - Summary metrics: results/de_results/cohort_phenotypic_summary.csv
  - 4x Individual Volcano Plots: results/de_results/figures/fig_volcano_<cohort>.png
  - LFC Correlation Heatmap: results/de_results/figures/fig_lfc_correlation_heatmap.png
  - Landmark Marker Effect Size Chart: results/de_results/figures/fig_marker_effect_sizes_by_cohort.png
  - DEG Set Overlap Chart: results/de_results/figures/fig_cohort_deg_overlap.png
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Cohort-Phenotyping")

COHORT_META = {
    "GSE107925": {"label": "Lifelong Germ-Free (GF)", "citation": "Thion et al., Cell 2018"},
    "GSE108045": {"label": "Acute Antibiotic Depletion (ABX)", "citation": "Thion et al., Cell 2018"},
    "GSE266602": {"label": "Microbiome Depletion Baseline", "citation": "Wang et al., 2024"},
    "GSE186210": {"label": "Dietary Fiber Starvation (Zero Fiber)", "citation": "Matt et al., J Neurosci 2023"}
}

KEY_MARKERS = [
    # Homeostatic / Anti-inflammatory
    "Tmem119", "Cx3cr1", "P2ry12", "Sall1", "Mertk", "Tsc22d3", "Ffar2",
    # Activation / Stress / Inflammatory
    "Tnf", "Il1b", "Il6", "Ccl2", "Nfkb1", "Ddit4", "Fos"
]

def render_volcano_plot(df: pd.DataFrame, cohort_id: str, out_path: str):
    """Render a publication-grade Volcano Plot for a single cohort."""
    plt.figure(figsize=(8.5, 6))
    sns.set_theme(style="whitegrid", font="DejaVu Sans")

    # Drop zero/missing p-values for plotting safety
    plot_df = df.dropna(subset=["log2FoldChange", "padj"]).copy()
    plot_df["neg_log10_padj"] = -np.log10(np.clip(plot_df["padj"], 1e-50, 1.0))
    plot_df["status"] = "Non-Significant"

    sig_up = (plot_df["padj"] < 0.05) & (plot_df["log2FoldChange"] >= 0.5)
    sig_down = (plot_df["padj"] < 0.05) & (plot_df["log2FoldChange"] <= -0.5)

    plot_df.loc[sig_up, "status"] = "Upregulated"
    plot_df.loc[sig_down, "status"] = "Downregulated"

    colors = {"Upregulated": "#d73027", "Downregulated": "#4575b4", "Non-Significant": "#cccccc"}

    # Scatter points
    for status, color in colors.items():
        subset = plot_df[plot_df["status"] == status]
        alpha = 0.4 if status == "Non-Significant" else 0.8
        size = 12 if status == "Non-Significant" else 24
        plt.scatter(subset["log2FoldChange"], subset["neg_log10_padj"], c=color, label=f"{status} ({len(subset):,})",
                    alpha=alpha, s=size, edgecolors="none")

    # Threshold guidelines
    plt.axvline(0.5, color="#555555", linestyle="--", linewidth=1.0, alpha=0.7)
    plt.axvline(-0.5, color="#555555", linestyle="--", linewidth=1.0, alpha=0.7)
    plt.axhline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=1.0, alpha=0.7)

    # Label top landmark genes and significant outliers
    top_candidates = plot_df[(plot_df["status"] != "Non-Significant")].sort_values("padj").head(10)["gene_symbol"].tolist()
    # Also add landmark markers if present and reasonably significant
    extra_markers = [g for g in KEY_MARKERS if g in plot_df["gene_symbol"].values and plot_df.loc[plot_df["gene_symbol"] == g, "padj"].values[0] < 0.15]
    label_genes = list(dict.fromkeys(top_candidates + extra_markers))[:15]

    max_y = plot_df["neg_log10_padj"].max()
    plt.ylim(-0.5, max_y * 1.25 + 1.5)

    for gene in label_genes:
        row = plot_df[plot_df["gene_symbol"] == gene]
        if not row.empty:
            x = row["log2FoldChange"].values[0]
            y = row["neg_log10_padj"].values[0]
            plt.annotate(
                gene, (x, y),
                xytext=(x + (0.35 if x >= 0 else -0.35), min(y + 0.6, max_y * 1.15)),
                fontsize=8.5, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#666666", alpha=0.85),
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.1", color="#333333", lw=0.7)
            )

    info = COHORT_META.get(cohort_id, {})
    plt.title(f"{cohort_id}: {info.get('label', '')}\nDifferential Expression Landscape ({info.get('citation', '')})",
              fontsize=12, fontweight="bold", pad=12)
    plt.xlabel(r"$\log_2$ Fold Change (Perturbed vs. Reference)", fontsize=11, fontweight="bold")
    plt.ylabel(r"$-\log_{10}(\mathrm{Adjusted\ } p\mathrm{-value})$", fontsize=11, fontweight="bold")
    plt.legend(loc="upper right", frameon=True, framealpha=0.9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    logger.info(f"Saved volcano plot: {out_path}")

def render_lfc_correlation_heatmap(deg_dict: dict, out_path: str):
    """Compute and plot Spearman rank correlation of log2FC across cohorts."""
    # Find common genes across all 4 cohorts
    common_genes = None
    for cid, df in deg_dict.items():
        genes = set(df["gene_symbol"].dropna())
        common_genes = genes if common_genes is None else common_genes.intersection(genes)

    logger.info(f"Common expressed genes for correlation: {len(common_genes):,}")
    common_list = sorted(list(common_genes))

    lfc_matrix = pd.DataFrame(index=common_list)
    for cid, df in deg_dict.items():
        sub = df.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")
        lfc_matrix[cid] = sub.loc[common_list, "log2FoldChange"]

    lfc_matrix = lfc_matrix.dropna()
    corr = lfc_matrix.corr(method="spearman")

    plt.figure(figsize=(7, 6))
    labels = [f"{cid}\n({COHORT_META[cid]['label'][:18]}...)" for cid in corr.columns]
    sns.heatmap(corr, annot=True, fmt=".3f", cmap="coolwarm", vmin=-0.2, vmax=0.8,
                xticklabels=labels, yticklabels=labels, cbar_kws={"label": "Spearman Correlation (ρ)"},
                linewidths=0.5, linecolor="white", annot_kws={"fontweight": "bold", "fontsize": 11})
    plt.title(f"Cross-Study Transcriptome Effect Size Correlation\n({len(lfc_matrix):,} Common Expressed Genes)",
              fontsize=12, fontweight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    logger.info(f"Saved correlation heatmap: {out_path}")

def render_marker_effect_sizes(deg_dict: dict, out_path: str):
    """Compare effect sizes of landmark homeostatic and inflammatory genes."""
    records = []
    for gene in KEY_MARKERS:
        for cid, df in deg_dict.items():
            row = df[df["gene_symbol"] == gene]
            if not row.empty:
                lfc = row["log2FoldChange"].values[0]
                padj = row["padj"].values[0]
                records.append({
                    "gene_symbol": gene,
                    "cohort": cid,
                    "log2FoldChange": lfc,
                    "padj": padj,
                    "significant": padj < 0.05
                })

    m_df = pd.DataFrame(records)
    if m_df.empty:
        return

    plt.figure(figsize=(12, 6))
    sns.set_theme(style="whitegrid", font="DejaVu Sans")
    palette = sns.color_palette("Set2", n_colors=len(deg_dict))

    ax = sns.barplot(x="gene_symbol", y="log2FoldChange", hue="cohort", data=m_df, palette=palette, edgecolor="black", linewidth=0.6)
    plt.axhline(0.0, color="black", linestyle="-", linewidth=1.0)
    plt.axhline(0.5, color="#d9534f", linestyle="--", linewidth=0.8, alpha=0.7)
    plt.axhline(-0.5, color="#4575b4", linestyle="--", linewidth=0.8, alpha=0.7)

    plt.title("Landmark Microglial Marker Effect Sizes Across Discovery Paradigms", fontsize=13, fontweight="bold", pad=12)
    plt.xlabel("Gene Symbol", fontsize=11, fontweight="bold")
    plt.ylabel(r"$\log_2$ Fold Change (Perturbed vs. Reference)", fontsize=11, fontweight="bold")
    plt.xticks(rotation=45, ha="right", fontweight="bold")
    plt.legend(title="Cohort", loc="lower right", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    logger.info(f"Saved marker effect sizes chart: {out_path}")

def render_deg_overlap_chart(deg_dict: dict, out_path: str):
    """Render summary DEG count comparisons across cohorts."""
    summary_data = []
    for cid, df in deg_dict.items():
        sig_up = ((df["padj"] < 0.05) & (df["log2FoldChange"] >= 0.5)).sum()
        sig_down = ((df["padj"] < 0.05) & (df["log2FoldChange"] <= -0.5)).sum()
        summary_data.append({"cohort": cid, "Direction": "Upregulated", "Count": sig_up})
        summary_data.append({"cohort": cid, "Direction": "Downregulated", "Count": sig_down})

    s_df = pd.DataFrame(summary_data)
    plt.figure(figsize=(8, 5))
    palette = {"Upregulated": "#d73027", "Downregulated": "#4575b4"}
    ax = sns.barplot(x="cohort", y="Count", hue="Direction", data=s_df, palette=palette, edgecolor="black", linewidth=0.8)

    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(f"{int(height)}",
                        (p.get_x() + p.get_width() / 2., height),
                        ha="center", va="bottom", fontsize=10, fontweight="bold", xytext=(0, 3),
                        textcoords="offset points")

    plt.title("Differentially Expressed Gene Counts Across Cohorts (padj < 0.05, |log2FC| >= 0.5)",
              fontsize=12, fontweight="bold", pad=12)
    plt.xlabel("Cohort ID", fontsize=11, fontweight="bold")
    plt.ylabel("Number of Significant DEGs", fontsize=11, fontweight="bold")
    plt.legend(title="Direction", loc="upper right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    logger.info(f"Saved DEG overlap chart: {out_path}")

def main():
    de_dir = "results/de_results"
    fig_dir = os.path.join(de_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    cohorts = ["GSE107925", "GSE108045", "GSE266602", "GSE186210"]
    deg_dict = {}

    summary_rows = []
    for cid in cohorts:
        deg_path = os.path.join(de_dir, f"{cid}_deg.csv")
        if not os.path.exists(deg_path):
            logger.warning(f"DEG file missing for {cid}: {deg_path}")
            continue
        df = pd.read_csv(deg_path)
        deg_dict[cid] = df

        sig_up = df[(df["padj"] < 0.05) & (df["log2FoldChange"] >= 0.5)]
        sig_down = df[(df["padj"] < 0.05) & (df["log2FoldChange"] <= -0.5)]

        top_up = ", ".join(sig_up.head(5)["gene_symbol"].tolist()) if not sig_up.empty else "None"
        top_down = ", ".join(sig_down.head(5)["gene_symbol"].tolist()) if not sig_down.empty else "None"

        summary_rows.append({
            "cohort": cid,
            "perturbation_paradigm": COHORT_META[cid]["label"],
            "total_expressed_genes": len(df),
            "sig_upregulated": len(sig_up),
            "sig_downregulated": len(sig_down),
            "total_significant_degs": len(sig_up) + len(sig_down),
            "top_upregulated_genes": top_up,
            "top_downregulated_genes": top_down
        })

        # Render Volcano Plot
        v_path = os.path.join(fig_dir, f"fig_volcano_{cid}.png")
        render_volcano_plot(df, cid, v_path)

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(de_dir, "cohort_phenotypic_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    logger.info(f"[SUCCESS] Saved cohort phenotypic summary: {summary_csv}")
    print("\n=== Cohort Differential Expression Phenotypic Summary ===")
    print(summary_df.to_string(index=False))

    # Cross-Cohort Comparative Figures
    render_lfc_correlation_heatmap(deg_dict, os.path.join(fig_dir, "fig_lfc_correlation_heatmap.png"))
    render_marker_effect_sizes(deg_dict, os.path.join(fig_dir, "fig_marker_effect_sizes_by_cohort.png"))
    render_deg_overlap_chart(deg_dict, os.path.join(fig_dir, "fig_cohort_deg_overlap.png"))

    logger.info("[SUCCESS] All Horizon 2 phenotypic diagnostic graphics generated successfully.")

if __name__ == "__main__":
    main()
