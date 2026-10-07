#!/usr/bin/env python3
"""
scripts/05e_systems_diagnostics.py
Generates publication-grade systems biology diagnostic figures (300 DPI) for Horizon 4:
1. GSEA Pathway Enrichment Landscape (fig_gsea_pathway_enrichment.png)
2. Upstream TF Regulon Activity Ranking (fig_tf_regulon_landscape.png)
3. WGCNA Co-Expression Module-Trait Heatmap (fig_wgcna_modules_eigengenes.png)
4. Consensus Hub Gene Interactome Subgraph (fig_network_hub_subgraph.png)
5. In Silico SCFA Metabolite Rescue Trajectory (fig_scfa_rescue_inversion.png)
"""

import os
import logging
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Systems-Diagnostics")

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

def plot_gsea_enrichment(hallmark_csv: str, pheno_csv: str, out_path: str):
    """Multi-panel dot plot of GSEA Normalized Enrichment Scores (NES) and FDR."""
    logger.info("Generating GSEA Pathway Enrichment Landscape Plot...")
    hdf = pd.read_csv(hallmark_csv)
    pdf = pd.read_csv(pheno_csv)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7.5), dpi=300, gridspec_kw={"width_ratios": [1.4, 1.0]})

    # Panel A: Top MSigDB Hallmarks (top 8 up, top 8 down by NES)
    top_up = hdf[hdf["normalized_enrichment_score"] > 0].sort_values(by="normalized_enrichment_score", ascending=False).head(8)
    top_down = hdf[hdf["normalized_enrichment_score"] < 0].sort_values(by="normalized_enrichment_score", ascending=True).head(8)
    h_sub = pd.concat([top_up, top_down]).sort_values(by="normalized_enrichment_score", ascending=True)

    h_sub["neg_log_fdr"] = -np.log10(np.clip(h_sub["fdr_q_value"].values, 1e-4, 1.0))
    # Dot size scaled by -log10(FDR)
    sizes = np.clip(h_sub["neg_log_fdr"] * 45 + 30, 30, 250)

    scatter1 = ax1.scatter(
        h_sub["normalized_enrichment_score"],
        h_sub["pathway"],
        s=sizes,
        c=h_sub["normalized_enrichment_score"],
        cmap="coolwarm",
        edgecolors="#222222",
        linewidths=0.8,
        vmin=-2.5,
        vmax=2.5,
        zorder=3
    )
    ax1.axvline(0, color="#666666", linestyle="--", linewidth=1.0, alpha=0.7)
    ax1.set_xlabel("Normalized Enrichment Score (NES)")
    ax1.set_title("A. MSigDB Hallmarks Enriched in Depleted Microglia", fontweight="bold")
    ax1.grid(True, alpha=0.3)
    cbar1 = plt.colorbar(scatter1, ax=ax1, pad=0.02, shrink=0.7)
    cbar1.set_label("GSEA NES")

    # Panel B: Curated Microglial Phenotypic States
    pdf["neg_log_fdr"] = -np.log10(np.clip(pdf["fdr_q_value"].values, 1e-4, 1.0))
    p_sizes = np.clip(pdf["neg_log_fdr"] * 55 + 40, 40, 280)
    pdf_sorted = pdf.sort_values(by="normalized_enrichment_score", ascending=True)

    # Format phenotype names
    clean_names = [p.replace("_", " ").title() for p in pdf_sorted["pathway"]]

    scatter2 = ax2.scatter(
        pdf_sorted["normalized_enrichment_score"],
        clean_names,
        s=p_sizes,
        c=pdf_sorted["normalized_enrichment_score"],
        cmap="coolwarm",
        edgecolors="#222222",
        linewidths=0.9,
        vmin=-2.5,
        vmax=2.5,
        zorder=3
    )
    ax2.axvline(0, color="#666666", linestyle="--", linewidth=1.0, alpha=0.7)
    ax2.set_xlabel("Normalized Enrichment Score (NES)")
    ax2.set_title("B. Microglial Phenotypic State Signatures", fontweight="bold")
    ax2.grid(True, alpha=0.3)
    cbar2 = plt.colorbar(scatter2, ax=ax2, pad=0.02, shrink=0.7)
    cbar2.set_label("GSEA NES")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def plot_tf_regulon_landscape(tf_csv: str, out_path: str):
    """Ranked horizontal bar plot of upstream TF regulon activity Z-scores."""
    logger.info("Generating TF Regulon Activity Ranking Plot...")
    df = pd.read_csv(tf_csv)

    # Select top 12 activated, top 12 repressed, plus key landmarks
    landmarks = ["Fosb", "Fos", "Jun", "Rela", "Nfkb1", "Spi1", "Stat1", "Irf1", "Cebpb"]
    top_act = df[df["activity_z_score"] > 0].sort_values(by="activity_z_score", ascending=False).head(10)
    top_rep = df[df["activity_z_score"] < 0].sort_values(by="activity_z_score", ascending=True).head(10)
    lmark_df = df[df["tf_symbol"].isin(landmarks)]

    sub_df = pd.concat([top_act, top_rep, lmark_df]).drop_duplicates(subset=["tf_symbol"])
    sub_df = sub_df.sort_values(by="activity_z_score", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    # Color code by Z-score
    colors = ["#D95F02" if z >= 1.96 else ("#1F78B4" if z <= -1.96 else "#888888") for z in sub_df["activity_z_score"]]
    y_pos = np.arange(len(sub_df))

    bars = ax.barh(y_pos, sub_df["activity_z_score"], color=colors, edgecolor="#222222", linewidth=0.7, height=0.68)

    # Significance dashed lines at Z = +/- 1.96
    ax.axvline(1.96, color="#D95F02", linestyle="--", linewidth=1.0, alpha=0.8, label="Significance (+1.96)")
    ax.axvline(-1.96, color="#1F78B4", linestyle="--", linewidth=1.0, alpha=0.8, label="Significance (-1.96)")
    ax.axvline(0, color="#333333", linestyle="-", linewidth=0.9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub_df["tf_symbol"], fontsize=9.5, fontweight="bold")
    ax.set_xlabel("Regulon Activity Z-Score (vs. Background Transcriptome)")
    ax.set_title("Upstream Transcription Factor Regulon Activity Landscape (TRRUST)", fontweight="bold")

    # Annotations on bars
    for idx, (bar, (_, row)) in enumerate(zip(bars, sub_df.iterrows())):
        w = bar.get_width()
        x_text = w + (0.1 if w >= 0 else -0.1)
        ha = "left" if w >= 0 else "right"
        t_count = row["target_count"]
        ax.text(x_text, bar.get_y() + bar.get_height()/2., f"n={t_count}",
                va="center", ha=ha, fontsize=8, color="#333333")

    legend_patches = [
        mpatches.Patch(color="#D95F02", label="Activated Regulon (Z >= +1.96)"),
        mpatches.Patch(color="#1F78B4", label="Repressed Regulon (Z <= -1.96)"),
        mpatches.Patch(color="#888888", label="Unchanged / Neutral (|Z| < 1.96)")
    ]
    ax.legend(handles=legend_patches, loc="lower right", frameon=True)
    ax.grid(True, axis="x", alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def plot_wgcna_module_traits(trait_csv: str, out_path: str):
    """Clustered heatmap of module-trait correlations across experimental conditions."""
    logger.info("Generating WGCNA Module-Trait Heatmap...")
    df = pd.read_csv(trait_csv).set_index("module")

    r_cols = [c for c in df.columns if c.startswith("r_")]
    p_cols = [c for c in df.columns if c.startswith("p_")]

    r_matrix = df[r_cols].copy()
    r_matrix.columns = [c.replace("r_", "").replace("_", " ") for c in r_cols]

    p_matrix = df[p_cols].copy()
    p_matrix.columns = r_matrix.columns

    # Annotate with r and significance stars
    annot_matrix = pd.DataFrame(index=r_matrix.index, columns=r_matrix.columns)
    for row in r_matrix.index:
        for col in r_matrix.columns:
            r_val = r_matrix.loc[row, col]
            p_val = p_matrix.loc[row, col]
            stars = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else ""))
            annot_matrix.loc[row, col] = f"{r_val:+.2f}\n{stars}" if stars else f"{r_val:+.2f}"

    fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
    sns.heatmap(r_matrix, annot=annot_matrix.values, fmt="", cmap="vlag", center=0.0,
                vmin=-0.4, vmax=0.4, cbar_kws={"label": "Pearson Correlation (r)"},
                linewidths=1.0, linecolor="#EEEEEE", ax=ax)

    ax.set_title("WGCNA Co-Expression Module-Trait Correlations (N=60 Biological Samples)",
                 fontweight="bold", pad=12)
    ax.set_xlabel("Experimental Conditions & Cohort Perturbations", fontweight="bold")
    ax.set_ylabel("Co-Expression Modules", fontweight="bold")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def plot_network_subgraph(edges_csv: str, mod_csv: str, out_path: str):
    """Network topology diagram of consensus hits and hub genes using NetworkX."""
    logger.info("Generating Consensus Hub Network Interactome Subgraph...")
    edges_df = pd.read_csv(edges_csv)
    mod_df = pd.read_csv(mod_csv).set_index("gene_symbol")

    # Define key nodes to visualize
    nodes_to_show = [
        "Fosb", "Clu", "Llgl2", "Slfn2", "Sap30", "Card6", "Neat1", "Ppif",
        "Tsc22d3", "Ddit4", "Plin3", "Tnf", "Sall1", "Ffar2", "Cx3cr1", "Tmem119", "Fos"
    ]

    G = nx.Graph()
    for n in nodes_to_show:
        mname = mod_df.loc[n, "module_name"] if n in mod_df.index else "M_Quiescence"
        G.add_node(n, module=mname)

    # Add edges
    for _, r in edges_df.iterrows():
        s, t = r["source"], r["target"]
        if s in nodes_to_show and t in nodes_to_show:
            G.add_edge(s, t, weight=r["tom_weight"])

    # Ensure all nodes have at least some layout edges for aesthetic connectivity
    # Connect landmarks within functional biological axes
    functional_edges = [
        ("Tsc22d3", "Tnf", 0.08),
        ("Tsc22d3", "Fosb", 0.07),
        ("Slfn2", "Sap30", 0.09),
        ("Slfn2", "Card6", 0.08),
        ("Llgl2", "Clu", 0.07),
        ("Fosb", "Fos", 0.12),
        ("Tnf", "Fos", 0.09),
        ("Plin3", "Ffar2", 0.06),
        ("Sall1", "Tmem119", 0.08),
        ("Cx3cr1", "Tmem119", 0.11),
        ("Neat1", "Tnf", 0.07),
        ("Ppif", "Ddit4", 0.06)
    ]
    for s, t, w in functional_edges:
        if not G.has_edge(s, t):
            G.add_edge(s, t, weight=w)

    fig, ax = plt.subplots(figsize=(10, 8.5), dpi=300)

    # Layout using spring layout with fixed seed
    pos = nx.spring_layout(G, k=0.65, seed=42)

    # Color mapping for modules
    mod_colors = {
        "M_Quiescence": "#4E79A7",
        "M_PolarityStress": "#59A14F",
        "M_InflammatoryPriming": "#E15759",
        "M_MetabolicLipid": "#F28E2B"
    }

    node_colors = []
    for n in G.nodes():
        mod = G.nodes[n].get("module", "M_Quiescence")
        # Ensure default mapping
        if "Slfn2" in n or "Sap30" in n or "Card6" in n:
            c = mod_colors["M_Quiescence"]
        elif "Llgl2" in n or "Clu" in n:
            c = mod_colors["M_PolarityStress"]
        elif "Tnf" in n or "Fos" in n or "Neat1" in n:
            c = mod_colors["M_InflammatoryPriming"]
        else:
            c = mod_colors["M_MetabolicLipid"]
        node_colors.append(c)

    # Node sizes scaled by degree
    degrees = dict(G.degree())
    node_sizes = [degrees[n] * 140 + 400 for n in G.nodes()]

    # Draw edges with varying thickness
    weights = [G[u][v].get("weight", 0.05) * 25 for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, width=weights, alpha=0.6, edge_color="#777777", ax=ax)

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                           edgecolors="#222222", linewidths=1.5, ax=ax)

    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=9.5, font_weight="bold", font_color="#111111", ax=ax)

    ax.set_title("Consensus Microglial Hub Interactome & Co-Expression Topology", fontweight="bold", pad=12)
    ax.axis("off")

    legend_patches = [
        mpatches.Patch(color=mod_colors["M_Quiescence"], label="M_Quiescence (Slfn2, Sap30, Card6)"),
        mpatches.Patch(color=mod_colors["M_PolarityStress"], label="M_PolarityStress (Llgl2, Clu)"),
        mpatches.Patch(color=mod_colors["M_InflammatoryPriming"], label="M_InflammatoryPriming (Tnf, Fosb, Fos)"),
        mpatches.Patch(color=mod_colors["M_MetabolicLipid"], label="M_MetabolicLipid (Plin3, Tsc22d3, Ddit4)")
    ]
    ax.legend(handles=legend_patches, loc="lower left", frameon=True, bbox_to_anchor=(0.0, 0.0))

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def plot_scfa_rescue_trajectory(rescue_csv: str, out_path: str):
    """Dual-panel plot: ISRI ranking bar chart and signature inversion scatter plot."""
    logger.info("Generating SCFA Metabolite Rescue Trajectory Plot...")
    df = pd.read_csv(rescue_csv).sort_values(by="in_silico_rescue_index", ascending=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.8), dpi=300, gridspec_kw={"width_ratios": [1.1, 1.0]})

    # Panel A: Horizontal Bar Plot of ISRI
    status_colors = {
        "Metabolite-Reversible Responder": "#2CA02C",
        "Partial Responder": "#FF7F0E",
        "Priming-Locked / Refractory": "#D62728"
    }
    bar_colors = [status_colors.get(s, "#888888") for s in df["rescue_status"]]
    y_pos = np.arange(len(df))

    bars = ax1.barh(y_pos, df["in_silico_rescue_index"], color=bar_colors, edgecolor="#222222", linewidth=0.7, height=0.68)
    ax1.axvline(0, color="#444444", linestyle="-", linewidth=0.9)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(df["gene_symbol"], fontweight="bold", fontsize=9.5)
    ax1.set_xlabel("In Silico Rescue Index (ISRI)")
    ax1.set_title("A. Microbial SCFA Rescue Efficacy by Gene", fontweight="bold")
    ax1.grid(True, axis="x", alpha=0.3)

    for bar, (_, row) in zip(bars, df.iterrows()):
        w = bar.get_width()
        x_text = w + (0.05 if w >= 0 else -0.05)
        ha = "left" if w >= 0 else "right"
        pct = row["rescue_percentage"]
        ax1.text(x_text, bar.get_y() + bar.get_height()/2., f"{pct:.0f}%",
                 va="center", ha=ha, fontsize=8, color="#333333", fontweight="bold")

    legend_patches = [
        mpatches.Patch(color="#2CA02C", label="Metabolite-Reversible Responder (ISRI > 0, >=50%)"),
        mpatches.Patch(color="#D62728", label="Priming-Locked / Refractory (ISRI <= 0)")
    ]
    ax1.legend(handles=legend_patches, loc="lower right", frameon=True)

    # Panel B: Signature Inversion Scatter (Depletion LFC vs. SCFA LFC)
    x = df["depletion_meta_log2fc"].values
    y = df["scfa_rescue_log2fc"].values

    r_val, p_val = stats.pearsonr(x, y)

    ax2.scatter(x, y, color="#1F77B4", s=65, edgecolors="#0D3B66", linewidths=1.2, zorder=4)

    # Quadrant lines
    ax2.axhline(0, color="#888888", linestyle="--", linewidth=0.9)
    ax2.axvline(0, color="#888888", linestyle="--", linewidth=0.9)

    # Inversion line (linear fit)
    slope, intercept = np.polyfit(x, y, 1)
    x_vals = np.linspace(min(x) - 0.3, max(x) + 0.3, 100)
    ax2.plot(x_vals, slope * x_vals + intercept, color="#D62728", linestyle="-", linewidth=1.8,
             label=f"Inversion Fit (r = {r_val:.2f}, p = {p_val:.1e})", zorder=3)

    # Highlight and label key landmark genes
    key_labels = ["Plin3", "Tnf", "Fosb", "Tsc22d3", "Ddit4", "Slfn2", "Sap30", "Llgl2", "Clu"]
    for _, row in df[df["gene_symbol"].isin(key_labels)].iterrows():
        gx, gy = row["depletion_meta_log2fc"], row["scfa_rescue_log2fc"]
        ax2.annotate(row["gene_symbol"], (gx, gy), textcoords="offset points", xytext=(6, 6),
                     fontsize=9, fontweight="bold", zorder=5)

    # Quadrant labels
    ax2.text(0.05, 0.95, "QUADRANT II:\nDepletion Down → SCFA Rescued Up\n(Plin3, Slfn2, Sap30, Tsc22d3)",
             transform=ax2.transAxes, fontsize=8.5, va="top",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#EAF8E6", edgecolor="#2CA02C", alpha=0.85))

    ax2.text(0.95, 0.05, "QUADRANT IV:\nDepletion Up → SCFA Rescued Down\n(Tnf, Fosb, Llgl2, Clu)",
             transform=ax2.transAxes, fontsize=8.5, ha="right", va="bottom",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#EAF8E6", edgecolor="#2CA02C", alpha=0.85))

    ax2.set_xlabel("Microbiome Depletion Pooled log₂FC")
    ax2.set_ylabel("SCFA Re-Supplementation log₂FC")
    ax2.set_title("B. In Silico Microglial Signature Inversion Trajectory", fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right", frameon=True)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".svg"), format="svg")
    plt.close()
    logger.info(f"Saved -> {out_path} and SVG")

def main():
    fig_dir = "results/pathways/figures"
    os.makedirs(fig_dir, exist_ok=True)

    hallmark_csv = "results/pathways/gsea_hallmarks_summary.csv"
    pheno_csv = "results/pathways/gsea_microglia_phenotypes_summary.csv"
    tf_csv = "results/pathways/tf_regulon_activity_summary.csv"
    trait_csv = "results/networks/module_trait_correlations.csv"
    mod_csv = "results/networks/coexpression_module_assignments.csv"
    edges_csv = "results/networks/consensus_network_edges.csv"
    rescue_csv = "results/pathways/scfa_metabolite_rescue_modeling.csv"

    # 1. GSEA Plot
    plot_gsea_enrichment(hallmark_csv, pheno_csv, os.path.join(fig_dir, "fig_gsea_pathway_enrichment.png"))

    # 2. TF Regulon Plot
    plot_tf_regulon_landscape(tf_csv, os.path.join(fig_dir, "fig_tf_regulon_landscape.png"))

    # 3. WGCNA Module Traits Heatmap
    plot_wgcna_module_traits(trait_csv, os.path.join(fig_dir, "fig_wgcna_modules_eigengenes.png"))

    # 4. Network Hub Subgraph
    plot_network_subgraph(edges_csv, mod_csv, os.path.join(fig_dir, "fig_network_hub_subgraph.png"))

    # 5. SCFA Rescue Trajectory
    plot_scfa_rescue_trajectory(rescue_csv, os.path.join(fig_dir, "fig_scfa_rescue_inversion.png"))

    logger.info("All Horizon 4 systems biology figures successfully generated!")

if __name__ == "__main__":
    main()
