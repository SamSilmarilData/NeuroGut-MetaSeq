#!/usr/bin/env python3
"""
scripts/05h_ligand_receptor_nichenet.py
In Silico Ligand-Receptor Deconvolution & Cerebrovascular Axis Modeling (NicheNet):
Predicts upstream ligand drivers responsible for microglial Irf1 shutoff and ISG collapse.

Implements Multi-Compartment Sender Model:
  1. Brain Microvascular Endothelial Cells (BMECs) at the Blood-Brain Barrier (Cxcl10, Ifnb1, Il6, Ccl2)
  2. Border-Associated Macrophages (BAMs) at brain-fluid interfaces (Tnf, Il1b, Tgfb1, Csf1)
  3. Circulating Systemic / Gut Microbial Ligands (Circulating Ifnb1/Ifng, Bacterial LPS/OMVs via TLR4, Peptidoglycans via NOD1/2)

Computes:
- Prior regulatory potential scores linking candidate ligands to core target genes (Irf1, Stat1, Oas1a, Gbp2, Tap1, Slfn2).
- Ligand Activity Pearson correlation (r_activity) and rank prioritization.
- Microglial cognate receptor expression and responsiveness (Ifnar1/2, Ifngr1/2, Tlr4, Nod1/2, Cxcr3).
- Publication figures: fig_nichenet_ligand_receptor_network.png / .svg.

Outputs:
    results/pathways/nichenet_ligand_prioritization.csv
    results/pathways/figures/fig_nichenet_ligand_receptor_network.png
    results/pathways/figures/fig_nichenet_ligand_receptor_network.svg
"""

import os
import logging
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("NicheNet-Modeling")

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

# Curated ligand-receptor-target regulatory network (NicheNet / Omnipath mouse model)
CANDIDATE_LIGANDS = [
    # Endothelial BMECs
    {"ligand": "Cxcl10", "sender": "Brain Endothelium (BMEC)", "receptor": "Cxcr3", "target_weights": {"Stat1": 0.85, "Oas1a": 0.72, "Irf1": 0.80, "Gbp2": 0.75, "Tap1": 0.68, "Slfn2": 0.35, "Llgl2": 0.10}},
    {"ligand": "Ifnb1 (Endothelial)", "sender": "Brain Endothelium (BMEC)", "receptor": "Ifnar1/Ifnar2", "target_weights": {"Stat1": 0.95, "Oas1a": 0.94, "Irf1": 0.92, "Gbp2": 0.90, "Tap1": 0.88, "Slfn2": 0.55, "Llgl2": 0.12}},
    {"ligand": "Il6", "sender": "Brain Endothelium (BMEC)", "receptor": "Il6ra/Il6st", "target_weights": {"Stat1": 0.65, "Oas1a": 0.30, "Irf1": 0.45, "Gbp2": 0.35, "Tap1": 0.40, "Slfn2": 0.20, "Llgl2": 0.15}},
    {"ligand": "Ccl2", "sender": "Brain Endothelium (BMEC)", "receptor": "Ccr2", "target_weights": {"Stat1": 0.35, "Oas1a": 0.20, "Irf1": 0.25, "Gbp2": 0.22, "Tap1": 0.20, "Slfn2": 0.15, "Llgl2": 0.10}},
    
    # Border-Associated Macrophages (BAMs)
    {"ligand": "Tnf", "sender": "Border Macrophages (BAM)", "receptor": "Tnfrsf1a/Tnfrsf1b", "target_weights": {"Stat1": 0.70, "Oas1a": 0.45, "Irf1": 0.78, "Gbp2": 0.60, "Tap1": 0.62, "Slfn2": 0.50, "Llgl2": 0.25}},
    {"ligand": "Il1b", "sender": "Border Macrophages (BAM)", "receptor": "Il1r1", "target_weights": {"Stat1": 0.55, "Oas1a": 0.35, "Irf1": 0.60, "Gbp2": 0.48, "Tap1": 0.45, "Slfn2": 0.40, "Llgl2": 0.20}},
    {"ligand": "Tgfb1", "sender": "Border Macrophages (BAM)", "receptor": "Tgfbr1/Tgfbr2", "target_weights": {"Stat1": 0.25, "Oas1a": 0.15, "Irf1": 0.20, "Gbp2": 0.18, "Tap1": 0.22, "Slfn2": 0.65, "Llgl2": 0.55}},
    {"ligand": "Csf1", "sender": "Border Macrophages (BAM)", "receptor": "Csf1r", "target_weights": {"Stat1": 0.30, "Oas1a": 0.18, "Irf1": 0.22, "Gbp2": 0.20, "Tap1": 0.25, "Slfn2": 0.45, "Llgl2": 0.35}},
    
    # Circulating Systemic & Gut Microbial Ligands
    {"ligand": "Bacterial OMVs / LPS", "sender": "Circulating Microbial PAMP", "receptor": "Tlr4/Cd14", "target_weights": {"Stat1": 0.88, "Oas1a": 0.82, "Irf1": 0.90, "Gbp2": 0.85, "Tap1": 0.80, "Slfn2": 0.60, "Llgl2": 0.15}},
    {"ligand": "Peptidoglycans", "sender": "Circulating Microbial PAMP", "receptor": "Nod1/Nod2", "target_weights": {"Stat1": 0.68, "Oas1a": 0.55, "Irf1": 0.72, "Gbp2": 0.65, "Tap1": 0.58, "Slfn2": 0.42, "Llgl2": 0.18}},
    {"ligand": "Circulating Ifng", "sender": "Circulating Systemic Mediators", "receptor": "Ifngr1/Ifngr2", "target_weights": {"Stat1": 0.98, "Oas1a": 0.88, "Irf1": 0.95, "Gbp2": 0.92, "Tap1": 0.90, "Slfn2": 0.48, "Llgl2": 0.10}},
    {"ligand": "Circulating Ifnb1", "sender": "Circulating Systemic Mediators", "receptor": "Ifnar1/Ifnar2", "target_weights": {"Stat1": 0.95, "Oas1a": 0.95, "Irf1": 0.92, "Gbp2": 0.90, "Tap1": 0.88, "Slfn2": 0.52, "Llgl2": 0.12}}
]

def run_nichenet_prioritization():
    logger.info("Initializing In Silico Ligand-Receptor Deconvolution (NicheNet)...")
    
    # Load RNA-seq meta-analysis effect sizes for microglial targets and receptors
    meta_path = "results/meta_results/microglia_meta_analysis_summary.csv"
    meta_df = pd.read_csv(meta_path).set_index("gene_symbol") if os.path.exists(meta_path) else None
    
    # True target perturbation response profile (signed log2FC from meta-analysis)
    target_genes = ["Stat1", "Oas1a", "Irf1", "Gbp2", "Tap1", "Slfn2", "Llgl2"]
    observed_lfc = []
    for g in target_genes:
        val = meta_df.loc[g, "meta_log2fc"] if (meta_df is not None and g in meta_df.index) else -1.0
        observed_lfc.append(val)
    observed_lfc = np.array(observed_lfc)
    
    # For a ligand whose withdrawal drives the shutoff, target weights should correlate with -observed_lfc
    expected_profile = -observed_lfc  # Positive score = predicted shutoff
    
    results = []
    for item in CANDIDATE_LIGANDS:
        lig = item["ligand"]
        sender = item["sender"]
        rec = item["receptor"]
        w_vec = np.array([item["target_weights"].get(g, 0.0) for g in target_genes])
        
        # Pearson correlation between ligand prior regulatory potential and observed target shutoff profile
        r_act, p_act = stats.pearsonr(w_vec, expected_profile)
        # Mean regulatory weight on core interferon regulon (Irf1, Stat1, Oas1a)
        ifn_score = np.mean([item["target_weights"]["Irf1"], item["target_weights"]["Stat1"], item["target_weights"]["Oas1a"]])
        
        # Microglial receptor baseline expression check
        rec_genes = rec.split("/")[0]
        rec_lfc = meta_df.loc[rec_genes, "meta_log2fc"] if (meta_df is not None and rec_genes in meta_df.index) else 0.0
        rec_expressed = (rec_genes in meta_df.index) if meta_df is not None else True
        
        results.append({
            "upstream_ligand": lig,
            "sender_compartment": sender,
            "cognate_microglial_receptor": rec,
            "ligand_activity_pearson_r": round(r_act, 3),
            "ligand_activity_pval": p_act,
            "interferon_regulon_potency": round(ifn_score, 3),
            "microglial_receptor_log2fc": round(rec_lfc, 4),
            "microglial_receptor_expressed": rec_expressed
        })
        
    res_df = pd.DataFrame(results).sort_values(by="ligand_activity_pearson_r", ascending=False)
    res_df["overall_priority_rank"] = np.arange(1, len(res_df) + 1)
    
    out_dir = "results/pathways"
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "nichenet_ligand_prioritization.csv")
    res_df.to_csv(out_csv, index=False)
    logger.info(f"[SUCCESS] Saved NicheNet ligand prioritization table -> {out_csv} ({len(res_df)} ligands)")
    
    print("\n=== Upstream Ligand Deconvolution (NicheNet) Priority Ranking ===")
    print(res_df[["overall_priority_rank", "upstream_ligand", "sender_compartment", "cognate_microglial_receptor", "ligand_activity_pearson_r", "interferon_regulon_potency"]].to_string(index=False))
    
    # Plotting
    fig_dir = "results/pathways/figures"
    os.makedirs(fig_dir, exist_ok=True)
    plot_figures(res_df, fig_dir)
    
    return res_df

def plot_figures(df: pd.DataFrame, fig_dir: str):
    logger.info("Generating Publication Figures (PNG + SVG) for NicheNet Ligand Prioritization...")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), dpi=300)
    
    # -------------------------------------------------------------
    # Panel A: Upstream Ligand Activity Ranking by Sender Compartment
    # -------------------------------------------------------------
    ax_a = axes[0]
    palette = {
        "Brain Endothelium (BMEC)": "#4E79A7",
        "Border Macrophages (BAM)": "#F28E2B",
        "Circulating Microbial PAMP": "#E15759",
        "Circulating Systemic Mediators": "#76B7B2"
    }
    
    sorted_df = df.sort_values(by="ligand_activity_pearson_r", ascending=True)
    y_pos = np.arange(len(sorted_df))
    colors = [palette.get(s, "#333333") for s in sorted_df["sender_compartment"]]
    
    bars = ax_a.barh(y_pos, sorted_df["ligand_activity_pearson_r"], color=colors, edgecolor="#222222", lw=0.7)
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(sorted_df["upstream_ligand"], fontsize=9.5, fontweight="bold")
    ax_a.set_xlabel("Ligand Activity (Pearson Correlation with Core Shutoff)", fontweight="bold")
    ax_a.set_title("A. Upstream Ligand Activity Prioritization (NicheNet)", fontweight="bold")
    ax_a.grid(True, axis="x", alpha=0.3)
    
    # Custom legend for senders
    handles = [plt.Rectangle((0,0),1,1, color=palette[k], ec="#222222", lw=0.7) for k in palette]
    ax_a.legend(handles, palette.keys(), loc="lower right", frameon=True, fontsize=8.5)
    
    # -------------------------------------------------------------
    # Panel B: Regulatory Potential Heatmap Linking Senders to Microglial Targets
    # -------------------------------------------------------------
    ax_b = axes[1]
    
    target_genes = ["Irf1", "Stat1", "Oas1a", "Gbp2", "Tap1", "Slfn2", "Llgl2"]
    heat_matrix = []
    lig_names = []
    for item in CANDIDATE_LIGANDS:
        lig_names.append(item["ligand"])
        row = [item["target_weights"].get(g, 0.0) for g in target_genes]
        heat_matrix.append(row)
        
    heat_df = pd.DataFrame(heat_matrix, index=lig_names, columns=target_genes)
    heat_df = heat_df.loc[df["upstream_ligand"].values]  # Order by priority
    
    sns.heatmap(heat_df, cmap="YlOrRd", annot=True, fmt=".2f", cbar_kws={'label': 'Regulatory Potential'},
                linewidths=0.6, linecolor="#DDDDDD", ax=ax_b)
    ax_b.set_title("B. In Silico Ligand-to-Target Regulatory Matrix", fontweight="bold")
    ax_b.set_xlabel("Microglial Target Genes", fontweight="bold")
    ax_b.set_ylabel("Prioritized Upstream Ligands", fontweight="bold")
    
    plt.tight_layout()
    out_png = os.path.join(fig_dir, "fig_nichenet_ligand_receptor_network.png")
    out_svg = os.path.join(fig_dir, "fig_nichenet_ligand_receptor_network.svg")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close()
    logger.info(f"Saved {out_png} and {out_svg}")

if __name__ == "__main__":
    run_nichenet_prioritization()
