#!/usr/bin/env python3
"""
scripts/04c_perturbation_subgroups.py
Subgroup Analysis & Heterogeneity Decomposition:
Disentangles Shared Microbial Tonic Surveillance from Model-Private Perturbation Axes
(Lifelong Germ-Free vs. Acute Antibiotic Cocktail vs. Dietary Fiber Starvation).

Outputs:
- results/meta_results/perturbation_subgroup_decomposition.csv
- results/meta_results/figures/fig_subgroup_perturbation_decomposition.png (and .svg)
"""

import os
import glob
import logging
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Subgroup-Analysis")

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

def load_data(meta_res_path="results/meta_results/microglia_meta_analysis_summary.csv",
              de_dir="results/de_results",
              counts_dir="data/processed",
              meta_dir="data/metadata"):
    """Loads meta-analysis summary, cohort DEG tables, and sample-level count matrices."""
    meta_df = pd.read_csv(meta_res_path)
    
    # Cohort DEG tables
    cohort_degs = {}
    for deg_file in sorted(glob.glob(os.path.join(de_dir, "*_deg.csv"))):
        cname = os.path.basename(deg_file).replace("_deg.csv", "")
        df = pd.read_csv(deg_file).drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")
        cohort_degs[cname] = df

    # Sample-level expression
    cohort_files = sorted(glob.glob(os.path.join(counts_dir, "*_counts.csv")))
    all_sample_dfs = []
    sample_meta_list = []

    for cfile in cohort_files:
        cname = os.path.basename(cfile).replace("_counts.csv", "")
        counts_df = pd.read_csv(cfile)
        if "gene_symbol" not in counts_df.columns:
            continue
        counts_df = counts_df.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")
        sample_cols = [c for c in counts_df.columns if c != "gene_id"]
        
        lib_sizes = counts_df[sample_cols].sum(axis=0)
        cpm = np.log2((counts_df[sample_cols].div(lib_sizes, axis=1) * 1e6) + 1.0)
        # Z-standardize within cohort to remove technical baseline shifts
        cpm_z = cpm.sub(cpm.mean(axis=1), axis=0).div(cpm.std(axis=1).replace(0, 1), axis=0)
        
        meta_f = os.path.join(meta_dir, f"{cname}_metadata.csv")
        mdf = pd.read_csv(meta_f).set_index("sample_id")
        
        for s in sample_cols:
            cond = mdf.loc[s, "condition"] if s in mdf.index else "unknown"
            sample_meta_list.append({
                "sample_id": s,
                "cohort": cname,
                "condition": cond,
                "model_type": "Germ-Free (Developmental)" if "107925" in cname else (
                              "Antibiotics (Acute Shock)" if "108045" in cname else (
                              "Fiber Starvation (Dietary)" if "186210" in cname else "Sham / Baseline"))
            })
        all_sample_dfs.append(cpm_z)

    # Common genes across expression matrices
    common_expr_genes = set(all_sample_dfs[0].index)
    for df in all_sample_dfs[1:]:
        common_expr_genes = common_expr_genes.intersection(df.index)
    common_expr_genes = sorted(list(common_expr_genes))
    
    expr_matrix = pd.concat([df.loc[common_expr_genes] for df in all_sample_dfs], axis=1)
    samples_df = pd.DataFrame(sample_meta_list).set_index("sample_id")
    expr_matrix = expr_matrix[samples_df.index]

    return meta_df, cohort_degs, expr_matrix, samples_df

def compute_subgroup_decomposition(meta_df: pd.DataFrame, cohort_degs: dict) -> pd.DataFrame:
    """
    Decomposes gene effects into:
    - Shared Microbial Core (concordant, low Q_between)
    - ABX-Specific Shock (disproportionate in GSE108045)
    - Fiber-Specific Metabolic (disproportionate in GSE186210)
    - Lifelong Germ-Free Developmental (dominant in GSE107925)
    """
    records = []
    
    c_gf = "GSE107925"
    c_abx = "GSE108045"
    c_fiber = "GSE186210"
    c_sham = "GSE266602"

    for _, row in meta_df.iterrows():
        gene = row["gene_symbol"]
        
        # Get cohort-specific estimates
        lfc_gf = cohort_degs[c_gf].loc[gene, "log2FoldChange"] if gene in cohort_degs[c_gf].index else np.nan
        se_gf = cohort_degs[c_gf].loc[gene, "lfcSE"] if gene in cohort_degs[c_gf].index else np.nan
        
        lfc_abx = cohort_degs[c_abx].loc[gene, "log2FoldChange"] if gene in cohort_degs[c_abx].index else np.nan
        se_abx = cohort_degs[c_abx].loc[gene, "lfcSE"] if gene in cohort_degs[c_abx].index else np.nan
        
        lfc_fiber = cohort_degs[c_fiber].loc[gene, "log2FoldChange"] if gene in cohort_degs[c_fiber].index else np.nan
        se_fiber = cohort_degs[c_fiber].loc[gene, "lfcSE"] if gene in cohort_degs[c_fiber].index else np.nan

        lfc_sham = cohort_degs[c_sham].loc[gene, "log2FoldChange"] if gene in cohort_degs[c_sham].index else np.nan
        se_sham = cohort_degs[c_sham].loc[gene, "lfcSE"] if gene in cohort_degs[c_sham].index else np.nan

        # Test between-subgroup heterogeneity across the 3 distinct perturbation classes (GF, ABX, Fiber)
        subgroup_thetas = []
        subgroup_vars = []
        for lfc, se in [(lfc_gf, se_gf), (lfc_abx, se_abx), (lfc_fiber, se_fiber)]:
            if not np.isnan(lfc) and not np.isnan(se) and se > 0:
                subgroup_thetas.append(lfc)
                subgroup_vars.append(se**2)

        if len(subgroup_thetas) >= 2:
            th_arr = np.array(subgroup_thetas)
            v_arr = np.array(subgroup_vars)
            w_arr = 1.0 / v_arr
            th_pooled = np.sum(w_arr * th_arr) / np.sum(w_arr)
            q_between = float(np.sum(w_arr * (th_arr - th_pooled)**2))
            df_between = len(subgroup_thetas) - 1
            p_q_between = float(1.0 - stats.chi2.cdf(q_between, df_between))
        else:
            q_between = 0.0
            p_q_between = 1.0

        i2 = row["i2_heterogeneity"]
        pooled_lfc = row["meta_log2fc"]
        models_lfc = [x for x in [lfc_gf, lfc_abx, lfc_fiber] if not np.isnan(x)]
        all_same_sign = (len(models_lfc) >= 2) and (all(x > 0 for x in models_lfc) or all(x < 0 for x in models_lfc))

        is_sig = (row["fdr_random_effects"] < 0.05 and abs(pooled_lfc) >= 0.4) or bool(row["significance_flag"])

        subgroup_axis = "Unassigned / Neutral"
        rationale = "Non-significant or low effect size"

        if is_sig and i2 < 35.0 and all_same_sign and p_q_between > 0.05:
            subgroup_axis = "Shared Microbial Core"
            rationale = f"Consensus invariant across models (I²={i2:.1f}%, Q_between p={p_q_between:.3f})"
        elif not np.isnan(lfc_abx) and abs(lfc_abx) >= 1.0 and i2 > 60.0:
            subgroup_axis = "ABX Mucosal Shock"
            rationale = f"Antibiotic-dominant perturbation effect (|LFC_ABX|={abs(lfc_abx):.2f}, I²={i2:.1f}%)"
        elif not np.isnan(lfc_fiber) and abs(lfc_fiber) >= 0.8 and i2 > 50.0:
            subgroup_axis = "Fiber Dietary Starvation"
            rationale = f"Fiber starvation-dominant metabolic effect (|LFC_Fiber|={abs(lfc_fiber):.2f}, I²={i2:.1f}%)"
        elif not np.isnan(lfc_gf) and abs(lfc_gf) >= 0.8 and i2 > 50.0:
            subgroup_axis = "Developmental Germ-Free"
            rationale = f"Lifelong embryonic microbial absence effect (|LFC_GF|={abs(lfc_gf):.2f}, I²={i2:.1f}%)"
        elif is_sig and all_same_sign:
            subgroup_axis = "Shared Microbial Core (Moderate Heterogeneity)"
            rationale = f"Concordant sign across models with moderate quantitative heterogeneity (I²={i2:.1f}%)"
        elif is_sig:
            subgroup_axis = "Model-Divergent Mixed Axis"
            rationale = f"Opposite directional regulation across perturbation models (I²={i2:.1f}%)"
        else:
            subgroup_axis = "Unassigned / Neutral"
            rationale = "Non-significant or low effect size"

        records.append({
            "gene_symbol": gene,
            "subgroup_axis": subgroup_axis,
            "classification_rationale": rationale,
            "meta_log2fc": pooled_lfc,
            "meta_se": row["meta_se"],
            "p_random_effects": row["p_random_effects"],
            "fdr_random_effects": row["fdr_random_effects"],
            "i2_heterogeneity": i2,
            "q_between_models": round(q_between, 3),
            "p_q_between": round(p_q_between, 4),
            "lfc_germ_free": round(lfc_gf, 4) if not np.isnan(lfc_gf) else np.nan,
            "lfc_antibiotics": round(lfc_abx, 4) if not np.isnan(lfc_abx) else np.nan,
            "lfc_fiber_starvation": round(lfc_fiber, 4) if not np.isnan(lfc_fiber) else np.nan,
            "lfc_sham_percoll": round(lfc_sham, 4) if not np.isnan(lfc_sham) else np.nan
        })

    sub_df = pd.DataFrame(records)
    return sub_df

def plot_subgroup_decomposition(sub_df: pd.DataFrame, expr_matrix: pd.DataFrame, samples_df: pd.DataFrame, out_png: str, out_svg: str):
    """Generates 4-panel multi-study decomposition figure."""
    logger.info("Generating Subgroup Decomposition Figures (PNG + SVG)...")
    
    fig = plt.figure(figsize=(16, 12), dpi=300)
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.28)
    
    # -------------------------------------------------------------
    # Panel A: Multi-Study Factor Analysis / SVD PCA (Factor 1 vs Factor 2)
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    
    # Run PCA on standardized expression of top variable genes across all 60 samples
    var_genes = expr_matrix.var(axis=1).sort_values(ascending=False).head(2000).index
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(expr_matrix.loc[var_genes].T)
    var_exp = pca.explained_variance_ratio_ * 100
    
    pca_df = pd.DataFrame(pca_coords, columns=["PC1", "PC2"], index=expr_matrix.columns)
    pca_df = pca_df.join(samples_df)
    
    cohort_palette = {
        "GSE107925": "#4E79A7",
        "GSE108045": "#F28E2B",
        "GSE186210": "#E15759",
        "GSE266602": "#76B7B2"
    }
    cond_markers = {"reference": "o", "perturbed": "s"}
    
    for cname, cgroup in pca_df.groupby("cohort"):
        for cond, sgroup in cgroup.groupby("condition"):
            ax_a.scatter(
                sgroup["PC1"], sgroup["PC2"],
                c=cohort_palette.get(cname, "#333333"),
                marker=cond_markers.get(cond, "o"),
                s=65, alpha=0.85, edgecolors="#222222", linewidth=0.8,
                label=f"{cname} ({cond.capitalize()})"
            )
            
    ax_a.set_xlabel(f"Factor 1: Microbial Tonic Depletion Axis ({var_exp[0]:.1f}% var)", fontweight="bold")
    ax_a.set_ylabel(f"Factor 2: Model Modality / Stress Axis ({var_exp[1]:.1f}% var)", fontweight="bold")
    ax_a.set_title("A. Multi-Study Factor Analysis (N=60 Biological Samples)", fontweight="bold")
    ax_a.axhline(0, color="#888888", linestyle=":", lw=0.8)
    ax_a.axvline(0, color="#888888", linestyle=":", lw=0.8)
    ax_a.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=True, fontsize=8)
    ax_a.grid(True, alpha=0.3)

    # -------------------------------------------------------------
    # Panel B: Subgroup Effect Size Concordance (GF vs ABX)
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    
    valid_mask = sub_df["lfc_germ_free"].notna() & sub_df["lfc_antibiotics"].notna()
    plot_sub = sub_df[valid_mask].copy()
    
    # Highlight categories
    shared_core = plot_sub[plot_sub["subgroup_axis"].str.contains("Shared Microbial Core")]
    abx_shock = plot_sub[plot_sub["subgroup_axis"] == "ABX Mucosal Shock"]
    fiber_shock = plot_sub[plot_sub["subgroup_axis"] == "Fiber Dietary Starvation"]
    neutral = plot_sub[plot_sub["subgroup_axis"].str.contains("Unassigned|Neutral|Divergent")]
    
    ax_b.scatter(neutral["lfc_germ_free"], neutral["lfc_antibiotics"],
                 c="#D0D4DC", alpha=0.25, s=12, label="Unassigned / Low Effect", rasterized=True)
    ax_b.scatter(shared_core["lfc_germ_free"], shared_core["lfc_antibiotics"],
                 c="#0072B2", alpha=0.85, s=35, edgecolors="#003B6F", lw=0.6,
                 label=f"Shared Microbial Core (n={len(shared_core)})", zorder=4)
    ax_b.scatter(abx_shock["lfc_germ_free"], abx_shock["lfc_antibiotics"],
                 c="#D55E00", alpha=0.85, s=40, edgecolors="#8C2D00", lw=0.7,
                 label=f"ABX Mucosal Shock (n={len(abx_shock)})", zorder=5)

    # Annotate landmark genes
    landmarks = [("Llgl2", 0.15, 0.2), ("Clu", -0.2, 0.2), ("Slfn2", -0.3, -0.3),
                 ("Tsc22d3", 0.2, -0.4), ("Ddit4", 0.2, -0.5), ("Plin3", -0.4, 0.2),
                 ("Fosb", 0.2, 0.2), ("Tnf", 0.2, 0.2)]
    
    for gene, dx, dy in landmarks:
        if gene in plot_sub["gene_symbol"].values:
            g_row = plot_sub[plot_sub["gene_symbol"] == gene].iloc[0]
            gx, gy = g_row["lfc_germ_free"], g_row["lfc_antibiotics"]
            ax_b.scatter(gx, gy, s=70, facecolors="none", edgecolors="#111111", lw=1.5, zorder=6)
            ax_b.annotate(gene, (gx, gy), (gx + dx, gy + dy),
                          fontweight="bold", fontsize=9,
                          arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8))

    ax_b.axline((0, 0), slope=1, color="#888888", linestyle="--", lw=1.0, alpha=0.7, label="Concordance Line (y = x)")
    ax_b.axhline(0, color="#CCCCCC", linestyle=":", lw=0.8)
    ax_b.axvline(0, color="#CCCCCC", linestyle=":", lw=0.8)
    ax_b.set_xlabel("Germ-Free Effect (log₂FC, GSE107925)", fontweight="bold")
    ax_b.set_ylabel("Antibiotic Cocktail Effect (log₂FC, GSE108045)", fontweight="bold")
    ax_b.set_title("B. Cross-Perturbation Concordance & Divergence", fontweight="bold")
    ax_b.legend(loc="upper left", frameon=True, fontsize=8)
    ax_b.grid(True, alpha=0.3)
    ax_b.set_xlim(-4, 4)
    ax_b.set_ylim(-5, 5)

    # -------------------------------------------------------------
    # Panel C: Subgroup Multi-Model Forest Profile for Key Signatures
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    
    key_genes_profile = [
        ("Llgl2", "Shared Invariant Core"),
        ("Clu", "Shared Invariant Core"),
        ("Slfn2", "Shared Invariant Core"),
        ("Tsc22d3", "ABX-Private Artifact"),
        ("Ddit4", "ABX-Private Artifact"),
        ("Plin3", "Fiber-Private Axis")
    ]
    
    y_positions = np.arange(len(key_genes_profile))
    bar_width = 0.22
    
    gf_vals = [sub_df[sub_df["gene_symbol"] == g]["lfc_germ_free"].values[0] for g, _ in key_genes_profile]
    abx_vals = [sub_df[sub_df["gene_symbol"] == g]["lfc_antibiotics"].values[0] for g, _ in key_genes_profile]
    fiber_vals = [sub_df[sub_df["gene_symbol"] == g]["lfc_fiber_starvation"].values[0] for g, _ in key_genes_profile]

    ax_c.barh(y_positions - bar_width, gf_vals, height=bar_width, color="#4E79A7", label="Germ-Free (GSE107925)")
    ax_c.barh(y_positions, abx_vals, height=bar_width, color="#F28E2B", label="Antibiotics (GSE108045)")
    ax_c.barh(y_positions + bar_width, fiber_vals, height=bar_width, color="#E15759", label="Fiber Starvation (GSE186210)")

    ax_c.set_yticks(y_positions)
    ax_c.set_yticklabels([f"{g}\n({role})" for g, role in key_genes_profile], fontsize=8.5, fontweight="bold")
    ax_c.axvline(0, color="#333333", linestyle="-", lw=0.9)
    ax_c.set_xlabel("log₂ Fold Change", fontweight="bold")
    ax_c.set_title("C. Perturbation-Specific Effect Size Decomposition", fontweight="bold")
    ax_c.legend(loc="lower right", frameon=True, fontsize=8)
    ax_c.grid(True, axis="x", alpha=0.3)

    # -------------------------------------------------------------
    # Panel D: Subgroup Architecture Breakdown Pie/Bar
    # -------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    
    sig_sub = sub_df[sub_df["fdr_random_effects"] < 0.05]
    axis_counts = sig_sub["subgroup_axis"].value_counts()
    
    colors_dict = {
        "Shared Microbial Core": "#0072B2",
        "Shared Microbial Core (Moderate Heterogeneity)": "#56B4E9",
        "ABX Mucosal Shock": "#D55E00",
        "Fiber Dietary Starvation": "#CC79A7",
        "Developmental Germ-Free": "#009E73",
        "Model-Divergent Mixed Axis": "#E69F00",
        "Unassigned / Neutral": "#999999"
    }
    
    ordered_axes = [k for k in colors_dict.keys() if k in axis_counts.index]
    counts_vals = [axis_counts[k] for k in ordered_axes]
    bar_colors = [colors_dict[k] for k in ordered_axes]
    
    bars = ax_d.barh(np.arange(len(ordered_axes)), counts_vals, color=bar_colors, edgecolor="#222222", lw=0.7)
    ax_d.set_yticks(np.arange(len(ordered_axes)))
    ax_d.set_yticklabels(ordered_axes, fontsize=8.5, fontweight="bold")
    ax_d.set_xlabel("Number of Meta-Significant Genes (FDR < 0.05)", fontweight="bold")
    ax_d.set_title("D. Global Transcriptome Partition Across Perturbation Axes", fontweight="bold")
    ax_d.grid(True, axis="x", alpha=0.3)
    
    for bar in bars:
        w = bar.get_width()
        pct = (w / len(sig_sub)) * 100
        ax_d.text(w + 5, bar.get_y() + bar.get_height()/2., f"{w} ({pct:.1f}%)",
                  va="center", ha="left", fontsize=8.5, fontweight="bold")
    ax_d.set_xlim(0, max(counts_vals) * 1.25)

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_svg, format="svg", bbox_inches="tight")
    plt.close()
    logger.info(f"Saved PNG -> {out_png}")
    logger.info(f"Saved SVG -> {out_svg}")

def main():
    logger.info("Initializing Subgroup Analysis & Heterogeneity Decomposition...")
    meta_res_path = "results/meta_results/microglia_meta_analysis_summary.csv"
    if not os.path.exists(meta_res_path):
        raise FileNotFoundError(f"{meta_res_path} not found.")

    meta_df, cohort_degs, expr_matrix, samples_df = load_data()
    sub_df = compute_subgroup_decomposition(meta_df, cohort_degs)
    
    out_csv = "results/meta_results/perturbation_subgroup_decomposition.csv"
    sub_df.to_csv(out_csv, index=False)
    logger.info(f"Saved subgroup decomposition -> {out_csv}")

    fig_dir = "results/meta_results/figures"
    os.makedirs(fig_dir, exist_ok=True)
    out_png = os.path.join(fig_dir, "fig_subgroup_perturbation_decomposition.png")
    out_svg = os.path.join(fig_dir, "fig_subgroup_perturbation_decomposition.svg")
    
    plot_subgroup_decomposition(sub_df, expr_matrix, samples_df, out_png, out_svg)
    logger.info("Subgroup analysis and figure generation complete!")

if __name__ == "__main__":
    main()
