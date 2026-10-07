#!/usr/bin/env python3
"""
scripts/01b_audit_data_landscape.py
Comprehensive diagnostic data landscape and quality control audit.
Evaluates:
  1. Sequencing depth, gene detection, and sparsity
  2. Microglial lineage purity vs. neuronal, astrocytic, oligodendrocytic, endothelial contaminants
  3. Ex vivo enzymatic dissociation stress scores (IEG activation) and treatment bias
Outputs:
  - Tabular audit: results/qc/horizon1_data_audit.csv
  - Diagnostic figures: results/qc/figures/fig_qc_*.png
"""

import os
import sys
import glob
import logging
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("QC-Audit")

MICROGLIAL_MARKERS = ["Tmem119", "Cx3cr1", "P2ry12", "Sall1", "Mertk", "Itgam", "Ptprc", "Tyrobp"]
CONTAMINANT_MARKERS = {
    "Neuronal": ["Rbfox3", "Snap25", "Map2"],
    "Astrocytic": ["Gfap", "Aldh1l1", "Aqp4"],
    "Oligodendrocytic": ["Mbp", "Mog", "Olig2"],
    "Endothelial": ["Cdh5", "Pecam1"]
}
ALL_CONTAMINANTS = [g for sub in CONTAMINANT_MARKERS.values() for g in sub]
STRESS_GENES = ["Fos", "Jun", "Atf3", "Egr1", "Nfkbia", "Dusp1"]

def compute_cpm(counts: np.ndarray) -> np.ndarray:
    """Compute Counts Per Million (CPM)."""
    total = np.sum(counts, axis=0)
    total[total == 0] = 1.0
    return (counts / total) * 1e6

def audit_cohort(cohort_id: str, counts_file: str, meta_file: str) -> pd.DataFrame:
    """Audit a single cohort and return sample-level QC metrics."""
    logger.info(f"Auditing cohort: {cohort_id}...")
    counts_df = pd.read_csv(counts_file)
    counts_df.columns = counts_df.columns.astype(str)
    meta_df = pd.read_csv(meta_file)
    meta_df["sample_id"] = meta_df["sample_id"].astype(str)

    sample_cols = meta_df["sample_id"].tolist()
    gene_symbols = counts_df["gene_symbol"].astype(str).tolist()
    mat = counts_df[sample_cols].values.astype(float)

    # 1. Depth & Sparsity
    depths = np.sum(mat, axis=0)
    detected_genes = np.sum(mat >= 5, axis=0)
    sparsity = np.mean(mat == 0, axis=0) * 100.0

    # CPM matrix for marker evaluation
    cpm = compute_cpm(mat)
    gene_dict = {g: i for i, g in enumerate(gene_symbols)}

    # 2. Lineage Purity Audit
    mg_indices = [gene_dict[g] for g in MICROGLIAL_MARKERS if g in gene_dict]
    contam_indices = [gene_dict[g] for g in ALL_CONTAMINANTS if g in gene_dict]

    mg_cpm = np.mean(cpm[mg_indices, :], axis=0) if mg_indices else np.zeros(len(sample_cols))
    contam_cpm = np.mean(cpm[contam_indices, :], axis=0) if contam_indices else np.zeros(len(sample_cols))

    purity_ratio = np.where(
        (mg_cpm + contam_cpm) > 0,
        (mg_cpm / (mg_cpm + contam_cpm)) * 100.0,
        100.0
    )

    # Specific contaminant CPMs
    contam_by_type = {}
    for ctype, clist in CONTAMINANT_MARKERS.items():
        c_idx = [gene_dict[g] for g in clist if g in gene_dict]
        contam_by_type[ctype] = np.mean(cpm[c_idx, :], axis=0) if c_idx else np.zeros(len(sample_cols))

    # 3. Ex Vivo Isolation Stress Audit
    stress_indices = [gene_dict[g] for g in STRESS_GENES if g in gene_dict]
    if stress_indices:
        stress_cpm = cpm[stress_indices, :]
        # Standardize each gene across samples
        z_stress = (stress_cpm - np.mean(stress_cpm, axis=1, keepdims=True)) / (np.std(stress_cpm, axis=1, keepdims=True) + 1e-6)
        composite_stress = np.mean(z_stress, axis=0)
    else:
        composite_stress = np.zeros(len(sample_cols))

    res_df = meta_df.copy()
    res_df["sequencing_depth"] = depths
    res_df["detected_genes_ge5"] = detected_genes
    res_df["sparsity_pct"] = np.round(sparsity, 2)
    res_df["microglial_marker_cpm"] = np.round(mg_cpm, 2)
    res_df["contaminant_cpm"] = np.round(contam_cpm, 2)
    res_df["purity_index_pct"] = np.round(purity_ratio, 2)
    res_df["purity_warning"] = res_df["purity_index_pct"] < 90.0
    res_df["isolation_stress_zscore"] = np.round(composite_stress, 3)

    for ctype, c_vals in contam_by_type.items():
        res_df[f"{ctype.lower()}_cpm"] = np.round(c_vals, 2)

    return res_df

def test_stress_confounding(audit_df: pd.DataFrame) -> pd.DataFrame:
    """Test if isolation stress scores differ significantly between reference and perturbed."""
    cohort_tests = []
    for cid, group in audit_df.groupby("cohort"):
        ref_stress = group[group["condition"] == "reference"]["isolation_stress_zscore"].values
        pert_stress = group[group["condition"] == "perturbed"]["isolation_stress_zscore"].values

        if len(ref_stress) > 1 and len(pert_stress) > 1:
            stat, p_val = stats.mannwhitneyu(ref_stress, pert_stress, alternative="two-sided")
            mean_diff = np.mean(pert_stress) - np.mean(ref_stress)
        else:
            stat, p_val, mean_diff = np.nan, 1.0, 0.0

        cohort_tests.append({
            "cohort": cid,
            "ref_n": len(ref_stress),
            "pert_n": len(pert_stress),
            "mean_stress_diff": np.round(mean_diff, 3),
            "mwu_stat": stat,
            "p_value": np.round(p_val, 4),
            "confounding_flag": p_val < 0.05
        })
    return pd.DataFrame(cohort_tests)

def plot_diagnostics(audit_df: pd.DataFrame, out_dir: str):
    """Render publication-grade diagnostic plots."""
    os.makedirs(out_dir, exist_ok=True)
    sns.set_theme(style="whitegrid", font="DejaVu Sans")

    # Plot 1: Library Depths
    plt.figure(figsize=(9, 5))
    palette = sns.color_palette("Set2", n_colors=audit_df["cohort"].nunique())
    ax = sns.boxplot(x="cohort", y="sequencing_depth", hue="cohort", data=audit_df, palette=palette, legend=False, showmeans=True, boxprops=dict(alpha=0.8))
    sns.stripplot(x="cohort", y="sequencing_depth", hue="condition", data=audit_df, dodge=True, palette="dark:gray", alpha=0.7, size=6, ax=ax)
    plt.title("Sequencing Depth Distribution Across Discovery Cohorts", fontsize=13, fontweight="bold", pad=12)
    plt.xlabel("Cohort ID", fontsize=11, fontweight="bold")
    plt.ylabel("Total Mapped Reads per Library", fontsize=11, fontweight="bold")
    plt.yscale("log")
    plt.tight_layout()
    fig1_path = os.path.join(out_dir, "fig_qc_library_depths.png")
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    logger.info(f"Saved: {fig1_path}")

    # Plot 2: Microglial Lineage Purity
    plt.figure(figsize=(9, 5))
    ax = sns.barplot(x="cohort", y="purity_index_pct", hue="condition", data=audit_df, palette="Blues_d", errorbar="sd", capsize=0.1)
    ax.axhline(90.0, color="#d9534f", linestyle="--", linewidth=1.5, label="Quality Threshold (90%)")
    plt.title("Microglial Lineage Purity Index Across Cohorts", fontsize=13, fontweight="bold", pad=12)
    plt.xlabel("Cohort ID", fontsize=11, fontweight="bold")
    plt.ylabel("Microglial Purity Index (%)", fontsize=11, fontweight="bold")
    plt.ylim(0, 105)
    plt.legend(loc="lower right")
    plt.tight_layout()
    fig2_path = os.path.join(out_dir, "fig_qc_microglial_purity.png")
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    logger.info(f"Saved: {fig2_path}")

    # Plot 3: Ex Vivo Isolation Stress
    plt.figure(figsize=(9, 5))
    palette_stress = {"reference": "#4575b4", "perturbed": "#d73027"}
    ax = sns.boxplot(x="cohort", y="isolation_stress_zscore", hue="condition", data=audit_df, palette=palette_stress, boxprops=dict(alpha=0.7))
    sns.stripplot(x="cohort", y="isolation_stress_zscore", hue="condition", data=audit_df, dodge=True, palette="dark:gray", alpha=0.6, size=5, ax=ax)
    plt.title("Ex Vivo Dissociation Stress Composite Score (Fos, Jun, Atf3, Egr1)", fontsize=13, fontweight="bold", pad=12)
    plt.xlabel("Cohort ID", fontsize=11, fontweight="bold")
    plt.ylabel("Standardized Stress Composite Z-Score", fontsize=11, fontweight="bold")
    # Clean up duplicate legends from stripplot
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[:2], labels[:2], title="Condition", loc="upper right")
    plt.tight_layout()
    fig3_path = os.path.join(out_dir, "fig_qc_isolation_stress.png")
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    logger.info(f"Saved: {fig3_path}")

def main():
    proc_dir = "data/processed"
    meta_dir = "data/metadata"
    qc_dir = "results/qc"
    fig_dir = os.path.join(qc_dir, "figures")
    os.makedirs(qc_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    cohorts = ["GSE107925", "GSE108045", "GSE266602", "GSE186210"]
    all_audits = []

    for cid in cohorts:
        counts_path = os.path.join(proc_dir, f"{cid}_counts.csv")
        meta_path = os.path.join(meta_dir, f"{cid}_metadata.csv")
        if not os.path.exists(counts_path) or not os.path.exists(meta_path):
            logger.warning(f"Skipping {cid}: Processed files not found.")
            continue
        audit_res = audit_cohort(cid, counts_path, meta_path)
        all_audits.append(audit_res)

    if not all_audits:
        logger.error("No cohorts audited. Ensure scripts/02_curate_metadata.py has run.")
        sys.exit(1)

    combined_audit = pd.concat(all_audits, ignore_index=True)
    audit_csv = os.path.join(qc_dir, "horizon1_data_audit.csv")
    combined_audit.to_csv(audit_csv, index=False)
    logger.info(f"[SUCCESS] Saved data audit table: {audit_csv} ({len(combined_audit)} samples)")

    # Confounding tests
    stress_confounding = test_stress_confounding(combined_audit)
    stress_csv = os.path.join(qc_dir, "horizon1_stress_confounding_test.csv")
    stress_confounding.to_csv(stress_csv, index=False)
    logger.info(f"[SUCCESS] Saved stress confounding test results: {stress_csv}")
    print("\n=== Ex Vivo Dissociation Stress Confounding Test ===")
    print(stress_confounding.to_string(index=False))

    # Summary table
    print("\n=== Microglial Purity & Library Depth Summary ===")
    summary = combined_audit.groupby("cohort").agg({
        "sample_id": "count",
        "sequencing_depth": ["mean", "median", "min", "max"],
        "purity_index_pct": ["mean", "min"],
        "purity_warning": "sum",
        "isolation_stress_zscore": "mean"
    })
    print(summary)

    # Diagnostic Plots
    plot_diagnostics(combined_audit, fig_dir)
    logger.info("[SUCCESS] Diagnostic QC plots generated successfully.")

if __name__ == "__main__":
    main()
