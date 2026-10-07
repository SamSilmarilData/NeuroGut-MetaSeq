"""
tests/test_de_results.py
Validates the outputs, column schemas, and statistical bounds of Horizon 2 differential expression results.
"""

import os
import pytest
import pandas as pd
import numpy as np

COHORTS = ["GSE107925", "GSE108045", "GSE266602", "GSE186210"]
EXPECTED_COLUMNS = ["gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
EXPECTED_FIGURES = [
    "results/de_results/figures/fig_volcano_GSE107925.png",
    "results/de_results/figures/fig_volcano_GSE108045.png",
    "results/de_results/figures/fig_volcano_GSE266602.png",
    "results/de_results/figures/fig_volcano_GSE186210.png",
    "results/de_results/figures/fig_lfc_correlation_heatmap.png",
    "results/de_results/figures/fig_marker_effect_sizes_by_cohort.png",
    "results/de_results/figures/fig_cohort_deg_overlap.png"
]

@pytest.mark.parametrize("cohort", COHORTS)
def test_de_result_file_exists_and_schema(cohort):
    csv_path = f"results/de_results/{cohort}_deg.csv"
    assert os.path.exists(csv_path), f"Missing DEG result file: {csv_path}"
    
    df = pd.read_csv(csv_path)
    assert len(df) >= 15000, f"Expected >= 15,000 tested genes in {cohort}, found {len(df)}"

    for col in EXPECTED_COLUMNS:
        assert col in df.columns, f"Missing required column '{col}' in {csv_path}"

    assert df["gene_symbol"].notnull().all(), f"Null gene symbols found in {cohort}"
    assert not np.isinf(df["log2FoldChange"].dropna()).any(), f"Infinite log2FC detected in {cohort}"
    
    # Statistical bounds
    valid_padj = df["padj"].dropna()
    assert ((valid_padj >= 0.0) & (valid_padj <= 1.0)).all(), f"padj out of bounds in {cohort}"
    valid_pval = df["pvalue"].dropna()
    assert ((valid_pval >= 0.0) & (valid_pval <= 1.0)).all(), f"pvalue out of bounds in {cohort}"

def test_cohort_phenotypic_summary_exists():
    summary_path = "results/de_results/cohort_phenotypic_summary.csv"
    assert os.path.exists(summary_path), f"Missing summary file: {summary_path}"
    df = pd.read_csv(summary_path)
    assert len(df) == 4, f"Expected 4 cohorts in phenotypic summary, found {len(df)}"
    assert set(df["cohort"].unique()) == set(COHORTS)
    assert (df["total_expressed_genes"] >= 15000).all()

def test_phenotypic_figures_exist_and_valid():
    for fig_path in EXPECTED_FIGURES:
        assert os.path.exists(fig_path), f"Missing figure: {fig_path}"
        assert os.path.getsize(fig_path) > 10000, f"Figure {fig_path} is too small / corrupt"
