"""
tests/test_meta_analysis_real.py
Automated validation of real cross-study meta-analysis results (Horizon 3).
Verifies:
1. Master summary, LOO sensitivity, and core consensus table existence and schemas.
2. Mathematical bounds: p-values, FDR, Higgins I^2, confidence intervals, and SE.
3. Leave-One-Out consistency and robustness score bounds.
4. Publication diagnostic figures existence and integrity.
"""

import os
import pytest
import pandas as pd
import numpy as np

META_DIR = "results/meta_results"
SUMMARY_FILE = os.path.join(META_DIR, "microglia_meta_analysis_summary.csv")
LOO_FILE = os.path.join(META_DIR, "microglia_meta_analysis_loo.csv")
CORE_FILE = os.path.join(META_DIR, "core_consensus_signature.csv")
FIG_DIR = os.path.join(META_DIR, "figures")

@pytest.fixture(scope="module")
def summary_df():
    assert os.path.exists(SUMMARY_FILE), f"Missing {SUMMARY_FILE}"
    df = pd.read_csv(SUMMARY_FILE)
    assert not df.empty, "Summary dataframe is empty"
    return df

@pytest.fixture(scope="module")
def loo_df():
    assert os.path.exists(LOO_FILE), f"Missing {LOO_FILE}"
    df = pd.read_csv(LOO_FILE)
    assert not df.empty, "LOO dataframe is empty"
    return df

@pytest.fixture(scope="module")
def core_df():
    assert os.path.exists(CORE_FILE), f"Missing {CORE_FILE}"
    df = pd.read_csv(CORE_FILE)
    assert not df.empty, "Core consensus dataframe is empty"
    return df

def test_summary_table_scale_and_schema(summary_df):
    """Verifies that meta-analysis covers >20,000 common genes and has required schema."""
    assert len(summary_df) >= 20000, f"Expected >= 20000 genes, found {len(summary_df)}"
    expected_cols = [
        "gene_symbol", "n_cohorts", "cohorts_detected", "meta_log2fc", "meta_se",
        "ci_lower", "ci_upper", "tau2", "cochran_q", "i2_heterogeneity",
        "heterogeneity_tier", "p_random_effects", "fdr_random_effects",
        "fisher_stat", "p_fisher", "fdr_fisher", "stouffer_z", "p_stouffer", "fdr_stouffer",
        "direction_concordance", "robustness_score", "max_lfc_shift", "loo_vulnerable_study",
        "significance_flag"
    ]
    for col in expected_cols:
        assert col in summary_df.columns, f"Missing column {col} in summary table"

def test_pvalue_and_fdr_bounds(summary_df):
    """Verifies all p-values and FDR values are strictly bounded in [0, 1]."""
    for col in ["p_random_effects", "fdr_random_effects", "p_fisher", "fdr_fisher", "p_stouffer", "fdr_stouffer"]:
        vals = summary_df[col].dropna()
        assert (vals >= 0.0).all(), f"Found negative values in {col}"
        assert (vals <= 1.0).all(), f"Found values > 1.0 in {col}"

def test_heterogeneity_metrics_bounds(summary_df):
    """Verifies Higgins I^2 is bounded in [0, 100] and tau^2, Cochran's Q are non-negative."""
    assert (summary_df["i2_heterogeneity"] >= 0.0).all()
    assert (summary_df["i2_heterogeneity"] <= 100.0).all()
    assert (summary_df["cochran_q"] >= 0.0).all()
    assert (summary_df["tau2"] >= 0.0).all()

def test_confidence_interval_and_se_integrity(summary_df):
    """Verifies CI lower <= meta_log2fc <= CI upper, and SE > 0."""
    assert (summary_df["meta_se"] > 0.0).all(), "Found non-positive meta_se"
    # ci_lower <= meta_log2fc <= ci_upper with small float tolerance
    assert (summary_df["ci_lower"] <= summary_df["meta_log2fc"] + 1e-4).all()
    assert (summary_df["ci_upper"] >= summary_df["meta_log2fc"] - 1e-4).all()

def test_loo_sensitivity_table_integrity(loo_df):
    """Verifies that LOO sensitivity table covers all 4 cohort omissions."""
    expected_cols = [
        "gene_symbol", "omitted_cohort", "k_remaining", "loo_log2fc",
        "loo_se", "loo_cochran_q", "loo_i2", "loo_p_value", "shift_from_full", "loo_fdr"
    ]
    for col in expected_cols:
        assert col in loo_df.columns, f"Missing column {col} in LOO table"

    omitted_cohorts = set(loo_df["omitted_cohort"].unique())
    assert {"GSE107925", "GSE108045", "GSE186210", "GSE266602"}.issubset(omitted_cohorts)
    assert len(loo_df) >= 30000, f"Expected >= 30000 LOO rows, found {len(loo_df)}"

def test_core_consensus_signature_properties(core_df):
    """Verifies that core consensus signature meets stringency filters."""
    assert (core_df["significance_flag"] == True).all()
    assert (core_df["n_cohorts"] >= 3).all()
    assert (core_df["i2_heterogeneity"] < 50.0).all()
    assert "consensus_tier" in core_df.columns
    assert set(core_df["consensus_tier"].unique()).issubset({
        "Tier 1: Omnipresent Core (k=4)",
        "Tier 2: Robust Broad Core (k=3)"
    })

def test_diagnostic_figures_exist():
    """Verifies all 5 publication-grade diagnostic figures exist and are non-empty."""
    figures = [
        "fig_meta_volcano.png",
        "fig_forest_plots_top.png",
        "fig_loo_stability.png",
        "fig_heterogeneity_distribution.png",
        "fig_consensus_heatmap.png"
    ]
    for fig in figures:
        p = os.path.join(FIG_DIR, fig)
        assert os.path.exists(p), f"Missing figure: {p}"
        assert os.path.getsize(p) > 10000, f"Figure {p} is suspiciously small ({os.path.getsize(p)} bytes)"
