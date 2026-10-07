"""
tests/test_data_audit.py
Validates the outputs and integrity of Horizon 1 data landscape and quality control audit.
"""

import os
import pytest
import pandas as pd

QC_AUDIT_CSV = "results/qc/horizon1_data_audit.csv"
QC_STRESS_CSV = "results/qc/horizon1_stress_confounding_test.csv"
EXPECTED_FIGURES = [
    "results/qc/figures/fig_qc_library_depths.png",
    "results/qc/figures/fig_qc_microglial_purity.png",
    "results/qc/figures/fig_qc_isolation_stress.png"
]

def test_qc_audit_csv_exists_and_valid():
    assert os.path.exists(QC_AUDIT_CSV), f"Missing {QC_AUDIT_CSV}"
    df = pd.read_csv(QC_AUDIT_CSV)
    assert len(df) == 60, f"Expected 60 curated samples across cohorts, found {len(df)}"
    
    expected_cohorts = {"GSE107925", "GSE108045", "GSE266602", "GSE186210"}
    assert set(df["cohort"].unique()) == expected_cohorts

    required_cols = [
        "sample_id", "cohort", "condition", "group_label",
        "sequencing_depth", "detected_genes_ge5", "sparsity_pct",
        "microglial_marker_cpm", "contaminant_cpm", "purity_index_pct",
        "isolation_stress_zscore"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing required column {col} in QC audit"

    # Depths should be positive
    assert (df["sequencing_depth"] > 1e5).all(), "Unrealistically low sequencing depth detected"
    # Purity should be between 0 and 100
    assert ((df["purity_index_pct"] >= 0) & (df["purity_index_pct"] <= 100)).all()

def test_stress_confounding_test_valid():
    assert os.path.exists(QC_STRESS_CSV), f"Missing {QC_STRESS_CSV}"
    df = pd.read_csv(QC_STRESS_CSV)
    assert len(df) == 4, f"Expected 4 cohorts in stress test, found {len(df)}"
    # Verify no systematic confounding (p > 0.05 for all cohorts)
    assert not df["confounding_flag"].any(), "Significant dissociation stress confounding detected!"

def test_qc_figures_exist_and_non_empty():
    for fig_path in EXPECTED_FIGURES:
        assert os.path.exists(fig_path), f"Missing figure: {fig_path}"
        assert os.path.getsize(fig_path) > 10000, f"Figure {fig_path} is too small / corrupt"
