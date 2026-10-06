"""
tests/test_count_matrix_integrity.py
Validates integrity of RNA-seq count matrices.
"""

import os
import glob
import pytest
import pandas as pd
import numpy as np

@pytest.mark.parametrize("cohort", ["GSE107925", "GSE108045", "GSE266602"])
def test_count_matrix_structure(cohort):
    counts_file = f"data/demo/{cohort}_counts.csv"
    meta_file = f"data/demo/{cohort}_metadata.csv"

    assert os.path.exists(counts_file), f"Missing count matrix {counts_file}"
    assert os.path.exists(meta_file), f"Missing metadata file {meta_file}"

    counts_df = pd.read_csv(counts_file)
    meta_df = pd.read_csv(meta_file)

    # Check gene identification columns
    assert "gene_id" in counts_df.columns, f"Missing gene_id in {counts_file}"
    assert "gene_symbol" in counts_df.columns, f"Missing gene_symbol in {counts_file}"
    assert counts_df["gene_symbol"].notnull().all(), f"Found null gene_symbol in {counts_file}"

    # Check sample columns match metadata
    sample_cols = [c for c in counts_df.columns if c not in ["gene_id", "gene_symbol"]]
    meta_samples = meta_df["sample_id"].tolist()
    assert set(sample_cols) == set(meta_samples), f"Sample column mismatch in {cohort}"

    # Verify counts are non-negative integers
    for s in sample_cols:
        col_vals = counts_df[s].values
        assert np.issubdtype(col_vals.dtype, np.integer) or np.all(col_vals == np.round(col_vals)), f"Non-integer counts in {s}"
        assert np.all(col_vals >= 0), f"Negative counts in {s}"
        assert not np.isnan(col_vals).any(), f"NaN values in {s}"
