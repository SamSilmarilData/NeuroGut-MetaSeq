"""
tests/test_metadata_schema.py
Validates metadata schemas for all cohorts.
"""

import os
import glob
import pytest
import pandas as pd

REQUIRED_COLUMNS = [
    "sample_id", "cohort", "condition", "group_label", "sex", "tissue", "sequencing_type"
]
ALLOWED_CONDITIONS = ["reference", "perturbed"]

def test_metadata_files_exist():
    meta_files = glob.glob("data/demo/*_metadata.csv")
    assert len(meta_files) >= 3, "Expected at least 3 demo metadata files"

@pytest.mark.parametrize("meta_file", glob.glob("data/demo/*_metadata.csv"))
def test_metadata_columns_and_values(meta_file):
    df = pd.read_csv(meta_file)
    # Check all required columns exist
    for col in REQUIRED_COLUMNS:
        assert col in df.columns, f"Missing required column '{col}' in {meta_file}"

    # Check non-empty sample IDs
    assert df["sample_id"].notnull().all(), f"Found null sample_id in {meta_file}"
    assert df["sample_id"].nunique() == len(df), f"Duplicate sample IDs in {meta_file}"

    # Check allowed conditions
    assert df["condition"].isin(ALLOWED_CONDITIONS).all(), f"Invalid condition values in {meta_file}"

    # Check tissue is Microglia
    assert (df["tissue"] == "Microglia").all(), f"Non-microglia tissue in {meta_file}"
