"""
tests/test_systems_biology_real.py
Automated validation of Systems Biology, Pathway Enrichment, TF Regulons,
WGCNA Co-expression Networks, and In Silico SCFA Metabolite Rescue (Horizon 4).

Verifies:
1. Reference gene set caches (MSigDB Hallmark, KEGG, TRRUST, microglia phenotypes).
2. Whole-transcriptome GSEA bounds, NES, FDR, and specific pathway enrichments.
3. Upstream TF regulon deconvolution activities and key TFs (IRF1, FOSB, RELA).
4. Co-expression network module assignments, eigengene-trait correlations, and hub genes.
5. In silico SCFA metabolite rescue indices (ISRI), rescue percentages, and inversion metrics.
6. Publication diagnostic figure generation and integrity.
"""

import os
import json
import pytest
import pandas as pd
import numpy as np

REF_DIR = "data/reference"
PATHWAYS_DIR = "results/pathways"
NETWORKS_DIR = "results/networks"
FIGURES_DIR = os.path.join(PATHWAYS_DIR, "figures")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def hallmark_gsea_df():
    path = os.path.join(PATHWAYS_DIR, "gsea_hallmarks_summary.csv")
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert not df.empty, "Hallmark GSEA summary is empty"
    return df


@pytest.fixture(scope="module")
def kegg_gsea_df():
    path = os.path.join(PATHWAYS_DIR, "gsea_kegg_summary.csv")
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert not df.empty, "KEGG GSEA summary is empty"
    return df


@pytest.fixture(scope="module")
def phenotypes_gsea_df():
    path = os.path.join(PATHWAYS_DIR, "gsea_microglia_phenotypes_summary.csv")
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert not df.empty, "Phenotypes GSEA summary is empty"
    return df


@pytest.fixture(scope="module")
def tf_regulon_df():
    path = os.path.join(PATHWAYS_DIR, "tf_regulon_activity_summary.csv")
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert not df.empty, "TF regulon activity summary is empty"
    return df


@pytest.fixture(scope="module")
def scfa_rescue_df():
    path = os.path.join(PATHWAYS_DIR, "scfa_metabolite_rescue_modeling.csv")
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert not df.empty, "SCFA metabolite rescue summary is empty"
    return df


@pytest.fixture(scope="module")
def coexp_modules_df():
    path = os.path.join(NETWORKS_DIR, "coexpression_module_assignments.csv")
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert not df.empty, "Coexpression module assignments table is empty"
    return df


@pytest.fixture(scope="module")
def module_traits_df():
    path = os.path.join(NETWORKS_DIR, "module_trait_correlations.csv")
    assert os.path.exists(path), f"Missing {path}"
    df = pd.read_csv(path)
    assert not df.empty, "Module-trait correlations table is empty"
    return df


# ---------------------------------------------------------------------------
# Test Reference Gene Set Caches
# ---------------------------------------------------------------------------

def test_gene_set_caches_integrity():
    """Validates that reference gene set JSON files are properly cached and populated."""
    expected_files = {
        "msigdb_hallmark_mouse.json": 50,
        "kegg_mouse.json": 200,          # At least 200 KEGG pathways
        "trrust_mouse.json": 300,        # At least 300 TRRUST TFs
        "microglia_phenotypes.json": 5   # Exactly 5 curated microglia phenotypes
    }
    for fname, min_keys in expected_files.items():
        fpath = os.path.join(REF_DIR, fname)
        assert os.path.exists(fpath), f"Missing reference file {fpath}"
        with open(fpath, "r") as f:
            data = json.load(f)
        assert isinstance(data, dict), f"{fname} is not a valid JSON dictionary"
        assert len(data) >= min_keys, f"{fname} has {len(data)} keys, expected >= {min_keys}"
        # Validate that each gene set is a non-empty list of gene symbols
        for key, genes in data.items():
            assert isinstance(genes, list), f"Value for {key} in {fname} is not a list"
            assert len(genes) > 0, f"Gene set {key} in {fname} is empty"


# ---------------------------------------------------------------------------
# Test Whole-Transcriptome GSEA & ORA
# ---------------------------------------------------------------------------

def test_hallmark_gsea_bounds_and_findings(hallmark_gsea_df):
    """Verifies Hallmark GSEA metrics and confirms Bloom 4.1 (Interferon Gamma repression)."""
    assert len(hallmark_gsea_df) == 50, f"Expected 50 hallmarks, got {len(hallmark_gsea_df)}"

    # Check bounds
    assert hallmark_gsea_df["nominal_p_value"].between(0.0, 1.0).all()
    assert hallmark_gsea_df["fdr_q_value"].between(0.0, 1.0).all()
    assert np.isfinite(hallmark_gsea_df["normalized_enrichment_score"]).all()

    # Bloom 4.1: Interferon Gamma Response must be significantly repressed
    ifn_gamma = hallmark_gsea_df[hallmark_gsea_df["pathway"] == "Interferon Gamma Response"]
    assert not ifn_gamma.empty, "Missing Interferon Gamma Response"
    nes = ifn_gamma["normalized_enrichment_score"].values[0]
    fdr = ifn_gamma["fdr_q_value"].values[0]
    assert nes < -1.5, f"Expected strong negative NES for IFN-gamma, got {nes}"
    assert fdr <= 0.05, f"Expected FDR <= 0.05 for IFN-gamma, got {fdr}"


def test_microglia_phenotype_gsea(phenotypes_gsea_df):
    """Verifies enrichment bounds for curated microglial phenotypes."""
    assert len(phenotypes_gsea_df) == 5
    assert phenotypes_gsea_df["fdr_q_value"].between(0.0, 1.0).all()

    # IRM (Interferon Responsive Microglia) should be strongly repressed
    irm = phenotypes_gsea_df[phenotypes_gsea_df["pathway"] == "Interferon_Responsive_Microglia_IRM"]
    assert not irm.empty, "Missing Interferon_Responsive_Microglia_IRM"
    assert irm["normalized_enrichment_score"].values[0] < -1.5
    assert irm["fdr_q_value"].values[0] <= 0.05


def test_ora_consensus_results():
    """Verifies that ORA analysis table exists and contains expected columns and valid p-values."""
    ora_path = os.path.join(PATHWAYS_DIR, "ora_consensus_pathways.csv")
    assert os.path.exists(ora_path), f"Missing {ora_path}"
    df = pd.read_csv(ora_path)
    assert not df.empty

    for p_col in ["core_p_value", "shock_p_value", "up_p_value", "core_fdr", "shock_fdr", "up_fdr"]:
        assert p_col in df.columns
        valid_vals = df[p_col].dropna()
        assert valid_vals.between(0.0, 1.0).all(), f"Values in {p_col} out of [0, 1]"


# ---------------------------------------------------------------------------
# Test Upstream TF Regulon Deconvolution
# ---------------------------------------------------------------------------

def test_tf_regulon_landscape(tf_regulon_df):
    """Verifies TF regulon activity deconvolution and confirms IRF1 repression."""
    assert len(tf_regulon_df) >= 300, f"Expected >= 300 TFs, got {len(tf_regulon_df)}"
    assert (tf_regulon_df["target_count"] >= 5).all(), "Found TF with fewer than 5 measured targets"

    # Verify p-value bounds
    for col in ["p_welch", "p_mann_whitney", "p_ks_test", "p_fisher_overlap", "fdr_welch", "fdr_mann_whitney", "fdr_ks", "fdr_fisher"]:
        assert col in tf_regulon_df.columns
        vals = tf_regulon_df[col].dropna()
        assert vals.between(0.0, 1.0).all(), f"Values in {col} out of [0, 1]"

    # Verify key microglia and immune TFs are present
    tfs = set(tf_regulon_df["tf_symbol"].str.lower())
    for key_tf in ["irf1", "fos", "rela", "stat1"]:
        assert key_tf in tfs, f"Key transcription factor {key_tf} missing from regulon table"

    # Bloom 4.1: IRF1 regulon activity must be significantly repressed (Z < 0, FDR < 0.05)
    irf1_row = tf_regulon_df[tf_regulon_df["tf_symbol"].str.lower() == "irf1"]
    assert not irf1_row.empty
    irf1_z = irf1_row["activity_z_score"].values[0]
    irf1_fdr = irf1_row["fdr_welch"].values[0]
    assert irf1_z < -1.5, f"Expected negative activity Z-score for Irf1, got {irf1_z}"
    assert irf1_fdr <= 0.05, f"Expected Irf1 FDR <= 0.05, got {irf1_fdr}"


# ---------------------------------------------------------------------------
# Test WGCNA Co-Expression Network Analysis
# ---------------------------------------------------------------------------

def test_coexpression_modules_and_traits(coexp_modules_df, module_traits_df):
    """Verifies co-expression clustering and module-trait correlations."""
    assert len(coexp_modules_df) >= 3000, f"Expected >= 3000 genes in network, got {len(coexp_modules_df)}"
    assert "module_name" in coexp_modules_df.columns
    assert "k_in" in coexp_modules_df.columns
    assert (coexp_modules_df["k_in"] >= 0).all(), "Intramodular connectivity k_in must be non-negative"

    # Trait correlation table
    r_cols = [c for c in module_traits_df.columns if c.startswith("r_")]
    p_cols = [c for c in module_traits_df.columns if c.startswith("p_")]
    assert len(r_cols) >= 3, "Expected at least 3 trait correlation columns"
    for r_col in r_cols:
        assert module_traits_df[r_col].between(-1.0, 1.0).all(), f"Correlations in {r_col} must be in [-1, 1]"
    for p_col in p_cols:
        assert module_traits_df[p_col].between(0.0, 1.0).all(), f"P-values in {p_col} must be in [0, 1]"


def test_hub_genes_and_edges():
    """Verifies hub genes summary and network edges."""
    hubs_path = os.path.join(NETWORKS_DIR, "hub_genes_summary.csv")
    edges_path = os.path.join(NETWORKS_DIR, "consensus_network_edges.csv")
    assert os.path.exists(hubs_path), f"Missing {hubs_path}"
    assert os.path.exists(edges_path), f"Missing {edges_path}"

    hubs_df = pd.read_csv(hubs_path)
    assert not hubs_df.empty
    assert "gene_symbol" in hubs_df.columns
    assert "intramodular_k_in" in hubs_df.columns

    edges_df = pd.read_csv(edges_path)
    assert not edges_df.empty
    assert "tom_weight" in edges_df.columns
    assert (edges_df["tom_weight"] >= 0).all(), "TOM weights must be non-negative"


# ---------------------------------------------------------------------------
# Test In Silico SCFA Metabolite Rescue Modeling
# ---------------------------------------------------------------------------

def test_scfa_metabolite_rescue(scfa_rescue_df):
    """Verifies Bloom 4.2: In silico SCFA signature inversion and rescue indexing."""
    assert len(scfa_rescue_df) >= 15, f"Expected >= 15 tested genes, got {len(scfa_rescue_df)}"
    required_cols = [
        "gene_symbol", "depletion_meta_log2fc", "scfa_rescue_log2fc",
        "net_post_rescue_log2fc", "in_silico_rescue_index", "rescue_percentage", "rescue_status"
    ]
    for col in required_cols:
        assert col in scfa_rescue_df.columns, f"Missing column {col} in SCFA rescue table"

    # Check key genes are included
    genes = set(scfa_rescue_df["gene_symbol"].values)
    for key_gene in ["Plin3", "Tnf", "Fosb", "Slfn2", "Sap30", "Tsc22d3"]:
        assert key_gene in genes, f"Key gene {key_gene} missing from SCFA rescue table"

    # Bloom 4.2 signature inversion: strong negative correlation between depletion and rescue log2FC
    r = np.corrcoef(scfa_rescue_df["depletion_meta_log2fc"], scfa_rescue_df["scfa_rescue_log2fc"])[0, 1]
    assert r < -0.5, f"Expected strong negative correlation (r < -0.5), got r = {r:.3f}"

    # Verify that the vast majority are identified as reversible responders
    responders = scfa_rescue_df[scfa_rescue_df["rescue_status"].str.contains("Reversible")]
    assert len(responders) >= 12, f"Expected at least 12 reversible responders, found {len(responders)}"

    # Specific gene checks: Plin3 (down in depletion, up in rescue)
    plin3 = scfa_rescue_df[scfa_rescue_df["gene_symbol"] == "Plin3"].iloc[0]
    assert plin3["depletion_meta_log2fc"] < 0
    assert plin3["scfa_rescue_log2fc"] > 0
    assert plin3["in_silico_rescue_index"] > 0

    # Specific gene checks: Tnf (up in depletion, down in rescue)
    tnf = scfa_rescue_df[scfa_rescue_df["gene_symbol"] == "Tnf"].iloc[0]
    assert tnf["depletion_meta_log2fc"] > 0
    assert tnf["scfa_rescue_log2fc"] < 0
    assert tnf["in_silico_rescue_index"] > 0


# ---------------------------------------------------------------------------
# Test Publication Diagnostic Figures
# ---------------------------------------------------------------------------

def test_systems_diagnostic_figures_exist():
    """Validates that all 5 Horizon 4 diagnostic figures exist and are non-empty high-res images."""
    expected_figures = [
        "fig_gsea_pathway_enrichment.png",
        "fig_tf_regulon_landscape.png",
        "fig_wgcna_modules_eigengenes.png",
        "fig_network_hub_subgraph.png",
        "fig_scfa_rescue_inversion.png"
    ]
    for fig_name in expected_figures:
        fig_path = os.path.join(FIGURES_DIR, fig_name)
        assert os.path.exists(fig_path), f"Figure missing: {fig_path}"
        file_size = os.path.getsize(fig_path)
        assert file_size > 50000, f"Figure {fig_name} is too small ({file_size} bytes), may be blank/corrupt"
