"""
tests/test_horizon5_narrative.py
Automated validation of Horizon 5: The Living Narrative.

Verifies:
1. Production web paper compilation, layout, dynamic SVG Forest Plot engine, and inlined gene database.
2. Asset synchronization in docs/assets/ (all 20 publication figures and key result CSVs).
3. Academic preprint manuscript (docs/MANUSCRIPT.md) structure, word count (>= 4,500 words), and sections.
4. Living Research Lab Notebook Entry 005 completeness and ADF adherence.
"""

import os
import pytest

DOCS_DIR = "docs"
ASSETS_DIR = os.path.join(DOCS_DIR, "assets")
INDEX_HTML = os.path.join(DOCS_DIR, "index.html")
MANUSCRIPT_MD = os.path.join(DOCS_DIR, "MANUSCRIPT.md")
NOTEBOOK_MD = os.path.join(DOCS_DIR, "LAB_NOTEBOOK.md")


def test_web_paper_production_build():
    """Validates that docs/index.html is compiled, high-capacity, and contains real consensus data."""
    assert os.path.exists(INDEX_HTML), f"Missing {INDEX_HTML}"
    file_size = os.path.getsize(INDEX_HTML)
    assert file_size >= 50000, f"docs/index.html too small ({file_size} bytes), expected >= 50 KB"

    with open(INDEX_HTML, "r", encoding="utf-8") as f:
        html = f.read()

    # Title & Metadata
    assert "Cross-Study Transcriptomic Meta-Analysis Reveals Basal Tonic Interferon Surveillance Collapse" in html
    assert "v1.0.0 Production Gold Release" in html
    assert "Dynamic Multi-Cohort Forest Plot" in html
    assert "GENE_DATABASE" in html

    # Landmark consensus and bloom genes inlined
    for marker in ["Llgl2", "Slfn2", "Clu", "Fosb", "Tsc22d3", "Ddit4", "Plin3", "Irf1", "Tnf", "Sap30"]:
        assert marker in html, f"Marker gene {marker} missing from index.html"

    # Tabbed gallery hooks
    assert "panel-qc-purity" in html
    assert "panel-glm-abx" in html
    assert "panel-meta-volcano" in html
    assert "panel-sys-gsea" in html


def test_web_assets_integrity():
    """Validates that all 20 publication figures and primary CSV results exist in docs/assets/."""
    assert os.path.exists(ASSETS_DIR), f"Missing {ASSETS_DIR}"

    expected_figures = [
        # QC (Horizon 1)
        "fig_qc_library_depths.png",
        "fig_qc_microglial_purity.png",
        "fig_qc_isolation_stress.png",
        # Cohort GLMs (Horizon 2)
        "fig_volcano_GSE107925.png",
        "fig_volcano_GSE108045.png",
        "fig_volcano_GSE186210.png",
        "fig_volcano_GSE266602.png",
        "fig_marker_effect_sizes_by_cohort.png",
        "fig_lfc_correlation_heatmap.png",
        "fig_cohort_deg_overlap.png",
        # Meta-analysis (Horizon 3)
        "fig_meta_volcano.png",
        "fig_forest_plots_top.png",
        "fig_loo_stability.png",
        "fig_heterogeneity_distribution.png",
        "fig_consensus_heatmap.png",
        # Systems biology & rescue (Horizon 4)
        "fig_gsea_pathway_enrichment.png",
        "fig_tf_regulon_landscape.png",
        "fig_wgcna_modules_eigengenes.png",
        "fig_network_hub_subgraph.png",
        "fig_scfa_rescue_inversion.png"
    ]

    for fig_name in expected_figures:
        fig_path = os.path.join(ASSETS_DIR, fig_name)
        assert os.path.exists(fig_path), f"Missing publication figure asset: {fig_path}"
        assert os.path.getsize(fig_path) > 10000, f"Figure asset {fig_name} is too small"

    expected_csvs = [
        "microglia_meta_analysis_summary.csv",
        "core_consensus_signature.csv",
        "microglia_meta_analysis_loo.csv",
        "gsea_hallmarks_summary.csv",
        "tf_regulon_activity_summary.csv",
        "scfa_metabolite_rescue_modeling.csv",
        "coexpression_module_assignments.csv",
        "hub_genes_summary.csv"
    ]

    for csv_name in expected_csvs:
        csv_path = os.path.join(ASSETS_DIR, csv_name)
        assert os.path.exists(csv_path), f"Missing CSV download asset: {csv_path}"
        assert os.path.getsize(csv_path) > 100, f"CSV asset {csv_name} is empty"


def test_manuscript_structure_and_completeness():
    """Validates that docs/MANUSCRIPT.md exists and adheres to Nature Neuroscience standards."""
    assert os.path.exists(MANUSCRIPT_MD), f"Missing {MANUSCRIPT_MD}"

    with open(MANUSCRIPT_MD, "r", encoding="utf-8") as f:
        text = f.read()

    words = text.split()
    assert len(words) >= 4500, f"Manuscript word count {len(words)} is below comprehensive threshold (4500 words)"

    required_sections = [
        "Abstract",
        "Significance Statement",
        "1. Introduction",
        "2. Results",
        "2.1 Multi-Cohort Harmonization",
        "2.3 Random-Effects Meta-Analysis",
        "2.5 Whole-Transcriptome GSEA",
        "2.6 Upstream Transcription Factor Regulon Deconvolution",
        "2.7 Weighted Gene Co-Expression Networks",
        "2.8 In Silico SCFA Metabolite Modeling",
        "3. Discussion",
        "4. Online Methods",
        "5. Tables & Figures",
        "6. Figure Legends",
        "7. References"
    ]

    for sec in required_sections:
        assert sec in text, f"Missing manuscript section: {sec}"

    # Verify key landmark biological findings discussed
    for bloom_concept in ["Llgl2", "Slfn2", "Clu", "Tsc22d3", "IRF1", "ISRI", "DerSimonian-Laird"]:
        assert bloom_concept in text, f"Key biological concept {bloom_concept} missing from manuscript"


def test_lab_notebook_entry_005():
    """Validates that Entry 005 exists in docs/LAB_NOTEBOOK.md concluding the ADF protocol."""
    assert os.path.exists(NOTEBOOK_MD), f"Missing {NOTEBOOK_MD}"

    with open(NOTEBOOK_MD, "r", encoding="utf-8") as f:
        text = f.read()

    assert "### Entry 005" in text, "Missing Entry 005 in LAB_NOTEBOOK.md"
    assert "Horizon 5: The Living Narrative" in text
    assert "PRODUCTION GOLD RELEASE (v1.0.0)" in text
    assert "1. Target Hypothesis & Dissemination Objective" in text
    assert "5. Project Culmination Retrospective" in text
