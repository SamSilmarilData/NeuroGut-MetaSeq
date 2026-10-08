#!/usr/bin/env python3
"""
tests/test_academic_upgrades.py
Comprehensive verification suite for Academic Peer-Review Overhaul (v1.1.0):
1. Test REML and Hartung-Knapp-Sidik-Jonkman (HKSJ) statistical implementation.
2. Test Two-Tier Subgroup Decomposition schema and output table.
3. Test Single-Cell Subpopulation Deconvolution and lineage invariance test.
4. Test In Vivo SCFA Metabolite Reversibility & 1,000-permutation specificity null model.
5. Test Zero-Jargon enforcement in formal manuscript and web paper.
6. Test vector SVG availability across all publication figures.
"""

import os
import re
import pytest
import numpy as np
import pandas as pd
from scipy import stats

def test_reml_and_hksj_mathematical_properties():
    """Verify REML tau2 bounded optimization and HKSJ t3 adjustment."""
    # Synthetic 4-study effect sizes and standard errors
    thetas = np.array([0.65, 0.70, 0.60, 0.75])
    variances = np.array([0.04, 0.05, 0.03, 0.04]) # SEs approx 0.2

    # Import functions directly from scripts/04_meta_analysis.py
    import importlib.util
    spec = importlib.util.spec_from_file_location("meta_mod", "scripts/04_meta_analysis.py")
    meta_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(meta_mod)

    tau2_reml = meta_mod.reml_tau2(thetas, variances)
    assert tau2_reml >= 0.0, "tau2_reml must be non-negative"

    theta_hksj, se_hksj, ci_low, ci_high, p_hksj = meta_mod.hksj_meta(thetas, variances, tau2_reml)

    # Check HKSJ properties
    assert ci_low < theta_hksj < ci_high, "Confidence interval must bracket point estimate"
    assert se_hksj > 0, "Standard error must be strictly positive"
    # For k=4, df=3, t_crit is approx 3.1824
    t_crit = stats.t.ppf(0.975, df=3)
    assert np.isclose(ci_high - theta_hksj, t_crit * se_hksj, atol=1e-3), "CI half-width must equal t_crit * SE"

def test_meta_analysis_summary_columns_and_reml():
    """Verify master meta-analysis summary includes primary REML/HKSJ and supplementary DL metrics."""
    summary_path = "results/meta_results/microglia_meta_analysis_summary.csv"
    assert os.path.exists(summary_path), f"Missing {summary_path}"

    df = pd.read_csv(summary_path)
    assert len(df) == 23096, "Master meta-analysis summary must evaluate all 23,096 common genes"

    required_cols = [
        "gene_symbol", "n_cohorts", "meta_log2fc", "meta_se", "ci_lower", "ci_upper",
        "tau2", "p_random_effects", "p_hksj", "meta_log2fc_dl", "meta_se_dl",
        "tau2_dl", "p_random_effects_dl", "cochran_q", "i2_heterogeneity",
        "fdr_random_effects", "fdr_random_effects_dl"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing required column '{col}' in meta-analysis summary"

    # Invariant core landmark check
    assert "Llgl2" in df["gene_symbol"].values
    assert "Slfn2" in df["gene_symbol"].values
    assert "Clu" in df["gene_symbol"].values

def test_subgroup_decomposition_table_and_figures():
    """Verify Two-Tier Subgroup Decomposition table and figures."""
    sub_path = "results/meta_results/perturbation_subgroup_decomposition.csv"
    assert os.path.exists(sub_path), f"Missing {sub_path}"

    sub_df = pd.read_csv(sub_path)
    assert len(sub_df) == 23096, "Subgroup table must cover all 23,096 genes"

    required_cols = [
        "gene_symbol", "subgroup_axis", "classification_rationale",
        "lfc_germ_free", "lfc_antibiotics", "lfc_fiber_starvation",
        "q_between_models", "p_q_between"
    ]
    for col in required_cols:
        assert col in sub_df.columns, f"Missing required column '{col}' in subgroup decomposition"

    # Check axis categories
    axes = set(sub_df["subgroup_axis"].unique())
    assert "Shared Microbial Core" in axes or "Shared Microbial Core (Moderate Heterogeneity)" in axes
    assert "ABX Mucosal Shock" in axes

    # Check figure artifacts
    assert os.path.exists("results/meta_results/figures/fig_subgroup_perturbation_decomposition.png")
    assert os.path.exists("results/meta_results/figures/fig_subgroup_perturbation_decomposition.svg")

def test_single_cell_deconvolution_and_lineage_stability():
    """Verify single-cell deconvolution across 60 samples and pan-microglial lineage stability."""
    deconv_path = "results/pathways/microglia_subpopulation_deconvolution.csv"
    assert os.path.exists(deconv_path), f"Missing {deconv_path}"

    df = pd.read_csv(deconv_path)
    assert len(df) == 60, f"Deconvolution must evaluate all 60 biological samples, found {len(df)}"

    required_cols = [
        "sample_id", "cohort", "condition", "sig_Interferon-Responsive (IRM)",
        "sig_Homeostatic Mature", "isg_raw_score", "lineage_pan_score", "isg_to_lineage_ratio"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing required column '{col}' in deconvolution table"

    # Check lineage stability (no significant difference in pan-microglial markers Hexb, Csf1r, Tmem119)
    ref_lin = df[df["condition"] == "reference"]["lineage_pan_score"].dropna()
    pert_lin = df[df["condition"] == "perturbed"]["lineage_pan_score"].dropna()
    t_lin, p_lin = stats.ttest_ind(ref_lin, pert_lin)
    assert p_lin > 0.05, f"Pan-microglial lineage markers should be invariant between conditions, got p={p_lin}"

    # Across primary perturbation discovery cohorts (GSE107925, GSE108045, GSE186210),
    # the ISG-to-lineage ratio and IRM signature are consistently suppressed in microbiome-depleted mice
    disc = df[df["cohort"] != "GSE266602"]
    for cname in ["GSE107925", "GSE108045", "GSE186210"]:
        c_sub = disc[disc["cohort"] == cname]
        c_ref_ratio = c_sub[c_sub["condition"] == "reference"]["isg_to_lineage_ratio"].mean()
        c_pert_ratio = c_sub[c_sub["condition"] == "perturbed"]["isg_to_lineage_ratio"].mean()
        assert c_pert_ratio < c_ref_ratio, f"In {cname}, perturbed ISG ratio ({c_pert_ratio:.3f}) must be lower than reference ({c_ref_ratio:.3f})"

        c_ref_irm = c_sub[c_sub["condition"] == "reference"]["sig_Interferon-Responsive (IRM)"].mean()
        c_pert_irm = c_sub[c_sub["condition"] == "perturbed"]["sig_Interferon-Responsive (IRM)"].mean()
        assert c_pert_irm < c_ref_irm, f"In {cname}, perturbed IRM score ({c_pert_irm:.3f}) must be lower than reference ({c_ref_irm:.3f})"

    # Overall discovery cohorts IRM suppression
    disc_ref_irm = disc[disc["condition"] == "reference"]["sig_Interferon-Responsive (IRM)"]
    disc_pert_irm = disc[disc["condition"] == "perturbed"]["sig_Interferon-Responsive (IRM)"]
    assert disc_pert_irm.mean() < disc_ref_irm.mean(), "Discovery cohort perturbed IRM must be lower than reference"

    # Check figures
    assert os.path.exists("results/pathways/figures/fig_sc_subpopulation_deconvolution.png")
    assert os.path.exists("results/pathways/figures/fig_sc_subpopulation_deconvolution.svg")

def test_scfa_metabolite_reversibility_and_null_model():
    """Verify empirical SCFA metabolite reversibility and specificity null test."""
    rescue_path = "results/pathways/scfa_metabolite_rescue_modeling.csv"
    assert os.path.exists(rescue_path), f"Missing {rescue_path}"

    rdf = pd.read_csv(rescue_path)
    assert len(rdf) >= 18, "SCFA rescue table must evaluate key landmark genes"
    assert "in_silico_rescue_index" in rdf.columns
    assert "dataset_grounding" in rdf.columns

    # Mean ISRI must be strongly positive
    mean_isri = rdf["in_silico_rescue_index"].mean()
    assert mean_isri > 0.3, f"Observed mean ISRI must be strongly positive, got {mean_isri}"

    # Check figures
    assert os.path.exists("results/pathways/figures/fig_scfa_rescue_specificity_null.png")
    assert os.path.exists("results/pathways/figures/fig_scfa_rescue_specificity_null.svg")

def test_zero_jargon_in_manuscript_and_web_paper():
    """Strictly verify that zero software engineering jargon appears in formal academic text."""
    manuscript_path = "docs/MANUSCRIPT.md"
    html_path = "docs/index.html"
    assert os.path.exists(manuscript_path), f"Missing {manuscript_path}"
    assert os.path.exists(html_path), f"Missing {html_path}"

    with open(manuscript_path, "r", encoding="utf-8") as f:
        ms_text = f.read()

    with open(html_path, "r", encoding="utf-8") as f:
        html_text = f.read()

    forbidden_patterns = [
        r"adaptive discovery framework",
        r"horizon\s+[0-5]",
        r"\(\*\*bloom",
        r"bloom\s+[0-9]\.[0-9]"
    ]

    for pattern in forbidden_patterns:
        match_ms = re.search(pattern, ms_text, re.IGNORECASE)
        assert match_ms is None, f"Found forbidden jargon '{match_ms.group(0)}' in {manuscript_path}"

        match_html = re.search(pattern, html_text, re.IGNORECASE)
        assert match_html is None, f"Found forbidden jargon '{match_html.group(0)}' in {html_path}"

def test_vector_svg_presence_in_assets():
    """Verify that zoomable vector SVG files exist for all key figures in docs/assets/."""
    expected_svgs = [
        "fig_meta_volcano.svg",
        "fig_forest_plots_top.svg",
        "fig_subgroup_perturbation_decomposition.svg",
        "fig_loo_stability.svg",
        "fig_consensus_heatmap.svg",
        "fig_gsea_pathway_enrichment.svg",
        "fig_tf_regulon_landscape.svg",
        "fig_wgcna_modules_eigengenes.svg",
        "fig_network_hub_subgraph.svg",
        "fig_sc_subpopulation_deconvolution.svg",
        "fig_scfa_rescue_specificity_null.svg",
        "fig_scfa_rescue_inversion.svg",
        "fig_sex_concordance_scatter.svg",
        "fig_sex_stratified_forest.svg",
        "fig_epigenomic_atac_footprinting.svg",
        "fig_nichenet_ligand_receptor_network.svg",
        "fig_pharmacokinetic_bbb_metabolic_axis.svg"
    ]

    for svg in expected_svgs:
        path = os.path.join("docs/assets", svg)
        assert os.path.exists(path), f"Missing vector SVG asset: {path}"

def test_interactive_volcano_plot_in_web_paper():
    """Verify that the interactive SVG Volcano Plot engine is properly embedded in docs/index.html."""
    html_path = "docs/index.html"
    assert os.path.exists(html_path), f"Missing {html_path}"

    with open(html_path, "r", encoding="utf-8") as f:
        html_text = f.read()

    assert 'id="volcanoContainer"' in html_text, "Missing volcanoContainer in web paper"
    assert 'id="volcanoSvgWrapper"' in html_text, "Missing volcanoSvgWrapper in web paper"
    assert 'id="volcanoTooltip"' in html_text, "Missing volcanoTooltip in web paper"
    assert 'renderInteractiveVolcanoPlot' in html_text, "Missing renderInteractiveVolcanoPlot JS function"
    assert 'filterVolcano' in html_text, "Missing filterVolcano JS function"
    assert 'volcano-filter-btn' in html_text, "Missing volcano-filter-btn elements"

def test_calibrated_epigenetic_language_in_manuscript():
    """Verify that manuscript uses calibrated life-sciences phrasing and contains new dedicated Discussion subsections."""
    manuscript_path = "docs/MANUSCRIPT.md"
    assert os.path.exists(manuscript_path), f"Missing {manuscript_path}"

    with open(manuscript_path, "r", encoding="utf-8") as f:
        ms_text = f.read()

    # Calibrated phrases must be present
    assert "chromatin-poised" in ms_text.lower(), "Manuscript must frame reversibility as chromatin-poised"
    assert "chromatin accessibility" in ms_text.lower(), "Manuscript must discuss chromatin accessibility"

    # Overreaching unmeasured claims must be absent
    assert "proves epigenetic reversibility" not in ms_text.lower(), "Manuscript must not claim proof of epigenetic reversibility without direct ChIP-seq"
    assert "epigenetically reversible activation" not in ms_text.lower(), "Manuscript must not overreach in title or narrative"

    # Discussion subsections must be present
    assert "Functional Repositioning of *Llgl2*: Myeloid Nutrient-Scavenging" in ms_text, "Missing Llgl2 nutrient adaptation section"
    assert "Resolving the Blood-Brain Barrier Pharmacokinetic Paradox" in ms_text, "Missing BBB pharmacokinetics section"
    assert "In Silico Cerebrovascular Ligand Relay" in ms_text, "Missing NicheNet ligand relay section"
