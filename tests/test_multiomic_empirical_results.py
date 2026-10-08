#!/usr/bin/env python3
"""
tests/test_multiomic_empirical_results.py
Verification suite for Multi-Omic & Mechanistic Empirical Upgrades (v1.2.0):
1. Test Sex-by-Condition interaction meta-regression and three-tier classification.
2. Test BayesPrism probabilistic single-cell deconvolution and per-cell ISG imputation.
3. Test Tripartite ATAC-seq chromatin accessibility & TOBIAS footprinting dynamics.
4. Test In Silico NicheNet upstream ligand-receptor deconvolution across sender compartments.
5. Test Llgl2/LAT1 metabolic co-expression and Three-Pillar pharmacokinetic BBB model.
"""

import os
import pytest
import numpy as np
import pandas as pd

def test_sex_dimorphism_meta_results():
    """Validates factorial sex-by-condition interaction meta-analysis across 51 samples."""
    csv_path = "results/meta_results/sex_dimorphism_meta_analysis.csv"
    assert os.path.exists(csv_path), f"Missing {csv_path}"
    
    df = pd.read_csv(csv_path)
    assert len(df) >= 20000, f"Expected >= 20,000 common genes, got {len(df)}"
    
    required_cols = [
        "gene_symbol", "pooled_male_log2fc", "pooled_male_se", "pooled_female_log2fc",
        "pooled_female_se", "interaction_log2fc", "interaction_pval", "interaction_fdr",
        "i2_sex_heterogeneity", "sex_dimorphism_tier"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing column {col} in sex dimorphism table"
        
    # Check classification tiers
    tier_counts = df["sex_dimorphism_tier"].value_counts()
    assert "Sex-Shared" in tier_counts, "Must have Sex-Shared tier"
    assert tier_counts["Sex-Shared"] > 10000, "Vast majority of genes must be sex-shared"
    
    # Check Irf1 and Stat1 invariance
    df_indexed = df.set_index("gene_symbol")
    assert df_indexed.loc["Irf1", "sex_dimorphism_tier"] == "Sex-Shared"
    assert df_indexed.loc["Irf1", "i2_sex_heterogeneity"] < 25.0
    assert df_indexed.loc["Stat1", "sex_dimorphism_tier"] == "Sex-Shared"
    assert df_indexed.loc["Stat1", "i2_sex_heterogeneity"] < 25.0
    assert df_indexed.loc["Llgl2", "sex_dimorphism_tier"] == "Sex-Shared"
    
    # Check figure generation
    fig_dir = "results/meta_results/figures"
    for fig in ["fig_sex_stratified_forest", "fig_sex_concordance_scatter"]:
        assert os.path.exists(os.path.join(fig_dir, f"{fig}.png"))
        assert os.path.exists(os.path.join(fig_dir, f"{fig}.svg"))

def test_bayesprism_single_cell_deconvolution():
    """Validates BayesPrism probabilistic deconvolution and per-cell expression imputation."""
    csv_path = "results/pathways/microglia_subpopulation_deconvolution.csv"
    assert os.path.exists(csv_path), f"Missing {csv_path}"
    
    df = pd.read_csv(csv_path)
    assert len(df) == 60, f"Expected 60 samples, got {len(df)}"
    
    prop_cols = [c for c in df.columns if c.startswith("prop_")]
    assert len(prop_cols) == 5, f"Expected 5 states (Homeostatic, IRM, DAM, Cycling, BAM), got {len(prop_cols)}"
    
    # Check sum-to-one constraint and non-negativity
    for _, row in df.iterrows():
        props = [row[c] for c in prop_cols]
        assert all(p >= 0.0 for p in props), "Cell proportions must be non-negative"
        assert np.isclose(sum(props), 1.0, atol=0.05), "Cell proportions must sum to approx 1.0"
        
    # Check per-cell homeostatic imputation columns
    for isg in ["Oas1a", "Stat1", "Gbp2", "Irf1"]:
        col = f"imputed_homeo_{isg}"
        assert col in df.columns, f"Missing imputed per-cell column {col}"
        assert (df[col] >= 0).all(), f"Imputed per-cell expression in {col} must be non-negative"
        
    # Check figures
    fig_dir = "results/pathways/figures"
    assert os.path.exists(os.path.join(fig_dir, "fig_sc_subpopulation_deconvolution.png"))
    assert os.path.exists(os.path.join(fig_dir, "fig_sc_subpopulation_deconvolution.svg"))

def test_epigenomic_atac_footprinting():
    """Validates tripartite microglial ATAC-seq accessibility and TOBIAS footprinting."""
    csv_path = "results/pathways/epigenomic_chromatin_footprinting.csv"
    assert os.path.exists(csv_path), f"Missing {csv_path}"
    
    df = pd.read_csv(csv_path)
    assert len(df) >= 10, f"Expected >= 10 loci, got {len(df)}"
    
    required_cols = [
        "gene_symbol", "transcription_factor_motif", "genomic_region",
        "atac_acc_spf", "atac_acc_depleted", "atac_acc_scfa_repleted",
        "delta_fp_depletion", "delta_fp_scfa_reversal", "chromatin_reversal_pct"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing column {col} in ATAC footprinting table"
        
    df_idx = df.set_index("gene_symbol")
    # Irf1, Stat1, Oas1a must show chromatin footprint collapse and SCFA reopening
    for g in ["Irf1", "Stat1", "Oas1a", "Gbp2"]:
        assert df_idx.loc[g, "delta_fp_depletion"] < -0.20, f"Expected significant footprint collapse for {g}"
        assert df_idx.loc[g, "delta_fp_scfa_reversal"] > 0.20, f"Expected significant SCFA reopening for {g}"
        assert df_idx.loc[g, "chromatin_reversal_pct"] >= 75.0, f"Expected >= 75% reversal for {g}"
        
    fig_dir = "results/pathways/figures"
    assert os.path.exists(os.path.join(fig_dir, "fig_epigenomic_atac_footprinting.png"))
    assert os.path.exists(os.path.join(fig_dir, "fig_epigenomic_atac_footprinting.svg"))

def test_nichenet_ligand_prioritization():
    """Validates in silico NicheNet upstream ligand-receptor deconvolution."""
    csv_path = "results/pathways/nichenet_ligand_prioritization.csv"
    assert os.path.exists(csv_path), f"Missing {csv_path}"
    
    df = pd.read_csv(csv_path)
    assert len(df) >= 10, f"Expected >= 10 evaluated ligands, got {len(df)}"
    
    required_cols = [
        "upstream_ligand", "sender_compartment", "cognate_microglial_receptor",
        "ligand_activity_pearson_r", "overall_priority_rank"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing column {col} in NicheNet table"
        
    # Check top prioritized ligands
    top_3 = df.head(3)["upstream_ligand"].tolist()
    assert any("OMVs" in l or "LPS" in l for l in top_3), "Microbial PAMPs must be prioritized"
    assert any("Ifnb1" in l for l in top_3), "Ifnb1 must be prioritized"
    
    fig_dir = "results/pathways/figures"
    assert os.path.exists(os.path.join(fig_dir, "fig_nichenet_ligand_receptor_network.png"))
    assert os.path.exists(os.path.join(fig_dir, "fig_nichenet_ligand_receptor_network.svg"))

def test_llgl2_metabolic_and_pk_axis():
    """Validates Llgl2-LAT1 co-expression and Three-Pillar pharmacokinetic BBB model."""
    csv_path = "results/pathways/llgl2_lat1_metabolic_coexpression.csv"
    assert os.path.exists(csv_path), f"Missing {csv_path}"
    
    df = pd.read_csv(csv_path)
    assert len(df) >= 5, f"Expected >= 5 evaluated genes, got {len(df)}"
    
    df_idx = df.set_index("target_gene")
    assert "Slc7a5" in df_idx.index, "LAT1 (Slc7a5) must be evaluated in metabolic axis"
    assert "Mtor" in df_idx.index, "mTOR must be evaluated in metabolic axis"
    
    # Slc7a5 must show positive upregulation in depletion (nutrient scavenging)
    assert df_idx.loc["Slc7a5", "depletion_log2fc"] > 0.3, "LAT1 must be upregulated in depletion"
    
    fig_dir = "results/pathways/figures"
    assert os.path.exists(os.path.join(fig_dir, "fig_pharmacokinetic_bbb_metabolic_axis.png"))
    assert os.path.exists(os.path.join(fig_dir, "fig_pharmacokinetic_bbb_metabolic_axis.svg"))
