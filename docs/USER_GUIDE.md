# NeuroGut-MetaSeq User Guide

This guide walks you through using, customizing, and reproducing the **NeuroGut-MetaSeq** analysis pipeline across both lightweight demo mode and full multi-cohort production execution.

---

## 1. Dual Execution Modes

NeuroGut-MetaSeq provides two synchronized execution tracks:

### Track 1: ⚡ 1-Minute Quickstart (Demo Mode)
Best for rapid smoke testing, local development, and CI environments without downloading gigabytes of raw data:

```bash
# 1. Initialize environment
make setup
source .venv/bin/activate

# 2. Run the demo analysis pipeline and build the web paper (~10 seconds)
make demo
```

Once complete, open `docs/index.html` in your web browser:
```bash
open docs/index.html    # macOS
# or
xdg-open docs/index.html # Linux
```

### Track 2: 🧬 Full Scientific Reproduction (Real 60-Sample Cohorts)
Best for reproducing the complete published meta-analysis across all 60 biological samples, 23,096 common genes, upstream regulons, and metabolite rescue models:

```bash
# Execute the full production pipeline end-to-end
make production
```

Or execute horizon-by-horizon:
```bash
# Horizon 1: Data Curation & QC Audit (Lineage purity >99%, dissociation stress test)
make horizon1

# Horizon 2: Differential Expression GLMs (PyDESeq2 across all 4 cohorts)
make horizon2

# Horizon 3: Statistical Meta-Analysis (Random Effects & LOO Sensitivity)
make horizon3

# Horizon 4: Systems Biology, Regulons, WGCNA & SCFA Rescue
make horizon4
```

---

## 2. Step-by-Step Command-Line Execution

For fine-grained command-line control, run the modular Python scripts directly:

### Step 2.1: Data Ingestion & Quality Control (Horizon 1)
```bash
# Harmonize count matrices and extract standardized metadata
python scripts/02_curate_metadata.py

# Execute microglial purity and ex vivo dissociation stress audit
python scripts/02b_qc_audit.py
```
*Outputs*: `data/metadata/*_metadata.csv`, `data/processed/*_counts.csv`, `results/qc/` diagnostic reports.

### Step 2.2: Cohort-Level Differential Expression (Horizon 2)
```bash
# Negative binomial generalized linear models controlling for biological sex
python scripts/03b_pydeseq2_analysis.py
```
*Outputs*: `results/de_results/*_deg.csv`, `results/figures/fig_volcano_*.png`.

### Step 2.3: Cross-Study Statistical Meta-Analysis & Subgroup Decomposition
```bash
# Restricted Maximum Likelihood (REML) with Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment,
# DerSimonian-Laird benchmark, Cochran's Q, Higgins I², and Leave-One-Out (LOO) sensitivity
python scripts/04_meta_analysis.py

# Subgroup meta-regression and two-tier biological decomposition (GF vs ABX vs Fiber Starvation)
python scripts/04c_perturbation_subgroups.py
```
*Outputs*: `results/meta_results/microglia_meta_analysis_summary.csv`, `results/meta_results/core_consensus_signature.csv`, `results/meta_results/microglia_meta_analysis_loo.csv`, `results/meta_results/perturbation_subgroup_decomposition.csv`, and subgroup forest plots.

### Step 2.4: Systems Biology, Networks, Deconvolution & Rescue Modeling
```bash
# 1. Cache reference gene set databases (MSigDB Hallmarks, KEGG, TRRUST, phenotypes)
python scripts/05a_cache_gene_sets.py

# 2. Run whole-transcriptome GSEA and hypergeometric ORA
python scripts/05_pathway_enrichment.py

# 3. Deconvolute upstream transcription factor regulons (357 TRRUST TFs)
python scripts/05b_tf_regulon_analysis.py

# 4. Construct WGCNA co-expression network and identify hub genes
python scripts/05c_coexpression_network.py

# 5. In vivo SCFA metabolite rescue modeling (GSE64977) & 1,000-permutation null model
python scripts/05d_metabolite_rescue.py

# 6. BayesPrism 5-state subpopulation deconvolution & ISG cell-intrinsic normalization
python scripts/05f_single_cell_deconvolution.py

# 7. Render publication diagnostic figures (PNG & SVG vector formats)
python scripts/05e_systems_diagnostics.py
```
*Outputs*: `results/pathways/`, `results/networks/`, and publication-quality raster (300 DPI PNG) and vector (SVG) figures in `results/pathways/figures/` and `results/figures/`.

### Step 2.5: Multi-Omic & Mechanistic Extensions (v1.2.0)
```bash
# 1. Factorial sex-interaction meta-regression across 51 sex-informative samples
python scripts/03d_sex_dimorphism_analysis.py

# 2. Tripartite microglial ATAC-seq peak accessibility and TOBIAS footprinting (Erny 2021)
python scripts/05g_epigenomic_footprinting.py

# 3. In silico NicheNet cerebrovascular ligand-receptor prioritization
python scripts/05h_ligand_receptor_nichenet.py

# 4. Myeloid Llgl2-LAT1 amino acid nutrient sensing and Three-Pillar BBB flux modeling
python scripts/05i_metabolic_llgl2_and_pharmacokinetics.py
```
*Outputs*: `results/meta_results/sex_dimorphism_meta_analysis.csv`, `results/pathways/epigenomic_chromatin_footprinting.csv`, `results/pathways/nichenet_ligand_prioritization.csv`, `results/pathways/llgl2_lat1_metabolic_coexpression.csv`, and publication graphics (PNG + SVG).

---

## 3. Interpreting Statistical & Multi-Omic Metrics

### A. Meta-Analysis Statistics (`results/meta_results/microglia_meta_analysis_summary.csv`)
| Metric | Interpretation |
|---|---|
| `meta_log2fc` | Primary pooled effect size across cohorts estimated via Restricted Maximum Likelihood (REML). Positive values indicate upregulation in microbiome-depleted microglia. |
| `[ci_lower, ci_upper]` | Hartung-Knapp-Sidik-Jonkman (HKSJ) adjusted 95% Confidence Interval ($t_{k-1}$ distribution). Protects against false positives under small cohort count ($k=4$). |
| `tau2` | REML between-study variance parameter ($\tau^2$). |
| `p_hksj` | Robust p-value under the HKSJ $t$-distribution adjustment. |
| `meta_log2fc_dl` | Benchmark DerSimonian-Laird pooled effect size, provided for historical and methodological comparison. |
| `i2_heterogeneity` | Higgins $I^2$ inconsistency metric ($0 - 100\%$). $I^2 < 25\%$ indicates low across-study heterogeneity; $I^2 > 50\%$ indicates substantial biological or technical heterogeneity across studies. |
| `fdr_random_effects` | Benjamini-Hochberg FDR under the primary REML random-effects model. |
| `fdr_fisher` | Combined p-value FDR via Fisher's $\chi^2$ method. Prioritizes genes with consistent evidence of differential expression across independent studies. |
| `direction_concordance` | Indicates whether the gene was consistently upregulated (`Concordant Up`), downregulated (`Concordant Down`), or had opposing directions (`Mixed`) across cohorts. |

### B. Factorial Sex-Dimorphism Statistics (`results/meta_results/sex_dimorphism_meta_analysis.csv`)
| Metric / Column | Interpretation |
|---|---|
| `interaction_log2fc` | Pooled sex-by-condition interaction effect ($\hat{\theta}_{\text{int}}$). Positive values indicate male-biased vulnerability to depletion. |
| `interaction_pval` | Two-tailed p-value testing for sex dimorphism. $p \ge 0.05$ indicates sex invariance. |
| `i2_sex_heterogeneity` | Higgins $I^2$ quantifying heterogeneity of the sex interaction across informative cohorts. |
| `sex_dimorphism_tier` | Classification: `Sex-Shared` (99.10% of genome, $I^2 = 0\%$), `Male-Biased Vulnerability` (0.32%, including *Slfn2* and *Oas1a*), or `Female-Biased Vulnerability` (0.08%). |

### C. BayesPrism Subpopulation Deconvolution (`results/pathways/microglia_subpopulation_deconvolution.csv`)
| Metric | Interpretation |
|---|---|
| `sig_Interferon-Responsive (IRM)` | Inferred proportion / score of the IRM subpopulation. Invariance ($17.5\%$ depleted vs $16.8\%$ ref) proves lack of physical IRM depletion. |
| `condition_index_kappa` | SVD condition index ($\kappa = 1.54 < 30$), confirming zero harmful collinearity between reference cell states. |
| `lineage_pan_score` | Mean CPM expression score of invariant pan-microglial lineage markers (*Hexb*, *Csf1r*, *Tmem119*, $p = 0.85$). |
| `isg_to_lineage_ratio` | Normalized ratio of ISG expression to microglial lineage markers ($p = 4.29 \times 10^{-6}$). Proves uniform per-cell ISG downregulation. |

### D. Epigenomic ATAC Footprinting (`results/pathways/epigenomic_chromatin_footprinting.csv`)
| Metric | Interpretation |
|---|---|
| `tobias_fp_depth_*` | Open chromatin transcription factor binding footprint depth in SPF, Depleted, and SCFA states. |
| `delta_fp_depletion` | Footprint collapse induced by microbiota depletion ($\text{FP}_{\text{dep}} - \text{FP}_{\text{SPF}}$). |
| `chromatin_reversal_pct` | Percentage restoration of open chromatin footprint upon SCFA treatment (e.g., 88.9% for *Irf1* and *Stat1*). |

### E. NicheNet Upstream Ligand Prioritization (`results/pathways/nichenet_ligand_prioritization.csv`)
| Metric | Interpretation |
|---|---|
| `ligand_activity_pearson_r` | Correlation between prior ligand regulatory potential and observed microglial target expression shifts. |
| `interferon_regulon_potency` | Specific regulatory potential score of candidate ligand against the 23 IRF1 downstream targets. |

### F. Myeloid Llgl2-LAT1 Nutrient Sensing (`results/pathways/llgl2_lat1_metabolic_coexpression.csv`)
| Metric | Interpretation |
|---|---|
| `pearson_r_with_llgl2` | Co-expression correlation with basolateral adapter *Llgl2* across the 60 samples (e.g., *Slc7a5* $r = +0.612$, *Mtor* $r = -0.784$). |
| `depletion_log2fc` | Meta-analytic log2 fold change demonstrating coordinated LAT1 upregulation and mTOR shutdown. |

---

## 4. Deploying the Interactive Paper to GitHub Pages

The web paper is pre-configured for automated deployment with **GitHub Pages**:

1. Push your repository to GitHub:
   ```bash
   git push origin main
   ```
2. In your repository on GitHub:
   - Navigate to **Settings** > **Pages**.
   - Under **Build and deployment**, select **Source: Deploy from a branch**.
   - Select the `main` branch and `/docs` folder.
3. Your interactive scientific paper will be live at:
   `https://<your-username>.github.io/NeuroGut-MetaSeq/`

---

## 5. Automated Verification & Testing

Execute the automated test suite covering metadata validation, count matrix integrity, random-effects mathematics, and multi-omic models:

```bash
# Run all 58 unit tests with verbose reporting
pytest tests/ -v
```

All 58 tests should report `PASSED` in ~1.3 seconds across 12 test modules:
- `test_multiomic_empirical_results.py` (Factorial sex interaction, BayesPrism $\kappa$, ATAC footprint %, NicheNet ranking, *Llgl2*-LAT1 co-expression)
- `test_academic_upgrades.py` (REML, HKSJ, subgroups, deconvolution, vector SVGs, zero-jargon)
- `test_count_matrix_integrity.py` (count matrix structure & types)
- `test_data_audit.py` (microglial purity & dissociation stress tests)
- `test_de_results.py` (DESeq2 tables and volcanic distributions)
- `test_framework_doc.py` (framework integrity)
- `test_horizon5_narrative.py` (web assets & manuscript structure)
- `test_meta_analysis_real.py` (full meta-analysis scale & bounds)
- `test_metadata_schema.py` (experimental sample attributes)
- `test_statistical_pipeline.py` (mathematical algorithm unit checks)
- `test_systems_biology_real.py` (pathways, regulons, WGCNA, rescue)
- `test_web_paper_build.py` (DOM integrity and literature review depth)
