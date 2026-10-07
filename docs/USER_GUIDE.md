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

### Step 2.3: Cross-Study Statistical Meta-Analysis (Horizon 3)
```bash
# DerSimonian-Laird Random Effects, Cochran's Q, Higgins I², and LOO sensitivity
python scripts/04_meta_analysis.py
```
*Outputs*: `results/meta_results/microglia_meta_analysis_summary.csv`, `results/meta_results/core_consensus_signature.csv`, `results/meta_results/microglia_meta_analysis_loo.csv`.

### Step 2.4: Systems Biology, Networks & Rescue Modeling (Horizon 4)
```bash
# 1. Cache reference gene set databases (MSigDB Hallmarks, KEGG, TRRUST, phenotypes)
python scripts/05a_cache_gene_sets.py

# 2. Run whole-transcriptome GSEA and hypergeometric ORA
python scripts/05_pathway_enrichment.py

# 3. Deconvolute upstream transcription factor regulons (357 TRRUST TFs)
python scripts/05b_tf_regulon_analysis.py

# 4. Construct WGCNA co-expression network and identify hub genes
python scripts/05c_coexpression_network.py

# 5. In silico SCFA metabolite rescue modeling & signature inversion
python scripts/05d_metabolite_rescue.py

# 6. Render 300 DPI publication diagnostic figures
python scripts/05e_systems_diagnostics.py
```
*Outputs*: `results/pathways/`, `results/networks/`, and 5 high-resolution figures in `results/pathways/figures/`.

---

## 3. Interpreting Statistical & Systems Biology Metrics

### A. Meta-Analysis Statistics (`results/meta_results/microglia_meta_analysis_summary.csv`)
| Metric | Interpretation |
|---|---|
| `meta_log2fc` | Pooled effect size across all cohorts using DerSimonian-Laird random effects. Positive values indicate upregulation in microbiome-depleted or perturbed microglia. |
| `[ci_lower, ci_upper]` | 95% Confidence Interval for the pooled effect size. If the interval excludes 0, the effect is statistically significant at $\alpha = 0.05$. |
| `i2_heterogeneity` | Higgins $I^2$ inconsistency metric ($0 - 100\%$). $I^2 < 25\%$ indicates low across-study heterogeneity; $I^2 > 50\%$ indicates substantial biological or technical heterogeneity across studies. |
| `fdr_random_effects` | Benjamini-Hochberg FDR under random-effects model. |
| `fdr_fisher` | Combined p-value FDR via Fisher's $\chi^2$ method. Prioritizes genes with consistent evidence of differential expression across independent studies. |
| `direction_concordance` | Indicates whether the gene was consistently upregulated (`Concordant Up`), downregulated (`Concordant Down`), or had opposing directions (`Mixed`) across cohorts. |

### B. Upstream TF Regulons (`results/pathways/tf_regulon_activity_summary.csv`)
| Metric | Interpretation |
|---|---|
| `activity_z_score` | Standardized shift in downstream target effect sizes relative to the background transcriptome. Negative $Z$ indicates repressed TF activity (e.g., `Irf1` $Z = -2.28$). |
| `fdr_welch` | Benjamini-Hochberg FDR of Welch's two-sample $t$-test. TFs with $\text{FDR} \le 0.05$ are considered master drivers of the perturbation state. |
| `target_count` | Number of empirically measured target genes evaluated in the regulon ($\ge 5$). |

### C. WGCNA Co-Expression Networks (`results/networks/`)
| Metric | Interpretation |
|---|---|
| `module_name` | Name of the discrete co-expression cluster (e.g., `M_Quiescence`). |
| `k_in` | Intramodular connectivity. High $k_{\text{in}}$ identifies module hub genes that coordinate cluster expression. |
| `module_trait_correlations` | Pearson correlation between the Module Eigengene (ME) and specific experimental traits (GF, ABX, Dietary Fiber Starvation). |

### D. In Silico SCFA Metabolite Rescue (`results/pathways/scfa_metabolite_rescue_modeling.csv`)
| Metric | Interpretation |
|---|---|
| `in_silico_rescue_index` | Direction-adjusted In Silico Rescue Index ($\text{ISRI} = -\operatorname{sign}(\hat{\theta}_{\text{depletion}}) \times \hat{\theta}_{\text{rescue}}$). Values $>0$ indicate successful restoration toward homeostatic baseline. |
| `rescue_percentage` | Proportion of the depletion-induced transcriptomic shift reversed by SCFA administration (clamped between $0\%$ and $100\%$). |
| `rescue_status` | Classification: `Metabolite-Reversible Responder` vs. `Irreversible/Non-responder`. |

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

Execute the automated test suite covering metadata validation, count matrix integrity, random-effects mathematics, and systems biology models:

```bash
# Run all 42 unit tests with verbose reporting
pytest tests/ -v
```

All 42 tests should report `PASSED` in ~1.1 seconds.
