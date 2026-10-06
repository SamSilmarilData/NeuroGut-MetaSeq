# NeuroGut-MetaSeq User Guide

This guide walks you through using, customizing, and deploying the **NeuroGut-MetaSeq** analysis pipeline.

---

## 1. Quickstart Reproduction

To run the complete pipeline and compile the scientific web paper locally:

```bash
# 1. Initialize environment
make setup

# 2. Run the demo analysis pipeline and build the web paper
make demo
```

Once complete, open `docs/index.html` in your web browser:
```bash
open docs/index.html   # On macOS
# or
xdg-open docs/index.html # On Linux
```

---

## 2. Deploying to GitHub Pages

The web paper is pre-configured for zero-configuration deployment with **GitHub Pages**:

1. Push your repository to GitHub:
   ```bash
   git push origin main
   ```
2. Navigate to your repository settings on GitHub:
   - Go to **Settings** > **Pages**.
   - Under **Build and deployment**, select **Source: GitHub Actions** (or select **Deploy from a branch** and choose `/docs` folder on `main`).
3. Your interactive scientific paper will be live at:
   `https://<your-username>.github.io/NeuroGut-MetaSeq/`

---

## 3. Adding a New RNA-Seq Dataset

You can integrate additional GEO accessions or private RNA-seq experiments by following these steps:

### Step 3.1: Update `config/datasets.yaml`
Add your dataset entry under `cohorts`:

```yaml
cohorts:
  GSE_CUSTOM:
    title: "Transcriptome analysis of microglia under novel metabolite treatment"
    publication: "Researcher et al., 2026"
    organism: "Mus musculus"
    tissue: "Microglia"
    sequencing_type: "Bulk RNA-seq"
    contrast_type: "treatment"
    reference_group: "Vehicle"
    treatment_group: "Metabolite_X"
    files:
      counts:
        url: "https://.../counts.txt.gz"
        local_path: "data/raw/GSE_CUSTOM/counts.txt.gz"
      series_matrix:
        url: "https://.../series_matrix.txt.gz"
        local_path: "data/raw/GSE_CUSTOM/series_matrix.txt.gz"
```

### Step 3.2: Curate Counts and Metadata
Add a parsing routine in `scripts/02_curate_metadata.py` to produce:
1. `data/processed/GSE_CUSTOM_counts.csv`
   - Must contain columns `gene_id`, `gene_symbol`, and one column per sample with integer counts.
2. `data/metadata/GSE_CUSTOM_metadata.csv`
   - Must contain columns: `sample_id`, `cohort`, `condition` (`reference` or `perturbed`), `group_label`, `sex`, `tissue`, `sequencing_type`.

### Step 3.3: Re-Run the Pipeline
```bash
python scripts/03b_pydeseq2_analysis.py
python scripts/04_meta_analysis.py
python scripts/05_pathway_enrichment.py
python scripts/06_generate_figures.py
python scripts/07_build_web_paper.py
```

---

## 4. Customizing Statistical Parameters

All statistical thresholds and database choices are configured in `config/analysis_params.yaml`:

- **Significance Cutoffs**:
  ```yaml
  differential_expression:
    alpha: 0.05
    log2fc_threshold: 0.585 # 1.5-fold change threshold
  ```
- **Meta-Analysis Models**:
  ```yaml
  meta_analysis:
    method: "random_effects" # DerSimonian-Laird inverse-variance
    significance_cutoffs:
      meta_fdr: 0.05
      meta_log2fc: 0.5
  ```
- **Pathway Databases**:
  Adjust target gene set databases under `pathway_enrichment.gene_set_databases`.

---

## 5. Interpreting Meta-Analysis Statistics

The primary output table is stored in `results/meta_results/microglia_meta_analysis_summary.csv`. Here is how to interpret key statistical columns:

| Metric | Interpretation |
|---|---|
| `meta_log2fc` | Pooled effect size across all cohorts using DerSimonian-Laird random effects. Positive values indicate upregulation in microbiome-depleted or perturbed microglia. |
| `[ci_lower, ci_upper]` | 95% Confidence Interval for the pooled effect size. If the interval excludes 0, the effect is statistically significant at $\alpha = 0.05$. |
| `i2_heterogeneity` | Higgins $I^2$ inconsistency metric ($0 - 100\%$). $I^2 < 25\%$ indicates low across-study heterogeneity; $I^2 > 50\%$ indicates substantial biological or technical heterogeneity across studies. |
| `fdr_fisher` | Benjamini-Hochberg adjusted p-value combining study significance via Fisher's $\chi^2$ method. Prioritizes genes with consistent evidence of differential expression across independent studies. |
| `direction_concordance` | Indicates whether the gene was consistently upregulated (`Concordant Up`), downregulated (`Concordant Down`), or had opposing directions (`Mixed`) across cohorts. |
