# Changelog

All notable changes to the **NeuroGut-MetaSeq** project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.0] - 2026-10-07

### Added
- **Comprehensive Academic Literature Review** (`docs/LITERATURE_REVIEW.md`):
  - 5,200+ word, 3-pillar scholarly review with 40 formal academic citations (PMIDs/DOIs).
  - Deep coverage of microglial ontogeny, homeostatic checkpoints (*Tmem119*, *Cx3cr1*, *P2ry12*), and microbiome depletion phenotypes (Erny et al. 2015, Thion et al. 2018).
  - Detailed biochemical signaling models: pan-HDAC inhibition, canonical NF-κB repression, GPCR signaling (*Ffar2*, *Ffar3*, *Hcar2*), and the "two-hit" neuroinflammatory priming paradox.
  - Comprehensive guide to computational RNA-seq meta-analysis: Negative Binomial GLMs, random-effects pooling, isolation stress auditing, and TF regulon inference.
- **Multi-Cohort RNA-Seq Ingestion Engine** (`scripts/01_download_geo.py`):
  - Automated downloader for public transcriptomic repositories using NCBI E-utilities and HTTPS FTP.
  - Curated initial dataset manifest (`config/datasets.yaml`) including landmark cohorts:
    - `GSE107925` (Thion et al., *Cell* 2018; Adult SPF vs GF microglia across sexes).
    - `GSE108045` (Thion et al., *Cell* 2018; Adult CTR vs ABX antibiotic-treated microglia).
    - `GSE266602` (Wang et al., 2024; Microbiome-depleted microglia under baseline and neuroinflammatory stress).
    - `GSE186210` (Matt et al., *J Neurosci* 2023; Dietary fiber and SCFA receptor perturbation).
- **Metadata Curation & Harmonization** (`scripts/02_curate_metadata.py`):
  - Standardized metadata schema (`sample_id`, `cohort`, `condition`, `group_label`, `sex`, `tissue`, `sequencing_type`).
  - Gene identifier mapping between Ensembl gene IDs and official MGI symbols.
  - Bundled synthetic Negative Binomial benchmark datasets in `data/demo/` for instantaneous CI and local dry-runs.
- **Dual-Engine Differential Gene Expression**:
  - Turnkey Python implementation (`scripts/03b_pydeseq2_analysis.py`) utilizing `pydeseq2` with negative binomial GLM and Wald testing.
  - Gold-standard R Bioconductor script (`scripts/03_deseq2_analysis.R`) using `DESeq2` and `apeglm`.
- **Hybrid Statistical Meta-Analysis Engine** (`scripts/04_meta_analysis.py`, `scripts/04_meta_analysis.R`):
  - DerSimonian-Laird inverse-variance random-effects pooling for $\log_2\text{FC} \pm 95\%\text{CI}$.
  - Study heterogeneity quantification using Cochran's $Q$ and Higgins $I^2$ metrics.
  - Non-parametric combined significance via Fisher's $\chi^2$ and sample-size weighted Stouffer's Z-test.
  - Benjamini-Hochberg False Discovery Rate (FDR) multiple testing adjustment.
  - Direction concordance evaluation across independent biological studies.
- **Functional Pathway Enrichment Engine** (`scripts/05_pathway_enrichment.py`, `scripts/05_pathway_enrichment.R`):
  - Hypergeometric over-representation analysis (ORA) and ranked gene set enrichment (GSEA).
  - Mapping to KEGG pathways (NF-κB, TNF, Cytokines, Chemokines, Phagosome) and Gene Ontology Biological Processes.
- **Publication Graphics Suite** (`scripts/06_generate_figures.py`):
  - 300 DPI PNG and vector SVG export for 5 publication figures:
    - Figure 1: Cross-cohort Principal Component Analysis.
    - Figure 2: Consensus Meta-Analysis Volcano Plot.
    - Figure 3: Multi-cohort Forest Plots for landmark genes (*Nfkb1*, *Tnf*, *Il1b*, *Cx3cr1*, *Tmem119*, *Trem2*, *Ffar2*, *Ccl2*).
    - Figure 4: Hierarchically clustered heatmap of top consensus genes.
    - Figure 5: Enriched inflammatory and microglial functional pathways.
- **Interactive Scientific Web Paper** (`docs/index.html`, `scripts/07_build_web_paper.py`):
  - Standalone, publication-ready academic article ready for immediate **GitHub Pages** hosting.
  - Embedded **Interactive Gene Explorer** allowing live search across all synthesized genes with dynamic effect sizes and statistical parameter reporting.
  - One-click CSV and vector graphics download hub.
- **FAIR Reproducibility & Automation**:
  - `Dockerfile` containerizing both R and Python environments.
  - `.github/workflows/reproducibility_ci.yml` for continuous testing and automated GitHub Pages deployment.
  - `Makefile` with single-command targets (`make setup`, `make demo`, `make paper`, `make test`, `make docker-build`).
  - Formal academic citation metadata in `CITATION.cff`.
  - Comprehensive technical and mathematical specifications in `SPECIFICATION.md`.
- **Automated Test Suite** (`tests/`):
  - 14 automated unit and integration tests covering metadata schemas, count matrix integrity, mathematical formulas, and web paper building.

### Fixed
- Resolved NumPy 2.0+ deprecation of `np.asfarray` by updating to `np.asarray(..., dtype=float)`.
- Escaped raw string formatting in web paper generator to prevent LaTeX `\url{}` unicode escape collisions.
- Dynamically imported numerical-prefixed script modules in test suites via `importlib.util`.
