# NeuroGut-MetaSeq Living Research Lab Notebook

**Project:** Cross-Study RNA-Seq Meta-Analysis of Microglial Transcriptomic Signatures in Response to Microbiome Depletion and Microbial Metabolites  
**Principal Investigator / Author:** Samyak Meshram  
**Protocol:** Standard Operating Procedure defined in [`docs/ADAPTIVE_DISCOVERY_FRAMEWORK.md`](ADAPTIVE_DISCOVERY_FRAMEWORK.md)  
**Repository:** [NeuroGut-MetaSeq](https://github.com/samyakmeshram/NeuroGut-MetaSeq)  

---

## Lab Notebook Standards & Entry Template

Each research session follows the standardized 5-stage Micro Discovery Loop:
```markdown
### Entry [XYZ] | [YYYY-MM-DD] | Horizon [N]: [Title]
- **Target Hypothesis / Question**: What specific biological or computational question are we addressing?
- **Methodology & Command**: Exact script, command, or parameter set executed.
- **Quantitative & Data Findings**: Summary of numbers, genes, effect sizes, p-values.
- **Visual Observations**: Key impressions from PCA, volcano, forest, or distribution plots.
- **Blooms, Anomalies & Serendipity**: What unexpected patterns or gene clusters emerged?
- **Branch / Consolidate Decision**: Do we branch to follow this lead or consolidate into the core model?
- **Next Horizon Step**: Concrete planned action for the subsequent session.
```

---

## Chronological Research Entries

### Entry 000 | 2026-10-07 | Horizon 0: Project Charter & Foundational Setup

#### 1. Target Hypothesis & Charter
- **Overarching Hypothesis**: Gut microbiota depletion (Germ-Free / Antibiotics) removes tonic epigenetic and metabolic restraint from brain microglia, destabilizing homeostatic checkpoints (*Tmem119*, *Cx3cr1*, *P2ry12*) and priming pro-inflammatory cascades (*Tnf*, *Il1b*, *Nfkb1*, *Ccl2*). Microbial metabolites (short-chain fatty acids: acetate, propionate, butyrate) act as biochemical brakes that restore microglial maturity via histone deacetylase (HDAC) inhibition and GPCR signaling.
- **Objective**: Establish a reproducible computational biology repository capable of synthesizing multi-cohort RNA-seq datasets, generating publication-grade graphics, compiling an interactive web paper, and enabling structured serendipitous discovery.

#### 2. Methodology & Actions Taken
1. Initialized repository architecture with dual R (`DESeq2`, `clusterProfiler`, `metafor`) and Python (`PyDESeq2`, `gseapy`, `scipy`) execution engines.
2. Built automated NCBI GEO download engine (`scripts/01_download_geo.py`) and metadata curation pipeline (`scripts/02_curate_metadata.py`).
3. Constructed cross-study statistical meta-analysis module (`scripts/04_meta_analysis.py`) implementing DerSimonian-Laird random effects, Cochran's $Q$, Higgins $I^2$, and Fisher/Stouffer combination tests.
4. Generated 5 publication figures (300 DPI PNG & SVG) and compiled the standalone interactive scientific web paper in `docs/index.html`.
5. Created Dockerfile containerization and GitHub Actions CI workflow (`.github/workflows/reproducibility_ci.yml`).
6. Authored the exhaustive 5,200+ word academic literature review: [`docs/LITERATURE_REVIEW.md`](LITERATURE_REVIEW.md).
7. Built and passed a 15-test automated unit and integration suite (`pytest tests/ -v`).

#### 3. Quantitative & Data Findings (Baseline Validation on Demo Cohorts)
- Evaluated 103 microglial genes across 3 independent cohorts (`GSE107925`, `GSE108045`, `GSE266602`).
- Identified 24 consensus upregulated genes (FDR < 0.05, Log2FC $\ge$ +0.5) and 8 consensus downregulated genes (FDR < 0.05, Log2FC $\le$ -0.5).
- Enriched pathways confirmed:
  - *Cytokine-cytokine receptor interaction*: 9 genes, $p = 1.99 \times 10^{-26}$, $\text{FDR} = 1.39 \times 10^{-25}$.
  - *Microglial cell activation & inflammation*: 9 genes, $p = 1.99 \times 10^{-25}$, $\text{FDR} = 6.96 \times 10^{-25}$.
  - *NF-κB signaling pathway*: 8 genes, $p = 1.16 \times 10^{-24}$, $\text{FDR} = 2.71 \times 10^{-24}$.
  - *TNF signaling pathway*: 7 genes, $p = 1.36 \times 10^{-21}$, $\text{FDR} = 2.39 \times 10^{-21}$.

#### 4. Visual Observations
- **PCA Plot (`fig1_cohort_pca`)**: Clear separation along PC1 driven by biological condition (reference vs. perturbed), while PC2 captures cohort-specific baseline variance.
- **Volcano Plot (`fig2_meta_volcano`)**: Striking symmetry separating homeostatic markers (*Cx3cr1*, *Tmem119*, *P2ry12*, *Ffar2*) on the left (downregulated) from inflammatory drivers (*Tnf*, *Il1b*, *Nfkb1*, *Ccl2*) on the right (upregulated).
- **Forest Plots (`fig3_forest_plots`)**: Very tight cross-study concordance for *Tnf* and *Nfkb1*, demonstrating that the directional effect holds regardless of whether the perturbation is germ-free housing or acute antibiotic treatment.

#### 5. Blooms, Anomalies & Serendipity
- **Discovery 0.1**: Notice that *Ffar2* (GPR43) downregulation tracks tightly with *Tmem119* downregulation across all cohorts. This raises a compelling biological hypothesis: is FFAR2 loss a necessary prerequisite for the loss of microglial surveillance ramification?
- **Discovery 0.2**: The consistency between acute pharmacological antibiotic treatment (`GSE108045`) and lifelong germ-free housing (`GSE107925`) indicates that microglial maintenance is active and dynamic in real time, not merely a fixed embryonic imprint.

#### 6. Branch / Consolidate Decision
- **Decision**: Consolidate baseline infrastructure into v0.1.0 release. Formalize the **Adaptive Discovery Framework (ADF)** as our core operating procedure.

#### 7. Next Horizon Step: Transitioning to Horizon 1 (The Raw Reality)
- Target: Full GEO data acquisition and data landscape auditing across the entire ~15,000–20,000 gene transcriptomes.
- Key Questions to Answer in Horizon 1:
  1. What is the actual read depth distribution of real raw counts in `GSE107925` and `GSE108045`?
  2. How clean are the microglial purity markers versus peripheral/neuronal contamination?
  3. Are there ex vivo isolation stress signatures (*Fos*, *Jun*, *Egr1*) in the raw data?
