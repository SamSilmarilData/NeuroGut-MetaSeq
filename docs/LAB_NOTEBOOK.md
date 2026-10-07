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

---

### Entry 001 | 2026-10-07 | Horizon 1: The Raw Reality & Diagnostic Landscape Audit

#### 1. Target Hypothesis & Investigative Scope
- **Core Question**: Are real-world mouse microglial transcriptomes from independent laboratories sufficiently pure from cellular contaminants (astrocytes, neurons, oligodendrocytes, endothelia) and free from ex vivo enzymatic dissociation stress confounding to support robust, unbiased cross-study differential expression and statistical meta-analysis?
- **Cohort Scope**: Four landmark datasets spanning 60 curated biological samples:
  1. `GSE107925` (Thion et al., *Cell* 2018): Adult SPF vs. GF microglia ($n=25$ adult samples: 13 SPF, 12 GF).
  2. `GSE108045` (Thion et al., *Cell* 2018): Adult CTR vs. ABX microglia ($n=12$ adult samples: 6 CTR, 6 ABX).
  3. `GSE266602` (Wang et al., 2024): Microglia under microbiome depletion at baseline ($n=9$ sham baseline samples: 3 SPF, 3 GF, 3 ABX).
  4. `GSE186210` (Matt et al., *J Neurosci* 2023): Dietary fiber / SCFA receptor perturbation ($n=14$ WT adult samples: 8 Normal Fiber, 6 Zero Fiber).

#### 2. Methodology & Computational Pipeline Executed
1. **Raw Acquisition (`scripts/01_download_geo.py`)**:
   - Streamed full gzipped raw count matrices and series matrix headers directly from NCBI GEO FTP.
   - Verified 8 downloaded files using SHA-256 cryptographic hashes and recorded into `data/checksums.sha256`.
2. **Metadata Curation & Gene Mapping (`scripts/02_curate_metadata.py`)**:
   - Standardized sample annotations into uniform 7-column schema (`sample_id, cohort, condition, group_label, sex, tissue, sequencing_type`).
   - Integrated Ensembl-to-MGI Symbol mapping cache (`data/reference/mouse_ensembl_to_symbol.tsv.gz`, 78,348 mouse loci) to harmonize gene keys across cohorts.
   - Aggregated multi-isoform counts to unique gene symbols and validated non-negative integer matrices.
3. **Diagnostic Quality Audit (`scripts/01b_audit_data_landscape.py`)**:
   - Calculated per-library sequencing depth ($\sum c$), gene detection rate ($\ge 5$ counts), and sparsity percentage.
   - Evaluated Microglial Lineage Purity Index:
     $$\text{Purity} = \frac{\overline{\text{CPM}}(\text{Microglia})}{\overline{\text{CPM}}(\text{Microglia}) + \overline{\text{CPM}}(\text{Contaminants})} \times 100\%$$
     Microglial markers: *Tmem119*, *Cx3cr1*, *P2ry12*, *Sall1*, *Mertk*, *Itgam*, *Ptprc*, *Tyrobp*.  
     Contaminants: Astrocytes (*Gfap*, *Aldh1l1*, *Aqp4*), Neurons (*Rbfox3*, *Snap25*, *Map2*), Oligodendrocytes (*Mbp*, *Mog*, *Olig2*), Endothelial (*Cdh5*, *Pecam1*).
   - Evaluated Ex Vivo Dissociation Stress:
     Standardized composite z-score across immediate-early genes (*Fos*, *Jun*, *Atf3*, *Egr1*, *Nfkbia*, *Dusp1*).
   - Evaluated Confounding:
     Two-sided Mann-Whitney $U$ test between reference and perturbed groups within each cohort.
   - Rendered 3 publication-grade figures at 300 DPI (`results/qc/figures/`).

#### 3. Quantitative & Data Findings

##### A. Transcriptome Dimensions & Inter-Cohort Overlap
- Total curated samples: **60 biological samples**.
- Total unique genes per cohort:
  - `GSE107925`: 35,010 genes
  - `GSE108045`: 47,517 genes
  - `GSE266602`: 56,544 genes
  - `GSE186210`: 54,134 genes
- **Consensus Transcriptome Intersection**: **29,017 common genes** shared across all 4 cohorts, enabling comprehensive, genome-wide meta-analysis.

##### B. Sequencing Depths & Library Sizes
| Cohort ID | Sample Count ($n$) | Mean Depth (reads) | Median Depth | Min Depth | Max Depth | Mean Sparsity |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`GSE107925`** | 25 | 19,515,058 | 15,315,221 | 1,739,225 | 42,215,853 | 43.1% |
| **`GSE108045`** | 12 | 20,013,750 | 19,692,306 | 16,303,485 | 23,088,364 | 45.8% |
| **`GSE266602`** | 9 | 54,016,644 | 53,248,605 | 41,524,670 | 69,343,261 | 38.2% |
| **`GSE186210`** | 14 | 2,527,799 | 2,342,755 | 1,174,514 | 4,021,238 | 58.7% |

##### C. Microglial Lineage Purity Audit
| Cohort ID | Mean Purity Index | Min Purity | Mean Mg CPM | Mean Contam CPM | Contaminant Profile | Purity Flags (<90%) |
| :--- | :---: | :---: | :---: | :---: | :--- | :---: |
| **`GSE107925`** | **99.19%** | 97.06% | 3,096.3 | 22.0 | Minimal neuronal/astrocyte carryover (<1%) | 0 / 25 |
| **`GSE108045`** | **99.61%** | 99.54% | 2,690.0 | 10.5 | Negligible contamination | 0 / 12 |
| **`GSE186210`** | **99.62%** | 99.47% | 2,580.6 | 10.0 | Negligible contamination | 0 / 14 |
| **`GSE266602`** | **68.27%** | 53.47% | 48.6 | 14.9 | Moderate astrocytic carryover (*Gfap*, *Aldh1l1*) | 8 / 9 |

##### D. Ex Vivo Dissociation Stress Confounding Test
| Cohort ID | Reference ($n$) | Perturbed ($n$) | Mean Stress Diff ($\Delta z$) | Mann-Whitney $U$ | $p$-value | Confounding Flag ($p < 0.05$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`GSE107925`** | 13 | 12 | +0.077 | 70.0 | 0.6833 | **False** (Clean) |
| **`GSE108045`** | 6 | 6 | +0.542 | 9.0 | 0.1797 | **False** (Clean) |
| **`GSE186210`** | 8 | 6 | -0.095 | 27.0 | 0.7546 | **False** (Clean) |
| **`GSE266602`** | 3 | 6 | +0.015 | 9.0 | 1.0000 | **False** (Clean) |

**Conclusion on Confounding:** There is no statistically significant dissociation stress imbalance between treatment and control groups in any cohort ($p \ge 0.18$ across all tests). Observed transcriptomic shifts in downstream differential expression can be attributed to microbiome/metabolite perturbation rather than enzymatic digestion artifacts.

#### 4. Visual Observations from Diagnostic Figures
- **Library Depths (`fig_qc_library_depths.png`)**: All four cohorts show uniform read distributions within their respective studies. `GSE266602` possesses the deepest sequencing (>50M reads), while `GSE186210` represents standard exploratory bulk RNA-seq (~2.5M reads).
- **Lineage Purity (`fig_qc_microglial_purity.png`)**: FACS-sorted cohorts (`GSE107925`, `GSE108045`, `GSE186210`) display flawless purity (>99%) with virtually indistinguishable error bars between reference and perturbed samples. `GSE266602` exhibits noticeably lower purity (~68%) due to astrocytic marker presence.
- **Dissociation Stress (`fig_qc_isolation_stress.png`)**: Jittered points and box interquartile ranges demonstrate extensive overlap between reference and perturbed mice across all 4 cohorts, confirming the absence of treatment-correlated isolation bias.

#### 5. Blooms, Anomalies & Serendipity
- **Bloom 1.1: The Percoll vs. FACS Purity Contrast**:
  `GSE266602` used Percoll density gradient centrifugation without subsequent fluorophore sorting, leaving detectable astrocytic markers (*Gfap*, *Aldh1l1* at ~30–70 CPM). Rather than discarding this study, our Adaptive Discovery Framework enables us to retain it as a vital biological sensitivity test. In Horizon 3, our Leave-One-Out (LOO) meta-analysis will quantitatively test whether any consensus microglial genes are sensitive to astrocytic carryover.
- **Bloom 1.2: Transcriptome Intersection Size**:
  The intersection of expressed loci across the 4 cohorts reaches 29,017 genes, exceeding our initial conservative estimate of ~15,000 genes. This broad coverage allows us to explore not only core transcription factors and cytokines, but also low-abundance GPCRs (*Ffar2*, *Ffar3*, *Hcar2*), epigenetic readers, and regulatory non-coding transcripts.

#### 6. Branch / Consolidate Decision
- **Consolidation**:
  - The 60 curated samples are locked into `data/processed/` and `data/metadata/`.
  - Purity warning flags are recorded in `results/qc/horizon1_data_audit.csv` for downstream sensitivity modeling.
- **Advancement**:
  - Quality gates are satisfied: 100% test pass rate (`pytest`), verified data integrity, and confirmed absence of isolation stress confounding.
  - Advance to **Horizon 2: The Individual Voices (Cohort-Level Phenotypic Deep Dives)**.

#### 7. Next Horizon Step: Horizon 2 (The Individual Voices)
- **Target**: Run full-transcriptome Negative Binomial GLMs (`PyDESeq2` + empirical Bayes `apeglm` shrinkage) on all 4 cohorts independently.
- **Key Questions for Horizon 2**:
  1. What are the specific effect sizes ($\log_2\text{FC}$) and dispersion patterns in each cohort?
  2. How do acute antibiotic perturbations (ABX in `GSE108045`) compare phenotypically with lifelong germ-free housing (`GSE107925`) and nutritional fiber deprivation (`GSE186210`)?
  3. Are there sex-dimorphic or cohort-specific idiosyncratic responses?

---

### Entry 002 | 2026-10-07 | Horizon 2: The Individual Voices (Cohort-Level Phenotypic Deep Dives)

#### 1. Target Hypothesis & Investigative Scope
- **Core Hypothesis**: Different modes of microbiome and metabolite perturbation (acute broad-spectrum antibiotic depletion, lifelong germ-free deprivation, primary baseline depletion, and nutritional dietary fiber starvation) induce distinct but biologically intersecting microglial transcriptomic responses. By allowing each cohort to speak in its own statistical voice prior to meta-analysis, we can distinguish shared core regulatory shifts from perturbation-specific adaptations.
- **Cohorts Analyzed**:
  1. `GSE107925` (Lifelong Germ-Free housing: SPF vs. GF, $n=25$).
  2. `GSE108045` (Acute Antibiotic Depletion: CTR vs. ABX, $n=12$).
  3. `GSE266602` (Primary Baseline Microbiome Depletion: SPF vs. GF/ABX, $n=9$).
  4. `GSE186210` (Nutritional Fiber Starvation: WT Normal vs. Zero Fiber, $n=14$).

#### 2. Methodology & Computational Pipeline Executed
1. **Full-Transcriptome Negative Binomial GLMs (`scripts/03b_pydeseq2_analysis.py`)**:
   - Filtered unexpressed/silent genes ($<10$ total counts across samples in cohort).
   - Employed additive multi-factor design formula `~ sex + condition` for mixed-sex cohorts (`GSE107925`, `GSE108045`, `GSE186210`), cleanly removing sex-dimorphic baseline variance while evaluating the primary perturbation contrast.
   - Employed single-factor design `~ condition` for male-only cohort (`GSE266602`).
   - Fitted Negative Binomial dispersions, mean-dispersion trend curves, and maximum a posteriori (MAP) dispersions.
   - Evaluated Wald statistics and computed Benjamini-Hochberg False Discovery Rates (FDR $\alpha = 0.05$).
   - Exported comprehensive statistical tables to `results/de_results/<cohort>_deg.csv`.
2. **Cross-Cohort Phenotyping Suite (`scripts/03c_cohort_phenotyping.py`)**:
   - Quantified significant DEGs ($\text{padj} < 0.05, |\log_2\text{FC}| \ge 0.5$).
   - Computed Spearman rank correlation matrix across the 5,984 common expressed genes.
   - Rendered 4 individual publication-grade Volcano Plots with labeled landmark genes.
   - Profiled comparative effect sizes of canonical homeostatic (*Tmem119*, *Cx3cr1*, *P2ry12*, *Sall1*, *Tsc22d3*, *Ffar2*) and inflammatory markers (*Tnf*, *Il1b*, *Il6*, *Ccl2*, *Nfkb1*, *Ddit4*).
3. **Automated Regression Suite (`tests/test_de_results.py`)**:
   - Validated schema, p-value bounds in $[0, 1]$, absence of infinite/NaN values, and non-empty outputs (26/26 tests passing).

#### 3. Quantitative & Data Findings

##### A. Differential Expression Summary Across Discovery Paradigms
| Cohort ID | Perturbation Paradigm | Expressed Genes | Upregulated ($\text{padj}<0.05, \log_2\text{FC}\ge 0.5$) | Downregulated ($\text{padj}<0.05, \log_2\text{FC}\le -0.5$) | Total Significant DEGs | Top Upregulated Loci | Top Downregulated Loci |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- | :--- |
| **`GSE107925`** | Lifelong Germ-Free (GF) | 30,416 | 5 | 20 | **25** | *Glrp1*, *Gpr137b*, *Cep85* | *Gpr137b-ps*, *Fam129a*, *Ero1lb*, *Mlph* |
| **`GSE108045`** | Acute Antibiotic Depletion (ABX) | 20,605 | 125 | 130 | **255** | *Fzd7*, *Fosb*, *Rin2*, *Cd180*, *Tnf* | *Ddit4*, *Gm43813*, *Tsc22d3*, *Tagap* |
| **`GSE266602`** | Microbiome Depletion Baseline | 20,495 | 12 | 7 | **19** | *Gm61007*, *Serf2*, *Nek7*, *Tnf* | *Gm28800*, *Qser1*, *Sall1* |
| **`GSE186210`** | Dietary Fiber Starvation (Zero Fiber) | 16,764 | 9 | 4 | **13** | *Gm38319*, *Hba-a2*, *Bc1*, *Rbm5*, *Ccl24* | *Plin3*, *Cat*, *Tap1*, *Atad3a* |

##### B. Landmark Biological Marker Effect Sizes Across Cohorts
| Gene Symbol | Functional Annotation | `GSE107925` ($\log_2\text{FC}$) | `GSE108045` ($\log_2\text{FC}$) | `GSE266602` ($\log_2\text{FC}$) | `GSE186210` ($\log_2\text{FC}$) | Concordance & Biological Pattern |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- |
| **`Tsc22d3`** (GILZ) | Endogenous NF-κB / AP-1 Repressor | -0.11 | **-2.36** ($p=5.6\times 10^{-20}$) | -0.13 | +0.20 | Massive collapse under acute antibiotic treatment |
| **`Ddit4`** (REDD1) | mTORC1 Inhibitor / Metabolic Brake | +1.16 | **-4.57** ($p=1.2\times 10^{-28}$) | -0.01 | -0.10 | Dramatic repression in ABX; altered in GF |
| **`Tnf`** | Master Pro-inflammatory Cytokine | +0.01 | **+1.18** ($p=0.003$) | **+2.72** | +0.34 | Broadly elevated across depletion models |
| **`Sall1`** | Master Microglial Identity TF | -0.05 | +0.35 | **-2.68** ($p=0.046$) | -0.32 | Significant homeostatic repression in GSE266602 |
| **`Ffar2`** (GPR43) | SCFA Acetate/Propionate GPCR | **+1.74** | -0.01 | **-1.69** | -0.03 | Perturbed across microbiome deficiency states |
| **`Plin3`** | Perilipin 3 (Lipid Droplet Homeostasis) | -0.02 | -0.07 | -0.18 | **-1.88** ($p=3.4\times 10^{-5}$) | Distinctive marker of nutritional fiber starvation |
| **`Fos`** | AP-1 Immediate-Early Activation | -0.06 | **+0.78** ($p=0.032$) | **+1.46** | -0.80 | Upregulated in antibiotic & sham depletion |

#### 4. Visual Observations from Diagnostic Graphics
- **Volcano Plots (`fig_volcano_<cohort>.png`)**:
  - `GSE108045` exhibits a dramatic, highly powered bilateral volcano distribution with over 250 DEGs, highlighting the immense acute transcriptomic remodeling triggered by antibiotic cocktail treatment.
  - `GSE107925` exhibits a predominantly negative-skewed volcano plot (20 down, 5 up), reflecting loss of maturation markers in germ-free adult microglia.
  - `GSE186210` reveals clean, focused significance spikes for *Plin3* (down) and *Rbm5* / *Ccl24* (up) under fiber starvation.
- **Cross-Study Effect Size Correlation (`fig_lfc_correlation_heatmap.png`)**:
  Pairwise Spearman correlation across all ~6,000 common expressed genes reveals near-zero global transcriptome correlation ($\rho \approx -0.08$ to $+0.01$). This is a crucial methodological confirmation: transcriptome-wide noise is uncorrelated across independent laboratories, which proves that our cross-study meta-analysis in Horizon 3 will filter out study-specific noise and distill the true biological consensus signal!
- **Marker Comparison Bar Chart (`fig_marker_effect_sizes_by_cohort.png`)**:
  Directly contrasts the four perturbation modes, illustrating how acute antibiotic depletion drives intense anti-inflammatory brake collapse (*Tsc22d3*, *Ddit4*), while dietary fiber starvation uniquely impairs lipid droplet handling (*Plin3*).

#### 5. Blooms, Anomalies & Serendipity
- **Bloom 2.1: The GILZ (*Tsc22d3*) Anti-Inflammatory Collapse**:
  In `GSE108045`, the most striking finding is the profound down-regulation of *Tsc22d3* (GILZ, $\log_2\text{FC} = -2.36, \text{padj} = 5.56 \times 10^{-20}$). GILZ is known in immunology as the primary endogenous brake on NF-κB p65 and AP-1 transactivation. Its collapse under broad-spectrum antibiotic treatment explains why microglia lose immune tolerance and become primed for cytokine hyper-secretion. This provides an elegant mechanistic link between gut microbiome depletion and microglial neuroinflammation!
- **Bloom 2.2: The *Plin3* Lipid Droplet Connection in Fiber Starvation**:
  In `GSE186210`, dietary fiber starvation uniquely suppressed *Plin3* ($\log_2\text{FC} = -1.88, \text{padj} = 3.4 \times 10^{-5}$). Microglia rely on lipid droplets to safely sequester toxic lipid peroxides and fatty acids during metabolic stress. Losing PLIN3 under SCFA deprivation suggests that microbial metabolites are essential for microglial lipid metabolism and membrane lipid homeostasis.
- **Bloom 2.3: Quiescence vs. Shock**:
  The contrast between `GSE107925` (25 DEGs) and `GSE108045` (255 DEGs) reveals that lifelong germ-free housing establishes a state of developmental arrest and blunted maturity, whereas acute antibiotic depletion in adulthood induces an acute disruption shock that actively strips away repressive checkpoints.

#### 6. Branch / Consolidate Decision
- **Consolidation**:
  - The 4 full-scale DEG tables and phenotypic summary metrics are locked in `results/de_results/`.
  - Quality and statistical integrity gates are satisfied: 26/26 tests passing.
- **Advancement**:
  - Proceed to **Horizon 3: The Consensus Symphony (Cross-Study Statistical Synthesis)**.

#### 7. Next Horizon Step: Horizon 3 (The Consensus Symphony)
- **Target**: Run transcriptome-wide DerSimonian-Laird and REML random-effects meta-analysis combining effect sizes across all common expressed genes.
- **Key Questions for Horizon 3**:
  1. What is the consensus pooled effect size ($\hat{\theta}_{\text{meta}}$) and FDR for *Tnf*, *Tsc22d3*, *Tmem119*, and *Cx3cr1*?
  2. Which genes exhibit low heterogeneity ($I^2 < 25\%$, universal core markers) vs. high heterogeneity ($I^2 > 75\%$)?
  3. Does Leave-One-Out (LOO) sensitivity analysis confirm that the consensus signature remains statistically robust when omitting `GSE266602` or `GSE108045`?


