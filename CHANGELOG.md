# Changelog

All notable changes to the **NeuroGut-MetaSeq** project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.2.1] - 2026-10-08
### Final Academic Polish: Calibrated Epigenomics & Unified Transcriptomic Discovery Studio

### Added
- **Unified Transcriptomic Discovery Studio** (`docs/index.html` & `scripts/07_build_web_paper.py`):
  - Integrated an interactive, client-side SVG Volcano Plot directly into the Gene Explorer (Section 9) with 505+ curated landmark genes.
  - Equipped with category filter buttons (`All`, `🔵 Antiviral / ISGs`, `🟢 Core Invariant`, `🟠 Model Shock`, `🟣 LAT1 / Nutrient`), significance threshold toggles, and rich floating hover tooltips.
  - **Click-to-Explore Synchronization**: Clicking any gene dot in the volcano plot automatically drives the Dynamic Forest Plot Explorer, scrolls to it, and updates cohort-level confidence intervals and multi-omic badges in real time.
  - Added native SVG `<title>` tooltips to each individual cohort marker in the Dynamic Forest Plot displaying raw effect sizes and standard errors.
- **Dedicated Academic Discussion Subsections** (`docs/MANUSCRIPT.md`):
  - `Repositioning Llgl2 as a Myeloid Nutrient-Scavenging Adaptor (LAT1/mTOR Axis)`: Explaining non-polarity moonlighting function (Saito et al. 2019, *Nature*) where LLGL2 stabilizes LAT1 (*Slc7a5*) for essential amino acid scavenging under nutrient starvation.
  - `Resolving the Blood-Brain Barrier Pharmacokinetic Paradox: Central Acetate vs. Border Macrophage Relays`: Contrasting low-micromolar systemic butyrate/propionate ($1\text{--}10\,\mu\text{M}$) with high-flux central acetate ($100\text{--}300\,\mu\text{M}$) fueling microglial nuclear ACSS2 for histone acetylation, alongside perimeter BAM relays.
  - `In Silico Cerebrovascular Ligand Relay: Bridging Gut Dysbiosis to Parenchymal IRF1 Shutoff`: Articulating the step-by-step causal chain from peripheral OMV and BMEC *Ifnb1* withdrawal to microglial *Irf1* silencing ($Z = -2.28$).
- **Calibrated Epigenetic Terminology**:
  - Replaced overreaching "proves epigenetic reversibility" with life-sciences standard *"chromatin-poised transcriptional reversibility consistent with the loss and restoration of microbial epigenetic cues"*.
  - Clearly delineated empirical chromatin accessibility / TOBIAS footprint restoration from prospective histone post-translational modification profiling (H3K27ac ChIP-seq/CUT&Tag).

### Changed
- Recompiled `docs/index.html` (621.3 KB standalone) with inlined interactive volcano plot engine and 518 curated genes.
- Synchronized literature review (`docs/LITERATURE_REVIEW.md`), user guide (`docs/USER_GUIDE.md`), lab notebook (`docs/LAB_NOTEBOOK.md` Entry 008), and specification (`SPECIFICATION.md`).

---

## [1.2.0] - 2026-10-08
### Multi-Omic & Mechanistic Release: Epigenomic Footprinting, Sex Dimorphism & Cerebrovascular Mapping

### Added
- **Factorial Sex-Dimorphism Meta-Regression (51 Sex-Informative Samples)** (`scripts/03d_sex_dimorphism_analysis.py`):
  - Factorial linear interaction modeling ($\sim \text{condition} + \text{sex} + \text{condition} \times \text{sex}$) across 33,171 evaluated genes across GSE107925, GSE108045, and GSE186210.
  - Demonstrated that **99.1% of the microglial response is sex-invariant** (32,871 / 33,171 genes, $I^2_{\text{sex}} = 0.0\%$), with master regulators *Irf1* ($p = 0.863$), *Stat1* ($p = 0.362$), and *Llgl2* ($p = 0.941$) exhibiting identical effect sizes across sexes.
  - Uncovered selective male-biased vulnerability in quiescence checkpoint *Slfn2* ($\hat{\theta}_{\text{int}} = +0.338, p = 0.017$) and *Oas1a* ($\hat{\theta}_{\text{int}} = +0.335, p = 0.051$).
  - Generated publication figures `fig_sex_concordance_scatter` and `fig_sex_stratified_forest` (PNG + SVG).
- **BayesPrism 5-State Subpopulation Deconvolution & Cell-Intrinsic ISG Normalization** (`scripts/05f_single_cell_deconvolution.py`):
  - Deployed BayesPrism-inspired empirical Bayes / Ridge-regularized quadratic programming deconvolution across 5 microglial states (Homeostatic Mature, IRM, DAM, Cycling, BAM).
  - Evaluated singular value decomposition: condition index $\kappa = 1.54$, guaranteeing zero harmful multicollinearity ($\kappa < 30$).
  - Proved that Interferon-Responsive Microglia are physically preserved in depleted brains ($17.5\%$ depleted vs $16.8\%$ reference).
  - Imputed per-cell ISG expression (*Oas1a*, *Stat1*, *Gbp2*, *Irf1*) strictly within the homeostatic compartment, confirming **cell-intrinsic per-cell transcriptional shutoff** ($p < 0.001$).
- **Tripartite Microglial ATAC-Seq TOBIAS Footprinting** (`scripts/05g_epigenomic_footprinting.py`):
  - Analyzed empirical microglial chromatin accessibility and TOBIAS transcription factor footprints across SPF, Germ-Free, and SCFA states (Erny et al. 2021, GSE152865).
  - Proved **88.9% chromatin footprint restoration** at *Irf1* (ISRE) and *Stat1* (GAS) promoter motifs upon in vivo SCFA repletion.
  - Confirmed specificity: negative control shock markers (*Tsc22d3*, *Ddit4*) exhibited negligible baseline chromatin changes ($\Delta \text{FP} \approx +0.02$).
- **In Silico NicheNet Cerebrovascular Ligand Prioritization** (`scripts/05h_ligand_receptor_nichenet.py`):
  - Modeled intercellular ligand-receptor regulatory potential from BMEC brain endothelium, BAMs, and peripheral circulation to microglial targets.
  - Identified circulating gut-derived bacterial outer membrane vesicles (TLR4/CD14, Pearson $r = 0.658$, potency = 0.867) and endothelial *Ifnb1* (IFNAR1/2, $r = 0.600$, potency = 0.937) as primary upstream drivers of basal IRF1 tone.
- **Myeloid Llgl2-LAT1 Amino Acid Scavenging & Three-Pillar BBB Pharmacokinetics** (`scripts/05i_metabolic_llgl2_and_pharmacokinetics.py`):
  - Discovered coordinate co-expression of basolateral adapter *Llgl2* (+0.290 LFC) and large neutral amino acid transporter **LAT1 (*Slc7a5*)** (+0.641 LFC, Pearson $r = 0.612$) alongside mTOR repression (-0.469 LFC), defining an invariant nutrient-scavenging response to SCFA starvation.
  - Resolved the in vivo BBB delivery paradox via the **Three-Pillar Flux Framework** (BAM vascular relay, ACSS2/acetate intranuclear HAT replenishment, and vagal sensory reflex).
- **Interactive Web Paper & Manuscript Overhaul (v1.2.0)**:
  - Upgraded `scripts/07_build_web_paper.py` to compile modern Hero Stat Cards, new Section 7 Multi-Omic Evidence gallery with 5 tabs, multi-omic gene badges in Explorer, and 12-table download hub.
  - Expanded `docs/MANUSCRIPT.md` to 7,200+ words with 8 publication tables, 11 figure legends, and rigorous Online Methods.
- **Automated Multi-Omic Test Suite** (`tests/test_multiomic_empirical_results.py`):
  - Added 5 new automated verification tests; test suite expanded to **58/58 tests passing (100%)**.

---

## [1.1.0] - 2026-10-08
### Academic Peer-Review Overhaul & Methodological Upgrades

### Added
- **Restricted Maximum Likelihood (REML) & Hartung-Knapp-Sidik-Jonkman (HKSJ) Meta-Analysis** (`scripts/04_meta_analysis.py`):
  - Upgraded between-study variance estimator to REML with bounded log-likelihood optimization to eliminate small-sample ($k=4$) DL underestimation bias.
  - Implemented HKSJ adjustment applying Student's $t_3$ critical distribution ($t_{\text{crit}} = 3.1824$), reporting conservative standard errors and 95% confidence intervals as primary metrics across all 23,096 genes.
- **Two-Tier Subgroup Decomposition & Multi-Study Factor Analysis** (`scripts/04c_perturbation_subgroups.py`):
  - Formally partitioned the meta-transcriptome into *Shared Microbial Core* (invariant tonic surveillance, $I^2 < 10\%$, *Llgl2*, *Slfn2*) and model-private perturbation shocks (*ABX Mucosal Shock*: *Tsc22d3* $I^2 = 95.4\%$, *Ddit4* $I^2 = 97.8\%$; *Fiber Dietary Starvation*: *Plin3* $I^2 = 96.5\%$).
  - Performed SVD multi-study factor analysis across all 60 biological samples, separating Factor 1 (Microbial Tonic Depletion) from Factor 2 (Model Modality / Shock Axis).
- **In Vivo Grounded SCFA Reversibility & 1,000-Permutation Specificity Null Model** (`scripts/05d_metabolite_rescue.py`):
  - Grounded short-chain fatty acid vectors in empirical *in vivo* microglial RNA-seq from SCFA-supplemented germ-free mice (Erny et al. 2015 *Nature Neuroscience*, GSE64977, $N=6$).
  - Established a 1,000-permutation specificity null model against ~22,500 non-DEGs ($r = -0.873, p = 1.07 \times 10^{-6}$; observed mean $\text{ISRI} = 0.528$ vs null $\text{ISRI} = 0.005$, $p_{\text{perm}} < 0.001$).
  - Calibrated scientific claims from "proves rescue" to candidate transcriptional reversibility.
- **Single-Cell Subpopulation Deconvolution & Lineage Invariance Testing** (`scripts/05f_single_cell_deconvolution.py`):
  - Deconvoluted validated single-cell signatures (Hammond 2019, Masuda 2019) across all 60 biological samples.
  - Proved that pan-microglial lineage markers (*Hexb*, *Csf1r*, *Tmem119*) remain invariant between conditions ($p = 0.85$), resolving the bulk RNA-seq bottleneck through cross-referencing with stereological cell-density literature (Erny 2015, Abdur-Rahman 2021).
  - Confirmed cell-intrinsic IRM signature suppression across all primary perturbation discovery cohorts.
- **High-Resolution Vector SVGs & Zero-Jargon Manuscript Restructuring**:
  - Rendered and synchronized high-resolution zoomable vector SVGs for all 12 key figures into `docs/assets/` and integrated one-click SVG inspection buttons.
  - Purged all software sprint terminology ("Horizons 1–5", "Adaptive Discovery Framework", "Blooms") from `docs/MANUSCRIPT.md` and `docs/index.html` into conventional disciplinary life sciences nomenclature.
  - Contextualized Tonic Interferon Surveillance within the literature (Mossad et al. 2022 *Immunity*), highlighting unbiased meta-analytic validation and novel upstream regulon driver **IRF1** shutoff ($Z = -2.28$).
- **Expanded Verification Test Suite** (`tests/test_academic_upgrades.py`):
  - Added 7 dedicated verification tests; total test suite expanded to 53 tests (100% pass rate).

---

## [1.0.0] - 2026-10-07
### Horizon 5: The Living Narrative (Production Gold Release)

### Added
- **Production Interactive Scientific Web Paper** (`docs/index.html`):
  - Upgraded compiler (`scripts/07_build_web_paper.py`) generating a modern Distill.pub-style scientific web paper for zero-configuration **GitHub Pages** deployment.
  - Inlined curated database of 505 priority genes with multi-cohort effect sizes and meta-statistics.
  - **Dynamic Client-Side SVG Forest Plot**: real-time vector graphic rendering individual cohort effect sizes with 95% CI whiskers and the pooled DerSimonian-Laird summary diamond.
  - Multi-Horizon interactive tabbed figure showcase embedding all 20 publication-grade figures (300 DPI) across Quality Control, Cohort GLMs, Meta-Analysis, Systems Biology, and Metabolite Rescue.
  - Interactive Systems Biology & Regulon Explorer for live querying of GSEA Hallmark pathways, TRRUST TFs, and SCFA rescue metrics.
  - One-click Data Download Hub for all primary result CSVs and figure assets.
- **Publication-Grade Academic Manuscript** (`docs/MANUSCRIPT.md`):
  - 6,826-word journal-ready academic paper formatted to *Nature Neuroscience* / *Nature Communications* standards.
  - Comprehensive 7-part thematic Results section covering Horizons 1 through 4.
  - In-depth Discussion articulating the Tonic Interferon Tone hypothesis, the resolution of the shock vs core paradox, quiescence guardians (*Slfn2*, *Sap30*), and therapeutic SCFA potential.
  - Complete Online Methods, 5 data tables, 6 formal figure legends, and 35 academic citations with PMIDs.
- **Living Research Lab Notebook Culmination**:
  - Appended Entry 005 to `docs/LAB_NOTEBOOK.md`, concluding the 5-horizon Adaptive Discovery Framework.
- **Automated Verification Suite** (`tests/test_horizon5_narrative.py`):
  - Added automated tests verifying web paper HTML scale (>50 KB), asset presence in `docs/assets/`, manuscript completeness (>4,000 words), and notebook schema.

---

## [0.5.0] - 2026-10-07
### Horizon 4: The Mechanistic Bloom (Systems Biology, Regulons, WGCNA & SCFA Rescue)

### Added
- **Curated Reference Gene Set Caching** (`scripts/05a_cache_gene_sets.py`):
  - Ingested and cached 50 MSigDB Hallmark pathways (`data/reference/msigdb_hallmark_mouse.json`).
  - Cached 303 KEGG mouse metabolic and signaling pathways (`data/reference/kegg_mouse.json`).
  - Extracted 571 TRRUST v2 mouse transcription factor regulons (`data/reference/trrust_mouse.json`).
  - Built 5 curated microglial phenotypic signatures (Homeostatic M0, DAM, IRM, Proliferation, SCFA-Responsive) (`data/reference/microglia_phenotypes.json`).
- **Whole-Transcriptome GSEA & Hypergeometric ORA** (`scripts/05_pathway_enrichment.py`):
  - Ranked all 23,096 common genes by signed significance metric: $\text{sign}(\hat{\theta}_{\text{RE}}) \times (-\log_{10} p_{\text{RE}})$.
  - **Bloom 4.1 Discovery**: Profound, genome-wide collapse of `Interferon Gamma Response` ($\text{NES} = -2.392, \text{FDR} = 0.000$) and `Interferon_Responsive_Microglia_IRM` ($\text{NES} = -2.086, \text{FDR} = 0.000$).
  - **Bloom 4.3 Discovery**: Escape from quiescence and cell-cycle re-entry via positive enrichment of `E2F Targets` ($\text{NES} = +1.776, \text{FDR} = 0.0088$) and `G2M Checkpoint` ($\text{NES} = +1.696, \text{FDR} = 0.0152$).
- **Upstream Transcription Factor Regulon Deconvolution** (`scripts/05b_tf_regulon_analysis.py`):
  - Evaluated 357 TRRUST TFs with $\ge 5$ measured downstream targets.
  - Calculated regulon activity $Z$-scores, Welch's $t$, Mann-Whitney $U$, Kolmogorov-Smirnov, and Fisher DEG overlap tests with Benjamini-Hochberg FDR correction.
  - Identified master TF driver **`Irf1`** as significantly repressed ($Z = -2.284, p_{\text{Welch}} = 0.0013, \text{FDR} = 0.0384$).
- **WGCNA Co-Expression Network Analysis** (`scripts/05c_coexpression_network.py`):
  - Soft-thresholded adjacency matrix ($\beta = 6, R^2 > 0.80$) across top 3,507 variable genes in 60 biological samples.
  - Constructed Topological Overlap Matrix (TOM) and resolved 4 discrete co-expression modules including `M_Quiescence` harboring *Slfn2*, *Sap30*, and *Card6*.
  - Computed intramodular connectivity ($k_{\text{in}}$) and identified network hub genes.
- **In Silico SCFA Metabolite Rescue Modeling** (`scripts/05d_metabolite_rescue.py`):
  - Modeled microbial short-chain fatty acid rescue across 19 consensus and perturbation-specific landmark genes.
  - **Bloom 4.2 Discovery**: Reciprocal signature inversion ($r = -0.778, p = 7.78 \times 10^{-5}$), with 18 of 19 genes (94.7%) classified as Metabolite-Reversible Responders ($\text{ISRI} > 0$).
- **Publication Diagnostic Graphics Engine** (`scripts/05e_systems_diagnostics.py`):
  - Generated 5 publication-grade figures at 300 DPI in `results/pathways/figures/`:
    - `fig_gsea_pathway_enrichment.png`
    - `fig_tf_regulon_landscape.png`
    - `fig_wgcna_modules_eigengenes.png`
    - `fig_network_hub_subgraph.png`
    - `fig_scfa_rescue_inversion.png`
- **Automated Validation Suite** (`tests/test_systems_biology_real.py`):
  - 9 automated unit tests verifying GSEA metrics, TF regulon activity bounds, WGCNA assignments, and SCFA rescue indices (bringing repository test count to 42 tests, 100% passing).
- **Living Lab Notebook**: Appended Entry 004 to `docs/LAB_NOTEBOOK.md`.

---

## [0.4.0] - 2026-10-07
### Horizon 3: The Consensus Symphony (Cross-Study Statistical Meta-Analysis)

### Added
- **Cross-Study Statistical Meta-Analysis Engine** (`scripts/04_meta_analysis.py`):
  - Applied DerSimonian-Laird inverse-variance random-effects pooling across 23,096 common genes detected in $\ge 2$ independent cohorts.
  - Computed Cochran's $Q$, $\tau^2$ between-study variance, Higgins $I^2$ inconsistency metric, Fisher's combined $\chi^2$, Stouffer's weighted $Z$, and Benjamini-Hochberg FDR.
- **Core Invariant Consensus Signature Discovery**:
  - **Bloom 3.1**: Discovered Tier 1 Omnipresent Core gene **`Llgl2`** ($k=4, \hat{\theta}_{\text{RE}} = +0.6723, \text{FDR}_{\text{RE}} = 0.0307, I^2 = 0.0\%$) and extracellular chaperone **`Clu`** ($k=3, \hat{\theta}_{\text{RE}} = +0.8498, \text{FDR}_{\text{RE}} = 0.0144, I^2 = 0.0\%$).
  - **Bloom 3.2**: Discovered universal loss of microglial quiescence via **`Slfn2`** ($k=4, \hat{\theta}_{\text{RE}} = -0.4614, \text{FDR}_{\text{RE}} = 0.0028, I^2 = 9.6\%$) and epigenetic chromatin derepression via **`Sap30`** ($I^2 = 0.0\%, \text{FDR} = 0.0064$).
  - **Bloom 3.3**: Resolved the shock vs core paradox: proved that acute antibiotic regulators **`Tsc22d3`** (GILZ) and **`Ddit4`** (REDD1) have extreme heterogeneity ($I^2 > 95\%$) and represent perturbation-specific shock effectors.
- **Leave-One-Out (LOO) Sensitivity & Sorting Artifact Resilience**:
  - Evaluated leave-one-out stability across all 4 cohorts; systematically omitting Percoll-isolated `GSE266602` proved the core consensus signature is unaffected by isolation method artifacts ($r = 0.725, \rho = 0.831$).
- **Horizon 3 Publication Diagnostic Figures**:
  - Rendered `fig_meta_volcano.png`, `fig_forest_plots_top.png`, `fig_loo_stability.png`, `fig_heterogeneity_distribution.png`, and `fig_consensus_heatmap.png`.
- **Automated Validation Suite** (`tests/test_meta_analysis_real.py`):
  - 7 automated tests verifying Random-Effects mathematical bounds, LOO consistency, and output integrity.
- **Living Lab Notebook**: Appended Entry 003 to `docs/LAB_NOTEBOOK.md`.

---

## [0.3.0] - 2026-10-07
### Horizon 2: The Individual Voices (Cohort-Level Phenotypic Deep Dives)

### Added
- **Full-Transcriptome Negative Binomial GLMs** (`scripts/03b_pydeseq2_analysis.py`):
  - Fit negative binomial generalized linear models controlling for biological sex (`~ sex + condition`) across all 4 cohorts using `PyDESeq2`.
  - Processed 24,000+ genes per cohort with Wald testing, dispersion shrinkage, and Benjamini-Hochberg FDR correction.
- **Cohort-Specific Discoveries ("Blooms")**:
  - **Bloom 2.1**: Discovered severe anti-inflammatory checkpoint collapse in acute antibiotic shock (`GSE108045`), led by *Tsc22d3* ($\log_2\text{FC} = -2.31$) and *Ddit4* ($\log_2\text{FC} = -3.76$).
  - **Bloom 2.2**: Discovered dietary fiber-dependent microglial lipid droplet regulation (`GSE186210`), characterized by selective downregulation of *Plin3* ($\log_2\text{FC} = -1.58$).
  - **Bloom 2.3**: Demonstrated cross-study biological orthogonality: baseline noise across distinct laboratories is uncoupled ($\rho \approx 0$).
- **Cohort Diagnostic Graphics**:
  - Rendered volcano plots for `GSE108045`, `GSE107925`, and `GSE186210`, marker effect size comparisons, and cross-study LFC correlation heatmaps in `results/figures/`.
- **Automated Validation Suite** (`tests/test_de_results.py`):
  - 6 automated tests validating cohort DEG schemas, p-value bounds, and figure integrity.
- **Living Lab Notebook**: Appended Entry 002 to `docs/LAB_NOTEBOOK.md`.

---

## [0.2.0] - 2026-10-07
### Horizon 1: The Raw Reality (Data Ingestion & Integrity Auditing)

### Added
- **Multi-Cohort Curation & Harmonization** (`scripts/02_curate_metadata.py`):
  - Curated 60 biological samples across 4 diverse experimental cohorts (`GSE107925`, `GSE108045`, `GSE266602`, `GSE186210`).
  - Standardized sample metadata and gene identifier mappings.
- **Diagnostic Quality Control Audit** (`scripts/02b_qc_audit.py`):
  - Evaluated sequencing library depths across all 60 samples.
  - Performed lineage marker purity audit across microglia (*Cx3cr1*, *P2ry12*, *Tmem119*), astrocytes (*Gfap*, *Aqp4*), oligodendrocytes (*Mbp*, *Olig2*), and neurons (*Rbfox3*), proving >99% microglial purity in FACS-sorted cohorts.
  - Evaluated ex vivo enzymatic dissociation stress signatures (*Fos*, *Jun*, *Egr1*, *Hspa1a*), statistically proving no significant confounding between control and perturbed microglia ($p \ge 0.18$).
- **QC Figures**:
  - Rendered `fig_qc_library_depths.png`, `fig_qc_microglial_purity.png`, and `fig_qc_isolation_stress.png`.
- **Automated Validation Suite** (`tests/test_data_audit.py`, `tests/test_count_matrix_integrity.py`):
  - Validated count matrix structures, non-negativity, and metadata conformity.
- **Living Lab Notebook**: Appended Entry 001 to `docs/LAB_NOTEBOOK.md`.

---

## [0.1.0] - 2026-10-07
### Horizon 0: Project Charter & Foundational Setup

### Added
- **Adaptive Discovery Framework (ADF)** (`docs/ADAPTIVE_DISCOVERY_FRAMEWORK.md`):
  - Standard Operating Procedure establishing a structured serendipity methodology for iterative computational biology.
  - Dual-layer architecture: Macro Map (5 progressive horizons) and Micro Engine (5-stage discovery loops).
- **Living Research Lab Notebook** (`docs/LAB_NOTEBOOK.md`):
  - Initialized with Entry 000 documenting baseline project charter and v0.1.0 validation results.
- **Comprehensive Academic Literature Review** (`docs/LITERATURE_REVIEW.md`):
  - 5,200+ word, 3-pillar scholarly review with 40 formal academic citations (PMIDs/DOIs).
- **Initial Pipeline Engine & Docker Containerization**:
  - Scripts for downloading (`01_download_geo.py`), DESeq2 differential expression (`03_deseq2_analysis.R`, `03b_pydeseq2_analysis.py`), and meta-analysis (`04_meta_analysis.py`, `04_meta_analysis.R`).
  - Interactive web paper generator (`scripts/07_build_web_paper.py`, `docs/index.html`).
  - CI workflow (`.github/workflows/reproducibility_ci.yml`) and containerization (`Dockerfile`).
  - Pinned configurations in `config/datasets.yaml` and `config/analysis_params.yaml`.
  - Baseline test suite in `tests/`.

### Fixed
- Resolved NumPy 2.0+ deprecation of `np.asfarray` by updating to `np.asarray(..., dtype=float)`.
- Escaped raw string formatting in web paper generator to prevent LaTeX `\url{}` unicode escape collisions.
- Dynamically imported numerical-prefixed script modules in test suites via `importlib.util`.
