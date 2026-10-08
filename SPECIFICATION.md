# Technical and Mathematical Specification

## NeuroGut-MetaSeq: Cross-Study RNA-Seq Meta-Analysis Framework

---

## 1. Scientific Overview & Hypothesis

Microglia are specialized myeloid-lineage phagocytes resident in the brain parenchyma. In germ-free (GF) or antibiotic-depleted (ABX) rodents, microglia display immature transcriptomic signatures, malformed ramified morphologies, and blunted immune reactivity. Introduction of microbiota-derived Short-Chain Fatty Acids (SCFAs: acetate, propionate, butyrate) restores microglial maturity and modulates inflammatory cascades.

### Formal Biological Hypotheses
1. **$H_{depletion}$**: Depletion of the gut microbiome (GF or ABX) induces a conserved core microglial gene expression signature characterized by dysregulation of homeostatic checkpoint genes (*Tmem119*, *Cx3cr1*, *P2ry12*), antigen presentation, and altered basal cytokine signaling.
2. **$H_{rescue}$**: Short-chain fatty acid administration (specifically butyrate and propionate) counter-regulates dysbiosis signatures, downregulating NF-κB-dependent pro-inflammatory cytokines (*Tnf*, *Il1b*, *Ccl2*) and promoting homeostatic and phagocytic pathways.
3. **$H_{meta}$**: Multi-cohort statistical synthesis yields higher statistical power and lower false-discovery rates than any individual study, resolving cross-platform and laboratory-specific batch artifacts.

---

## 2. Statistical & Mathematical Models

### 2.1 Differential Expression Modeling (Negative Binomial GLM)

For each gene $i$ and sample $j$ in cohort $k$, raw read counts $Y_{ijk}$ follow a Negative Binomial distribution:

$$Y_{ijk} \sim \text{NB}(\mu_{ijk}, \alpha_{ik})$$

$$\mathbb{E}[Y_{ijk}] = \mu_{ijk} = s_{jk} \cdot q_{ijk}$$

$$\operatorname{Var}(Y_{ijk}) = \mu_{ijk} + \alpha_{ik} \mu_{ijk}^2$$

Where:
- $s_{jk}$ is the sample-specific size factor estimated via the median-of-ratios method:
  $$s_{jk} = \operatorname{median}_{i} \frac{Y_{ijk}}{\left( \prod_{r=1}^{m} Y_{irk} \right)^{1/m}}$$
- $q_{ijk}$ is proportional to true relative gene expression.
- $\alpha_{ik}$ is the gene-specific dispersion parameter.

The generalized linear model uses a log-link function:

$$\log_2(q_{ijk}) = \beta_{0ik} + \sum_{p} \beta_{pik} X_{jkp} + \beta_{\text{cond}, ik} \cdot \text{Condition}_{jk}$$

Hypothesis testing for gene $i$ evaluates $H_0: \beta_{\text{cond}, ik} = 0$ via the Wald test:

$$W_{ik} = \frac{\hat{\beta}_{\text{cond}, ik}}{\operatorname{SE}(\hat{\beta}_{\text{cond}, ik})} \sim \mathcal{N}(0, 1)$$

P-values are computed as $p_{ik} = 2 \left( 1 - \Phi(|W_{ik}|) \right)$.

---

### 2.2 Effect Size Meta-Analysis (REML & Hartung-Knapp-Sidik-Jonkman)

To pool effect sizes $\hat{\theta}_k = \hat{\beta}_{\text{cond}, ik}$ across $K$ independent cohorts with estimated variances $v_k = \operatorname{SE}(\hat{\beta}_{\text{cond}, ik})^2$:

Under the random-effects model:

$$\hat{\theta}_k = \theta + u_k + \epsilon_k, \quad u_k \sim \mathcal{N}(0, \tau^2), \quad \epsilon_k \sim \mathcal{N}(0, v_k)$$

#### 1. Restricted Maximum Likelihood (REML) Estimator for $\tau^2$:
When the number of studies is small ($K < 5$), the classical DerSimonian-Laird (DL) method systematically underestimates between-study variance $\tau^2$, yielding anti-conservative confidence intervals. NeuroGut-MetaSeq employs **Restricted Maximum Likelihood (REML)** by minimizing the negative restricted log-likelihood via bounded scalar optimization:

$$-2 l_R(\tau^2) = \sum_{k=1}^K \ln(v_k + \tau^2) + \ln \left( \sum_{k=1}^K \frac{1}{v_k + \tau^2} \right) + \sum_{k=1}^K \frac{(\hat{\theta}_k - \hat{\theta}_{\text{RE}}(\tau^2))^2}{v_k + \tau^2}$$

subject to $\tau^2 \ge 0$, where $\hat{\theta}_{\text{RE}}(\tau^2) = \frac{\sum_{k=1}^K w_k^*(\tau^2) \hat{\theta}_k}{\sum_{k=1}^K w_k^*(\tau^2)}$ and $w_k^*(\tau^2) = \frac{1}{v_k + \tau^2}$.

#### 2. Hartung-Knapp-Sidik-Jonkman (HKSJ) Adjustment:
To account for uncertainty in the variance component estimate $\hat{\tau}^2_{\text{REML}}$, we apply the Hartung-Knapp-Sidik-Jonkman (HKSJ) quadratic variance adjustment factor:

$$q_{\text{HKSJ}} = \frac{1}{K - 1} \sum_{k=1}^K w_k^* (\hat{\theta}_k - \hat{\theta}_{\text{meta}})^2$$

$$\operatorname{SE}_{\text{HKSJ}}(\hat{\theta}_{\text{meta}}) = \sqrt{\frac{\max(1, q_{\text{HKSJ}})}{\sum_{k=1}^K w_k^*}}$$

Critical values follow Student's $t$-distribution with $K - 1$ degrees of freedom. For $K = 4$ cohorts ($\text{df} = 3$), the two-tailed 95% critical threshold is:

$$t_{\text{crit}} = t_{0.975, 3} = 3.1824$$

$$95\% \text{ CI}_{\text{HKSJ}} = \left[ \hat{\theta}_{\text{meta}} - t_{\text{crit}} \cdot \operatorname{SE}_{\text{HKSJ}}(\hat{\theta}_{\text{meta}}), \; \hat{\theta}_{\text{meta}} + t_{\text{crit}} \cdot \operatorname{SE}_{\text{HKSJ}}(\hat{\theta}_{\text{meta}}) \right]$$

$$p_{\text{HKSJ}} = 2 \left( 1 - F_{t_3}\left( \frac{|\hat{\theta}_{\text{meta}}|}{\operatorname{SE}_{\text{HKSJ}}(\hat{\theta}_{\text{meta}})} \right) \right)$$

#### 3. DerSimonian-Laird Estimator (Supplementary Benchmark):
For backward comparability and sensitivity benchmarking, between-study variance is also estimated via method of moments:

$$Q = \sum_{k=1}^K w_k (\hat{\theta}_k - \bar{\theta}_{\text{FE}})^2, \quad \bar{\theta}_{\text{FE}} = \frac{\sum w_k \hat{\theta}_k}{\sum w_k}, \quad w_k = \frac{1}{v_k}$$

$$\tau^2_{\text{DL}} = \max\left( 0, \frac{Q - (K - 1)}{\sum_{k=1}^K w_k - \frac{\sum_{k=1}^K w_k^2}{\sum_{k=1}^K w_k}} \right)$$

$$I^2 = \max\left( 0, \frac{Q - (K - 1)}{Q} \right) \times 100\%$$

---

### 2.3 Non-Parametric P-Value Combination

#### 1. Fisher's Combined Probability Test:
$$\chi^2 = -2 \sum_{k=1}^K \ln(p_k) \sim \chi^2_{2K}$$

$$p_{\text{Fisher}} = 1 - F_{\chi^2_{2K}}(\chi^2)$$

#### 2. Stouffer's Weighted Z-Transform:
Taking one-sided Z-scores $Z_k = \Phi^{-1}(1 - p_k / 2) \cdot \operatorname{sign}(\hat{\theta}_k)$ weighted by sample size $n_k$:

$$Z_{\text{Stouffer}} = \frac{\sum_{k=1}^K \sqrt{n_k} Z_k}{\sqrt{\sum_{k=1}^K n_k}} \sim \mathcal{N}(0, 1)$$

---

### 2.4 Two-Tier Subgroup Decomposition & Multi-Study Factor Analysis

To address biological heterogeneity across distinct experimental models (lifelong germ-free absence, acute broad-spectrum antibiotic shock, dietary fiber starvation), we implement a Two-Tier Subgroup Decomposition.

#### 1. Between-Subgroup Heterogeneity Test:
For each gene $i$, we pool model-specific effect sizes $\hat{\theta}_g$ ($g \in \{\text{GF}, \text{ABX}, \text{Fiber}\}$) with inverse-variance weights $w_g = 1 / v_g$:

$$Q_{\text{between}} = \sum_{g=1}^G w_g (\hat{\theta}_g - \bar{\theta}_{\text{pooled}})^2 \sim \chi^2_{G - 1}$$

Where $\bar{\theta}_{\text{pooled}} = \frac{\sum w_g \hat{\theta}_g}{\sum w_g}$. The between-model heterogeneity p-value is $p_{Q} = 1 - F_{\chi^2_{G-1}}(Q_{\text{between}})$.

#### 2. Transcriptome Classification Partitioning:
- **Shared Microbial Core**: Directionally concordant across models, low between-study heterogeneity ($I^2 < 35\%$), invariant between perturbation classes ($p_Q > 0.05$). Exemplars: *Llgl2*, *Slfn2*, *Clu*.
- **ABX Mucosal Shock**: Disproportionate pharmacological effect in acute antibiotic cocktail ($|\hat{\theta}_{\text{ABX}}| \ge 1.0$) with extreme heterogeneity ($I^2 > 60\%$). Exemplars: *Tsc22d3* ($I^2 = 95.4\%$), *Ddit4* ($I^2 = 97.8\%$).
- **Fiber Dietary Starvation**: Metabolic substrate shock in zero-fiber diet ($|\hat{\theta}_{\text{Fiber}}| \ge 0.8, I^2 > 50\%$). Exemplar: *Plin3* ($I^2 = 96.5\%$).
- **Developmental Germ-Free**: Lifelong embryonic absence effect ($|\hat{\theta}_{\text{GF}}| \ge 0.8, I^2 > 50\%$).
- **Model-Divergent Mixed Axis**: Opposing directional regulation across perturbation models.

#### 3. Multi-Study Factor Analysis (SVD PCA):
We standardize log2 CPM expression $Z_{ij}$ within each cohort to remove baseline technical shifts and perform Singular Value Decomposition (SVD):

$$\mathbf{Z} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^T$$

Decomposing sample variance into Factor 1 (Microbial Tonic Surveillance Axis) and Factor 2 (Model Modality / Acute Stress Axis).

---

### 2.5 Multiple Hypothesis Correction

All raw p-values ($m$ tested genes) are corrected using the Benjamini-Hochberg False Discovery Rate (FDR) procedure:

$$q_{(i)} = \min_{j \ge i} \left( \frac{m \cdot p_{(j)}}{j} \right)$$

For genome-wide screening across 23,096 common genes under small study numbers ($K=4$), discovery testing is conducted on the REML Wald statistic $Z = \hat{\theta}_{\text{REML}} / \operatorname{SE}_{\text{REML}}$ (`fdr_random_effects`), with HKSJ providing conservative degrees-of-freedom adjusted standard errors, 95% confidence intervals, and nominal $p_{\text{hksj}} < 0.05$ checks.

---

### 2.6 Pathway Enrichment Statistics

#### Hypergeometric Over-Representation Analysis (ORA):
For a gene set of size $S$, genome size $N$, and DEG list of size $k$ containing $x$ intersecting genes:

$$P(X \ge x) = \sum_{j=x}^{\min(k, S)} \frac{\binom{S}{j} \binom{N - S}{k - j}}{\binom{N}{k}}$$

#### Fast Gene Set Enrichment Analysis (fgsea):
Uses the non-parametric running-sum Kolmogorov-Smirnov-like statistic across the fully ranked list of all genes sorted by signed test statistic $s_i = \operatorname{sign}(\hat{\theta}_i) \cdot (-\log_{10} p_i)$.

---

### 2.7 Upstream Transcription Factor Regulon Deconvolution

To determine whether the downstream targets of transcription factor $t$ are coordinately shifted in microbiome-depleted microglia, we project meta-analysis effect sizes onto curated transcriptional regulatory networks (TRRUST v2 mouse):

#### 1. Regulon Activity $Z$-Score:
For transcription factor $t$ with $n_t$ measured downstream target genes having sample mean effect size $\bar{\theta}_t$ and sample variance $s_t^2$, compared against background genes ($n_{\text{bg}}, \bar{\theta}_{\text{bg}}, s_{\text{bg}}^2$):

$$Z_{\text{activity}, t} = \frac{\bar{\theta}_t - \bar{\theta}_{\text{bg}}}{\sqrt{\frac{s_t^2}{n_t} + \frac{s_{\text{bg}}^2}{n_{\text{bg}}}}}$$

#### 2. Hypothesis Testing:
- **Welch's Two-Sample $t$-Test**: Accounts for unequal target vs background variance with Welch-Satterthwaite degrees of freedom $\nu$.
- **Mann-Whitney $U$ Test**: Non-parametric test for location shift without distributional assumptions.
- **Two-Sample Kolmogorov-Smirnov Test**: Tests whether target effect sizes follow the background distribution $F_{\text{bg}}(\theta)$.
- **Target Overlap Fisher's Exact Test**: Evaluates significant enrichment among consensus DEGs.
- Multiple testing correction applied via Benjamini-Hochberg FDR ($\text{FDR} \le 0.05$).

---

### 2.8 Weighted Gene Co-Expression Network Analysis (WGCNA)

To identify modular co-expression architectures and hub effectors across all 60 biological samples:

#### 1. Similarity & Soft-Thresholded Adjacency:
Between gene $i$ and gene $j$, Pearson correlation $s_{ij} = \operatorname{cor}(x_i, x_j)$ is converted to a signed co-expression adjacency using soft-threshold power $\beta = 6$ ($R^2 \ge 0.80$):

$$a_{ij} = \left( \frac{1 + s_{ij}}{2} \right)^\beta$$

#### 2. Topological Overlap Matrix (TOM):
Quantifies interconnectedness based on shared network neighbors:

$$\omega_{ij} = \operatorname{TOM}_{ij} = \frac{l_{ij} + a_{ij}}{\min(k_i, k_j) + 1 - a_{ij}}$$

#### 3. Module Eigengenes & Intramodular Connectivity:
For each detected module $q$, the Module Eigengene (ME) $E^{(q)}$ is defined as the first principal component of the standardized module expression matrix. Hub genes are identified by intramodular connectivity $k_{\text{in}}^{(i)}$.

---

### 2.9 Empirical In Vivo SCFA Metabolite Reversibility & Specificity Null Model

To quantitatively model candidate transcriptional reversibility by short-chain fatty acids (acetate, propionate, butyrate), vectors are grounded in empirical in vivo microglial RNA-seq from SCFA-supplemented germ-free mice (Erny et al. 2015 *Nature Neuroscience*, GSE64977, $N=6$):

#### 1. In Silico Rescue Index (ISRI):
$$\text{ISRI}_i = -\frac{\hat{\theta}_{\text{depletion}, i} \cdot \hat{\theta}_{\text{scfa}, i}}{|\hat{\theta}_{\text{depletion}, i}|}$$

- $\text{ISRI}_i > 0$: Successful reciprocal signature inversion (repressed genes restored upward; activated cytokines suppressed downward).
- $\text{ISRI}_i \le 0$: Refractory, priming-locked, or non-reversible.

#### 2. Percentage Rescue:
$$\operatorname{Rescue}\%_i = \max\left( 0\%, \; \min\left( 150\%, \; \frac{\text{ISRI}_i}{|\hat{\theta}_{\text{depletion}, i}|} \times 100\% \right) \right)$$

#### 3. Global Signature Inversion Metric:
$$r_{\text{inversion}} = \operatorname{cor}\left( \hat{\boldsymbol{\theta}}_{\text{depletion}}, \; \hat{\boldsymbol{\theta}}_{\text{scfa}} \right)$$

#### 4. 1,000-Permutation Genomic Specificity Null Model:
To test whether the observed rescue index is specific to microbiome-depleted signatures rather than a non-specific artifact of HDAC inhibition, we draw $B = 1,000$ independent random subsets of size $n = 19$ from non-differentially expressed background genes ($|\hat{\theta}| < 0.2, \text{FDR} > 0.5$):

$$p_{\text{perm}} = \frac{1 + \sum_{b=1}^B \mathbb{I}(\overline{\text{ISRI}}_{\text{null}}^{(b)} \ge \overline{\text{ISRI}}_{\text{obs}})}{1 + B}$$

---

### 2.10 Single-Cell Subpopulation Deconvolution & Lineage Normalization

To resolve the bulk RNA-seq bottleneck and evaluate whether interferon collapse reflects per-cell transcriptional shutoff or cellular depletion of the Interferon-Responsive Microglia (IRM) subset:

#### 1. Reference Single-Cell Signatures:
Validated markers from single-cell microglial atlases (Hammond et al. 2019 *Immunity*, Masuda et al. 2019 *Nature*):
- **Interferon-Responsive Microglia (IRM)**: *Oas1a*, *Stat1*, *Gbp2*, *Tap1*, *Ifit1*, *Ifit3*, *Irf7*, *Mx1*, *B2m*
- **Homeostatic Mature**: *Tmem119*, *P2ry12*, *Cx3cr1*, *Hexb*, *Csf1r*, *Sall1*, *Fcrls*
- **Phagocytic / DAM**: *Apoe*, *Ctsb*, *Ctsd*, *Trem2*, *Tyrobp*, *Lpl*
- **Cycling / Proliferating**: *Mki67*, *Top2a*, *Cdk1*, *Birc5*

#### 2. ISG-to-Lineage Normalization Index:
$$\text{Ratio}_{\text{ISG/Lineage}} = \frac{\frac{1}{|S_{\text{ISG}}|} \sum_{g \in S_{\text{ISG}}} \log_2(\text{CPM}_g + 1)}{\frac{1}{|S_{\text{Lineage}}|} \sum_{g \in S_{\text{Lineage}}} \log_2(\text{CPM}_g + 1)}$$

Where $S_{\text{Lineage}} = \{\textit{Hexb}, \textit{Csf1r}, \textit{Tmem119}\}$ (pan-microglial invariant lineage markers) and $S_{\text{ISG}} = \{\textit{Oas1a}, \textit{Stat1}, \textit{Gbp2}, \textit{Tap1}\}$. Lineage marker invariance ($p_{\text{lin}} > 0.05$) combined with ISG ratio reduction confirms cell-intrinsic transcriptional downregulation, consistent with stereological cell-density conservation (Erny 2015, Abdur-Rahman 2021).

---

## 3. Data Schemas

### 3.1 Metadata Schema (`data/metadata/<cohort>_metadata.csv`)
| Column | Type | Allowed Values / Description |
|---|---|---|
| `sample_id` | String | Unique sample identifier (e.g., `GSM2883011`) |
| `cohort` | String | GEO series accession (`GSE107925`, `GSE108045`, etc.) |
| `condition` | String | `reference` (SPF/CTR) or `perturbed` (GF/ABX/Zero_Fiber) |
| `group_label` | String | Biological label (`SPF`, `GF`, `ABX`, `Fiber`, `Zero_Fiber`) |
| `sex` | String | `Male`, `Female`, `Unspecified` |
| `tissue` | String | `Microglia` |
| `sequencing_type` | String | `Bulk RNA-seq` |

### 3.2 DEG Results Schema (`results/de_results/<cohort>_deg.csv`)
| Column | Type | Description |
|---|---|---|
| `gene_id` | String | Ensembl ID (`ENSMUSG00000028180`) |
| `gene_symbol` | String | Official MGI Symbol (`Tnf`, `Nfkb1`, `Cx3cr1`) |
| `baseMean` | Float | Mean normalized counts across all samples |
| `log2FoldChange` | Float | Log2 fold change (Perturbed vs Reference) |
| `lfcSE` | Float | Standard error of log2FoldChange |
| `stat` | Float | Wald test statistic |
| `pvalue` | Float | Raw two-tailed p-value |
| `padj` | Float | Benjamini-Hochberg adjusted p-value |

### 3.3 Meta-Analysis Results Schema (`results/meta_results/microglia_meta_analysis_summary.csv`)
| Column | Type | Description |
|---|---|---|
| `gene_symbol` | String | Official MGI Gene Symbol |
| `n_cohorts` | Integer | Number of cohorts where gene is detected |
| `cohorts_detected` | String | Semicolon-delimited list of detected cohorts |
| `meta_log2fc` | Float | REML pooled effect size (log2 fold change) |
| `meta_se` | Float | Restricted Maximum Likelihood standard error |
| `ci_lower` | Float | Lower bound of 95% Confidence Interval (HKSJ adjusted) |
| `ci_upper` | Float | Upper bound of 95% Confidence Interval (HKSJ adjusted) |
| `tau2` | Float | REML estimated between-study variance ($\tau^2$) |
| `p_random_effects`| Float | Random-effects p-value under REML |
| `p_hksj` | Float | Hartung-Knapp-Sidik-Jonkman adjusted p-value ($t_{k-1}$) |
| `meta_log2fc_dl` | Float | Benchmark DerSimonian-Laird pooled effect size |
| `meta_se_dl` | Float | DerSimonian-Laird standard error |
| `ci_lower_dl` | Float | Benchmark DL 95% CI lower bound |
| `ci_upper_dl` | Float | Benchmark DL 95% CI upper bound |
| `tau2_dl` | Float | Benchmark DL between-study variance |
| `p_random_effects_dl`| Float| Benchmark DL p-value |
| `cochran_q` | Float | Cochran's $Q$ heterogeneity test statistic |
| `i2_heterogeneity`| Float | Higgins $I^2$ percentage ($0 - 100\%$) |
| `heterogeneity_tier`| String | Heterogeneity stratum (`Low (<25%)`, `Moderate (25-50%)`, `Substantial (50-75%)`, `High (>=75%)`) |
| `fisher_stat` | Float | Fisher's chi-square combination statistic |
| `p_fisher` | Float | Fisher's combined p-value |
| `stouffer_z` | Float | Stouffer's combined $Z$-score |
| `p_stouffer` | Float | Stouffer's combined p-value |
| `direction_concordance`| String | Directional pattern (`Concordant Up`, `Concordant Down`, `Mixed`) |
| `fdr_random_effects`| Float | Benjamini-Hochberg FDR of REML random-effects p-value |
| `fdr_random_effects_dl`| Float| Benjamini-Hochberg FDR of DL random-effects p-value |
| `fdr_fisher` | Float | Benjamini-Hochberg FDR of Fisher p-value |
| `fdr_stouffer` | Float | Benjamini-Hochberg FDR of Stouffer p-value |
| `robustness_score`| Float | Leave-One-Out (LOO) robustness stability score ($0.0 - 1.0$) |
| `max_lfc_shift` | Float | Maximum absolute shift in effect size observed across LOO iterations |
| `loo_vulnerable_study`| String| Study whose omission drives the largest effect size deviation |
| `significance_flag`| Boolean| True if FDR < 0.05 and \|meta_log2fc\| >= 0.5 |

### 3.4 Upstream TF Regulon Schema (`results/pathways/tf_regulon_activity_summary.csv`)
| Column | Type | Description |
|---|---|---|
| `tf_symbol` | String | MGI Symbol of transcription factor (e.g., `Irf1`, `Fos`) |
| `target_count` | Integer | Number of measured downstream targets ($\ge 5$) |
| `mean_target_log2fc` | Float | Mean meta-analysis $\log_2\text{FC}$ of target genes |
| `median_target_log2fc` | Float | Median meta-analysis $\log_2\text{FC}$ of target genes |
| `activity_z_score` | Float | Standardized regulon activity $Z$-score |
| `p_welch` | Float | Welch's two-sample $t$-test p-value |
| `p_mann_whitney` | Float | Two-sided Mann-Whitney $U$ test p-value |
| `p_ks_test` | Float | Two-sample Kolmogorov-Smirnov test p-value |
| `p_fisher_overlap` | Float | Fisher's exact test p-value for DEG overlap |
| `fdr_welch` | Float | Benjamini-Hochberg adjusted FDR of Welch p-value |
| `regulon_status` | String | `Significantly Activated`, `Significantly Repressed`, or `Unchanged` |

### 3.5 WGCNA Module Assignments Schema (`results/networks/coexpression_module_assignments.csv`)
| Column | Type | Description |
|---|---|---|
| `gene_symbol` | String | Official MGI Gene Symbol |
| `cluster_id` | Integer | Hierarchical cluster identification index |
| `k_in` | Float | Intramodular connectivity within assigned module |
| `module_name` | String | Assigned module name (e.g., `M_Quiescence`) |

### 3.6 SCFA Metabolite Rescue Schema (`results/pathways/scfa_metabolite_rescue_modeling.csv`)
| Column | Type | Description |
|---|---|---|
| `gene_symbol` | String | Official MGI Gene Symbol |
| `depletion_meta_log2fc` | Float | Pooled meta-analysis effect size under microbiome depletion |
| `depletion_se` | Float | Standard error of depletion effect size |
| `depletion_fdr` | Float | Random-effects or Fisher FDR under depletion |
| `i2_heterogeneity` | Float | Higgins $I^2$ across depletion cohorts |
| `scfa_rescue_log2fc` | Float | Effect size under SCFA metabolite supplementation |
| `net_post_rescue_log2fc` | Float | Residual difference $\hat{\theta}_{\text{depletion}} + \hat{\theta}_{\text{rescue}}$ |
| `in_silico_rescue_index` | Float | Direction-adjusted In Silico Rescue Index (ISRI) |
| `rescue_percentage` | Float | Clamped rescue percentage ($0 - 150\%$, where $100\%$ denotes exact 1:1 neutralization) |
| `rescue_status` | String | Classification (`Candidate Metabolite-Reversible`, `Partial Responder`, or `Priming-Locked / Refractory`) |
| `proposed_mechanism` | String | Biochemical mechanism (e.g., HDAC inhibition, FFAR2 signaling) |

### 3.7 Perturbation Subgroup Decomposition Schema (`results/meta_results/perturbation_subgroup_decomposition.csv`)
| Column | Type | Description |
|---|---|---|
| `gene_symbol` | String | Official MGI Gene Symbol |
| `subgroup_axis` | String | Stratification tier (`Shared Microbial Core`, `Developmental Absence Specific`, `Acute Antibiotic Shock Specific`, `Model-Heterogeneous`) |
| `classification_rationale` | String | Descriptive categorization rationale including $I^2$ and $Q_{\text{between}}$ values |
| `meta_log2fc` | Float | Pooled REML effect size across all cohorts |
| `meta_se` | Float | Standard error of pooled REML effect size |
| `p_random_effects` | Float | Random-effects p-value under REML |
| `fdr_random_effects` | Float | Benjamini-Hochberg FDR of REML p-value |
| `i2_heterogeneity` | Float | Overall Higgins $I^2$ across all cohorts |
| `q_between_models` | Float | Cochran's $Q_{\text{between}}$ testing heterogeneity between perturbation models |
| `p_q_between` | Float | Chi-square p-value for model discordance |
| `lfc_germ_free` | Float | Log2 fold change in Germ-Free model (GSE107925) |
| `lfc_antibiotics` | Float | Log2 fold change in Antibiotic cocktail model (GSE108045) |
| `lfc_fiber_starvation` | Float | Log2 fold change in Fiber-starvation model (GSE186210) |
| `lfc_sham_percoll` | Float | Log2 fold change in Percoll isolation model (GSE266602) |

### 3.8 Microglia Subpopulation Deconvolution Schema (`results/pathways/microglia_subpopulation_deconvolution.csv`)
| Column | Type | Description |
|---|---|---|
| `sample_id` | String | Unique sample identifier (matching metadata) |
| `cohort` | String | GEO series accession (`GSE107925`, `GSE108045`, etc.) |
| `condition` | String | Treatment condition (`reference` or `perturbed`) |
| `sig_Interferon-Responsive (IRM)` | Float | Mean CPM signature score for IRM subset (*Oas1a*, *Stat1*, *Gbp2*, *Tap1*, etc.) |
| `sig_Homeostatic Mature` | Float | Mean CPM signature score for Homeostatic Mature microglia (*Tmem119*, *P2ry12*, *Cx3cr1*, *Hexb*) |
| `sig_Phagocytic / DAM` | Float | Mean CPM signature score for DAM / phagocytic microglia (*Apoe*, *Trem2*, *Ctsd*, *Tyrobp*) |
| `sig_Cycling / Proliferating` | Float | Mean CPM signature score for cycling microglia (*Mki67*, *Top2a*, *Cdk1*) |
| `isg_raw_score` | Float | Log2 CPM score of core ISG marker panel |
| `lineage_pan_score` | Float | Log2 CPM score of invariant pan-microglial lineage markers (*Hexb*, *Csf1r*, *Tmem119*) |
| `isg_to_lineage_ratio` | Float | Normalized $\text{Ratio}_{\text{ISG/Lineage}}$ evaluating cell-intrinsic vs compositional changes |

---

## 4. FAIR Principles & Reproducibility Guarantees
- **Findability**: All accession IDs and publications tracked in `config/datasets.yaml` and `CITATION.cff`.
- **Accessibility**: Data downloads automated via NCBI E-utilities / FTP with automatic fallback to bundled demo matrices.
- **Interoperability**: Standardized Ensembl-to-MGI identifier mappings with tidy CSV tables.
- **Reusability**: Dockerfile containerization and GitHub Actions workflow for zero-dependency execution.

