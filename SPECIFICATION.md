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

### 2.2 Effect Size Meta-Analysis (Inverse-Variance Random-Effects)

To pool effect sizes $\hat{\theta}_k = \hat{\beta}_{\text{cond}, ik}$ across $K$ independent cohorts with estimated variances $v_k = \operatorname{SE}(\hat{\beta}_{\text{cond}, ik})^2$:

Under the random-effects assumption:

$$\hat{\theta}_k = \theta + u_k + \epsilon_k, \quad u_k \sim \mathcal{N}(0, \tau^2), \quad \epsilon_k \sim \mathcal{N}(0, v_k)$$

#### DerSimonian-Laird Estimator for Between-Study Variance ($\tau^2$):
First compute the fixed-effects weights $w_k = \frac{1}{v_k}$ and Cochran's heterogeneity statistic $Q$:

$$Q = \sum_{k=1}^K w_k (\hat{\theta}_k - \bar{\theta}_{\text{FE}})^2, \quad \bar{\theta}_{\text{FE}} = \frac{\sum_{k=1}^K w_k \hat{\theta}_k}{\sum_{k=1}^K w_k}$$

The between-study variance $\tau^2$ is estimated as:

$$\tau^2 = \max\left( 0, \frac{Q - (K - 1)}{\sum_{k=1}^K w_k - \frac{\sum_{k=1}^K w_k^2}{\sum_{k=1}^K w_k}} \right)$$

Higgins & Thompson's $I^2$ inconsistency metric:

$$I^2 = \max\left( 0, \frac{Q - (K - 1)}{Q} \right) \times 100\%$$

#### Random-Effects Weights and Pooled Effect:
$$w_k^* = \frac{1}{v_k + \tau^2}$$

$$\hat{\theta}_{\text{meta}} = \frac{\sum_{k=1}^K w_k^* \hat{\theta}_k}{\sum_{k=1}^K w_k^*}, \quad \operatorname{SE}(\hat{\theta}_{\text{meta}}) = \sqrt{\frac{1}{\sum_{k=1}^K w_k^*}}$$

$$95\% \text{ CI} = \left[ \hat{\theta}_{\text{meta}} - 1.96 \cdot \operatorname{SE}(\hat{\theta}_{\text{meta}}), \; \hat{\theta}_{\text{meta}} + 1.96 \cdot \operatorname{SE}(\hat{\theta}_{\text{meta}}) \right]$$

$$Z_{\text{meta}} = \frac{\hat{\theta}_{\text{meta}}}{\operatorname{SE}(\hat{\theta}_{\text{meta}})}, \quad p_{\text{meta, RE}} = 2 \left( 1 - \Phi(|Z_{\text{meta}}|) \right)$$

---

### 2.3 Non-Parametric P-Value Combination

#### 1. Fisher's Combined Probability Test:
$$\chi^2 = -2 \sum_{k=1}^K \ln(p_k) \sim \chi^2_{2K}$$

$$p_{\text{Fisher}} = 1 - F_{\chi^2_{2K}}(\chi^2)$$

#### 2. Stouffer's Weighted Z-Transform:
Taking one-sided Z-scores $Z_k = \Phi^{-1}(1 - p_k / 2) \cdot \operatorname{sign}(\hat{\theta}_k)$ weighted by sample size $n_k$:

$$Z_{\text{Stouffer}} = \frac{\sum_{k=1}^K \sqrt{n_k} Z_k}{\sqrt{\sum_{k=1}^K n_k}} \sim \mathcal{N}(0, 1)$$

---

### 2.4 Multiple Hypothesis Correction

All raw p-values ($m$ tested genes) are corrected using the Benjamini-Hochberg False Discovery Rate (FDR) procedure:

$$q_{(i)} = \min_{j \ge i} \left( \frac{m \cdot p_{(j)}}{j} \right)$$

Genes with $q \le 0.05$ and $|\hat{\theta}_{\text{meta}}| \ge 0.5$ are designated as statistically significant meta-DEGs.

---

### 2.5 Pathway Enrichment Statistics

#### Hypergeometric Over-Representation Analysis (ORA):
For a gene set of size $S$, genome size $N$, and DEG list of size $k$ containing $x$ intersecting genes:

$$P(X \ge x) = \sum_{j=x}^{\min(k, S)} \frac{\binom{S}{j} \binom{N - S}{k - j}}{\binom{N}{k}}$$

#### Fast Gene Set Enrichment Analysis (fgsea):
Uses the non-parametric running-sum Kolmogorov-Smirnov-like statistic across the fully ranked list of all genes sorted by signed test statistic $s_i = \operatorname{sign}(\hat{\theta}_i) \cdot (-\log_{10} p_i)$.

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
| `gene_symbol` | String | MGI Gene Symbol |
| `n_cohorts` | Integer | Number of cohorts where gene is detected |
| `meta_log2fc` | Float | Pooled DerSimonian-Laird effect size |
| `meta_se` | Float | Standard error of pooled effect size |
| `ci_lower` | Float | Lower bound of 95% Confidence Interval |
| `ci_upper` | Float | Upper bound of 95% Confidence Interval |
| `cochran_q` | Float | Cochran's Q test statistic |
| `i2_heterogeneity` | Float | Higgins I² percentage (0 - 100%) |
| `fisher_stat` | Float | Fisher's chi-square test statistic |
| `p_fisher` | Float | Fisher's combined p-value |
| `fdr_fisher` | Float | Benjamini-Hochberg FDR of Fisher p-value |
| `p_stouffer` | Float | Stouffer's combined p-value |
| `fdr_stouffer` | Float | Benjamini-Hochberg FDR of Stouffer p-value |
| `direction_concordance`| String | `Concordant Up`, `Concordant Down`, or `Mixed` |
| `significance_flag` | Boolean | True if FDR < 0.05 and |meta_log2fc| >= 0.5 |

---

## 4. FAIR Principles & Reproducibility Guarantees
- **Findability**: All accession IDs and publications tracked in `config/datasets.yaml` and `CITATION.cff`.
- **Accessibility**: Data downloads automated via NCBI E-utilities / FTP with automatic fallback to bundled demo matrices.
- **Interoperability**: Standardized Ensembl-to-MGI identifier mappings with tidy CSV tables.
- **Reusability**: Dockerfile containerization and GitHub Actions workflow for zero-dependency execution.
