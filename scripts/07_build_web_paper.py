#!/usr/bin/env python3
"""
scripts/07_build_web_paper.py
Compiles an interactive, publication-ready scientific web paper into docs/index.html.
Embeds interactive Gene Explorer, Pathway Browser, high-res figures, and data downloads.
"""

import os
import sys
import json
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("Build-Web-Paper")

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NeuroGut-MetaSeq | Microglial Transcriptomic Meta-Analysis</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Crimson+Pro:ital,wght@0,400;0,600;0,700;1,400&family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
  <style>
    :root {
      --primary: #1a365d;
      --secondary: #2b6cb0;
      --accent: #c53030;
      --accent-blue: #2b5c8f;
      --bg: #fdfdfd;
      --text: #2d3748;
      --text-muted: #718096;
      --border: #e2e8f0;
      --card-bg: #ffffff;
      --code-bg: #f7fafc;
      --highlight: #feebc8;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: 'Crimson Pro', Georgia, serif;
      font-size: 19px;
      line-height: 1.7;
      color: var(--text);
      background-color: var(--bg);
      text-rendering: optimizeLegibility;
      -webkit-font-smoothing: antialiased;
    }
    header, main, footer {
      max-width: 900px;
      margin: 0 auto;
      padding: 0 24px;
    }
    header {
      padding-top: 50px;
      padding-bottom: 40px;
      border-bottom: 2px solid var(--border);
      text-align: center;
    }
    .badge {
      display: inline-block;
      padding: 4px 12px;
      font-family: 'Inter', sans-serif;
      font-size: 12px;
      font-weight: 600;
      letter-spacing: 0.5px;
      text-transform: uppercase;
      color: #ffffff;
      background: var(--primary);
      border-radius: 999px;
      margin-bottom: 16px;
    }
    h1 {
      font-family: 'Inter', sans-serif;
      font-size: 38px;
      font-weight: 700;
      line-height: 1.25;
      color: var(--primary);
      margin-bottom: 18px;
    }
    .subtitle {
      font-size: 22px;
      color: var(--text-muted);
      font-style: italic;
      margin-bottom: 24px;
    }
    .authors {
      font-family: 'Inter', sans-serif;
      font-size: 15px;
      color: var(--text);
      font-weight: 500;
      margin-bottom: 16px;
    }
    .meta-nav {
      display: flex;
      justify-content: center;
      gap: 12px;
      margin-top: 20px;
      flex-wrap: wrap;
    }
    .btn {
      font-family: 'Inter', sans-serif;
      font-size: 13px;
      font-weight: 600;
      padding: 8px 16px;
      border-radius: 6px;
      text-decoration: none;
      display: inline-flex;
      align-items: center;
      gap: 6px;
      transition: all 0.2s ease;
      cursor: pointer;
    }
    .btn-primary { background: var(--secondary); color: white; border: 1px solid var(--secondary); }
    .btn-primary:hover { background: var(--primary); }
    .btn-outline { background: white; color: var(--text); border: 1px solid var(--border); }
    .btn-outline:hover { background: var(--code-bg); border-color: var(--text-muted); }
    
    .abstract-box {
      background: #f7fafc;
      border-left: 4px solid var(--secondary);
      padding: 24px 28px;
      margin: 40px 0;
      border-radius: 0 8px 8px 0;
    }
    .abstract-title {
      font-family: 'Inter', sans-serif;
      font-size: 14px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.8px;
      color: var(--secondary);
      margin-bottom: 10px;
    }
    
    h2 {
      font-family: 'Inter', sans-serif;
      font-size: 26px;
      font-weight: 700;
      color: var(--primary);
      margin-top: 48px;
      margin-bottom: 18px;
      border-bottom: 1px solid var(--border);
      padding-bottom: 8px;
    }
    h3 {
      font-family: 'Inter', sans-serif;
      font-size: 20px;
      font-weight: 600;
      color: var(--text);
      margin-top: 28px;
      margin-bottom: 12px;
    }
    p { margin-bottom: 20px; }
    
    .figure-container {
      margin: 36px 0;
      background: white;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 20px;
      box-shadow: 0 2px 4px rgba(0,0,0,0.03);
    }
    .figure-img {
      width: 100%;
      height: auto;
      border-radius: 4px;
      display: block;
    }
    .figure-caption {
      margin-top: 14px;
      font-size: 15px;
      line-height: 1.5;
      color: var(--text);
      border-top: 1px solid var(--border);
      padding-top: 10px;
    }
    .figure-caption strong {
      font-family: 'Inter', sans-serif;
      font-size: 14px;
      color: var(--primary);
    }
    
    /* Interactive Explorer */
    .interactive-card {
      background: #f8fafc;
      border: 2px solid #cbd5e1;
      border-radius: 10px;
      padding: 24px;
      margin: 36px 0;
    }
    .interactive-title {
      font-family: 'Inter', sans-serif;
      font-size: 18px;
      font-weight: 700;
      color: var(--primary);
      margin-bottom: 12px;
    }
    .search-bar {
      display: flex;
      gap: 10px;
      margin-bottom: 16px;
    }
    .search-input {
      flex: 1;
      font-family: 'Inter', sans-serif;
      font-size: 16px;
      padding: 10px 14px;
      border: 1px solid var(--border);
      border-radius: 6px;
      outline: none;
    }
    .search-input:focus { border-color: var(--secondary); }
    .quick-tags {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-bottom: 20px;
    }
    .tag-btn {
      font-family: 'JetBrains Mono', monospace;
      font-size: 12px;
      padding: 4px 10px;
      background: white;
      border: 1px solid var(--border);
      border-radius: 4px;
      cursor: pointer;
      transition: all 0.15s ease;
    }
    .tag-btn:hover { background: var(--secondary); color: white; border-color: var(--secondary); }
    .gene-result-card {
      background: white;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 18px;
      font-family: 'Inter', sans-serif;
      font-size: 14px;
    }
    .stats-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin: 14px 0;
    }
    .stat-pill {
      background: #f1f5f9;
      padding: 10px;
      border-radius: 6px;
    }
    .stat-label { font-size: 11px; text-transform: uppercase; color: var(--text-muted); font-weight: 600; }
    .stat-value { font-size: 16px; font-weight: 700; color: var(--primary); margin-top: 2px; }

    /* Tables */
    table {
      width: 100%;
      border-collapse: collapse;
      font-family: 'Inter', sans-serif;
      font-size: 14px;
      margin: 24px 0;
    }
    th, td {
      padding: 10px 14px;
      border: 1px solid var(--border);
      text-align: left;
    }
    th { background: #f8fafc; font-weight: 600; color: var(--primary); }
    tr:nth-child(even) { background: #fcfcfc; }

    /* Code & Citation */
    pre, code {
      font-family: 'JetBrains Mono', monospace;
      font-size: 13px;
      background: var(--code-bg);
      border-radius: 4px;
    }
    pre {
      padding: 16px;
      overflow-x: auto;
      border: 1px solid var(--border);
      margin: 20px 0;
    }
    code { padding: 2px 6px; }

    footer {
      border-top: 2px solid var(--border);
      padding: 40px 24px;
      margin-top: 60px;
      font-family: 'Inter', sans-serif;
      font-size: 14px;
      color: var(--text-muted);
      text-align: center;
    }
  </style>
</head>
<body>

<header>
  <span class="badge">Open-Science Transcriptomics</span>
  <h1>Cross-Study Meta-Analysis of Microglial Gene Expression in Response to Microbiome Metabolites</h1>
  <div class="subtitle">Systematic integration of bulk RNA-seq cohorts reveals gut microbiota control over microglial neuroinflammatory checkpoints</div>
  <div class="authors">
    <strong>NeuroGut-MetaSeq Consortium</strong> &bull; Samyak Meshram
  </div>
  <div class="meta-nav">
    <a href="#interactive-explorer" class="btn btn-primary">🔍 Interactive Gene Explorer</a>
    <a href="#results" class="btn btn-outline">📊 Explore Figures</a>
    <a href="https://github.com/samyakmeshram/NeuroGut-MetaSeq" target="_blank" class="btn btn-outline">💻 GitHub Repository</a>
    <a href="assets/microglia_meta_analysis_summary.csv" download class="btn btn-outline">⬇️ Download Results CSV</a>
  </div>
</header>

<main>

  <div class="abstract-box">
    <div class="abstract-title">Abstract</div>
    <p>
      Microglia are brain-resident myeloid sentinels whose physiological maturation and immune alertness are continuously calibrated by gut microbiota-derived metabolites, particularly Short-Chain Fatty Acids (SCFAs). However, individual RNA-sequencing studies are frequently limited by statistical power, protocol differences, and cohort-specific artifacts. Here, we present <strong>NeuroGut-MetaSeq</strong>, an open-science computational meta-analysis framework synthesizing independent transcriptomic cohorts (spanning germ-free, antibiotic-depleted, and metabolite-perturbed rodent models). By applying negative binomial generalized linear models and inverse-variance DerSimonian-Laird random-effects pooling alongside Fisher and Stouffer combination tests, we define a high-confidence consensus regulon of microbiome-dependent microglial genes. Microbiome depletion consistently blunts homeostatic markers (<em>Tmem119</em>, <em>Cx3cr1</em>, <em>P2ry12</em>, <em>Ffar2</em>) while priming pro-inflammatory cascades (<em>Tnf</em>, <em>Il1b</em>, <em>Nfkb1</em>, <em>Ccl2</em>). Pathway enrichment uncovers prominent modulation of NF-κB, chemokine signaling, and phagocytic receptors. We provide an interactive web platform and FAIR-compliant containerized pipeline to guarantee complete academic reproducibility.
    </p>
  </div>

  <section id="introduction">
    <h2>1. Introduction</h2>
    <p>
      Over the past decade, the microbiota-gut-brain axis has emerged as a fundamental regulator of neurological health and central nervous system (CNS) immunology. Microglia, originating early in embryogenesis from yolk sac progenitors, persist throughout adulthood through self-renewal and actively survey the brain parenchyma. Landmark studies have established that germ-free (GF) mice or mice chronically treated with broad-spectrum antibiotics (ABX) exhibit microglia with immature, hyper-ramified morphologies and compromised responses to viral and bacterial challenges.
    </p>
    <p>
      Crucially, oral administration of short-chain fatty acids (SCFAs)—principally acetate, propionate, and butyrate produced by bacterial fermentation of non-digestible dietary fibers—is sufficient to rescue microglial maturity and suppress maladaptive neuroinflammation. Despite these findings, independent transcriptomic datasets often yield diverging lists of statistically significant differentially expressed genes (DEGs), attributable to small cohort sizes, disparate isolation methodologies (e.g. CD11b FACS vs magnetic beads), and varying sequencing platforms. Computational meta-analysis across multiple independent cohorts provides a robust mathematical solution to separate universal biological signals from study-specific noise.
    </p>
  </section>

  <section id="landscape">
    <h2>2. Multi-Cohort Transcriptomic Landscape</h2>
    <p>
      We systematically curated public RNA-seq datasets from the NCBI Gene Expression Omnibus (GEO) conforming to strict biological and computational standards: (1) isolated mouse brain microglia, (2) high-throughput bulk RNA sequencing, (3) raw integer read count matrices available, and (4) well-annotated experimental conditions contrasting microbiome depletion or metabolite exposure against colonized controls.
    </p>

    <table>
      <thead>
        <tr>
          <th>GEO Accession</th>
          <th>Study Reference</th>
          <th>Experimental Comparison</th>
          <th>Samples ($n$)</th>
          <th>Tissue / Sorting</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td><strong>GSE107925</strong></td>
          <td>Thion et al., <em>Cell</em> 2018</td>
          <td>Germ-Free (GF) vs SPF Control (Adult Males & Females)</td>
          <td>12</td>
          <td>FACS CD11b+ CD45low adult microglia</td>
        </tr>
        <tr>
          <td><strong>GSE108045</strong></td>
          <td>Thion et al., <em>Cell</em> 2018</td>
          <td>Antibiotic (ABX) Gut Depletion vs SPF Control</td>
          <td>12</td>
          <td>Adult brain microglia</td>
        </tr>
        <tr>
          <td><strong>GSE266602</strong></td>
          <td>Wang et al., <em>Exp Neurol</em> 2024</td>
          <td>Germ-Free (GF_Sham) vs Colonized (SPF_Sham)</td>
          <td>6</td>
          <td>Adult primary microglia</td>
        </tr>
        <tr>
          <td><strong>GSE186210</strong></td>
          <td>Matt et al., <em>J Neurosci</em> 2023</td>
          <td>Zero-Fiber (SCFA-Deficient) vs High-Fiber Diet</td>
          <td>46</td>
          <td>Adult brain microglia</td>
        </tr>
      </tbody>
    </table>
  </section>

  <section id="methods">
    <h2>3. Statistical & Mathematical Methodology</h2>
    <p>
      Our analytical workflow couples individual negative binomial generalized linear models (GLMs) with multi-cohort random-effects pooling:
    </p>
    <div style="background: #f8fafc; padding: 18px 22px; border-left: 4px solid var(--primary); margin: 24px 0; border-radius: 4px;">
      <p style="margin-bottom: 8px;"><strong>1. Within-Study Differential Expression:</strong> Modeled with Negative Binomial GLMs $\operatorname{NB}(\mu_{ij}, \alpha_i)$ with Wald test statistics $W_i = \hat{\beta}_i / \operatorname{SE}(\hat{\beta}_i)$.</p>
      <p style="margin-bottom: 8px;"><strong>2. Random-Effects Effect Size Pooling:</strong> DerSimonian-Laird inverse-variance weighting estimating pooled $\hat{\theta}_{\text{meta}}$ and between-study variance $\tau^2$ alongside Cochran's $Q$ and Higgins $I^2$ inconsistency metrics.</p>
      <p style="margin-bottom: 8px;"><strong>3. Non-Parametric P-Value Integration:</strong> Fisher's combination ($\chi^2 = -2 \sum \ln p_i$) and sample-size weighted Stouffer's Z-transform.</p>
      <p style="margin-bottom: 0;"><strong>4. Multiple Hypothesis Adjustment:</strong> Benjamini-Hochberg False Discovery Rate (FDR $\le 0.05$) across all synthesized genes.</p>
    </div>
  </section>

  <section id="results">
    <h2>4. Results & Meta-Analysis Findings</h2>

    <h3>Figure 1: Cross-Cohort Principal Component Analysis</h3>
    <div class="figure-container">
      <img src="assets/fig1_cohort_pca.png" alt="Figure 1: PCA Plot" class="figure-img">
      <div class="figure-caption">
        <strong>Figure 1 | Cross-cohort principal component analysis of microglial transcriptomes.</strong>
        Samples from independent cohorts cluster distinctly along PC1 and PC2, reflecting both biological condition (reference vs microbiome-perturbed) and study-specific sequencing batch characteristics, underscoring the necessity of statistical meta-analytic synthesis over naive sample pooling.
      </div>
    </div>

    <h3>Figure 2: Consensus Meta-Analysis Volcano Plot</h3>
    <div class="figure-container">
      <img src="assets/fig2_meta_volcano.png" alt="Figure 2: Volcano Plot" class="figure-img">
      <div class="figure-caption">
        <strong>Figure 2 | Multi-cohort meta-analysis volcano plot.</strong>
        Consensus effect size ($\log_2\text{FC}$) versus $-\log_{10}(\text{Fisher FDR adjusted } p\text{-value})$. Red points denote consensus significantly upregulated genes (including <em>Tnf</em>, <em>Il1b</em>, <em>Nfkb1</em>, <em>Ccl2</em>, <em>Trem2</em>); blue points denote consensus downregulated genes (including <em>Tmem119</em>, <em>Cx3cr1</em>, <em>P2ry12</em>, <em>Ffar2</em>).
      </div>
    </div>

    <h3>Figure 3: Forest Plots of Landmark Genes</h3>
    <div class="figure-container">
      <img src="assets/fig3_forest_plots.png" alt="Figure 3: Forest Plots" class="figure-img">
      <div class="figure-caption">
        <strong>Figure 3 | Forest plots of landmark microglial and inflammatory markers.</strong>
        Cohort-specific effect sizes ($\log_2\text{FC} \pm 95\%\text{CI}$) across GSE107925, GSE108045, and GSE266602, alongside the pooled DerSimonian-Laird random-effects summary diamond (navy). Notice the high cross-study concordance for core inflammatory regulators (<em>Tnf</em>, <em>Il1b</em>, <em>Nfkb1</em>) and homeostatic checkpoints (<em>Tmem119</em>, <em>Cx3cr1</em>).
      </div>
    </div>

    <h3>Figure 4: Clustered Heatmap of Top Consensus Genes</h3>
    <div class="figure-container">
      <img src="assets/fig4_consensus_heatmap.png" alt="Figure 4: Heatmap" class="figure-img">
      <div class="figure-caption">
        <strong>Figure 4 | Hierarchical clustering of top consensus meta-DEGs across cohorts.</strong>
        Heatmap displays standardized effect sizes ($\log_2\text{FC}$) across independent studies alongside the pooled meta-estimate, highlighting core gene modules uniformly altered by gut microbiome state.
      </div>
    </div>

    <h3>Figure 5: Functional Pathway Enrichment</h3>
    <div class="figure-container">
      <img src="assets/fig5_pathway_enrichment.png" alt="Figure 5: Pathways" class="figure-img">
      <div class="figure-caption">
        <strong>Figure 5 | Over-representation and pathway enrichment of consensus microglial regulon.</strong>
        Significant enrichment ($-\log_{10}\text{FDR}$) of NF-κB signaling, cytokine-cytokine receptor interaction, chemokine signaling, and phagosomal maturation pathways.
      </div>
    </div>
  </section>

  <!-- Interactive Explorer Section -->
  <section id="interactive-explorer">
    <h2>5. Interactive Gene Explorer</h2>
    <p>
      Search any microglial, inflammatory, or metabolic gene in our meta-analysis database to inspect its cross-cohort statistics and pooled random-effects parameters:
    </p>

    <div class="interactive-card">
      <div class="interactive-title">Microglial Meta-Transcriptome Search</div>
      <div class="search-bar">
        <input type="text" id="geneSearch" class="search-input" placeholder="Type a gene symbol (e.g. Tnf, Cx3cr1, Nfkb1, Trem2, Ffar2)..." value="Tnf">
        <button id="searchBtn" class="btn btn-primary">Search Gene</button>
      </div>

      <div class="quick-tags">
        <span style="font-size: 13px; font-weight: 600; color: var(--text-muted); align-self: center;">Quick Picks:</span>
        <button class="tag-btn" onclick="selectGene('Tnf')">Tnf</button>
        <button class="tag-btn" onclick="selectGene('Il1b')">Il1b</button>
        <button class="tag-btn" onclick="selectGene('Nfkb1')">Nfkb1</button>
        <button class="tag-btn" onclick="selectGene('Ccl2')">Ccl2</button>
        <button class="tag-btn" onclick="selectGene('Cx3cr1')">Cx3cr1</button>
        <button class="tag-btn" onclick="selectGene('Tmem119')">Tmem119</button>
        <button class="tag-btn" onclick="selectGene('P2ry12')">P2ry12</button>
        <button class="tag-btn" onclick="selectGene('Ffar2')">Ffar2</button>
        <button class="tag-btn" onclick="selectGene('Trem2')">Trem2</button>
        <button class="tag-btn" onclick="selectGene('Apoe')">Apoe</button>
        <button class="tag-btn" onclick="selectGene('Stat1')">Stat1</button>
      </div>

      <div id="geneResultBox" class="gene-result-card">
        <!-- Injected via JavaScript -->
      </div>
    </div>
  </section>

  <section id="discussion">
    <h2>6. Discussion & Biological Implications</h2>
    <p>
      Our cross-study meta-analysis establishes that microbiome depletion does not merely silence microglial activity; rather, it destabilizes their homeostatic transcriptional identity. Microglia in germ-free or antibiotic-treated brains experience significant downregulation of resting checkpoints including <em>Tmem119</em>, <em>Cx3cr1</em>, and <em>P2ry12</em>, alongside SCFA receptors such as <em>Ffar2</em> (GPR43). Simultaneously, these cells display baseline hyperactivity of transcriptional programs driven by <em>Nfkb1</em>, <em>Rela</em>, and <em>Stat1</em>, with elevated transcripts for classic inflammatory mediators (<em>Tnf</em>, <em>Il1b</em>, <em>Ccl2</em>).
    </p>
    <p>
      This apparent paradox—immature structural phenotype paired with hypersensitive pro-inflammatory gene readiness—explains why microbiome-deficient animals exhibit exaggerated and dysregulated neuroinflammation upon subsequent injury or pathogenic challenge. Gut microbial metabolites, particularly butyrate and propionate, act as vital epigenetic and signaling brakes that keep microglial immune cascades properly tuned.
    </p>
  </section>

  <section id="reproducibility">
    <h2>7. Reproducibility & FAIR Data Access</h2>
    <p>
      In accordance with open-science principles, all data matrices, analysis code, Docker configurations, and results tables are publicly accessible and fully reproducible in one command:
    </p>

    <pre><code># 1. Clone repository
git clone https://github.com/samyakmeshram/NeuroGut-MetaSeq.git
cd NeuroGut-MetaSeq

# 2. Run complete end-to-end pipeline and compile this web paper
make all</code></pre>

    <div style="margin: 24px 0;">
      <a href="assets/microglia_meta_analysis_summary.csv" download class="btn btn-primary">⬇️ Download Consensus Meta-Analysis CSV</a>
      <a href="assets/pathway_summary.csv" download class="btn btn-outline">⬇️ Download Pathway Enrichment CSV</a>
      <a href="https://github.com/samyakmeshram/NeuroGut-MetaSeq/blob/main/Dockerfile" target="_blank" class="btn btn-outline">🐳 View Dockerfile</a>
    </div>
  </section>

  <section id="citation">
    <h2>8. How to Cite</h2>
    <pre><code>@software{meshram2026neurogut,
  author       = {Samyak Meshram},
  title        = {NeuroGut-MetaSeq: Cross-Study RNA-Seq Meta-Analysis of Microglial Transcriptomic Signatures in Response to Microbiome Depletion and Microbial Metabolites},
  year         = {2026},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{https://github.com/samyakmeshram/NeuroGut-MetaSeq}}
}</code></pre>
  </section>

</main>

<footer>
  <p>&copy; 2026 NeuroGut-MetaSeq Open-Science Initiative &bull; Built with reproducibility and FAIR data principles.</p>
</footer>

<script>
  // Embedded meta-analysis dataset
  const META_DATA = __META_DATA_JSON__;

  function renderGeneCard(gene) {
    const box = document.getElementById("geneResultBox");
    const entry = META_DATA.find(d => d.gene_symbol.toLowerCase() === gene.toLowerCase());

    if (!entry) {
      box.innerHTML = `<div style="color: var(--accent); font-weight: 600;">Gene '<strong>${gene}</strong>' not found in consensus meta-analysis dataset. Please verify the gene symbol.</div>`;
      return;
    }

    const isSig = entry.significance_flag;
    const isUp = entry.meta_log2fc > 0;
    const color = isUp ? "#c53030" : "#2b6cb0";
    const statusBadge = isSig ? (isUp ? "<span style='color:#c53030; font-weight:bold;'>Consensus Upregulated (FDR < 0.05)</span>" : "<span style='color:#2b6cb0; font-weight:bold;'>Consensus Downregulated (FDR < 0.05)</span>") : "<span style='color:#718096;'>Not Stably Differential</span>";

    box.innerHTML = `
      <div style="display: flex; justify-content: space-between; align-items: baseline; border-bottom: 1px solid var(--border); padding-bottom: 8px;">
        <h3 style="margin: 0; font-size: 22px; color: var(--primary);">Gene: <em>${entry.gene_symbol}</em></h3>
        <span style="font-size: 13px;">${statusBadge}</span>
      </div>
      <div class="stats-grid">
        <div class="stat-pill">
          <div class="stat-label">Pooled Effect (Log2FC)</div>
          <div class="stat-value" style="color: ${color};">${entry.meta_log2fc > 0 ? "+" : ""}${entry.meta_log2fc.toFixed(3)}</div>
        </div>
        <div class="stat-pill">
          <div class="stat-label">95% Confidence Interval</div>
          <div class="stat-value" style="font-size: 14px;">[${entry.ci_lower.toFixed(3)}, ${entry.ci_upper.toFixed(3)}]</div>
        </div>
        <div class="stat-pill">
          <div class="stat-label">Fisher Combined FDR</div>
          <div class="stat-value">${entry.fdr_fisher < 0.001 ? entry.fdr_fisher.toExponential(2) : entry.fdr_fisher.toFixed(4)}</div>
        </div>
        <div class="stat-pill">
          <div class="stat-label">Heterogeneity (I²)</div>
          <div class="stat-value">${entry.i2_heterogeneity.toFixed(1)}% (Q=${entry.cochran_q.toFixed(2)})</div>
        </div>
        <div class="stat-pill">
          <div class="stat-label">Cohorts Detected</div>
          <div class="stat-value">${entry.n_cohorts} (${entry.cohorts_detected})</div>
        </div>
        <div class="stat-pill">
          <div class="stat-label">Direction Concordance</div>
          <div class="stat-value" style="font-size: 14px;">${entry.direction_concordance}</div>
        </div>
      </div>
    `;
  }

  function selectGene(g) {
    document.getElementById("geneSearch").value = g;
    renderGeneCard(g);
  }

  document.getElementById("searchBtn").addEventListener("click", () => {
    const val = document.getElementById("geneSearch").value.trim();
    if (val) renderGeneCard(val);
  });

  document.getElementById("geneSearch").addEventListener("keyup", (e) => {
    if (e.key === "Enter") {
      const val = e.target.value.trim();
      if (val) renderGeneCard(val);
    }
  });

  // Initial render
  renderGeneCard("Tnf");
</script>

</body>
</html>
"""

def main():
    meta_path = "results/meta_results/microglia_meta_analysis_summary.csv"
    if not os.path.exists(meta_path):
        logger.error(f"Meta-analysis summary not found at {meta_path}. Run step 04 first.")
        sys.exit(1)

    meta_df = pd.read_csv(meta_path)
    records = meta_df.to_dict(orient="records")
    meta_json = json.dumps(records)

    html_content = HTML_TEMPLATE.replace("__META_DATA_JSON__", meta_json)

    os.makedirs("docs/assets", exist_ok=True)
    out_html = "docs/index.html"
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html_content)

    # Copy CSV files to docs/assets/ for public download
    meta_df.to_csv("docs/assets/microglia_meta_analysis_summary.csv", index=False)
    if os.path.exists("results/pathways/pathway_summary.csv"):
        import shutil
        shutil.copy("results/pathways/pathway_summary.csv", "docs/assets/pathway_summary.csv")

    logger.info("=" * 60)
    logger.info(f"[SUCCESS] Standalone publication web paper compiled -> {out_html}")
    logger.info(f"File size: {os.path.getsize(out_html):,} bytes")
    logger.info(f"Embedded {len(records)} genes in interactive explorer.")
    logger.info("Ready for GitHub Pages deployment!")
    logger.info("=" * 60)

if __name__ == "__main__":
    main()
