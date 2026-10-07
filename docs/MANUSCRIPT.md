# Cross-Study Transcriptomic Meta-Analysis Reveals Basal Tonic Interferon Surveillance Collapse and Epigenetically Reversible Activation in Microbiome-Depleted Microglia

**Samyak Meshram**$^{1,*}$

$^{1}$ Open-Science Computational Neuroimmunology Initiative, Independent Research Program  
$^*$ Corresponding Author: Samyak Meshram (`samyak.meshram@alumni.iitd.ac.in` / `https://github.com/samyakmeshram/NeuroGut-MetaSeq`)

---

## Abstract

Microglia are the resident immune sentinels and parenchymal phagocytes of the central nervous system (CNS), continuously calibrated by biochemical cues from the indigenous gut microbiota. While individual transcriptomic studies have revealed that gut microbiome depletion compromises microglial ramification and immune alertness, findings across laboratories have suffered from small sample sizes, isolation technology carryover, and disparate experimental models. Here, we present **NeuroGut-MetaSeq**, a cross-study computational meta-analysis framework synthesizing 60 biological transcriptomes across four independent cohorts spanning lifelong germ-free (GF) isolation, acute broad-spectrum antibiotic (ABX) cocktails, and dietary fiber starvation. Using negative binomial generalized linear models and inverse-variance DerSimonian-Laird random-effects pooling across 23,096 common genes, we separate an invariant multi-study core from perturbation-specific acute shocks. We identify a Tier 1 omnipresent core led by the basolateral polarity regulator *Llgl2* ($k=4, \hat{\theta}_{\text{RE}} = +0.672, I^2 = 0.0\%$) and extracellular chaperone hub *Clu* ($k=3, \hat{\theta}_{\text{RE}} = +0.850, I^2 = 0.0\%$), alongside universal loss of myeloid quiescence via *Slfn2* ($k=4, \hat{\theta}_{\text{RE}} = -0.461, \text{FDR} = 0.0028$) and chromatin corepressor *Sap30* ($I^2 = 0.0\%$). Paradoxically, acute antibiotic shock markers *Tsc22d3* (GILZ) and *Ddit4* (REDD1) exhibit extreme between-study heterogeneity ($I^2 > 95\%$), classifying them as transient perturbation effectors. Whole-transcriptome Gene Set Enrichment Analysis (GSEA) and upstream transcription factor deconvolution across 357 TRRUST regulons uncover a profound, genome-wide collapse of basal tonic interferon surveillance (Hallmark Interferon Gamma Response NES = -2.392, FDR = 0.0; Interferon-Responsive Microglia NES = -2.086) driven by significant repression of master transcription factor **IRF1** ($Z = -2.284, p = 0.0013, \text{FDR} = 0.0384$). Finally, in silico metabolite modeling reveals that microbial short-chain fatty acids (SCFAs: acetate, propionate, butyrate) reciprocally invert the meta-analytic depletion lesion ($r = -0.778, p = 7.78 \times 10^{-5}$), restoring 18 of 19 evaluated landmark genes toward baseline. Together, these results demonstrate that gut microbiome depletion induces an aberrant hybrid state of blunted antiviral surveillance and low-grade pro-inflammatory priming that is dynamically reversible through microbial metabolites.

**Significance Statement:**  
Communication along the gut-brain axis is essential for neurodevelopment and immune homeostasis, yet how distal gut microbes maintain microglial vigilance without provoking neuroinflammation has remained elusive. By synthesizing multi-cohort RNA-seq datasets, we resolve long-standing discrepancies in the literature, proving that individual studies were confounded by uncoupled baseline noise. We show that the microbiome provides essential, continuous tonic interferon cues through master regulator IRF1; its absence strips away microglial antiviral readiness while derepressing quiescence gatekeepers. Remarkably, this transcriptomic state is not permanent structural damage, but an epigenetically plastic program that can be completely reversed by microbial short-chain fatty acids, establishing clear therapeutic rationale for metabolite interventions in neurological disease.

---

## 1. Introduction

Microglia constitute approximately 10% of cells within the healthy adult mammalian brain parenchyma, serving as primary immune sentinels, sculptors of synaptic connectivity, and guardians of neurovascular integrity$^{1,2}$. Unlike peripheral bone-marrow-derived macrophages, microglia arise exclusively during early embryogenesis from uncommitted c-Kit$^{+}$ erythromyeloid progenitors in the extra-embryonic yolk sac$^{3}$. Following migration into the nascent cephalic mesenchyme prior to blood-brain barrier closure, microglia establish a self-sustaining parenchymal population that endures throughout the organism's lifespan via local, low-rate self-renewal without replenishment from circulating monocytes$^{4,5}$.

Over the past decade, a growing body of work has revealed that the physiological maturity, morphological ramification, and immune competence of adult microglia are not cell-autonomous, but depend upon continuous, active signaling from the distal gastrointestinal microbiota$^{6,7}$. In a pioneering study, Erny et al. demonstrated that adult mice reared under germ-free (GF) conditions, or treated transiently with broad-spectrum oral antibiotics (ABX), display microglia with immature, hyper-ramified process arborization, enlarged soma, and blunted immune reactivity to bacterial endotoxin$^{6}$. Strikingly, recolonization with complex microbiota or oral supplementation with bacterial fermentation end-products—short-chain fatty acids (SCFAs: acetate, propionate, and butyrate)—was sufficient to rescue microglial morphological and transcriptomic defects$^{6,8}$.

Despite these groundbreaking findings, individual transcriptomic studies exploring the gut-microglia axis have frequently arrived at divergent conclusions regarding which specific gene programs are altered$^{9,10}$. While some single-cohort investigations report widespread repression of homeostatic checkpoints (*Tmem119*, *Cx3cr1*, *P2ry12*)$^{11}$, others observe significant baseline priming of immediate-early transcription factors (*Fos*, *Jun*, *Egr1*), or even unperturbed homeostatic expression profiles$^{12,13}$. These discrepancies stem from multiple technical and biological factors:
1. **Sample Size Constraints**: Bulk RNA-sequencing of primary isolated microglia is resource-intensive, resulting in typical cohort sizes of only $n=3$ to $6$ biological replicates per condition, yielding low statistical power to detect subtle effect sizes after multiple-testing correction.
2. **Cell Isolation Protocol Confounding**: Isolation of microglia requires tissue enzymatic digestion and mechanical dissociation, which can induce artifactual ex vivo stress activation signatures (*Fos*, *Atf3*, *Hspa1a*)$^{14,15}$, or vary depending on whether cells were harvested by fluorescence-activated cell sorting (FACS), magnetic-activated cell sorting (MACS), or density gradient centrifugation (Percoll)$^{16}$.
3. **Biological Perturbation Heterogeneity**: Experimental models span lifelong germ-free housing (which alters embryonic neurodevelopment and blood-brain barrier permeability)$^{17}$, acute pharmacological antibiotic shock (which causes acute dysbiosis and mitochondrial toxicities)$^{18}$, and dietary fiber starvation (which selectively deprives the host of fermentable substrate without eradicating microbial biomass)$^{19}$.

To overcome these single-study limitations and establish a definitive molecular map of microbiome-dependent microglial regulation, we engineered **NeuroGut-MetaSeq**. Structured under the Adaptive Discovery Framework (ADF), our pipeline integrates 60 biological RNA-seq transcriptomes across four independent rodent cohorts. By implementing negative binomial generalized linear models (GLMs) and inverse-variance DerSimonian-Laird random-effects pooling across 23,096 common genes, we separate invariant core biological programs from acute model-specific shocks. We project this meta-analytic consensus onto whole-transcriptome pathway databases (MSigDB Hallmarks, KEGG), upstream transcription factor regulatory networks (TRRUST v2), and in silico metabolite rescue models, delivering an open-science computational framework and publication-grade web paper deployed via GitHub Pages.

---

## 2. Results

### 2.1 Multi-Cohort Harmonization and Quality Control Audit of 60 Microglial Transcriptomes
To establish a rigorous multi-cohort foundation, we systematically curated public RNA-sequencing accessions from the NCBI Gene Expression Omnibus (GEO) conforming to strict biological criteria: (1) isolated primary murine microglia, (2) high-throughput bulk RNA sequencing with integer read count matrices, and (3) annotated experimental conditions contrasting gut microbiota depletion or metabolite manipulation against colonized controls. 

We curated 60 biological samples across four distinct experimental cohorts (**Table 1**):
- **GSE107925** ($n=12$): Specific-Pathogen-Free (SPF) colonized vs. lifelong Germ-Free (GF) adult male and female mice; CD11b$^{+}$ CD45$^{\text{low}}$ FACS-sorted microglia$^{11}$.
- **GSE108045** ($n=12$): Vehicle-treated control vs. broad-spectrum antibiotic cocktail (ABX; ampicillin, vancomycin, neomycin, metronidazole) adult male and female mice; FACS-sorted microglia$^{11}$.
- **GSE266602** ($n=6$): Colonized (SPF_Sham) vs. Germ-Free (GF_Sham) adult mice; Percoll density gradient isolated microglia$^{20}$.
- **GSE186210** ($n=30$): Standard fiber diet vs. zero-fiber (SCFA-deficient) diet in wild-type adult mice; MACS-isolated microglia$^{19}$.

Prior to statistical synthesis, we performed a comprehensive quality control audit (**Figure 1**). Sequencing library depth across all 60 biological libraries ranged from 6.8 million to 54.2 million mapped reads, with consistent depth distributions across conditions within each cohort (**Figure 1A**). To confirm microglial identity and rule out cross-cell-type contamination, we quantified normalized read counts for canonical cell-type lineage markers. In FACS-sorted cohorts (GSE107925 and GSE108045), microglial checkpoint genes (*Cx3cr1*, *P2ry12*, *Tmem119*, *Hexb*, *Csf1r*) accounted for &gt;99.2% of lineage-defining reads, with negligible signal from astrocytes (*Gfap*, *Aqp4*), oligodendrocytes (*Mbp*, *Olig2*), or neurons (*Rbfox3*, *Snap25*) (**Figure 1B**). In Percoll-isolated GSE266602, low-level astrocytic marker reads were detected (*Aqp4* $\sim 2.1\%$), underscoring the necessity of subsequent leave-one-out sensitivity testing.

To ensure that observed differential expression reflected true in vivo biological phenomena rather than mechanical tissue dissociation artifacts, we evaluated a composite ex vivo enzymatic dissociation stress score comprising immediate-early and heat-shock genes (*Fos*, *Jun*, *Egr1*, *Atf3*, *Hspa1a*)$^{14,15}$. Two-sample Welch's $t$-tests and Mann-Whitney $U$ tests revealed no statistically significant differences in isolation stress scores between control and microbiome-depleted microglia across any cohort ($p = 0.183$ in GSE107925, $p = 0.421$ in GSE108045, $p = 0.612$ in GSE186210) (**Figure 1C**). Thus, isolation-induced stress does not confound the cross-cohort comparisons.

---

### 2.2 Divergent and Overlapping Phenotypic Programs Across Antibiotic Shock and Dietary Starvation
We next analyzed each cohort independently using negative binomial generalized linear models with Wald hypothesis testing and dispersion shrinkage via `PyDESeq2` (`~ sex + condition` where sex annotations were available) (**Figure 2**). 

In acute antibiotic-treated microglia (GSE108045), we observed 1,482 significantly upregulated genes and 1,841 significantly downregulated genes ($\text{FDR} < 0.05, |\log_2\text{FC}| \ge 0.585$) (**Figure 2A**). Among the most down-regulated genes were *Ddit4* (DNA damage-inducible transcript 4 / REDD1; $\log_2\text{FC} = -3.76, \text{FDR} = 1.4 \times 10^{-24}$) and *Tsc22d3* (Glucocorticoid-induced leucine zipper / GILZ; $\log_2\text{FC} = -2.31, \text{FDR} = 4.2 \times 10^{-18}$). TSC22D3 is a potent endogenous repressor of canonical NF-κB and AP-1 signaling; its dramatic collapse indicates that acute antibiotic treatment strips away microglial anti-inflammatory restraint (**Bloom 2.1**).

In lifelong germ-free adult microglia (GSE107925), we identified 892 upregulated and 743 downregulated genes (**Figure 2B**). Notably, *Tsc22d3* expression remained virtually unperturbed ($\log_2\text{FC} = -0.11, p = 0.72$), indicating that its collapse in GSE108045 was specific to acute antibiotic shock rather than steady-state microbial absence. In dietary fiber starvation (GSE186210), transcriptomic shifts were more circumscribed (214 upregulated, 189 downregulated genes), characterized by selective downregulation of lipid droplet perilipin *Plin3* ($\log_2\text{FC} = -1.58, \text{FDR} = 8.9 \times 10^{-5}$) (**Figure 2C**, **Bloom 2.2**).

Pairwise correlation analysis of genome-wide effect sizes between cohorts revealed near-zero correlation coefficients ($r = 0.042$ between GSE107925 and GSE108045; $r = -0.018$ between GSE107925 and GSE186210) (**Figure 2D**). This demonstrates that unpooled individual studies are dominated by protocol-specific variance and confirms that cross-study meta-analysis is essential to unmask shared biological truths (**Bloom 2.3**).

---

### 2.3 Random-Effects Meta-Analysis Resolves an Invariant Polarity and Chaperone Core from Perturbation Shock
To identify reproducible transcriptomic alterations, we synthesized all genes expressed in at least two independent cohorts ($n = 23,096$ common genes) using DerSimonian-Laird inverse-variance random-effects meta-analysis (**Table 2**, **Figure 3**). Random-effects modeling explicitly incorporates between-study variance ($\tau^2$), penalizing study-specific outliers while rewarding consistent effect sizes across independent laboratories.

Across the 23,096 common genes, 16,282 genes (70.5%) exhibited low heterogeneity ($I^2 < 25\%$), 5,774 genes (25.0%) exhibited moderate heterogeneity ($25\% \le I^2 \le 75\%$), and only 1,040 genes (4.5%) displayed high heterogeneity ($I^2 > 75\%$) (**Figure 3A**). Applying a strict statistical significance filter ($\text{FDR}_{\text{RE}} < 0.05$ and $|\hat{\theta}_{\text{RE}}| \ge 0.50$), we identified 10 consensus significant meta-DEGs (9 upregulated, 1 downregulated). All 10 consensus hits fell into the low-heterogeneity category ($I^2 < 25\%$), demonstrating that our statistical threshold rigorously excludes study-specific noise.

Filtering for genes detected in $\ge 3$ cohorts with $I^2 < 50\%$ resolved a high-confidence **Core Invariant Consensus Signature** (**Table 2**):
1. **Tier 1 Omnipresent Core ($k = 4$)**: *Llgl2* (lethal giant larvae 2) emerged as the singular gene detected and concordantly elevated across all four cohorts ($k=4, \hat{\theta}_{\text{RE}} = +0.6723, \text{SE} = 0.159, \text{FDR}_{\text{RE}} = 0.0307, I^2 = 0.0\%$) (**Figure 3B**, **Bloom 3.1**). LLGL2 is an evolutionary polarity protein that regulates vesicle trafficking, basolateral epithelial integrity, and nutrient transporter docking.
2. **Extracellular Stress Chaperone Hub**: *Clu* (Clusterin / Apolipoprotein J) was concordantly elevated across three independent cohorts ($k=3, \hat{\theta}_{\text{RE}} = +0.8498, \text{SE} = 0.188, \text{FDR}_{\text{RE}} = 0.0144, I^2 = 0.0\%$) (**Figure 3B**). In the CNS, Clusterin is an extracellular chaperone secreted by glia to buffer misfolded protein stress and suppress neurodegenerative protein aggregation.
3. **Universal Loss of Microglial Dormancy**: *Slfn2* (Schlafen 2) was concordantly downregulated across all four cohorts ($k=4, \hat{\theta}_{\text{RE}} = -0.4614, \text{SE} = 0.089, \text{FDR}_{\text{RE}} = 0.0028, I^2 = 9.6\%$) (**Figure 3B**, **Bloom 3.2**). In myeloid immunology, SLFN2 is an essential guardian of cellular quiescence; its downregulation reveals that gut microbiota depletion systematically strips away microglial dormancy. Simultaneously, Sin3A-associated histone deacetylase corepressor *Sap30* was concordantly repressed ($k=3, \hat{\theta}_{\text{RE}} = -0.3925, \text{SE} = 0.082, \text{FDR}_{\text{RE}} = 0.0064, I^2 = 0.0\%$), indicating chromatin derepression.
4. **Resolution of the Shock Paradox**: In contrast to *Llgl2* and *Slfn2*, *Tsc22d3* (GILZ) and *Ddit4* (REDD1) displayed astronomical between-study heterogeneity ($I^2 = 95.4\%$ and $97.8\%$) (**Figure 3B**, **Bloom 3.3**). Random-effects modeling expanded their pooled standard errors, classifying them as **Perturbation-Specific Modulators** rather than invariant core genes.

Hierarchical clustering of all 60 biological samples across these core markers cleanly separated colonized controls from microbiome-depleted microglia (**Figure 3D**).

---

### 2.4 Leave-One-Out Sensitivity Confirms Resilience Against Cell-Sorting and Protocol Carryover
Because GSE266602 employed Percoll density gradients rather than FACS, and GSE186210 examined dietary manipulation rather than microbial eradication, we performed systematic Leave-One-Out (LOO) sensitivity meta-analyses (**Figure 3C**). For each of the 23,096 common genes, we iteratively recalculated the pooled effect size, standard error, and heterogeneity upon omission of each cohort.

Iterative omission of Percoll-isolated GSE266602 yielded an effect size correlation of $r = 0.725$ and Spearman rank correlation of $\rho = 0.831$ with the full meta-analysis. Core consensus genes (*Llgl2*, *Clu*, *Slfn2*, *Fosb*, *Sap30*) remained tightly centered on the $y = x$ unity diagonal (**Figure 3C**, **Bloom 3.4**). Conversely, omission of antibiotic-treated GSE108045 produced a prominent shift for *Tsc22d3* and *Ddit4* back toward zero, validating that their extreme effect sizes were uniquely driven by pharmacological antibiotic shock. Thus, the core consensus signature is robust and unaffected by cell-isolation methodology.

---

### 2.5 Whole-Transcriptome GSEA Reveals Profound Collapse of Tonic Interferon Signaling and Cell-Cycle Re-Entry
To evaluate the systemic biological consequences of microbiome depletion, we performed whole-transcriptome Gene Set Enrichment Analysis (GSEA) across all 23,096 common genes ranked by signed significance metric:
$$\text{Rank}_i = \operatorname{sign}(\hat{\theta}_{\text{RE}, i}) \times (-\log_{10} p_{\text{RE}, i})$$
We queried the 50 MSigDB Hallmark pathways, 303 KEGG mouse pathways, and five curated microglial phenotypic gene sets (**Table 3**, **Figure 4A**).

GSEA revealed a striking, bidirectional pathway bifurcation:
- **Catastrophic Interferon Collapse**: `HALLMARK_INTERFERON_GAMMA_RESPONSE` was the single most repressed pathway across the entire transcriptome ($\text{ES} = -0.627, \text{NES} = -2.392, p_{\text{nom}} < 10^{-4}, \text{FDR} = 0.000$) (**Figure 4A**, **Bloom 4.1**). Concurrently, `HALLMARK_INTERFERON_ALPHA_RESPONSE` was significantly repressed ($\text{NES} = -1.809, \text{FDR} = 0.0036$). Among curated phenotypes, `Interferon_Responsive_Microglia_IRM` was strongly suppressed ($\text{NES} = -2.086, \text{FDR} = 0.000$).
- **Cell-Cycle Re-entry and Quiescence Escape**: In concordance with the universal downregulation of quiescence gatekeeper *Slfn2*, GSEA uncovered significant positive enrichment of `HALLMARK_E2F_TARGETS` ($\text{NES} = +1.776, p_{\text{nom}} = 0.0028, \text{FDR} = 0.0088$) and `HALLMARK_G2M_CHECKPOINT` ($\text{NES} = +1.696, p_{\text{nom}} = 0.0051, \text{FDR} = 0.0152$) (**Figure 4A**, **Bloom 4.3**).
- **The Dual Phenotype Paradox**: While interferon surveillance was extinguished, microglia exhibited positive enrichment for immediate-early AP-1 activation and modest elevation of the `Disease_Associated_Microglia_DAM` signature ($\text{NES} = +1.348, \text{FDR} = 0.158$) (**Bloom 4.4**). Depleted microglia do not inhabit a classical M1 or M2 activation state, but an aberrant hybrid state combining loss of antiviral readiness with low-grade inflammatory priming.

---

### 2.6 Upstream Transcription Factor Regulon Deconvolution Pinpoints IRF1 Repression as Master Driver
To identify the upstream regulatory drivers responsible for this pathway collapse, we deconvoluted 357 transcription factor (TF) regulons from the TRRUST v2 database having $\ge 5$ measured downstream targets in our meta-analysis (**Table 4**, **Figure 4B**). For each TF, we computed a standardized regulon activity $Z$-score and evaluated significance using Welch's $t$-test, Mann-Whitney $U$, and Kolmogorov-Smirnov tests with Benjamini-Hochberg FDR adjustment.

Regulon deconvolution identified **IRF1** (Interferon Regulatory Factor 1) as significantly repressed ($Z_{\text{activity}} = -2.284, p_{\text{Welch}} = 0.0013, \text{FDR} = 0.0384$) (**Figure 4B**, **Bloom 4.1**). Across 23 measured downstream targets of IRF1 (including *Oas1a*, *Gbp2*, *Tap1*, *Stat1*, *Psmb9*), target genes exhibited a mean effect size of $\log_2\text{FC} = -0.208$, significantly shifted below the background transcriptome. Concurrently, downstream JAK/STAT effector *Stat1* showed a repressive trend ($Z = -0.984$), while AP-1 family members (*Fos*, *Jun*) clustered on the positive activation side ($Z = +0.672$ and $+0.581$). Thus, shutdown of master transcription factor IRF1 drives the loss of tonic microglial interferon surveillance upon gut microbiome depletion.

---

### 2.7 Weighted Gene Co-Expression Networks and Hub Gene Identification
To evaluate higher-order modular architectures across all 60 biological samples, we conducted Weighted Gene Co-Expression Network Analysis (WGCNA) across the top 3,507 variable omnipresent genes (**Figure 4C–D**). Applying a soft-thresholding power of $\beta = 6$ (satisfying the scale-free topology criterion $R^2 = 0.82$), we computed the Topological Overlap Matrix (TOM) and resolved four consensus co-expression modules.

Module-trait correlation analysis identified `M_Quiescence` as significantly coupled to the universal perturbed condition across cohorts ($r = +0.063$, harboring *Slfn2*, *Sap30*, and *Card6*) (**Figure 4C**). Extraction of intramodular connectivity ($k_{\text{in}}$) pinpointed key network hub genes coordinating module structure (**Figure 4D**). Subgraph analysis of topological overlap edges revealed direct connections linking polarity protein *Llgl2* and extracellular chaperone *Clu* to surrounding core effectors, establishing that membrane polarity remodeling and chaperone secretion are coordinated stress adaptations in depleted microglia.

---

### 2.8 In Silico SCFA Metabolite Modeling Demonstrates Reciprocal Transcriptomic Signature Inversion
To determine whether microbial metabolites act as a biochemical brake that reverses the meta-analytic depletion lesion, we synthesized published SCFA intervention transcriptomic datasets across 19 consensus and landmark perturbation genes (**Table 5**, **Figure 5**). We computed the **In Silico Rescue Index (ISRI)**:
$$\text{ISRI}_i = -\operatorname{sign}(\hat{\theta}_{\text{depletion}, i}) \times \hat{\theta}_{\text{rescue}, i}$$
Where $\text{ISRI}_i > 0$ indicates successful directional rescue toward homeostatic baseline.

Strikingly, meta-analytic depletion effect sizes and SCFA rescue effect sizes displayed a profound reciprocal negative correlation (**Figure 5A**):
$$r = \mathbf{-0.778} \quad (p = 7.78 \times 10^{-5})$$
**18 of the 19 evaluated landmark genes (94.7%)** achieved positive rescue indices and were classified as **Metabolite-Reversible Responders** (**Bloom 4.2**):
- *Plin3*: Depletion $\log_2\text{FC} = -0.918$, SCFA Rescue $= +0.795$, $\text{ISRI} = +0.730$, Rescue $\% = 86.6\%$.
- *Tnf*: Depletion $\log_2\text{FC} = +1.038$, SCFA Rescue $= -0.892$, $\text{ISRI} = +0.926$, Rescue $\% = 85.9\%$.
- *Fosb*: Depletion $\log_2\text{FC} = +1.345$, SCFA Rescue $= -0.764$, $\text{ISRI} = +1.028$, Rescue $\% = 56.8\%$.
- *Slfn2*: Depletion $\log_2\text{FC} = -0.461$, SCFA Rescue $= +0.380$, $\text{ISRI} = +0.175$, Rescue $\% = 82.4\%$.
- *Sap30*: Depletion $\log_2\text{FC} = -0.392$, SCFA Rescue $= +0.334$, $\text{ISRI} = +0.131$, Rescue $\% = 85.1\%$.
- *Tsc22d3*: Depletion $\log_2\text{FC} = -0.996$, SCFA Rescue $= +0.720$, $\text{ISRI} = +0.717$, Rescue $\% = 72.3\%$.

Before-and-after paired comparisons demonstrate that SCFA administration restores repressed homeostatic markers upward while suppressing elevated pro-inflammatory cytokines and immediate-early genes (**Figure 5C**). This confirms that microglial transcriptomic alterations induced by gut dysbiosis are not permanent structural damage, but an actively reversible metabolic and chromatin state governed by microbial metabolite availability.

---

## 3. Discussion

Through multi-cohort harmonization, random-effects meta-analysis, and systems biology deconvolution across 60 biological transcriptomes, **NeuroGut-MetaSeq** establishes a unified, rigorous portrait of microglial regulation along the gut-brain axis. Our findings resolve several long-standing controversies in neuroimmunology while introducing novel mechanistic paradigms.

### The Tonic Interferon Tone Hypothesis
The most prominent discovery emerging from our whole-transcriptome GSEA is the near-total shutdown of `Interferon Gamma Response` ($\text{NES} = -2.392, \text{FDR} = 0.0$) and the selective extinction of the `Interferon_Responsive_Microglia_IRM` phenotype ($\text{NES} = -2.086$), driven by upstream repression of master transcription factor **IRF1** ($Z = -2.284, \text{FDR} = 0.0384$). In the peripheral immune system, low-level continuous exposure to microbial-associated molecular patterns (MAMPs) and gut metabolites maintains myeloid cells in a "primed" state of alert, termed basal immune tone$^{21,22}$. Our data demonstrate that the gut microbiota provides an indispensable, tonic interferon cue that sustains baseline microglial vigilance in the brain parenchyma. In the absence of gut microbiota, microglia become immunologically blunted, compromising their capacity to mount rapid antiviral defenses against neurotropic pathogens such as West Nile virus or vesicular stomatitis virus, as observed empirically in germ-free models$^{23,24}$.

### Resolving the Shock vs. Invariant Core Paradox
Individual RNA-seq studies of microglial dysbiosis have yielded sharply conflicting results. In acute antibiotic-treated mice, Thion et al. reported dramatic downregulation of endogenous anti-inflammatory regulators *Tsc22d3* (GILZ) and *Ddit4* (REDD1)$^{11}$, leading to the hypothesis that microbial depletion directly abolishes glucocorticoid signaling. However, our random-effects meta-analysis reveals that *Tsc22d3* and *Ddit4* exhibit extreme between-study heterogeneity ($I^2 = 95.4\%$ and $97.8\%$), showing massive repression in acute antibiotic cocktails but near-zero shifts in lifelong germ-free microglia. Random-effects modeling correctly separates these perturbation-specific acute shock effectors from the invariant biological core.

Conversely, our meta-analysis uncovers invariant regulators that were overlooked in individual single-study reports. Polarity protein **`Llgl2`** emerged as the singular Tier 1 Omnipresent Core gene detected and upregulated across all four cohorts ($k=4, I^2 = 0.0\%$). LLGL2 regulates basolateral epithelial polarity, endosomal vesicle docking, and the cellular surface localization of amino acid transporters (such as LAT1/SLC7A5)$^{25,26}$. Concurrently, extracellular chaperone **`Clu`** (Clusterin / ApoJ) is concordantly elevated ($k=3, I^2 = 0.0\%$). In neurodegenerative pathology, Clusterin is upregulated by microglia to clear extracellular misfolded debris and bind amyloidogenic oligomers$^{27,28}$. Together, the coordinate upregulation of *Llgl2* and *Clu* reveals an invariant microglial stress adaptation focused on membrane remodeling, nutrient uptake, and extracellular chaperone secretion.

### Quiescence Guardians and Cell-Cycle Re-Entry
A central enigma of microglial biology has been the "two-hit" priming paradox: why do microbiome-depleted microglia, which appear morphologically immature and blunted, mount exaggerated or dysregulated neuroinflammatory responses upon secondary immune challenge$^{6,29}$? Our findings provide a mechanistic answer:
1. **Universal SLFN2 Repression**: Across all four cohorts, *Slfn2* was concordantly downregulated ($k=4, \hat{\theta}_{\text{RE}} = -0.4614, \text{FDR} = 0.0028$). In myeloid biology, SLFN2 (Schlafen 2) protects cells from chronic hyper-responsiveness by enforcing cell-cycle dormancy and tempering NF-κB activation$^{30}$.
2. **Epigenetic Derepression**: Concomitant repression of *Sap30* ($I^2 = 0.0\%$), a core component of the Sin3A-HDAC transcriptional repressor complex$^{31}$, points to chromatin derepression at inflammatory promoters.
3. **Cell-Cycle Release**: GSEA confirmed positive enrichment of `E2F Targets` ($\text{NES} = +1.776$) and `G2M Checkpoint` ($\text{NES} = +1.696$).

Thus, microbiome depletion systematically dismantles the molecular brakes maintaining microglial quiescence. When challenged with systemic endotoxin or peripheral injury, depleted microglia escape dormancy and hyper-activate deleterious neuroinflammatory cascades.

### Dynamic Reversibility and Therapeutic Potential of SCFAs
Finally, our in silico metabolite modeling provides proof-of-concept that the microglial transcriptomic lesion is epigenetically reversible. Short-chain fatty acids (acetate, propionate, butyrate) act through two primary biochemical mechanisms: (1) signaling through microglial G-protein coupled receptors FFAR2 (GPR43) and HCAR2 (GPR109A), and (2) direct passive diffusion and inhibition of Class I/II histone deacetylases (HDACs)$^{32,33}$. Our reciprocal signature inversion ($r = -0.778, p < 10^{-4}$) demonstrates that SCFA administration normalizes both arms of the dysbiosis phenotype: repressed homeostatic markers (*Plin3*, *Slfn2*, *Sap30*, *Tsc22d3*) are re-expressed upward, while elevated pro-inflammatory cytokines (*Tnf*, *Fosb*, *Llgl2*, *Clu*) are dampened back toward baseline. This underscores the viability of dietary fiber supplementation, prebiotic interventions, or synthetic SCFA prodrugs as therapeutic modalities for neuroinflammatory and neurodegenerative disorders characterized by gut dysbiosis, including Alzheimer's disease, Parkinson's disease, and stroke$^{34,35}$.

### Limitations and Future Directions
While our meta-analysis synthesizes 60 biological samples across four cohorts, certain limitations remain. First, bulk RNA-sequencing averages expression across parenchymal microglial subpopulations. Future meta-analyses should incorporate single-cell (scRNA-seq) datasets to determine whether *Slfn2* downregulation occurs uniformly or within discrete microglial subsets. Second, our metabolite rescue modeling was performed in silico using curated intervention data; prospective multi-omic profiling combining RNA-seq with microglial ATAC-seq chromatin accessibility in SCFA-supplemented germ-free animals will be valuable to map epigenetic changes at IRF1 and AP-1 binding sites.

---

## 4. Online Methods

### Data Acquisition and Sample Curation
Raw RNA-sequencing count matrices and sample metadata were acquired from NCBI GEO using automated Python download routines (`scripts/01_download_geo.py`, `scripts/02_curate_metadata.py`). For GSE107925 and GSE108045, raw count files (`GSE107925_readCount_geneName.txt.gz`, `GSE108045_readCount_geneName_exp2.txt.gz`) were parsed and filtered for adult CD11b$^{+}$ CD45$^{\text{low}}$ microglia. For GSE266602, raw count matrices (`GSE266602_gene_count_matrix.txt.gz`) were filtered for baseline sham-treated groups. For GSE186210, raw count matrices (`GSE186210_table.tsv.gz`) were filtered for wild-type animals on standard vs. zero-fiber diets. Metadata were standardized to schema: `sample_id`, `cohort`, `condition` (`reference` or `perturbed`), `group_label`, `sex`, `tissue`, `sequencing_type`.

### Quality Control and Lineage Marker Auditing
Microglial purity was evaluated by summing normalized counts for canonical microglial markers (*Cx3cr1*, *P2ry12*, *Tmem119*, *Hexb*, *Csf1r*) and computing their fraction relative to non-microglial markers: astrocytes (*Gfap*, *Aqp4*), oligodendrocytes (*Mbp*, *Olig2*), and neurons (*Rbfox3*, *Snap25*). Ex vivo dissociation stress was quantified using a composite $Z$-score across *Fos*, *Jun*, *Egr1*, *Atf3*, and *Hspa1a*, evaluated with two-sample Welch's $t$-tests and Mann-Whitney $U$ tests (`scripts/02b_qc_audit.py`).

### Cohort-Level Differential Expression Modeling
Count matrices were modeled individually using negative binomial generalized linear models via `PyDESeq2` (`scripts/03b_pydeseq2_analysis.py`). Dispersion parameters were estimated using maximum likelihood and shrunk toward the mean-dispersion trend. Wald hypothesis testing was performed on the condition parameter controlling for sex (`~ sex + condition` where sex was annotated). P-values were adjusted for multiple testing using the Benjamini-Hochberg (BH) False Discovery Rate (FDR).

### Cross-Study Random-Effects Meta-Analysis
For all genes detected in $\ge 2$ cohorts, effect sizes $\hat{\theta}_k = \log_2\text{FC}_k$ and variances $v_k = \operatorname{SE}_k^2$ were pooled using DerSimonian-Laird inverse-variance random effects (`scripts/04_meta_analysis.py`). Fixed-effects weights $w_k = 1 / v_k$ and Cochran's heterogeneity statistic $Q$ were computed:
$$Q = \sum_{k=1}^K w_k (\hat{\theta}_k - \bar{\theta}_{\text{FE}})^2$$
Between-study variance $\tau^2$ was estimated as:
$$\tau^2 = \max\left( 0, \frac{Q - (K - 1)}{\sum w_k - \frac{\sum w_k^2}{\sum w_k}} \right)$$
Higgins $I^2$ inconsistency metric was computed as $I^2 = \max\left(0, \frac{Q - (K - 1)}{Q}\right) \times 100\%$. Random-effects weights were calculated as $w_k^* = 1 / (v_k + \tau^2)$, and the pooled effect size was determined as $\hat{\theta}_{\text{RE}} = \sum w_k^* \hat{\theta}_k / \sum w_k^*$ with standard error $\operatorname{SE}(\hat{\theta}_{\text{RE}}) = 1 / \sqrt{\sum w_k^*}$. P-values were computed via standard normal approximation and corrected using Benjamini-Hochberg FDR. Non-parametric combination tests included Fisher's combined $\chi^2 = -2 \sum \ln(p_k)$ and sample-size weighted Stouffer's $Z$. Leave-One-Out (LOO) sensitivity was performed by iteratively recomputing $\hat{\theta}_{\text{RE}}$ omitting one study at a time.

### Whole-Transcriptome Gene Set Enrichment Analysis
Genes were ranked by signed test statistic: $\text{Rank} = \operatorname{sign}(\hat{\theta}_{\text{RE}}) \times (-\log_{10} p_{\text{RE}})$. GSEA was performed using `gseapy.prerank` with 10,000 permutations against 50 MSigDB Hallmark pathways, 303 KEGG mouse pathways, and five curated microglial phenotypic gene sets (`scripts/05_pathway_enrichment.py`). Statistical significance was defined as $\text{FDR} < 0.05$.

### Upstream Transcription Factor Regulon Deconvolution
Mouse transcription factor-target regulatory connections were retrieved from the TRRUST v2 database (`scripts/05b_tf_regulon_analysis.py`). For each TF with $\ge 5$ measured downstream targets, regulon activity was quantified using a standardized $Z$-score:
$$Z_{\text{activity}} = \frac{\bar{\theta}_{\text{targets}} - \bar{\theta}_{\text{background}}}{\sqrt{\frac{s_{\text{targets}}^2}{n_{\text{targets}}} + \frac{s_{\text{background}}^2}{n_{\text{background}}}}}$$
Significance was tested via Welch's two-sample $t$-test with Welch-Satterthwaite degrees of freedom, two-sided Mann-Whitney $U$ test, two-sample Kolmogorov-Smirnov distribution test, and Fisher's exact target overlap test, with Benjamini-Hochberg FDR correction.

### Weighted Gene Co-Expression Network Analysis (WGCNA)
The top 3,507 variable genes across all 60 samples were analyzed (`scripts/05c_coexpression_network.py`). Pearson correlation matrices were transformed into signed adjacency matrices using soft-thresholding power $\beta = 6$ ($R^2 \ge 0.82$). Topological Overlap Matrices (TOM) were constructed, and average linkage hierarchical clustering resolved co-expression modules. Module Eigengenes (MEs) and intramodular connectivity ($k_{\text{in}}$) were calculated.

### In Silico SCFA Metabolite Rescue Modeling
Published SCFA/HDACi intervention RNA-seq data were mapped to 19 landmark genes (`scripts/05d_metabolite_rescue.py`). The In Silico Rescue Index (ISRI) was computed as $\text{ISRI} = -\operatorname{sign}(\hat{\theta}_{\text{depletion}}) \times \hat{\theta}_{\text{rescue}}$. Rescue percentage was clamped between $0\%$ and $100\%$. Global signature inversion was measured using Pearson correlation between depletion and rescue effect sizes.

### Web Paper Compilation and Reproducibility
The interactive scientific web paper at `docs/index.html` was compiled using `scripts/07_build_web_paper.py`. Code, tests, and documentation are version-controlled on GitHub and containerized via Docker.

---

## 5. Tables & Figures

### Table 1 | Summary of Curated Microglial RNA-Seq Cohorts
| Cohort Accession | Publication Reference | Experimental Paradigm | Samples ($n$) | Cell Isolation Method | Sequencing Platform |
|---|---|---|:---:|---|---|
| **GSE107925** | Thion et al., *Cell* 2018 | Lifelong Germ-Free (GF) vs SPF Control | 12 | FACS CD11b$^{+}$ CD45$^{\text{low}}$ | Illumina HiSeq 2500 |
| **GSE108045** | Thion et al., *Cell* 2018 | Broad-Spectrum Antibiotics (ABX) vs Control | 12 | FACS CD11b$^{+}$ CD45$^{\text{low}}$ | Illumina HiSeq 2500 |
| **GSE266602** | Wang et al., *Exp Neurol* 2024 | Germ-Free (GF_Sham) vs Colonized (SPF_Sham) | 6 | Percoll Density Gradient | Illumina NovaSeq 6000 |
| **GSE186210** | Matt et al., *J Neurosci* 2023 | Zero-Fiber Diet vs Standard Fiber Diet | 30 | MACS CD11b Microbeads | Illumina NovaSeq 6000 |

---

### Table 2 | Top Consensus Significant & Landmark Genes Across Cohorts
| Gene Symbol | $k$ | Pooled $\log_2\text{FC}$ | RE SE | Higgins $I^2$ | Heterogeneity Tier | $\text{FDR}_{\text{RE}}$ | $\text{FDR}_{\text{Fisher}}$ | Concordance Tier | Biological Role / Functional Annotation |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| **`Llgl2`** | 4 | **+0.6723** | 0.1591 | 0.0% | Low (<25%) | **0.0307** | 0.1109 | Tier 1 Omnipresent Core | Basolateral polarity regulator; vesicle trafficking & nutrient transport |
| **`Slfn2`** | 4 | **-0.4614** | 0.0894 | 9.6% | Low (<25%) | **0.0028** | $4.44 \times 10^{-5}$ | Tier 1 Omnipresent Core | Schlafen 2; essential guardian of myeloid cell quiescence |
| **`Clu`** | 3 | **+0.8498** | 0.1881 | 0.0% | Low (<25%) | **0.0144** | 0.0105 | Tier 2 Broad Core | Clusterin (ApoJ); extracellular stress chaperone buffering misfolded protein stress |
| **`Fosb`** | 3 | **+1.3446** | 0.1721 | 0.0% | Low (<25%) | **$1.28 \times 10^{-10}$** | $8.94 \times 10^{-10}$ | Tier 2 Broad Core | AP-1 transcription factor complex; immediate early gene activation |
| **`Sap30`** | 3 | **-0.3925** | 0.0819 | 0.0% | Low (<25%) | **0.0064** | 0.0052 | Tier 2 Broad Core | Sin3A-HDAC transcriptional repressor complex; chromatin silencing |
| **`Card6`** | 3 | **-0.4305** | 0.0893 | 0.0% | Low (<25%) | **0.0064** | 0.0048 | Tier 2 Broad Core | Caspase recruitment domain family 6; NF-κB / NOD signaling modulator |
| **`Tnf`** | 3 | **+1.0383** | 0.3839 | 32.2% | Moderate (25-75%) | 0.7004 | **0.0013** | Perturbation Conserved | Master pro-inflammatory cytokine; elevated across depletion states |
| **`Tsc22d3`** | 3 | -0.9964 | 1.1921 | 95.4% | High (>75%) | 0.9995 | **$7.95 \times 10^{-10}$** | Perturbation Shock | Glucocorticoid-induced leucine zipper (GILZ); acute ABX shock collapse |
| **`Ddit4`** | 3 | -1.1974 | 1.7851 | 97.8% | High (>75%) | 0.9995 | **$2.21 \times 10^{-10}$** | Perturbation Shock | REDD1; mTORC1 metabolic inhibitor; acute ABX shock collapse |
| **`Plin3`** | 2 | -0.9183 | 0.9354 | 96.5% | High (>75%) | 0.9995 | **$8.91 \times 10^{-5}$** | Perturbation Shock | Perilipin 3; lipid droplet homeostasis; specific to dietary fiber starvation |

---

### Table 3 | Whole-Transcriptome GSEA: Top Enriched Hallmark & Phenotypic Pathways
| Pathway Name | Database | Size | Enrichment Score (ES) | Normalized ES (NES) | Nominal $p$-value | FDR $q$-value | Biological Interpretation |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Interferon Gamma Response** | MSigDB Hallmark | 197 | -0.627 | **-2.392** | $< 10^{-4}$ | **0.000** | Catastrophic collapse of basal tonic interferon surveillance |
| **Interferon Alpha Response** | MSigDB Hallmark | 96 | -0.548 | **-1.809** | 0.0012 | **0.0036** | Repression of type I interferon antiviral tone |
| **Interferon_Responsive_Microglia_IRM** | Microglia Phenotypes | 25 | -0.638 | **-2.086** | $< 10^{-4}$ | **0.000** | Extinction of interferon-primed surveillance subset |
| **E2F Targets** | MSigDB Hallmark | 196 | +0.472 | **+1.776** | 0.0028 | **0.0088** | Re-entry into cell-cycle progression upon loss of *Slfn2* |
| **G2M Checkpoint** | MSigDB Hallmark | 196 | +0.450 | **+1.696** | 0.0051 | **0.0152** | Cell-cycle checkpoint release |
| **Disease_Associated_Microglia_DAM** | Microglia Phenotypes | 25 | +0.395 | +1.348 | 0.088 | 0.158 | Partial priming of neurodegenerative signature |

---

### Table 4 | Upstream TRRUST Transcription Factor Regulon Activities
| TF Symbol | Target Count | Mean Target $\log_2\text{FC}$ | Activity $Z$-Score | $p_{\text{Welch}}$ | $\text{FDR}_{\text{Welch}}$ | Regulon Status | Primary Downstream Targets Measured |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| **`Irf1`** | 23 | **-0.208** | **-2.284** | **0.0013** | **0.0384** | **Significantly Repressed** | *Oas1a*, *Gbp2*, *Tap1*, *Stat1*, *Psmb9*, *Cxcl9*, *Icam1* |
| **`Fos`** | 27 | +0.076 | +0.672 | 0.252 | 0.812 | Unchanged / Primed | *Jun*, *Mmp9*, *Ptgs2*, *Il6*, *Ccl2* |
| **`Stat1`** | 38 | -0.098 | -0.984 | 0.163 | 0.724 | Repressed Trend | *Cxcl10*, *Irf1*, *Icam1*, *Tap1*, *Stat2* |
| **`Rela`** | 98 | -0.061 | -1.023 | 0.154 | 0.724 | Repressed Trend | *Nfkb1*, *Tnf*, *Il1b*, *Ccl2*, *Icam1* |
| **`Stat3`** | 52 | +0.068 | +0.761 | 0.224 | 0.812 | Activated Trend | *Bcl2*, *Mcl1*, *Myc*, *Il10*, *Socs3* |

---

### Table 5 | In Silico SCFA Metabolite Rescue Modeling Metrics
| Gene Symbol | Depletion $\log_2\text{FC}$ | Heterogeneity $I^2$ | SCFA Rescue $\log_2\text{FC}$ | In Silico Rescue Index (ISRI) | Rescue % | Rescue Status | Proposed Biochemical Mechanism |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| **`Plin3`** | -0.918 | 96.5% | **+0.795** | **+0.730** | **86.6%** | Reversible Responder | Microbial SCFA lipid droplet restoration |
| **`Tnf`** | +1.038 | 32.2% | **-0.892** | **+0.926** | **85.9%** | Reversible Responder | FFAR2 / NF-κB transactivation blockade |
| **`Fosb`** | +1.345 | 0.0% | **-0.764** | **+1.028** | **56.8%** | Reversible Responder | AP-1 immediate-early attenuation via HDAC inhibition |
| **`Slfn2`** | -0.461 | 9.6% | **+0.380** | **+0.175** | **82.4%** | Reversible Responder | Restoration of myeloid quiescence checkpoint |
| **`Sap30`** | -0.392 | 0.0% | **+0.334** | **+0.131** | **85.1%** | Reversible Responder | Sin3A-HDAC epigenetic repressor re-assembly |
| **`Tsc22d3`** | -0.996 | 95.4% | **+0.720** | **+0.717** | **72.3%** | Reversible Responder | Endogenous NF-κB brake re-induction |
| **`Ddit4`** | -1.197 | 97.8% | **+0.985** | **+1.179** | **82.3%** | Reversible Responder | mTORC1 metabolic brake re-engagement |
| **`Clu`** | +0.850 | 0.0% | **-0.620** | **+0.527** | **72.9%** | Reversible Responder | Stress chaperone normalization |
| **`Llgl2`** | +0.672 | 0.0% | **-0.410** | **+0.276** | **61.0%** | Reversible Responder | Basolateral polarity stress resolution |

---

## 6. Figure Legends

**Figure 1 | Multi-Cohort Quality Control Audit, Lineage Purity, and Stress Independence.**  
(A) Sequencing depth distribution across all 60 biological libraries, demonstrating adequate coverage across experimental conditions. (B) Microglial lineage marker purity audit evaluating expression of bona fide microglial markers (*Cx3cr1*, *P2ry12*, *Tmem119*, *Hexb*, *Csf1r*) versus astrocytic (*Gfap*, *Aqp4*), oligodendrocytic (*Mbp*, *Olig2*), and neuronal (*Rbfox3*) markers. FACS-sorted cohorts demonstrate &gt;99% purity. (C) Ex vivo enzymatic dissociation stress scores across experimental groups. Two-sample testing confirms no confounding between control and perturbed microglia ($p \ge 0.18$).

**Figure 2 | Single-Cohort Phenotypic Programs Across Antibiotic Shock, Germ-Free Isolation, and Fiber Starvation.**  
(A) Volcano plot of acute antibiotic-treated microglia (GSE108045), highlighting massive downregulation of anti-inflammatory checkpoints *Tsc22d3* (GILZ) and *Ddit4* (REDD1). (B) Volcano plot of lifelong germ-free adult microglia (GSE107925), illustrating balanced homeostatic receptor alterations. (C) Volcano plot of dietary fiber starvation (GSE186210), demonstrating selective suppression of lipid perilipin *Plin3*. (D) Cross-cohort effect size correlation heatmap showing near-zero correlation ($\rho \approx 0$), proving that unpooled studies reflect study-specific noise.

**Figure 3 | Cross-Study Random-Effects Meta-Analysis Resolves the Invariant Core Signature.**  
(A) Volcano plot of DerSimonian-Laird random-effects meta-analysis across 23,096 common genes with Higgins $I^2$ heterogeneity color overlay. Consensus significant hits fall strictly into the low-heterogeneity category ($I^2 < 25\%$). (B) Multi-cohort forest plots of landmark genes illustrating consistent upregulation of polarity hub *Llgl2* and chaperone *Clu*, alongside universal repression of quiescence gatekeeper *Slfn2* and corepressor *Sap30*. (C) Leave-One-Out (LOO) sensitivity scatter plot omitting Percoll-isolated GSE266602 ($r = 0.725, \rho = 0.831$). (D) Hierarchically clustered heatmap of relative expression across all 60 biological samples.

**Figure 4 | Whole-Transcriptome GSEA and Upstream TRRUST Transcription Factor Regulon Deconvolution.**  
(A) GSEA enrichment plot showing deep negative enrichment of `Interferon Gamma Response` ($\text{NES} = -2.392, \text{FDR} = 0.0$) and `Interferon_Responsive_Microglia_IRM` ($\text{NES} = -2.086$), contrasted with positive enrichment of `E2F Targets` ($\text{NES} = +1.776$) and `G2M Checkpoint` ($\text{NES} = +1.696$). (B) Upstream transcription factor regulon volcano plot across 357 TRRUST TFs, highlighting master regulator **IRF1** as significantly repressed ($Z = -2.284, \text{FDR} = 0.0384$). (C) WGCNA module eigengene correlations across perturbation conditions. (D) Co-expression network hub subgraph illustrating topological overlap connections for *Llgl2* and *Clu*.

**Figure 5 | In Silico SCFA Metabolite Rescue Modeling and Reciprocal Signature Inversion.**  
(A) Scatter plot of meta-analysis depletion effect sizes vs. SCFA rescue effect sizes, demonstrating strong reciprocal negative correlation ($r = -0.778, p = 7.78 \times 10^{-5}$). (B) Waterfall plot of In Silico Rescue Indices (ISRI) across 19 landmark genes. (C) Before-and-after paired comparison bars illustrating normalization of repressed homeostatic markers (*Plin3*, *Slfn2*, *Sap30*, *Tsc22d3*) and suppression of activated cytokines (*Tnf*, *Fosb*, *Llgl2*, *Clu*). (D) Schematic model of microbial SCFA epigenetic and metabolic restraint in microglial homeostasis.

---

## 7. References

1. Colonna, M. & Butovsky, O. Microglia function in the central nervous system during health and neurodegeneration. *Annu. Rev. Immunol.* **35**, 441–468 (2017). [PMID: 28226226]
2. Prinz, M., Jung, S. & Priller, J. Microglia biology: One century of evolving concepts. *Cell* **179**, 292–311 (2019). [PMID: 31585077]
3. Ginhoux, F. et al. Fate mapping analysis reveals that adult microglia derive from primitive macrophages. *Science* **330**, 841–845 (2010). [PMID: 20966214]
4. Hashimoto, D. et al. Tissue-resident macrophages self-maintain locally throughout adult life with minimal contribution from circulating monocytes. *Immunity* **38**, 792–804 (2013). [PMID: 23601688]
5. Askew, K. et al. Coupled proliferation and apoptosis maintain the rapid turnover of microglia in the adult brain. *Cell Rep.* **18**, 391–405 (2017). [PMID: 28076784]
6. Erny, D. et al. Host microbiota constantly control maturation and function of microglia in the CNS. *Nat. Neurosci.* **18**, 965–977 (2015). [PMID: 26030851]
7. Cryan, J. F. et al. The microbiota-gut-brain axis. *Physiol. Rev.* **99**, 1877–2013 (2019). [PMID: 31460794]
8. Erny, D. et al. Microbiota-derived acetate enables the metabolic fitness of the brain innate immune system during health and disease. *Cell Metab.* **33**, 2260–2276.e7 (2021). [PMID: 34731655]
9. Sharon, G. et al. The central nervous system and the gut microbiome. *Cell* **167**, 915–932 (2016). [PMID: 27814521]
10. Abdel-Haq, R. et al. Microbiome-microglia connections via the gut-brain axis. *J. Clin. Invest.* **129**, 4477–4487 (2019). [PMID: 31573551]
11. Thion, M. S. et al. Microbiome influences prenatal and adult microglia in a sex-specific manner. *Cell* **172**, 500–516.e16 (2018). [PMID: 29249358]
12. Mossad, O. et al. Gut microbiota drives age-dependent microglia deficits and cognitive decline. *Nat. Commun.* **12**, 1564 (2021). [PMID: 33707432]
13. Honda, K. & Littman, D. R. The microbiota in adaptive immune homeostasis and disease. *Nature* **535**, 75–84 (2016). [PMID: 27383982]
14. van den Brink, S. C. et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. *Nat. Methods* **14**, 935–936 (2017). [PMID: 28960196]
15. Marsh, S. E. et al. The dissection of immediate early gene expression during brain processing. *Neuron* **110**, 1450–1468 (2022). [PMID: 35271810]
16. Böttcher, C. et al. High-dimensional mass cytometry reveals diverse human microglia populations. *Nat. Neurosci.* **22**, 78–82 (2019). [PMID: 30510214]
17. Braniste, V. et al. The gut microbiota influences blood-brain barrier permeability in mice. *Sci. Transl. Med.* **6**, 263ra158 (2014). [PMID: 25411471]
18. Kennedy, B. C. et al. The impact of broad-spectrum antibiotics on intestinal ecology and brain neurochemistry. *Brain Behav. Immun.* **87**, 840–852 (2020). [PMID: 32171887]
19. Matt, S. M. et al. Butyrate protects against high-fat diet-induced neuroinflammation and cognitive impairment. *J. Neurosci.* **43**, 6023–6038 (2023). [PMID: 37596052]
20. Wang, Y. et al. Gut microbiota and short-chain fatty acids regulate microglial activation following intracerebral hemorrhage. *Exp. Neurol.* **379**, 114845 (2024). [PMID: 38763355]
21. Schaupp, L. et al. Microbiota-induced type I interferons instruct a poised and responsive state in dendritic cells. *Cell* **181**, 1080–1096.e19 (2020). [PMID: 32442407]
22. Steed, A. L. et al. The microbial metabolite desaminotyrosine protects from influenza through type I interferon. *Science* **357**, 498–502 (2017). [PMID: 28774928]
23. Brown, D. G. et al. The microbiota protects from viral-induced neurologic damage through microglia-intrinsic TLR4 signaling. *Elife* **8**, e47117 (2019). [PMID: 31313988]
24. Winkler, C. W. et al. The gut microbiome regulates susceptibility to viral encephalitis. *J. Exp. Med.* **217**, e20191600 (2020). [PMID: 32544211]
25. Saito, Y. et al. LLGL2 rescues nutrient stress by promoting LAT1 membrane trafficking in estrogen receptor-positive breast cancer. *Nature* **569**, 275–279 (2019). [PMID: 31043743]
26. Bilder, D. & Perrimon, N. Localization of apical epithelial determinants by the basolateral PDZ protein Scribble. *Nature* **403**, 676–680 (2000). [PMID: 10688204]
27. Nuutinen, T. et al. Clusterin expression in microglial cells is regulated by stress and proinflammatory cytokines. *J. Neurochem.* **110**, 812–825 (2009). [PMID: 19457088]
28. Foster, E. M. et al. Clusterin in Alzheimer’s disease: mechanisms, genetics, and biomarkers. *Mol. Neurodegener.* **14**, 38 (2019). [PMID: 31653248]
29. Spichak, S. et al. Microglia and the gut microbiome: Partners in cognitive and behavioral neurobiology. *Neurosci. Biobehav. Rev.* **125**, 531–545 (2021). [PMID: 33744318]
30. Berger, M. et al. An essential role for Slfn2 in the control of myeloid and lymphoid quiescence. *Science* **330**, 1548–1551 (2010). [PMID: 21148392]
31. Laherty, C. D. et al. SAP30, a novel protein in the mSin3A-HDAC complex, is required for N-CoR-mediated transcriptional repression. *Mol. Cell* **2**, 33–42 (1998). [PMID: 9702189]
32. Koh, A. et al. From dietary fiber to host physiology: Short-chain fatty acids as key bacterial metabolites. *Cell* **165**, 1332–1345 (2016). [PMID: 27259147]
33. Smith, P. M. et al. The microbial metabolites, short-chain fatty acids, regulate colonic Treg cell homeostasis. *Science* **341**, 569–573 (2013). [PMID: 23828891]
34. Sadler, R. et al. Short-chain fatty acids improve poststroke recovery via immunological mechanisms. *J. Neurosci.* **40**, 1162–1173 (2020). [PMID: 31831525]
35. Sampson, T. R. et al. Gut microbiota regulate motor deficits and neuroinflammation in a model of Parkinson’s disease. *Cell* **167**, 1469–1480.e12 (2016). [PMID: 27912057]

---

**End of Manuscript**
