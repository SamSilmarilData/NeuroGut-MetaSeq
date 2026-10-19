# bioRxiv Submission Package Dossier: NeuroGut-MetaSeq (v1.3.0)

This dossier contains the complete metadata, author information, compliance statements, and upload files required for submission to [bioRxiv](https://submit.biorxiv.org/).

---

## 1. Article Metadata

- **Article Title**: Cross-Study Meta-Analysis of the Gut-Microbiota-Microglia Axis Uncovers Cell-Intrinsic Interferon Shutoff, Invariant Nutrient-Sensing Adapters, and Multi-Omic Reversibility
- **Short Title / Running Head**: Multi-Omic Meta-Analysis of the Gut-Microglia Axis
- **Primary Subject Category**: **Bioinformatics** (Recommended)
- **Secondary Subject Category**: **Neuroscience**
- **Manuscript Type**: New Results / Regular Article
- **License Selected**: Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) or CC BY 4.0

---

## 2. Authors and Affiliations

1. **Samyak Meshram** (Corresponding Author)
   - **Affiliation**: Independent Researcher, Open-Science Initiative
   - **Email**: `sam.meshram18@gmail.com`
   - **Role**: Conceived and designed the study, developed the analytical meta-analysis pipelines, performed multi-omic integrations, engineered the software, drafted the manuscript.
2. **Dr. Soumya Dhokey**
   - **Affiliation**: Independent Researcher, Open-Science Initiative
   - **Role**: Conceived and designed the study, evaluated the biological findings, neuroimmunological models, and clinical interpretations, refined and edited the manuscript.

---

## 3. bioRxiv-Compliant Submission Abstract

*Conforms to bioRxiv $\le 2,500$ characters including spaces requirement (Words: 236, Characters: 2,022):*

> Microglia, the resident immune sentinels of the central nervous system, are continuously calibrated by biochemical cues from the gut microbiota. However, transcriptomic findings across laboratories have suffered from small sample sizes, isolation artifacts, and disparate models. Here, we present NeuroGut-MetaSeq (v1.3.0), synthesizing 60 biological transcriptomes across four independent rodent cohorts (lifelong germ-free, acute antibiotics, dietary fiber starvation) using Restricted Maximum Likelihood (REML) variance estimation and Hartung-Knapp-Sidik-Jonkman (HKSJ) random-effects pooling across 23,096 common genes.
>
> We separate an invariant multi-study core led by basolateral polarity adapter Llgl2 (k=4, theta=+0.672, 95% CI [0.166, 1.178], I^2=0%) and extracellular chaperone Clu (k=3, theta=+0.873, I^2=0%), alongside universal loss of myeloid quiescence via Slfn2 (theta=-0.465, p=0.0197) and corepressor Sap30 (I^2=0%). Acute antibiotic shock markers (Tsc22d3, Ddit4) exhibit extreme heterogeneity (I^2 > 95%), identifying them as model-private artifacts. Factorial sex-interaction modeling (N=51) reveals that 99.1% of microglial dysbiosis is sex-invariant, with selective male-biased vulnerability in Slfn2. Whole-transcriptome GSEA and TRRUST deconvolution pinpoint master regulator IRF1 (Z=-2.453, FDR=0.0380) driving tonic interferon surveillance collapse. BayesPrism deconvolution confirms interferon-responsive microglia (IRM) are physically retained (17.5% vs 16.8%), proving an authentic cell-intrinsic transcriptional shutoff validated against stereological cell counts. Tripartite ATAC-seq footprinting demonstrates 88.9% open-chromatin footprint restoration at Irf1 and Stat1 motifs upon short-chain fatty acid (SCFA) repletion. Together, these multi-omic findings prove that gut microbiome depletion induces an aberrant hybrid state of blunted antiviral surveillance and metabolic scavenging characterized by chromatin-poised transcriptional reversibility upon microbial metabolite restoration.

---

## 4. Mandatory bioRxiv Compliance Statements

### Competing Interest Statement
> "The authors declare no competing interests."

### Data Availability Statement
> "All analyzed raw datasets are publicly available in the NCBI Gene Expression Omnibus (GEO) under accession numbers **GSE107925**, **GSE108045**, **GSE186210**, **GSE266602**, **GSE152865**, and **GSE67858**. All curated sample metadata, normalized read count matrices, and meta-analytic summary tables are deposited and accessible in the study repository at https://github.com/samyakmeshram/NeuroGut-MetaSeq."

### Code Availability Statement
> "All analytical pipelines, statistical modeling scripts, and visualization code are publicly available under the open-source MIT license at https://github.com/samyakmeshram/NeuroGut-MetaSeq. The workflows are fully reproducible using Python 3.9+ and R 4.3+."

---

## 5. Submission Artifacts
- **Primary Manuscript PDF**: `paper/manuscript.pdf` (Single-column layout, line-numbered via `lineno`, 8 `booktabs` tables, 11 figures with legends).
- **Supplementary Information PDF**: `paper/supplementary.tex` -> `paper/supplementary.pdf` (Supplementary Figures S1--S6, Supplementary Tables S1--S8, Supplementary Methods).
- **Permanent Archive (Zenodo DOI)**: Release tag `v1.3.0-preprint` linked via Zenodo.
