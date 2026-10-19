# bioRxiv Submission Pre-Flight Checklist

Use this checklist prior to and during submission to the [bioRxiv Submission Portal](https://submit.biorxiv.org/).

---

## Pre-Flight Verification Checklist

### 1. Document & File Preparation
- [x] **Primary Manuscript PDF** (`paper/manuscript.pdf`):
  - Formatted in clean single-column preprint format.
  - Line numbers activated continuously (`\linenumbers`) for referee citation.
  - All 11 figures rendered at high resolution (300+ DPI) with clear subpanel lettering (A, B, C...).
  - All 8 tables typeset with `booktabs` (`\toprule`, `\midrule`, `\bottomrule`) without overflowing margins.
  - Citations ordered sequentially with superscript numbers matching `references.bib`.
- [x] **Supplementary Information PDF** (`paper/supplementary.pdf`):
  - Supplementary Figures S1--S6 and Supplementary Tables S1--S8.
  - Independent page numbering (S1, S2...).
- [x] **Plain-Text Abstract** (`biorxiv_submission/ABSTRACT_SUBMISSION.txt`):
  - Length: 236 words, 2,022 characters with spaces (strictly below bioRxiv's 2,500-character ceiling).

### 2. Portal Form Fields Copy-Paste Guide
- [ ] **Manuscript Title**:
  ```
  Cross-Study Meta-Analysis of the Gut-Microbiota-Microglia Axis Uncovers Cell-Intrinsic Interferon Shutoff, Invariant Nutrient-Sensing Adapters, and Multi-Omic Reversibility
  ```
- [ ] **Subject Area**: Select **Bioinformatics** (Primary) or **Neuroscience** (Secondary).
- [ ] **Authors & Order**:
  1. Samyak Meshram (Corresponding: `sam.meshram18@gmail.com`, Affiliation: Independent Researcher, Open-Science Initiative)
  2. Soumya Dhokey (Co-author, Affiliation: Independent Researcher, Open-Science Initiative)
- [ ] **Abstract**: Copy and paste exactly from `biorxiv_submission/ABSTRACT_SUBMISSION.txt`.
- [ ] **Competing Interest Statement**: Select "No" / Enter:
  ```
  The authors declare no competing interests.
  ```
- [ ] **Data Availability Statement**: Enter:
  ```
  All analyzed raw datasets are publicly available in the NCBI Gene Expression Omnibus (GEO) under accession numbers GSE107925, GSE108045, GSE186210, GSE266602, GSE152865, and GSE64977. All curated metadata, normalized count matrices, and meta-analytic summary tables are deposited in the study GitHub repository (https://github.com/samyakmeshram/NeuroGut-MetaSeq).
  ```
- [ ] **License**: Select **CC BY-NC-ND 4.0** (standard for preprints) or **CC BY 4.0** (open-access re-use).

### 3. File Upload Step
- [ ] Upload `paper/manuscript.pdf` as **Manuscript PDF** (or Main Article).
- [ ] Upload `paper/supplementary.pdf` as **Supplementary File** / **Data Supplement**.

### 4. Post-Conversion Inspection
- [ ] Review bioRxiv auto-generated consolidated PDF proof.
- [ ] Verify that figures and equations render cleanly in the web previewer.
- [ ] Confirm corresponding author receives the submission confirmation email.
