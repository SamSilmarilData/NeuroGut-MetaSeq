# Zenodo Code & Data Freeze Guide: Tag `v1.4.0`

bioRxiv strongly recommends linking preprint submissions to a citable, permanently archived version of all code and processed data on [Zenodo](https://zenodo.org/). Follow these steps to execute the freeze and obtain a citable DOI.

---

## Step 1: Create and Push Git Release Tag

From the project repository root:

```bash
# 1. Verify working directory is clean
git status

# 2. Create annotated tag
git tag -a v1.4.0 -m "NeuroGut-MetaSeq v1.4.0: Production Academic Release & bioRxiv Preprint Code/Data Freeze"

# 3. Push tag to GitHub
git push origin v1.4.0
```

---

## Step 2: Connect GitHub Repository to Zenodo

1. Navigate to [https://zenodo.org/account/settings/github/](https://zenodo.org/account/settings/github/).
2. Log in using your GitHub account (`samyakmeshram`).
3. Locate `samyakmeshram/NeuroGut-MetaSeq` in the list of repositories.
4. Toggle the switch to **ON** (Enabled).

---

## Step 3: Trigger the Archive & Mint the DOI

1. Go to GitHub: `https://github.com/samyakmeshram/NeuroGut-MetaSeq/releases/new`.
2. Select tag: **`v1.4.0`**.
3. Release Title: `NeuroGut-MetaSeq v1.4.0: Production Academic Release`.
4. Release Notes:
   ```markdown
   Official release freeze for the bioRxiv preprint:
   "Cross-Study Meta-Analysis of the Gut-Microbiota-Microglia Axis Uncovers Cell-Intrinsic Interferon Shutoff, Invariant Nutrient-Sensing Adapters, and Multi-Omic Reversibility"

   Includes:
   - Consolidated 6 Biological Acts results architecture and dedicated full-page display plates
   - Master figure engine (scripts/build_unified_figures.py) with edge-to-edge panoramic schematics (Figs 1A, 6A) and Fig 2D Consensus Concordance Matrix
   - Complete multi-cohort count matrices and curated clinical metadata (60 samples, 4 cohorts)
   - PyDESeq2, REML-HKSJ random-effects meta-analysis pipelines
   - Factorial sex-dimorphism interaction models
   - BayesPrism single-cell deconvolution and cell-intrinsic imputation
   - Microglial ATAC-seq TOBIAS footprinting and NicheNet intercellular signaling
   - Publication-grade LaTeX source and compiled manuscript (21 pages) & supplementary (13 pages) PDFs
   ```
5. Click **Publish Release**.
6. Within 1--2 minutes, Zenodo will automatically archive the release tarball and assign a permanent DOI (e.g., `10.5281/zenodo.XXXXXXX`).

---

## Step 4: Add DOI Badge to README and Paper

Once assigned:
1. Update `README.md` with:
   ```markdown
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.<YOUR-DOI>.svg)](https://doi.org/10.5281/zenodo.<YOUR-DOI>)
   ```
2. Update `CITATION.cff` and `paper/manuscript.tex` under Code Availability.
