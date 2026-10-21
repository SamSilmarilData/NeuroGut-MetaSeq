# NeuroGut-MetaSeq: Publication-Grade LaTeX PDF Compilation Pipeline

This directory contains the complete, reproducible LaTeX source code and build automation scripts for generating publication-grade PDF documents for **NeuroGut-MetaSeq (v1.4.0)**.

---

## 📌 Document Status: Academic Release (v1.4.0)

> [!NOTE]
> This represents the **production academic release** (`v1.4.0`) of the preprint. It introduces the unified **6 Biological Acts** narrative, edge-to-edge panoramic schematics (Figures 1A and 6A), single-panel Table-Forest visualizations, Bland-Altman MA concordance transforms, and the Cross-Cohort Consensus Concordance Matrix (Figure 2D). All statistical results, multi-omic analyses, tables, figures, and literature citations compile cleanly with zero errors.

---

## 📂 File Architecture

| File | Purpose | Key Details |
| :--- | :--- | :--- |
| [`manuscript.tex`](manuscript.tex) | Primary Manuscript LaTeX Source | Line-numbered (`lineno`), `booktabs` Tables 1--8, Figures 1--6 inline in Results, unnumbered title header. |
| [`supplementary.tex`](supplementary.tex) | Standalone Supplementary Information | Independent pagination (`S1, S2...`), full mathematical derivations, Tables S1--S4, Figures S1--S9. |
| [`references.bib`](references.bib) | BibTeX Literature Database | 46 fully resolved entries with normalized author names, journal titles, DOIs, and PMIDs. |
| [`ref_key_map.json`](ref_key_map.json) | Citation Mapping Dictionary | Maps bracketed citation numbers from `docs/MANUSCRIPT.md` to BibTeX keys. |
| [`build.sh`](build.sh) | Standalone Shell Compilation Script | Executable script running 3 passes of `pdflatex` and `bibtex` for both documents. |
| [`Makefile`](Makefile) | Standard Unix Build Target | Provides `make all` and `make clean` targets. |
| [`figures/`](figures/) | High-Resolution Figure Assets | Contains all 300+ DPI PNG and PDF assets embedded directly into the LaTeX document. |
| [`manuscript.pdf`](manuscript.pdf) | Compiled Primary Manuscript PDF | 21 pages, ~10 MB, publication-grade bioRxiv layout with zero undefined references. |
| [`supplementary.pdf`](supplementary.pdf) | Compiled Supplementary Information PDF | 13 pages, ~6 MB, complete with table of contents, derivations, tables S1--S4, and figures S1--S9. |

---

## 🚀 How to Compile

Compilation requires a standard TeX Live distribution (such as MacTeX or TeX Live 2025) with `pdflatex` and `bibtex`.

### Method 1: Using the Build Script
```bash
cd paper
./build.sh
```

### Method 2: Using Make
```bash
cd paper
make clean
make all
```

### Method 3: Manual Step-by-Step
```bash
cd paper

# Compile Primary Manuscript
pdflatex -interaction=nonstopmode manuscript.tex
bibtex manuscript
pdflatex -interaction=nonstopmode manuscript.tex
pdflatex -interaction=nonstopmode manuscript.tex

# Compile Standalone Supplementary Information
pdflatex -interaction=nonstopmode supplementary.tex
bibtex supplementary || true
pdflatex -interaction=nonstopmode supplementary.tex
pdflatex -interaction=nonstopmode supplementary.tex
```

---

## 🎨 Typesetting & Design Standards

1. **Continuous Narrative Flow**:
   - Section headings flow naturally with balanced page-filling.
   - Dedicated full-page display plates host 4-panel composite figures cleanly within publication float constraints.
2. **Sequential Figure & Table Chronology**:
   - Every composite figure (Figures 1 through 6) and table (Tables 1 through 8) is embedded inline within its respective Results section in strictly monotonic narrative order.
3. **Peer-Review Ready Line Numbering**:
   - `\linenumbers` is activated starting at `\begin{abstract}`, keeping author affiliations clean while ensuring line-by-line reviewer commentary across the text.
4. **Professional Tabular Formatting**:
   - All tables utilize standard publication `booktabs` horizontal rules (`\toprule`, `\midrule`, `\bottomrule`) wrapped in `\resizebox{\textwidth}{!}{...}` with explicit `p{...}` column widths to prevent hyphenation ribboning.
5. **Full Literature Verification**:
   - All 46 citations from the literature review are indexed in `references.bib` with active DOIs and PMIDs, rendered through `\bibliographystyle{unsrtnat}` and `\nocite{*}`.

---

## 🗺️ Editorial & Visual Enhancements Implemented (v1.4.0)

- [x] **Conceptual Cartoon Schematics**:
  - Engineered panoramic widescreen overview schematic (`Fig 1A`) and Three-Pillar in vivo BBB flux model (`Fig 6A`), enlarged by +27% vertically and decoupled to span edge-to-edge (0 px / 1 px gap).
- [x] **Visual Hierarchy & Color Harmony**:
  - Replaced crowded 60-sample heatmap with Cross-Cohort Consensus Concordance Matrix (`Fig 2D`) featuring exact numerical LFCs, significance asterisks, and REML synthesis box.
  - Implemented Single-Panel Table-Forest plots (`Fig 2B`) and Bland-Altman MA concordance scatter plots (`Fig 3A`).
- [x] **Storytelling & Narrative Pacing**:
  - Restructured Results into **6 Biological Acts** (`docs/MANUSCRIPT_SKELETON_DESIGN.md`), optimizing page budget and reading momentum.
- [x] **Interactive Web Synchronization**:
  - Synchronized [`docs/index.html`](../docs/index.html) and [`docs/MANUSCRIPT.md`](../docs/MANUSCRIPT.md) continuously with LaTeX sources and figures.
