# NeuroGut-MetaSeq: Publication-Grade LaTeX PDF Compilation Pipeline

This directory contains the complete, reproducible LaTeX source code and build automation scripts for generating publication-grade PDF documents for **NeuroGut-MetaSeq (v1.3.0-draft1)**.

---

## 📌 Document Status: Preprint First Draft (v1.3.0-draft1)

> [!NOTE]
> This represents the **first compiled draft** of the preprint. While all statistical results, multi-omic analyses, tables, figures, and literature citations are fully generated and compiled with zero errors, iterative editorial refinements—including enhanced conceptual cartoon diagrams, visual aesthetics, narrative pacing, and typography adjustments—will be implemented in subsequent iterations.

---

## 📂 File Architecture

| File | Purpose | Key Details |
| :--- | :--- | :--- |
| [`manuscript.tex`](manuscript.tex) | Primary Manuscript LaTeX Source | Line-numbered (`lineno`), `booktabs` Tables 1--8, Figures 1--11 inline in Results, unnumbered title header. |
| [`supplementary.tex`](supplementary.tex) | Standalone Supplementary Information | Independent pagination (`S1, S2...`), full mathematical derivations, Tables S1--S4, Figures S1--S6. |
| [`references.bib`](references.bib) | BibTeX Literature Database | 46 fully resolved entries with normalized author names, journal titles, DOIs, and PMIDs. |
| [`ref_key_map.json`](ref_key_map.json) | Citation Mapping Dictionary | Maps bracketed citation numbers from `docs/MANUSCRIPT.md` to BibTeX keys. |
| [`build.sh`](build.sh) | Standalone Shell Compilation Script | Executable script running 3 passes of `pdflatex` and `bibtex` for both documents. |
| [`Makefile`](Makefile) | Standard Unix Build Target | Provides `make all` and `make clean` targets. |
| [`figures/`](figures/) | High-Resolution Figure Assets | Contains all 300+ DPI PNG assets embedded directly into the LaTeX document. |
| [`manuscript.pdf`](manuscript.pdf) | Compiled Primary Manuscript PDF | 23 pages, 9.0 MB, publication-grade bioRxiv layout with zero undefined references. |
| [`supplementary.pdf`](supplementary.pdf) | Compiled Supplementary Information PDF | 8 pages, 3.6 MB, complete with table of contents, derivations, tables, and figures. |

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
   - Artificial `\clearpage` breaks are restricted strictly to major structural transitions (e.g. before Discussion and Supplementary sections), preventing sparse or half-empty pages.
2. **Sequential Figure & Table Chronology**:
   - Every figure (Figures 1 through 11) and table (Tables 1 through 8) is embedded inline within its respective Results subsection in strictly monotonic narrative order.
3. **Peer-Review Ready Line Numbering**:
   - `\linenumbers` is activated starting at `\begin{abstract}`, keeping author affiliations clean while ensuring line-by-line reviewer commentary across the text.
4. **Professional Tabular Formatting**:
   - All tables utilize standard publication `booktabs` horizontal rules (`\toprule`, `\midrule`, `\bottomrule`) wrapped in `\resizebox{\textwidth}{!}{...}` with explicit `p{...}` column widths to prevent hyphenation ribboning.
5. **Full Literature Verification**:
   - All 46 citations from the literature review are indexed in `references.bib` with active DOIs and PMIDs, rendered through `\bibliographystyle{unsrtnat}` and `\nocite{*}`.

---

## 🗺️ Editorial & Visual Iteration Roadmap

While this First Draft is statistically complete and reviewer-proofed, upcoming iterations will focus on:

- [ ] **Conceptual Cartoon Schematics**:
  - Developing multi-panel visual cartoons (e.g. gut-vascular-meningeal interface, Three-Pillar BBB flux relay, and *Llgl2*-LAT1 nutrient scavenging axis) to complement empirical data plots.
- [ ] **Visual Hierarchy & Color Harmony**:
  - Refining color palettes across multi-cohort forest plots and scatter plots for maximum accessibility and visual contrast.
- [ ] **Storytelling & Narrative Pacing**:
  - Polishing prose transitions between Results movements and sharpening Discussion emphasis to optimize reader engagement.
- [ ] **Interactive Web Synchronization**:
  - Keeping [`docs/index.html`](../docs/index.html) and [`docs/MANUSCRIPT.md`](../docs/MANUSCRIPT.md) continuously synchronized with LaTeX edits.
