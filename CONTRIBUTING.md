# Contributing to NeuroGut-MetaSeq

We welcome contributions from computational biologists, bioinformaticians, and neuroimmunology researchers!

---

## 🔬 How You Can Contribute

1. **Adding New RNA-Seq Datasets**:
   - Suggest or integrate newly published microglial RNA-seq datasets examining microbiome perturbations, specific microbial metabolites (e.g., indole-3-propionic acid, deoxycholic acid, acetate), or diet-induced dysbiosis.
2. **Improving Statistical Meta-Analysis**:
   - Enhancing random-effects models (e.g., meta-regression on sex or age, robust variance estimation with small-sample corrections, empirical Bayes shrinkage).
3. **Multi-Omics & Spatial Transcriptomics Expansion**:
   - Integrating microglial ATAC-seq chromatin accessibility data to validate TF regulon predictions (e.g., chromatin openness at *Irf1* and AP-1 binding sites).
   - Spatial transcriptomics (e.g., 10x Visium, MERFISH, Stereo-seq) to map consensus hits (*Llgl2*, *Clu*, *Slfn2*) across anatomical brain regions.
   - Extending single-cell subpopulation deconvolution (implemented in v1.1.0) with single-cell trajectory and RNA-velocity modeling.
4. **Cross-Species Translation**:
   - Extending the cross-study meta-analysis framework to human post-mortem brain transcriptomics in gut dysbiosis cohorts (e.g., Parkinson's, Alzheimer's, ASD).
5. **Enhancing the Interactive Web Paper**:
   - Adding interactive volcano plot tooltips or dynamic pathway network graphs using D3.js or Plotly for GitHub Pages.

---

## 🛠️ Development Workflow

### 1. Set Up Your Environment
```bash
git clone https://github.com/samyakmeshram/NeuroGut-MetaSeq.git
cd NeuroGut-MetaSeq

make setup
source .venv/bin/activate
```

### 2. Run the Test Suite
Ensure all existing tests pass before making modifications:
```bash
pytest tests/ -v  # Verify 53/53 tests passing
```

### 3. Adding a New Dataset
To integrate a new GEO series accession:
1. Add the accession and sample metadata details to `config/datasets.yaml`.
2. Add the download URLs for raw counts and series matrices.
3. Update `scripts/02_curate_metadata.py` to extract sample condition labels (reference vs perturbed) and map gene IDs.
4. Execute `make all` to verify that differential expression and meta-analysis incorporate the new cohort without errors.

### 4. Code Style & Standards
- Python code should adhere to PEP 8 standards and include type hints where practical.
- Keep data tables tidy with standardized column names matching `SPECIFICATION.md`.
- Avoid committing large raw genomic files (>50 MB); raw files should be downloaded on-demand and excluded via `.gitignore`.

---

## 📝 Commit Convention
Please use conventional commit prefixes:
- `feat:` for new features or datasets.
- `fix:` for bug fixes.
- `docs:` for documentation updates.
- `test:` for adding or updating unit tests.
- `refactor:` for code restructuring without behavioral changes.
