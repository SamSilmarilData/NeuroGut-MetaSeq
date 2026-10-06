# Contributing to NeuroGut-MetaSeq

We welcome contributions from computational biologists, bioinformaticians, and neuroimmunology researchers!

---

## 🔬 How You Can Contribute

1. **Adding New RNA-Seq Datasets**:
   - Suggest or integrate newly published bulk or single-cell microglial datasets examining microbiome perturbations, specific microbial metabolites (e.g. indole-3-propionic acid, deoxycholic acid, acetate), or diet-induced dysbiosis.
2. **Improving Statistical Meta-Analysis**:
   - Enhancing random-effects models (e.g., meta-regression, robust variance estimation with small-sample corrections, empirical Bayes shrinkage).
3. **Extending Pathway & Network Analysis**:
   - Integrating weighted gene co-expression network analysis (WGCNA) or transcription factor regulon inferencing (e.g., SCENIC / pySCENIC).
4. **Enhancing the Interactive Web Paper**:
   - Adding interactive volcano plot tooltips or dynamic pathway network graphs using D3.js or Plotly.

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
make test
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
