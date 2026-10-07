# NeuroGut-MetaSeq Makefile
# Publication & Reproducibility Orchestration

.PHONY: help setup demo data deseq meta pathways figures paper test clean docker-build docker-run all horizon1 horizon2 horizon3 horizon4 systems production

PYTHON ?= .venv/bin/python
PIP ?= .venv/bin/pip
PYTEST ?= .venv/bin/pytest

help:
	@echo "======================================================================"
	@echo "NeuroGut-MetaSeq: Computational Meta-Analysis Pipeline"
	@echo "======================================================================"
	@echo "Quickstart (Demo):"
	@echo "  make demo         : Run end-to-end pipeline on bundled demo data (10s)"
	@echo "  make paper        : Compile publication web paper in docs/index.html"
	@echo ""
	@echo "Scientific Reproduction (Real 60-Sample Cohorts):"
	@echo "  make horizon1     : Ingest & audit 60 samples across 4 cohorts (QC)"
	@echo "  make horizon2     : Run full-transcriptome GLMs (~ sex + condition)"
	@echo "  make horizon3     : Run Random-Effects meta-analysis across 23,096 genes"
	@echo "  make horizon4     : Run Systems Biology (GSEA, Regulons, WGCNA, Rescue)"
	@echo "  make systems      : Alias for make horizon4"
	@echo "  make production   : Run full real-cohort pipeline end-to-end"
	@echo ""
	@echo "Validation & Containerization:"
	@echo "  make test         : Run full 42-test automated unit test suite"
	@echo "  make docker-build : Build reproducible Docker container"
	@echo "  make docker-run   : Run pipeline inside Docker container"
	@echo "  make clean        : Remove intermediate generated results"
	@echo "======================================================================"

setup:
	python3 -m venv .venv
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

demo:
	@echo "[*] Step 1: Curating metadata for demo cohorts..."
	$(PYTHON) scripts/02_curate_metadata.py --demo
	@echo "[*] Step 2: Running differential gene expression..."
	$(PYTHON) scripts/03b_pydeseq2_analysis.py --demo
	@echo "[*] Step 3: Running cross-study meta-analysis..."
	$(PYTHON) scripts/04_meta_analysis.py --demo
	@echo "[*] Step 4: Running pathway enrichment..."
	$(PYTHON) scripts/05_pathway_enrichment.py --demo
	@echo "[*] Step 5: Generating publication figures..."
	$(PYTHON) scripts/06_generate_figures.py --demo
	@echo "[*] Step 6: Building interactive scientific web paper..."
	$(PYTHON) scripts/07_build_web_paper.py
	@echo "[SUCCESS] Demo workflow and web paper complete! View at docs/index.html"

horizon1:
	@echo "[*] Horizon 1: Curating 60 biological samples and running lineage QC..."
	$(PYTHON) scripts/02_curate_metadata.py
	$(PYTHON) scripts/02b_qc_audit.py

horizon2:
	@echo "[*] Horizon 2: Running negative binomial GLMs across 4 cohorts..."
	$(PYTHON) scripts/03b_pydeseq2_analysis.py

horizon3:
	@echo "[*] Horizon 3: Running DerSimonian-Laird Random-Effects meta-analysis..."
	$(PYTHON) scripts/04_meta_analysis.py

horizon4:
	@echo "[*] Horizon 4: Executing Systems Biology & Regulon Networks..."
	$(PYTHON) scripts/05a_cache_gene_sets.py
	$(PYTHON) scripts/05_pathway_enrichment.py
	$(PYTHON) scripts/05b_tf_regulon_analysis.py
	$(PYTHON) scripts/05c_coexpression_network.py
	$(PYTHON) scripts/05d_metabolite_rescue.py
	$(PYTHON) scripts/05e_systems_diagnostics.py
	@echo "[SUCCESS] Horizon 4 complete! Figures in results/pathways/figures/"

systems: horizon4

production: horizon1 horizon2 horizon3 horizon4 test

data:
	$(PYTHON) scripts/01_download_geo.py

deseq: horizon2

meta: horizon3

pathways:
	$(PYTHON) scripts/05_pathway_enrichment.py

figures:
	$(PYTHON) scripts/06_generate_figures.py

paper: demo

test:
	$(PYTEST) tests/ -v

clean:
	rm -rf data/processed/* results/de_results/* results/meta_results/* results/pathways/* results/figures/* docs/index.html
	touch data/processed/.gitkeep

docker-build:
	docker build -t neurogut-metaseq:latest .

docker-run:
	docker run --rm -v $(PWD)/docs:/workspace/docs -v $(PWD)/results:/workspace/results neurogut-metaseq:latest

all: demo test

