# NeuroGut-MetaSeq Makefile
# Publication & Reproducibility Orchestration

.PHONY: help setup demo data deseq meta pathways figures paper test clean docker-build docker-run all

PYTHON ?= .venv/bin/python
PIP ?= .venv/bin/pip
PYTEST ?= .venv/bin/pytest

help:
	@echo "======================================================================"
	@echo "NeuroGut-MetaSeq: Computational Meta-Analysis Pipeline"
	@echo "======================================================================"
	@echo "Available commands:"
	@echo "  make setup        : Create Python venv and install dependencies"
	@echo "  make demo         : Run full end-to-end pipeline on bundled demo data"
	@echo "  make paper        : Compile publication web paper in docs/index.html"
	@echo "  make data         : Download real RNA-seq cohorts from NCBI GEO"
	@echo "  make deseq        : Run cohort-level differential expression analysis"
	@echo "  make meta         : Run cross-study statistical meta-analysis"
	@echo "  make pathways     : Run GO and KEGG pathway enrichment"
	@echo "  make figures      : Generate publication-ready figures in results/figures"
	@echo "  make test         : Run automated test suite"
	@echo "  make docker-build : Build reproducible Docker container"
	@echo "  make docker-run   : Run pipeline inside Docker container"
	@echo "  make all          : Full pipeline execution and web paper generation"
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

data:
	$(PYTHON) scripts/01_download_geo.py

deseq:
	$(PYTHON) scripts/03b_pydeseq2_analysis.py

meta:
	$(PYTHON) scripts/04_meta_analysis.py

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
