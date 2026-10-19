#!/usr/bin/env python3
"""
tests/test_biorxiv_and_latex_pipeline.py
Unit tests for the bioRxiv preprint submission package, LaTeX PDF compilation,
BibTeX database, and scientific hardening statements.
"""

import os
import json
import re
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def test_biorxiv_abstract_compliance():
    abstract_path = os.path.join(REPO_ROOT, "biorxiv_submission", "ABSTRACT_SUBMISSION.txt")
    assert os.path.exists(abstract_path), "ABSTRACT_SUBMISSION.txt does not exist"
    
    with open(abstract_path) as f:
        text = f.read().strip()
        
    chars_with_spaces = len(text)
    words = text.split()
    
    # bioRxiv character requirement: <= 2,500 characters including spaces
    assert chars_with_spaces <= 2500, f"Abstract exceeds bioRxiv limit: {chars_with_spaces} > 2500"
    assert 200 <= len(words) <= 300, f"Abstract word count outside expected 200-300: {len(words)}"

def test_biorxiv_metadata_json_integrity():
    meta_path = os.path.join(REPO_ROOT, "biorxiv_submission", "METADATA.json")
    assert os.path.exists(meta_path), "METADATA.json does not exist"
    
    with open(meta_path) as f:
        meta = json.load(f)
        
    required_keys = [
        "article_title", "primary_subject_area", "secondary_subject_area",
        "authors", "competing_interests", "data_availability", "code_availability"
    ]
    for k in required_keys:
        assert k in meta, f"Missing key in METADATA.json: {k}"
        
    assert meta["primary_subject_area"] == "Bioinformatics"
    assert meta["competing_interests"] == "The authors declare no competing interests."
    
    # Verify all 6 GEO accessions in data availability
    for acc in ["GSE107925", "GSE108045", "GSE186210", "GSE266602", "GSE152865", "GSE67858"]:
        assert acc in meta["data_availability"], f"Missing accession {acc} in data availability"

def test_bibtex_database_integrity():
    bib_path = os.path.join(REPO_ROOT, "paper", "references.bib")
    assert os.path.exists(bib_path), "references.bib does not exist"
    
    with open(bib_path) as f:
        bib_text = f.read()
        
    entries = re.findall(r'@article\{([^,]+),', bib_text)
    assert len(entries) == 46, f"Expected exactly 46 BibTeX entries, found {len(entries)}"

def test_compiled_pdfs_exist_and_non_empty():
    manuscript_pdf = os.path.join(REPO_ROOT, "paper", "manuscript.pdf")
    supplementary_pdf = os.path.join(REPO_ROOT, "paper", "supplementary.pdf")
    
    assert os.path.exists(manuscript_pdf), "manuscript.pdf does not exist"
    assert os.path.exists(supplementary_pdf), "supplementary.pdf does not exist"
    
    # Check sizes (manuscript should be > 2MB with figures, supplementary > 1MB)
    assert os.path.getsize(manuscript_pdf) > 2_000_000, "manuscript.pdf is suspiciously small"
    assert os.path.getsize(supplementary_pdf) > 1_000_000, "supplementary.pdf is suspiciously small"

def test_scientific_hardening_statements():
    doc_path = os.path.join(REPO_ROOT, "docs", "MANUSCRIPT.md")
    tex_path = os.path.join(REPO_ROOT, "paper", "manuscript.tex")
    
    with open(doc_path) as f:
        md = f.read()
    with open(tex_path) as f:
        tex = f.read()
        
    for text_content in [md, tex]:
        # Calibrated reversibility phrasing
        assert "chromatin-poised transcriptional reversibility" in text_content.lower()
        # Stereological validation for deconvolution
        assert "stereological" in text_content.lower()
        # Antibiotic caveat and leave-one-out robustness
        assert "leave-one-out" in text_content.lower()
        assert "mitochondrial" in text_content.lower()
