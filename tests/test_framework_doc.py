"""
tests/test_framework_doc.py
Verifies that the Adaptive Discovery Framework and Living Lab Notebook exist
and contain mandatory structural sections.
"""

import os
import pytest

def test_adaptive_discovery_framework_exists():
    path = "docs/ADAPTIVE_DISCOVERY_FRAMEWORK.md"
    assert os.path.exists(path), f"Missing {path}"
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()

    assert "Structured Serendipity" in content
    assert "Macro Map" in content
    assert "Micro Engine" in content
    assert "5 Macro Horizons" in content
    assert "Micro Discovery Loop" in content
    assert len(content.split()) >= 800, "Framework document seems incomplete"

def test_lab_notebook_exists():
    path = "docs/LAB_NOTEBOOK.md"
    assert os.path.exists(path), f"Missing {path}"
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()

    assert "Living Research Lab Notebook" in content
    assert "Entry Template" in content
    assert "Entry 000" in content
    assert "Horizon 0: Project Charter" in content
    assert "Target Hypothesis" in content
    assert "Blooms, Anomalies & Serendipity" in content
