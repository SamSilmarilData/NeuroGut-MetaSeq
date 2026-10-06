"""
tests/test_web_paper_build.py
Verifies that the scientific web paper builds and contains necessary assets.
"""

import os
import pytest

def test_web_paper_files_exist():
    assert os.path.exists("docs/index.html"), "Missing docs/index.html web paper"

def test_web_paper_content():
    with open("docs/index.html", "r", encoding="utf-8") as f:
        html = f.read()

    assert "NeuroGut-MetaSeq" in html
    assert "Abstract" in html
    assert "Interactive Gene Explorer" in html
    assert "Figure 1 |" in html
    assert "Figure 2 |" in html
    assert "Figure 3 |" in html
    assert "META_DATA" in html
    assert len(html) > 5000, "Web paper HTML seems incomplete"
