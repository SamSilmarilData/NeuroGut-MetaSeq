"""
tests/test_statistical_pipeline.py
Unit tests verifying mathematical accuracy of meta-analysis algorithms.
"""

import os
import sys
import importlib.util
import pytest
import numpy as np
from scipy import stats

# Dynamically import scripts/04_meta_analysis.py
spec = importlib.util.spec_from_file_location("meta_analysis", "scripts/04_meta_analysis.py")
meta_mod = importlib.util.module_from_spec(spec)
sys.modules["meta_analysis"] = meta_mod
spec.loader.exec_module(meta_mod)

benjamini_hochberg = meta_mod.benjamini_hochberg
dersimonian_laird = meta_mod.dersimonian_laird
combine_pvalues_fisher = meta_mod.combine_pvalues_fisher
combine_pvalues_stouffer = meta_mod.combine_pvalues_stouffer

def test_benjamini_hochberg_accuracy():
    pvals = np.array([0.01, 0.04, 0.03, 0.20])
    fdr = benjamini_hochberg(pvals)
    assert len(fdr) == 4
    assert fdr[0] == pytest.approx(0.04, rel=1e-3)
    assert fdr[1] <= 0.06
    assert fdr[2] <= 0.06
    assert fdr[3] == pytest.approx(0.20, rel=1e-3)

def test_dersimonian_laird_homogeneity():
    thetas = np.array([1.5, 1.5, 1.5])
    variances = np.array([0.1, 0.1, 0.1])
    theta_re, se_re, ci_low, ci_high, tau2, Q, I2, p_re = dersimonian_laird(thetas, variances)

    assert theta_re == pytest.approx(1.5, rel=1e-3)
    assert Q == pytest.approx(0.0, abs=1e-6)
    assert I2 == pytest.approx(0.0, abs=1e-6)
    assert tau2 == pytest.approx(0.0, abs=1e-6)
    assert ci_low < 1.5 < ci_high

def test_dersimonian_laird_heterogeneity():
    thetas = np.array([2.0, -1.0, 3.0])
    variances = np.array([0.05, 0.05, 0.05])
    theta_re, se_re, ci_low, ci_high, tau2, Q, I2, p_re = dersimonian_laird(thetas, variances)

    assert Q > 10.0
    assert I2 > 50.0
    assert tau2 > 0.0

def test_fisher_pvalue_combination():
    pvals = np.array([0.01, 0.02])
    chi2_stat, p_comb = combine_pvalues_fisher(pvals)

    expected_chi2 = -2.0 * (np.log(0.01) + np.log(0.02))
    expected_p = 1.0 - stats.chi2.cdf(expected_chi2, df=4)

    assert chi2_stat == pytest.approx(expected_chi2, rel=1e-4)
    assert p_comb == pytest.approx(expected_p, rel=1e-4)

def test_stouffer_z_combination():
    pvals = np.array([0.01, 0.01])
    thetas = np.array([1.0, 1.0])
    ns = np.array([10, 10])

    z, p = combine_pvalues_stouffer(pvals, thetas, ns)
    assert z > 0
    assert p < 0.01
