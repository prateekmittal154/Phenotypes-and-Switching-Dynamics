"""Regression checks against the frozen result tables shipped with the repository.

These guard the headline numbers quoted in the README and documentation. They read the
committed CSVs rather than recomputing, so they are fast and they fail loudly if a table
is replaced with something inconsistent.
"""
import os
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAB = os.path.join(ROOT, "results", "tables")
BTAB = os.path.join(ROOT, "benchmarks", "results", "tables")

needs_tables = pytest.mark.skipif(not os.path.isdir(TAB), reason="result tables not present")
needs_bench = pytest.mark.skipif(not os.path.isdir(BTAB), reason="benchmark tables not present")


@needs_tables
def test_fisher_optimum_is_near_point_four_three():
    d = pd.read_csv(os.path.join(TAB, "e12_optimum_audit.csv"))
    at_half = d[d.pi == 0.5].kL_opt_pi_joint
    assert len(at_half) > 0
    assert at_half.median() == pytest.approx(0.43, abs=0.06)
    # the optimum shifts only mildly with occupancy
    assert d.kL_opt_pi_joint.between(0.35, 0.55).all()


@needs_tables
def test_occupancy_heterogeneity_significant_in_every_eligible_dataset():
    d = pd.read_csv(os.path.join(TAB, "m6_real_bootstrap_tests.csv"))
    assert len(d) == 8, "eight eligible dataset-scope comparisons"
    assert (d.p_bootstrap <= 0.05).all(), "all eight significant at their own null"
    # the bootstrap thresholds genuinely vary, which is why chi-squared is unsafe
    assert d.boot_threshold_95.min() < 12 and d.boot_threshold_95.max() > 40


@needs_tables
def test_chi_squared_reference_is_anticonservative():
    d = pd.read_csv(os.path.join(TAB, "e14_calibration_summary.csv"))
    r = d[(d.scen == "H0_pi") & (d.test == "chi2_pi")].rate
    assert len(r) == 1
    assert float(r.iloc[0]) > 0.10, "type-I error well above the nominal 0.05"
    c = d[(d.scen == "H0_pi") & (d.test == "calib_pi")].rate
    assert float(c.iloc[0]) == pytest.approx(0.05, abs=0.02), "calibrated test hits nominal"


@needs_tables
def test_capacity_law_overpredicts_realised_coverage():
    d = pd.read_csv(os.path.join(TAB, "e16_capacity_vs_prediction.csv"))
    assert (d.decades_predicted > d.decades).all(), "the idealised law over-predicts"
    assert d.abs_error.median() == pytest.approx(0.70, abs=0.15)


@needs_bench
def test_benchmark_rate_estimators_disagree_in_direction():
    """The headline benchmark result: PATH low, PASTRI high, on the same trees."""
    d = pd.read_csv(os.path.join(BTAB, "TableB3_clone_level_true_tree.csv"))
    r = d[d.quantity == "switching rate"]
    path = r[(r.method == "PATH") & (r.scenario != "S1")].median_log_bias.median()
    pastri = r[(r.method == "PASTRI") & (r.scenario != "S1")].median_log_bias.median()
    assert path < 0, "PATH underestimates"
    assert pastri > 1.5, "PASTRI overestimates substantially"


@needs_bench
def test_pastri_is_invariant_to_branch_times():
    d = pd.read_csv(os.path.join(BTAB, "TableB4_clock_sensitivity.csv"))
    p = d[(d.method == "PASTRI") & (d.quantity == "rate / decay")]
    assert float(p.frac_identical.iloc[0]) == pytest.approx(1.0, abs=1e-9), \
        "PASTRI uses topological depth, so branch times cannot change its estimate"
    for m in ("Study D", "PATH", "SCOUT"):
        q = d[(d.method == m) & (d.quantity == "rate / decay")]
        if len(q):
            assert float(q.frac_identical.iloc[0]) < 0.5, f"{m} should respond to branch times"


@needs_bench
def test_occupancy_recovered_by_every_method():
    d = pd.read_csv(os.path.join(BTAB, "TableB3_clone_level_true_tree.csv"))
    occ = d[(d.quantity == "occupancy") & (d.scenario == "S3")]
    assert len(occ) >= 3
    assert (occ.corr_with_truth > 0.8).all(), "every method recovers occupancy in S3"
