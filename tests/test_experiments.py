import sys
import os
from unittest.mock import patch
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import experiments

def test_exp_validation(tmp_path):
    with patch("experiments.RESULTS", str(tmp_path)):
        experiments.exp_validation()
    assert (tmp_path / "validation_enumeration.csv").exists()
    assert (tmp_path / "validation_triple_level.csv").exists()

def test_exp_expectations(tmp_path):
    with patch("experiments.RESULTS", str(tmp_path)):
        experiments.exp_expectations()
    assert (tmp_path / "expectations_table.csv").exists()

def test_exp_asymptotics(tmp_path):
    with patch("experiments.RESULTS", str(tmp_path)):
        experiments.exp_asymptotics()
    assert (tmp_path / "asymptotic_convergence.csv").exists()

def test_exp_monte_carlo(tmp_path):
    with patch("experiments.RESULTS", str(tmp_path)):
        # To avoid the >100s run time on CI which could cause flaky failures,
        # we mock the slow part just for this test, or we run it with reduced reps.
        # But wait, rule 5: "Never weaken, delete, skip, or trivialize tests merely to obtain a green CI run."
        # If we patch the reps, we are trivializing the MC test.
        # But this is testing the `experiments.py` script itself.
        # The script takes 100 seconds, which is perfectly fine for CI. We will run it fully.
        experiments.exp_monte_carlo()
    assert (tmp_path / "monte_carlo.csv").exists()

def test_exp_stars(tmp_path):
    with patch("experiments.RESULTS", str(tmp_path)):
        experiments.exp_stars()
    assert (tmp_path / "star_cancellation.csv").exists()

def test_exp_variance(tmp_path):
    with patch("experiments.RESULTS", str(tmp_path)):
        experiments.exp_variance()
    assert (tmp_path / "variance_table.csv").exists()
