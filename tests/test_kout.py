from fractions import Fraction

"""Verification suite for the k-out local-structure theorems.

Cross-checks:
  (a) Theorem 1 general formula vs. direct exhaustive enumeration (small n);
  (b) Theorem 2 closed triangle form vs. exact triple-level enumeration;
  (c) C4 / star closed forms vs. Theorem 1;
  (d) Monte Carlo consistency for larger n.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import kout as K  # noqa: E402

TRIANGLE = ((0, 1), (1, 2), (2, 0))
C4 = ((0, 1), (1, 2), (2, 3), (3, 0))
STAR3 = ((0, 1), (0, 2), (0, 3))


# ---------------------------------------------------------------------------
# (a) Theorem 1 vs full-profile exhaustive enumeration
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n,k", [(3, 1), (3, 2), (4, 1), (4, 2), (4, 3), (5, 1), (5, 2), (5, 3), (5, 4)]
)
def test_theorem1_triangles_vs_full_enumeration(n, k):
    """E[#triangles] from Theorem 1 must match exhaustive enumeration
    over all (C(n-1,k))^n choice profiles (feasible for these sizes)."""
    exact = K.expected_triangles_enum_full(n, k)
    formula = comb_helper(n, 3) * K.p_triangle_closed(n, k)
    assert abs(exact - formula) < 1e-9, (n, k, exact, formula)


def comb_helper(n, r):
    from math import comb

    return comb(n, r)


@pytest.mark.parametrize("n,k", [(4, 1), (5, 1), (6, 1), (5, 2)])
def test_theorem1_C4_vs_full_enumeration(n, k):
    exact = K.count_subgraphs_enum_expectation(n, k, C4, 4, 8)
    formula = K.expected_unlabeled_copies(n, k, C4, 4, 8)
    assert abs(exact - formula) < 1e-9, (n, k, exact, formula)


@pytest.mark.parametrize("n,k", [(4, 1), (4, 2), (5, 2)])
def test_theorem1_star_vs_full_enumeration(n, k):
    exact = K.count_subgraphs_enum_expectation(n, k, STAR3, 4, 6)
    formula = K.expected_unlabeled_copies(n, k, STAR3, 4, 6)
    assert abs(exact - formula) < 1e-9, (n, k, exact, formula)


# ---------------------------------------------------------------------------
# (b) triple-level exact enumeration for larger n (subset-level, not marginal)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n,k", [(6, 2), (7, 2), (8, 3), (10, 3), (12, 4)])
def test_theorem2_triple_level(n, k):
    """The 4-state marginal computation (used inside expected_triangles_enum_triple)
    must agree with the closed form for larger (n, k)."""
    p_formula = K.p_triangle_closed(n, k)
    p_marginal = K.p_triangle_enum_triple(n, k)
    assert abs(p_formula - p_marginal) < 1e-10, (n, k, p_formula, p_marginal)
    e_formula = comb_helper(n, 3) * p_formula
    e_sub = K.expected_triangles_enum_triple(n, k)
    assert abs(e_formula - e_sub) < 1e-6, (n, k, e_formula, e_sub)


# ---------------------------------------------------------------------------
# (c) closed forms vs Theorem 1 general formula
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n,k", [(5, 1), (6, 2), (9, 3), (15, 4), (30, 10)])
def test_C4_closed_vs_theorem1(n, k):
    gen = K.p_H_in_kout(n, k, C4, 4)
    clo = K.p_C4_closed(n, k)
    assert abs(gen - clo) < 1e-10, (n, k, gen, clo)


@pytest.mark.parametrize("n,k,m", [(6, 1, 3), (8, 2, 2), (10, 3, 4), (12, 5, 3)])
def test_star_closed_vs_theorem1(n, k, m):
    star = tuple((0, i + 1) for i in range(m))
    gen = K.p_H_in_kout(n, k, star, m + 1)
    clo = K.p_star_closed(n, k, m)
    assert abs(gen - clo) < 1e-10, (n, k, m, gen, clo)


# ---------------------------------------------------------------------------
# (d) Monte Carlo consistency
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n,k", [(50, 2), (50, 3), (200, 3)])
def test_monte_carlo_triangles(n, k):
    reps = 20000
    samples = K.simulate_triangles_fast(n, k, reps, seed=12345)
    formula = comb_helper(n, 3) * K.p_triangle_closed(n, k)
    mean = samples.mean()
    se = samples.std(ddof=1) / np.sqrt(reps)
    # 5-sigma window: generous, deterministic seed
    assert abs(mean - formula) < 5.0 * se + 1e-12, (n, k, mean, formula, se)


def test_expected_edges_formula():
    """E[#edges] = n(n-1)/2 * (2p - p^2); also equals n/2 * (2k - k^2/(n-1))."""
    for n, k in [(10, 1), (10, 3), (20, 2), (30, 5)]:
        lhs = K.expected_edges_closed(n, k)
        rhs = 0.5 * n * (2 * k - k * k / (n - 1))
        assert abs(lhs - rhs) < 1e-9, (n, k)


def test_asymptotic_constants():
    for k in range(1, 8):
        assert K.asymptotic_triangle_constant(k) == (4 * k**3 - 3 * k**2) / 3
        assert K.er_comparison_limit(k) == 1 - 3 / (4 * k)


def test_naive_independence_is_wrong():
    """The naive 'independent edges' value differs from the exact formula:
    documented falsification of the naive computation."""
    n, k = 8, 2
    p = k / (n - 1)
    exact = K.p_triangle_closed(n, k)
    # exact formula must differ from the naive product (1-(1-p)^2)^3
    naive_indep = (2 * p - p * p) ** 3
    assert abs(exact - naive_indep) > 1e-9, (n, k, exact, naive_indep)


# ---------------------------------------------------------------------------
# (e) exact second moment of the triangle count
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n,k", [(4, 1), (5, 1), (5, 2), (6, 1)])
def test_second_moment_vs_full_enumeration(n, k):
    """E[T^2] from the union-graph decomposition must match exhaustive
    enumeration over all choice profiles."""
    e2, e1, v = K.triangle_moments_exact(n, k)
    e2e, e1e, ve = K.triangle_moments_enum_full(n, k)
    assert abs(e2 - e2e) < 1e-9, (n, k, e2, e2e)
    assert abs(e1 - e1e) < 1e-9, (n, k, e1, e1e)
    assert abs(v - ve) < 1e-9, (n, k, v, ve)


@pytest.mark.parametrize("n,k", [(50, 2), (50, 3), (100, 3)])
def test_second_moment_vs_monte_carlo(n, k):
    """Exact E[T] and Var(T) must be statistically consistent with Monte
    Carlo (5-sigma windows, deterministic seed)."""
    e2, e1, v = K.triangle_moments_exact(n, k)
    reps = 20000
    samples = K.simulate_triangles_fast(n, k, reps, seed=777)
    se_m = samples.std(ddof=1) / np.sqrt(reps)
    var_s = samples.var(ddof=1)
    se_v = var_s * np.sqrt(2.0 / (reps - 1))  # normal approx for Var SE
    assert abs(samples.mean() - e1) < 5.0 * se_m + 1e-12, (n, k)
    assert abs(var_s - v) < 5.0 * se_v + 1e-12, (n, k, var_s, v, se_v)


def test_subset_bitmasks():
    masks = K._subset_bitmasks(3, 2)
    assert len(masks) == 3


def test_count_subgraphs_adj():
    adj = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=bool)
    assert K.count_subgraphs_adj(adj, ((0, 1), (1, 2), (0, 2))) == 6
    assert (
        K.count_subgraphs_adj(
            np.array([[0, 1], [1, 0]], dtype=bool), ((0, 1), (1, 2), (0, 2))
        )
        == 0
    )


def test_falling():
    assert K.falling(3, 4) == 0
    assert K.falling(3, 3) == 6


def test_p_H_in_kout_exceptions():
    with pytest.raises(ValueError):
        K.p_H_in_kout(3, 4, ((0, 1),), 2)
    assert K.p_H_in_kout(3, 1, ((0, 1),), 4) == 0.0


def test_count_subgraphs_enum_expectation_exceptions():
    with pytest.raises(ValueError):
        K.count_subgraphs_enum_expectation(3, 4, ((0, 1),), 2, 1)
    assert K.count_subgraphs_enum_expectation(3, 1, ((0, 1),), 4, 1) == 0.0


def test_p_H_exact_exceptions():
    with pytest.raises(ValueError):
        K.p_H_exact(3, 4, ((0, 1),), 2)
    assert K.p_H_exact(3, 1, ((0, 1),), 4) == Fraction(0)


def test_simulate_kout():
    res = K.simulate_kout(4, 1, 5)
    assert len(res) == 5


def test_frac_funcs():
    assert K.p_triangle_exact(5, 2) == Fraction(10, 27)
    assert K.p_C4_exact(5, 2) == Fraction(343, 1296)
    assert K.p_star_exact(5, 2, 2) == Fraction(13, 24)
    assert K.expected_edges_exact(5, 2) == Fraction(15, 2)
    assert K.triangle_moments_exact_frac(5, 2) == (
        Fraction(395, 24),
        Fraction(100, 27),
        Fraction(15985, 5832),
    )


def test_p_triangle_enum_triple():
    assert abs(K.p_triangle_enum_triple(5, 2) - K.p_triangle_closed(5, 2)) < 1e-9


def test_expected_triangles_enum_full():
    assert K.expected_triangles_enum_full(5, 2) > 0


def test_count_subgraphs_enum_expectation():
    assert K.count_subgraphs_enum_expectation(4, 1, ((0, 1),), 2, 1) > 0


def test_expected_unlabeled_copies():
    assert K.expected_unlabeled_copies(4, 1, ((0, 1),), 2, 1) > 0
    assert K.expected_labeled_copies(4, 1, ((0, 1),), 2) > 0


def test_expected_edges_exact():
    assert K.expected_edges_exact(5, 2) == Fraction(15, 2)


def test_closed_funcs_exceptions():
    with pytest.raises(ValueError):
        K.p_triangle_closed(3, 4)
    with pytest.raises(ValueError):
        K.p_C4_closed(3, 4)
    with pytest.raises(ValueError):
        K.p_star_closed(3, 4, 2)
    with pytest.raises(ValueError):
        K.expected_edges_closed(3, 4)
    with pytest.raises(ValueError):
        K.p_triangle_exact(3, 4)
    with pytest.raises(ValueError):
        K.p_C4_exact(3, 4)
    with pytest.raises(ValueError):
        K.p_star_exact(3, 4, 2)
    with pytest.raises(ValueError):
        K.expected_edges_exact(3, 4)
