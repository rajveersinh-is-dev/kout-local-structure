"""Experiments for the k-out local-structure paper.

All exact quantities are computed with Fraction arithmetic (integer
binomials); Monte Carlo uses a fixed seed.  Every table below is written to
results/ as CSV, and results/summary.md summarizes the headline numbers.

Run:  python src/experiments.py        (about 3-6 minutes on a laptop)
"""
import csv
import os
import sys
import time
from fractions import Fraction
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np  # noqa: E402

import kout as K  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "results")
os.makedirs(RESULTS, exist_ok=True)

TRIANGLE = ((0, 1), (1, 2), (2, 0))
C4 = ((0, 1), (1, 2), (2, 3), (3, 0))
STAR_M3 = ((0, 1), (0, 2), (0, 3))


def write_csv(name, header, rows):
    """Save csv to file.
    
    Args:
        name:
        header:
        rows (list):
    
    """
    path = os.path.join(RESULTS, name)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {name} ({len(rows)} rows)")


def t0():
    """T0.
    
    Returns:
        The computed result
    
    """
    return time.perf_counter()


# ---------------------------------------------------------------------------
# 1. Exhaustive-enumeration validation (Theorem 1 vs brute force)
# ---------------------------------------------------------------------------


def exp_validation():
    """Exp validation.
    
    """
    print("[1] exhaustive-enumeration validation")
    start = t0()
    rows = []

    def add(case, n, k, formula, enumv):
        """Add.
        
        Args:
            case:
            n:
            k:
            formula:
            enumv:
        
        """
        rows.append([case, n, k, repr(formula), f"{float(formula):.12f}",
                     f"{float(enumv):.12f}", f"{abs(float(formula)-float(enumv)):.3e}"])

    for n, k in [(3, 1), (3, 2), (4, 1), (4, 2), (4, 3), (5, 1), (5, 2),
                 (5, 3), (5, 4), (6, 1)]:
        formula = Fraction(comb(n, 3), 1) * K.p_triangle_exact(n, k)
        enumv = K.expected_triangles_enum_full(n, k)
        add("triangles", n, k, formula, enumv)
    for n, k in [(4, 1), (5, 1), (6, 1), (5, 2)]:
        # E[unlabeled C4] = 3*C(n,4) * P(C4): (n)_4/8 = 3*C(n,4), |Aut(C4)|=8
        formula = Fraction(3 * comb(n, 4), 1) * K.p_C4_exact(n, k)
        enumv = K.count_subgraphs_enum_expectation(n, k, C4, 4, 8)
        add("C4", n, k, formula, enumv)
    for n, k in [(4, 1), (4, 2), (5, 2)]:
        # E[unlabeled K_{1,3}] = (n)_4/6 * P_star, |Aut(K_{1,3})| = 3! = 6
        formula = Fraction(K.falling(n, 4), 6) * K.p_star_exact(n, k, 3)
        enumv = K.count_subgraphs_enum_expectation(n, k, STAR_M3, 4, 6)
        add("star_K13", n, k, formula, enumv)
    for n, k in [(4, 1), (5, 1), (5, 2), (6, 1)]:
        e2, e1, v = K.triangle_moments_exact_frac(n, k)
        e2e, e1e, ve = K.triangle_moments_enum_full(n, k)
        add("E[T^2]", n, k, e2, e2e)
        add("Var(T)", n, k, v, ve)
    write_csv("validation_enumeration.csv",
              ["case", "n", "k", "formula_exact", "formula_float",
               "enumeration", "abs_diff"], rows)
    print(f"  done in {t0()-start:.1f}s")

    # triple-level exact enumeration (larger n)
    rows = []
    start = t0()
    for n, k in [(6, 2), (7, 2), (8, 3), (10, 3), (12, 4), (15, 4)]:
        p_exact = K.p_triangle_exact(n, k)
        p_sub = K.expected_triangles_enum_triple(n, k) / comb(n, 3)
        rows.append([n, k, f"{float(p_exact):.15f}", f"{float(p_sub):.15f}",
                     f"{abs(float(p_exact)-float(p_sub)):.3e}"])
    write_csv("validation_triple_level.csv",
              ["n", "k", "P_closed_exact", "P_triple_enum", "abs_diff"], rows)
    print(f"  triple-level done in {t0()-start:.1f}s")


# ---------------------------------------------------------------------------
# 2. Exact expectation tables
# ---------------------------------------------------------------------------


def exp_expectations():
    """Exp expectations.
    
    """
    print("[2] exact expectation tables")
    start = t0()
    rows = []
    for n in [10, 20, 50, 100, 200, 500, 1000]:
        for k in range(1, 7):
            Etri = float(Fraction(comb(n, 3), 1) * K.p_triangle_exact(n, k))
            Ec4 = float(Fraction(3 * comb(n, 4), 1) * K.p_C4_exact(n, k))
            Eedges = float(K.expected_edges_exact(n, k))
            Estar2 = float(Fraction(K.falling(n, 3), 2) * K.p_star_exact(n, k, 2))
            Estar3 = float(Fraction(K.falling(n, 4), 6) * K.p_star_exact(n, k, 3))
            q = 2 * k / (n - 1) - (k / (n - 1)) ** 2
            Enaive = comb(n, 3) * q ** 3
            rows.append([n, k, f"{Etri:.10f}", f"{Ec4:.10f}", f"{Eedges:.10f}",
                         f"{Estar2:.10f}", f"{Estar3:.10f}",
                         f"{Enaive:.10f}", f"{(Enaive-Etri)/Etri:.6f}"])
    write_csv("expectations_table.csv",
              ["n", "k", "E_triangles", "E_C4", "E_edges", "E_star_K12",
               "E_star_K13", "E_triangles_naive", "rel_naive_error"], rows)
    print(f"  done in {t0()-start:.1f}s")


# ---------------------------------------------------------------------------
# 3. Asymptotic convergence (exact arithmetic)
# ---------------------------------------------------------------------------
# machine-exact series coefficients (verified symbolically in
# src/verify_symbolic.py)
E1_COEF = lambda k: -k ** 2 * (k - 1) * (4 * k - 3) / 2       # e1
E2_COEF = lambda k: k ** 2 * (6 * k ** 3 - 39 * k ** 2 + 46 * k - 15) / 6
R1_COEF = lambda k: 3 * (4 * k - 3) / (8 * k)                 # ratio correction
V1_COEF = lambda k: -k ** 2 * (16 * k ** 2 - 21 * k + 7) / 2  # Var correction


def exp_asymptotics():
    """Exp asymptotics.
    
    """
    print("[3] asymptotic convergence")
    start = t0()
    rows = []
    for k in range(1, 6):
        lim = (4 * k ** 3 - 3 * k ** 2) / 3
        for n in [50, 100, 200, 400, 800, 1600]:
            E = float(Fraction(comb(n, 3), 1) * K.p_triangle_exact(n, k))
            q = 2 * k / (n - 1) - (k / (n - 1)) ** 2
            ratio = E / (comb(n, 3) * q ** 3)
            rows.append([k, n, f"{E:.10f}", f"{lim:.10f}",
                         f"{(E-lim)*n:.8f}", f"{E1_COEF(k):.8f}",
                         f"{ratio:.10f}", f"{1-3/(4*k):.10f}",
                         f"{ratio - (1-3/(4*k) + R1_COEF(k)/n):+.3e}"])
    write_csv("asymptotic_convergence.csv",
              ["k", "n", "E_tri", "limit", "(E-limit)*n", "e1",
               "ratio_to_ER", "ratio_limit", "ratio_minus_2term"], rows)
    print(f"  done in {t0()-start:.1f}s")


# ---------------------------------------------------------------------------
# 4. Monte Carlo validation
# ---------------------------------------------------------------------------


def exp_monte_carlo():
    """Exp monte carlo.
    
    """
    print("[4] Monte Carlo validation")
    start = t0()
    rows = []
    cases = [(50, 2, 20000), (50, 3, 20000), (100, 2, 20000), (100, 3, 20000),
             (200, 3, 20000), (200, 4, 20000), (500, 3, 4000)]
    for n, k, reps in cases:
        st = t0()
        samples = K.simulate_triangles_fast(n, k, reps, seed=20261018)
        e2, e1, v = K.triangle_moments_exact_frac(n, k)
        e1, v = float(e1), float(v)
        mean = float(samples.mean())
        var = float(samples.var(ddof=1))
        se_m = float(samples.std(ddof=1)) / np.sqrt(reps)
        se_v = var * np.sqrt(2.0 / (reps - 1))
        rows.append([n, k, reps, f"{mean:.6f}", f"{e1:.6f}",
                     f"{(mean-e1)/se_m:+.2f}", f"{var:.6f}", f"{v:.6f}",
                     f"{(var-v)/se_v:+.2f}", f"{t0()-st:.1f}"])
    write_csv("monte_carlo.csv",
              ["n", "k", "reps", "MC_mean", "exact_E", "z_mean", "MC_var",
               "exact_Var", "z_var", "seconds"], rows)
    print(f"  done in {t0()-start:.1f}s")


# ---------------------------------------------------------------------------
# 5. Star cancellation: n^m * P(K_{1,m}) -> c_m(k)
# ---------------------------------------------------------------------------
STAR_C = {   # first nonzero series coefficient of P(K_{1,m}), from sympy
    2: lambda k: k * (4 * k - 1),
    3: lambda k: 2 * k * (4 * k ** 2 - 3 * k + 1),
    4: lambda k: k * (16 * k ** 3 - 24 * k ** 2 + 19 * k - 6),
    5: lambda k: 2 * k * (16 * k ** 4 - 40 * k ** 3 + 55 * k ** 2 - 40 * k + 12),
    6: lambda k: k * (64 * k ** 5 - 240 * k ** 4 + 500 * k ** 3
                      - 615 * k ** 2 + 418 * k - 120),
}


def exp_stars():
    """Exp stars.
    
    """
    print("[5] star cancellation")
    start = t0()
    import math as _m
    rows = []
    for m in range(2, 7):
        for k in (2, 3):
            cm = STAR_C[m](k)
            for n in [100, 200, 400, 800, 1600]:
                P = float(K.p_star_exact(n, k, m))
                E = float(Fraction(K.falling(n, m + 1), 1)
                          * K.p_star_exact(n, k, m) / _m.factorial(m))
                rows.append([m, k, n, f"{n**m * P:.8f}", f"{cm:.8f}",
                             f"{n**m * P - cm:+.3e}", f"{E:.6f}"])
    write_csv("star_cancellation.csv",
              ["m", "k", "n", "n^m*P", "c_m(k)", "diff", "E_unlabeled"], rows)
    print(f"  done in {t0()-start:.1f}s")


# ---------------------------------------------------------------------------
# 6. Variance table
# ---------------------------------------------------------------------------


def exp_variance():
    """Exp variance.
    
    """
    print("[6] variance table")
    start = t0()
    rows = []
    for k in range(1, 6):
        lim = (4 * k ** 3 - 3 * k ** 2) / 3
        for n in [50, 100, 200, 500, 1000]:
            e2, e1, v = K.triangle_moments_exact_frac(n, k)
            e1, v = float(e1), float(v)
            rows.append([k, n, f"{e1:.10f}", f"{v:.10f}", f"{v/e1:.8f}",
                         f"{(v-lim)*n:.6f}", f"{V1_COEF(k):.6f}"])
    write_csv("variance_table.csv",
              ["k", "n", "E_T", "Var_T", "Var_over_E", "(Var-limit)*n",
               "v1"], rows)
    print(f"  done in {t0()-start:.1f}s")


if __name__ == "__main__":
    exp_validation()
    exp_expectations()
    exp_asymptotics()
    exp_monte_carlo()
    exp_stars()
    exp_variance()
    print("all experiments completed")
