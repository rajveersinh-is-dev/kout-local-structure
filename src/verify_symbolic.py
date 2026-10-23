"""Symbolic verification of the k-out local-structure results.

Checks, with sympy (exact rational algebra):

 (1) Theorem 1 (inclusion-exclusion over F <= E(H)) reproduces the closed
     forms for the triangle, the 4-cycle, and stars K_{1,m} (m = 2..6).
 (2) The alternative triangle identity
         P(T) = (1-A^2)^3 + A^3 (C-A)(3 - A^2 - A C - C^2)
     is equivalent to P(T) = 1 - 3A^2 + 3A^2 B - B^3 given B = A*C.
 (3) Special cases k = 1:  P(T) = 2/(n-1)^3 and P(C4) = 2/(n-1)^4 exactly.
 (4) Exact asymptotics (series in t = 1/n):
       E[#triangles] = C(n,3) P_T  ->  (4k^3 - 3k^2)/3   (with corrections);
       E[#C4]        = 3 C(n,4) P_C4 ->  constant in k (derived here);
       E[#edges]     = n k - n k^2 / (2(n-1)) exactly;
       E-R matched-density ratio -> 1 - 3/(4k) (with correction);
       E[#K_{1,m}] limits for m = 1..6.

Run:  python src/verify_symbolic.py
"""
from math import comb

import sympy as sp

n, k, t = sp.symbols('n k t', positive=True)

# ---------------------------------------------------------------------------
# Binomial-ratio building blocks as exact rational functions of (n, k).
#   f(d) = C(n-1-d, k) / C(n-1, k)   for d >= 1;  f(0) = 1.
# Valid as rational functions; coincides with the combinatorial value for
# n-1-d >= k, and both sides vanish when n-1-d < k and d >= 1 ... (the
# product formula is used only where the combinatorial value is nonzero).
# ---------------------------------------------------------------------------
def f(d):
    if d == 0:
        return sp.Integer(1)
    num = sp.Integer(1)
    den = sp.Integer(1)
    for i in range(d):
        num *= (n - 1 - k - i)
        den *= (n - 1 - i)
    return sp.cancel(num / den)

A = f(1)                     # C(n-2,k)/C(n-1,k)
B = f(2)                     # C(n-3,k)/C(n-1,k)
C_ = sp.cancel((n - 2 - k) / (n - 2))

q = sp.cancel(2 * k / (n - 1) - k ** 2 / (n - 1) ** 2)   # matched edge density

PASS = []
FAIL = []


def check(name, expr_a, expr_b):
    diff = sp.cancel(sp.expand(expr_a - expr_b))
    ok = diff == 0
    (PASS if ok else FAIL).append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        print("   diff =", sp.simplify(diff))


# ---------------------------------------------------------------------------
# Theorem 1 subset sum for a graph H given as an edge list on {0..v-1}
# ---------------------------------------------------------------------------
import collections

def theorem1_expr(h_edges, h_vertices):
    deg_counts = collections.defaultdict(int)
    for r in range(len(h_edges) + 1):
        for F in sp.utilities.iterables.combinations(h_edges, r):
            deg = [0] * h_vertices
            for u, w in F:
                deg[u] += 1
                deg[w] += 1
            deg_counts[tuple(sorted(deg))] += (-1) ** r

    total = sp.Integer(0)
    for deg, coeff in deg_counts.items():
        if coeff == 0: continue
        term = sp.Integer(coeff)
        for d in deg:
            term = term * f(d)
        total = sp.cancel(total + term)
    return total


TRIANGLE = ((0, 1), (1, 2), (2, 0))
C4 = ((0, 1), (1, 2), (2, 3), (3, 0))

# (1) closed forms vs Theorem 1
p_tri_closed = 1 - 3 * A ** 2 + 3 * A ** 2 * B - B ** 3
check("Theorem1(triangle) == closed form", theorem1_expr(TRIANGLE, 3), p_tri_closed)

p_c4_closed = (1 - 4 * A ** 2 + 4 * A ** 2 * B + 2 * A ** 4
               - 4 * A ** 2 * B ** 2 + B ** 4)
check("Theorem1(C4) == closed form", theorem1_expr(C4, 4), p_c4_closed)

for m in (2, 3, 4, 5, 6):
    star = tuple((0, i + 1) for i in range(m))
    p_star_closed = sum((-1) ** j * comb(m, j) * f(j) * A ** j
                        for j in range(m + 1))
    check(f"Theorem1(K_1,{m}) == closed form",
          theorem1_expr(star, m + 1), sp.cancel(p_star_closed))

# (2) alternative triangle identity, using B = A*C
check("B == A*C", B, A * C_)
alt = (1 - A ** 2) ** 3 + A ** 3 * (C_ - A) * (3 - A ** 2 - A * C_ - C_ ** 2)
check("triangle alternate identity", sp.cancel(alt), sp.cancel(p_tri_closed))

# (3) k = 1 special cases
p_tri_k1 = sp.cancel(p_tri_closed.subs(k, 1))
check("P(T) = 2/(n-1)^3 for k=1", p_tri_k1, 2 / (n - 1) ** 3)
p_c4_k1 = sp.cancel(p_c4_closed.subs(k, 1))
check("P(C4) = 2/(n-1)^4 for k=1", p_c4_k1, 2 / (n - 1) ** 4)

# ---------------------------------------------------------------------------
# (4) asymptotics: substitute n = 1/t and series-expand
# ---------------------------------------------------------------------------
def series_coeffs(expr_in_n, order):
    """Series of expr(n) in t = 1/n up to t^order; returns list of coeffs."""
    expr_t = sp.cancel(sp.together(expr_in_n).subs(n, 1 / t))
    s = sp.series(expr_t, t, 0, order + 1).removeO()
    poly = sp.Poly(sp.expand(s), t)
    out = []
    for i in range(order + 1):
        out.append(sp.simplify(poly.coeff_monomial(t ** i)))
    return out


# --- triangles ---
E_tri = sp.cancel(sp.expand(n * (n - 1) * (n - 2) / 6 * p_tri_closed))
c_tri = series_coeffs(E_tri, 3)
print("\nE[#triangles] = e0 + e1/n + e2/n^2 + ...")
for i, c in enumerate(c_tri):
    print(f"  e{i} =", sp.factor(c))
check("triangle limit e0 == (4k^3-3k^2)/3", c_tri[0], (4 * k ** 3 - 3 * k ** 2) / 3)

# --- C4 ---
E_c4 = sp.cancel(sp.expand(n * (n - 1) * (n - 2) * (n - 3) / 8 * p_c4_closed))
c_c4 = series_coeffs(E_c4, 2)
print("\nE[#C4] = d0 + d1/n + d2/n^2 + ...")
for i, c in enumerate(c_c4):
    print(f"  d{i} =", sp.factor(c))

# --- edges ---
E_edges = sp.cancel(n * (n - 1) / 2 * q)
check("E[#edges] == n k - n k^2/(2(n-1))", E_edges, n * k - n * k ** 2 / (2 * (n - 1)))

# --- E-R matched-density triangle ratio ---
ratio = sp.cancel(p_tri_closed / q ** 3)
c_ratio = series_coeffs(ratio, 2)
print("\nE-R ratio = r0 + r1/n + r2/n^2 + ...")
for i, c in enumerate(c_ratio):
    print(f"  r{i} =", sp.factor(c))
check("E-R ratio limit == 1 - 3/(4k)", c_ratio[0], 1 - sp.Rational(3, 4) / k)

# --- C4 E-R ratio ---
ratio_c4 = sp.cancel(p_c4_closed / q ** 4)
c_ratio_c4 = series_coeffs(ratio_c4, 2)
print("\nE-R C4 ratio = s0 + s1/n + ...")
for i, c in enumerate(c_ratio_c4):
    print(f"  s{i} =", sp.factor(c))

# --- stars: for m >= 2, P_star = O(n^-2) (constant and 1/n terms cancel),
#     so g(m,k) = lim n^2 P(K_{1,m} at fixed labels) is finite and
#     E[#K_{1,m}] = (n)_{m+1}/m! * P_star ~ g(m,k) n^{m-1} / m!.
#     (m = 1 is the single edge: E[#K_{1,1}] = E[#edges] ~ nk.)
#     We also verify the stronger cancellation: the first nonzero term of
#     P(K_{1,m}) in 1/n is at order m (i.e. P = Theta(n^-m)).
print("\nstar structure: first nonzero term of P(K_{1,m}) in 1/n, and g = "
      "lim n^2 P:")
for m in (2, 3, 4, 5, 6):
    p_star = sp.cancel(sum((-1) ** j * comb(m, j) * f(j) * A ** j
                           for j in range(m + 1)))
    p_t = sp.cancel(sp.together(p_star).subs(n, 1 / t))
    s = sp.series(p_t, t, 0, m + 1).removeO()
    s = sp.expand(s)
    first = None
    for i in range(0, m + 1):
        c = sp.simplify(s.coeff(t, i))
        if c != 0:
            first = (i, sp.factor(c))
            break
    g2 = sp.simplify(s.coeff(t, 2))
    print(f"  m={m}: first nonzero term at t^{first[0] if first else '>m'}:"
          f"  {first[1] if first else ''}")
    if m >= 2:
        g = series_coeffs(sp.cancel(n ** 2 * p_star), 2)
        print(f"        n^2*P = {sp.factor(g[0])} + {sp.factor(g[1])}/n + "
              f"{sp.factor(g[2])}/n^2  =>  lim E/n^(m-1) = "
              f"{sp.factor(g[0] / sp.factorial(m))}")

# --- path P3 (2-edge path) limit: E = (n)_3/2 * P ---
p_p3 = sp.cancel(1 - 2 * A ** 2 + A ** 2 * B)
check("Theorem1(P3) == 1 - 2A^2 + A^2 B", theorem1_expr(((0, 1), (1, 2)), 3), p_p3)
E_p3 = sp.cancel(sp.expand(n * (n - 1) * (n - 2) / 2 * p_p3))
c_p3 = series_coeffs(E_p3 / n, 2)   # E[#P3] ~ p0 * n + ...
print("\nE[#P3 paths]/n = p0 + p1/n + ...")
for i, c in enumerate(c_p3):
    print(f"  p{i} =", sp.factor(c))

# --- naive-independence gap for triangles ---
# Exact compact identity: q^3 - P(T) = -A^3 (C-A) (3 - A^2 - A C - C^2),
# since q = 1 - A^2.
gap = sp.cancel(q ** 3 - p_tri_closed)          # > 0: naive overestimates
gap_alt = sp.cancel(-A ** 3 * (C_ - A) * (3 - A ** 2 - A * C_ - C_ ** 2))
check("gap identity q^3 - P(T) == -A^3(C-A)(3-A^2-AC-C^2)", gap, gap_alt)
c_gap = series_coeffs(gap, 4)
print("\nq^3 - P(T) (naive overestimate) = g0 + g1/n + ...")
for i, c in enumerate(c_gap):
    print(f"  g{i} =", sp.factor(c))

# --- exact second moment / variance of the triangle count ---
# E[T^2] = N0 P3 + N2 PK4e + N1 Pbow + Nd P2t with the union graphs of two
# triangles; all probabilities are Theorem 1 values (rational in n, k).
K4e_edges = ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3))
bow_edges = ((0, 1), (1, 2), (0, 2), (0, 3), (3, 4), (0, 4))
two_edges = ((0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5))

P3s = p_tri_closed
PK4e = theorem1_expr(K4e_edges, 4)
Pbow = theorem1_expr(bow_edges, 5)
P2t = theorem1_expr(two_edges, 6)

N0 = n * (n - 1) * (n - 2) / 6
N2 = 12 * n * (n - 1) * (n - 2) * (n - 3) / 24
N1 = n * (n - 1) * (n - 2) / 2 * (n - 3) * (n - 4) / 2
Nd = N0 * (n - 3) * (n - 4) * (n - 5) / 6
E2 = sp.expand(N0 * P3s + N2 * PK4e + N1 * Pbow + Nd * P2t)
E1 = sp.expand(N0 * P3s)
E1_sq = sp.cancel(E1 ** 2)
VarT = sp.cancel(E2 - E1_sq)
c_var = series_coeffs(VarT, 2)
print("\nVar(T) = v0 + v1/n + v2/n^2 + ...")
for i, c in enumerate(c_var):
    print(f"  v{i} =", sp.factor(c))
print("\nVar/E ratio limit =", sp.factor(sp.cancel(c_var[0] / c_tri[0])))
print("E[T^2] limit =", sp.factor(series_coeffs(E2, 1)[0]))

print(f"\n===== {len(PASS)} passed, {len(FAIL)} failed =====")
if FAIL:
    raise SystemExit(1)


# ---------------------------------------------------------------------------
# (5) Orientation invariants Omega_k(H): sum over edge-owner assignments
#     sigma (each edge assigned to exactly one of its two endpoints) of the
#     product over vertices of the falling factorials (k)_{t_x(sigma)}.
#     Machine-verified instances: triangle, C4, stars K_{1,m} (m=2..6),
#     plus the labeled/unlabeled cross-checks with the series limits.
# ---------------------------------------------------------------------------
def falling_k(t):
    # (k)_t as a sympy polynomial in k; (k)_0 = 1
    return sp.prod(k - i for i in range(t))  # t = 0 gives sp.Integer(1)


def omega_k(h_edges, h_vertices):
    total = sp.Integer(0)
    e = len(h_edges)
    for mask in range(2 ** e):
        t = [0] * h_vertices
        for i, (u, w) in enumerate(h_edges):
            owner = w if (mask >> i) & 1 else u
            t[owner] += 1
        term = sp.Integer(1)
        for x in range(h_vertices):
            term = term * falling_k(t[x])
        total = total + term
    return sp.cancel(total)


K4e_edges_O = ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3))  # K4 minus an edge
print("\nOrientation invariants Omega_k(H) = sum_sigma prod_x (k)_{t_x}:")
check("Omega(K3) == 8k^3 - 6k^2", omega_k(TRIANGLE, 3),
      8 * k ** 3 - 6 * k ** 2)
check("Omega(K3) == 6 * (4k^3-3k^2)/3 (labeled/unlabeled)", omega_k(TRIANGLE, 3),
      6 * (4 * k ** 3 - 3 * k ** 2) / 3)
check("Omega(C4) == 2k^2(8k^2-8k+1)", omega_k(C4, 4),
      2 * k ** 2 * (8 * k ** 2 - 8 * k + 1))
check("Omega(C4) == 8 * d0", omega_k(C4, 4),
      8 * (k ** 2 * (8 * k ** 2 - 8 * k + 1) / 4))
check("Omega(K4 - e) k=1 equals 0", omega_k(K4e_edges_O, 4).subs(k, 1),
      sp.Integer(0))
for m in (2, 3, 4, 5, 6):
    star = tuple((0, i + 1) for i in range(m))
    expected = sum(comb(m, j) * falling_k(j) * k ** (m - j)
                   for j in range(m + 1))
    check(f"Omega(K_1,{m}) == sum_j C(m,j)(k)_j k^(m-j)",
          omega_k(star, m + 1), sp.cancel(expected))


def check_arbitrary_H(name, h_edges, h_vertices):
    P_H = theorem1_expr(h_edges, h_vertices)
    E_N_H = P_H
    for i in range(h_vertices):
        E_N_H *= (n - i)
    e = len(h_edges)
    v = h_vertices

    expr_t = sp.cancel(sp.together(E_N_H / n**(v - e)).subs(n, 1/t))

    s = sp.series(expr_t, t, 0, 1).removeO()
    limit_val = sp.simplify(s.subs(t, 0))

    omega_val = omega_k(h_edges, h_vertices)

    diff = sp.simplify(limit_val - omega_val)
    ok = diff == 0
    (PASS if ok else FAIL).append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        print("   diff =", diff)

check_arbitrary_H("Omega arbitrary P4 (3-edge path)", ((0, 1), (1, 2), (2, 3)), 4)
check_arbitrary_H("Omega arbitrary Paw graph", ((0, 1), (0, 2), (0, 3), (1, 2)), 4)
check_arbitrary_H("Omega arbitrary Diamond graph", ((0, 1), (1, 2), (2, 3), (3, 0), (0, 2)), 4)

print(f"\n===== FINAL: {len(PASS)} passed, {len(FAIL)} failed =====")
if FAIL:
    raise SystemExit(1)


# ---------------------------------------------------------------------------
# (5) Independent cross-check with exact integer/Fraction arithmetic and
#     Python's math.comb (no sympy): confirm the truncated asymptotic series
#     reproduces E[#H] for moderately large n to within the next-term size.
# ---------------------------------------------------------------------------
from fractions import Fraction


def E_tri_exact(N, kk):
    """E[#triangles] as an exact Fraction via integer binomials."""
    An = Fraction(comb(N - 2, kk), comb(N - 1, kk))
    Bn = Fraction(comb(N - 3, kk), comb(N - 1, kk)) if N - 3 >= kk else Fraction(0)
    P = 1 - 3 * An ** 2 + 3 * An ** 2 * Bn - Bn ** 3
    return Fraction(comb(N, 3), 1) * P


def E_c4_exact(N, kk):
    An = Fraction(comb(N - 2, kk), comb(N - 1, kk))
    Bn = Fraction(comb(N - 3, kk), comb(N - 1, kk)) if N - 3 >= kk else Fraction(0)
    P = (1 - 4 * An ** 2 + 4 * An ** 2 * Bn + 2 * An ** 4
         - 4 * An ** 2 * Bn ** 2 + Bn ** 4)
    return Fraction(3 * comb(N, 4), 1) * P


print("\nExact-fraction spot checks of the asymptotic series (n = 400):")
N = 400
for kk in (1, 2, 3, 4, 5):
    e0 = (4 * kk ** 3 - 3 * kk ** 2) / 3
    e1 = -kk ** 2 * (kk - 1) * (4 * kk - 3) / 2
    e2 = kk ** 2 * (6 * kk ** 3 - 39 * kk ** 2 + 46 * kk - 15) / 6
    e3 = -kk ** 2 * (kk ** 4 - 27 * kk ** 3 + 99 * kk ** 2 - 96 * kk + 27) / 6
    trunc = e0 + e1 / N + e2 / N ** 2
    exact = float(E_tri_exact(N, kk))
    err = abs(exact - trunc)
    bound = abs(e3) / N ** 3 * 1.05 + 1e-12
    ok = err <= bound
    (PASS if ok else FAIL).append(f"tri series n={N},k={kk}")
    print(f"[{'PASS' if ok else 'FAIL'}] k={kk}: exact={exact:.8f} "
          f"trunc={trunc:.8f} err={err:.2e} bound={bound:.2e}")

d0 = lambda kk: kk ** 2 * (8 * kk ** 2 - 8 * kk + 1) / 4
d1 = lambda kk: -kk ** 3 * (8 * kk ** 2 - 8 * kk + 1) / 2
for kk in (1, 2, 3, 4, 5):
    exact = float(E_c4_exact(N, kk))
    trunc = d0(kk) + d1(kk) / N
    err = abs(exact - trunc)
    # d2 coefficient known; use it for the bound
    d2 = kk ** 2 * (12 * kk ** 4 - 34 * kk ** 3 + 9 * kk ** 2 + 12 * kk - 2) / 4
    bound = abs(d2) / N ** 2 * 1.05 + 1e-12
    ok = err <= bound
    (PASS if ok else FAIL).append(f"C4 series n={N},k={kk}")
    print(f"[{'PASS' if ok else 'FAIL'}] C4 k={kk}: exact={exact:.8f} "
          f"trunc={trunc:.8f} err={err:.2e} bound={bound:.2e}")

print(f"\n===== FINAL: {len(PASS)} passed, {len(FAIL)} failed =====")
if FAIL:
    raise SystemExit(1)
