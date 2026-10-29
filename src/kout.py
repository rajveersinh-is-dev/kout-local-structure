"""Exact local-structure analysis of the random k-out graph.

Model (standard symmetrized k-out graph, cf. Frieze & Karonski, Introduction to
Random Graphs, Ch. 22): every vertex v in {0,...,n-1} independently and uniformly
selects a k-subset S_v of the other n-1 vertices; the undirected simple graph
G(n,k-out) has edge {u,v} iff u chose v or v chose u.

Main object: exact expected number of copies of a fixed finite simple graph H.
"""

from __future__ import annotations

import itertools
from fractions import Fraction
from math import comb, factorial

import numpy as np

# ---------------------------------------------------------------------------
# Model primitives
# ---------------------------------------------------------------------------


def _subset_bitmasks(m: int, k: int) -> np.ndarray:
    """All k-subsets of {0,...,m-1} encoded as integer bitmasks (ascending)."""
    out = []
    for combo in itertools.combinations(range(m), k):
        mask = 0
        for x in combo:
            mask |= 1 << x
        out.append(mask)
    return np.array(out, dtype=np.int64)


def _ksubset_masks_excluding(n: int, k: int, v: int) -> np.ndarray:
    """All k-subsets of {0,...,n-1} \\ {v}, as n-bit integer bitmasks.

    These are the possible choice sets of vertex v in the k-out model.
    """
    out = []
    for combo in itertools.combinations([x for x in range(n) if x != v], k):
        mask = 0
        for x in combo:
            mask |= 1 << x
        out.append(mask)
    return np.array(out, dtype=np.int64)


def _drop_self_targets(perm: np.ndarray, k: int) -> np.ndarray:
    """From a batch of permutations of {0,...,n-1} (rows), take, per row, the
    first k entries that differ from that row's vertex label.

    Rows are grouped per replicate: `perm` has shape (R, n) where consecutive
    blocks of n rows belong to distinct replicates; the row's vertex label is
    its position within the block.  Returns the (R, k) target array.  The
    resulting k-set is uniform over k-subsets of [n] \\ {v} by symmetry of the
    random permutation.
    """
    R, n = perm.shape
    part = perm[:, : k + 1]  # (R, k+1)
    v = np.tile(np.arange(n), R // n)[:, None]  # row vertex labels
    cols = np.arange(k + 1)[None, :]
    rank = np.where(part != v, cols, k + 1)  # column rank if kept
    order = np.sort(rank, axis=1)[:, :k]  # columns of the k kept
    return np.take_along_axis(part, order, axis=1)


def choices_to_adjacency(choices: np.ndarray, n: int) -> np.ndarray:
    """Symmetrized adjacency from a choice profile.

    choices: (n,) array of integer bitmasks; bit j set in choices[v] means
    "vertex v chose vertex j". Returns (n, n) boolean adjacency of the
    underlying simple graph.
    """
    bits = ((np.asarray(choices)[:, None] >> np.arange(n)[None, :]) & 1).astype(bool)
    return bits | bits.T


def count_triangles_adj(adj: np.ndarray) -> int:
    """Number of triangles of a symmetric boolean adjacency matrix."""
    a = adj.astype(np.int64)
    return int(np.trace(a @ a @ a) // 6)


def count_subgraphs_adj(adj: np.ndarray, h_edges: tuple[tuple[int, int], ...]) -> int:
    """Number of (labeled) copies of H in adj via brute force over vertex maps."""
    n = adj.shape[0]
    v = max(max(e) for e in h_edges) + 1
    if n < v:
        return 0
    count = 0
    for phi in itertools.permutations(range(n), v):
        ok = all(adj[phi[u], phi[w]] for u, w in h_edges)
        if ok:
            count += 1
    return count


def count_subgraphs_enum_expectation(
    n: int, k: int, h_edges: tuple[tuple[int, int], ...], h_vertices: int, aut_size: int
) -> float:
    """E[# unlabeled copies of H] by EXHAUSTIVE enumeration of choice profiles.

    Independent of Theorem 1: for every choice profile (there are
    C(n-1,k)^n of them) the adjacency matrix is built and the labeled copy
    count of H is computed directly from it (no inclusion-exclusion, no
    algebra).  The average is divided by |Aut(H)|.
    Only feasible for small n (product of binomials).
    """
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    if h_vertices > n:
        return 0.0
    e0 = np.array([u for u, _ in h_edges], dtype=np.int64)
    e1 = np.array([w for _, w in h_edges], dtype=np.int64)
    embs = np.array(list(itertools.permutations(range(n), h_vertices)), dtype=np.int64)
    a_idx = embs[:, e0]  # (M, |E(H)|) source vertices
    b_idx = embs[:, e1]  # (M, |E(H)|) target vertices
    mask_lists = [_ksubset_masks_excluding(n, k, v) for v in range(n)]
    total_labeled = 0
    n_profiles = 0
    for profile in itertools.product(*[range(len(ml)) for ml in mask_lists]):
        choices = np.array(
            [mask_lists[v][j] for v, j in enumerate(profile)], dtype=np.int64
        )
        adj = choices_to_adjacency(choices, n)
        ok = adj[a_idx, b_idx].all(axis=1)
        total_labeled += int(ok.sum())
        n_profiles += 1
    return total_labeled / n_profiles / aut_size


# ---------------------------------------------------------------------------
# Theorem 1: exact expected number of copies of a fixed graph H
# ---------------------------------------------------------------------------


def falling(n: int, m: int) -> int:
    """Falling factorial (n)_m = n!/(n-m)!; 0 when m > n."""
    if m > n:
        return 0
    return factorial(n) // factorial(n - m)


def p_H_in_kout(
    n: int, k: int, h_edges: tuple[tuple[int, int], ...], h_vertices: int
) -> float:
    """P(a fixed labeled embedding of H occurs in G(n,k-out)).

    Theorem 1:  P = sum_{F subset of E(H)} (-1)^{|F|} * prod_u C(n-1-d_u(F), k)
                / C(n-1, k)^{v}.
    """
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    if h_vertices > n:
        return 0.0
    denom = comb(n - 1, k) if n - 1 >= 0 else 0
    total = 0.0
    for r in range(len(h_edges) + 1):
        for F in itertools.combinations(h_edges, r):
            deg = [0] * h_vertices
            for u, w in F:
                deg[u] += 1
                deg[w] += 1
            term = (-1) ** r
            for u in range(h_vertices):
                m = n - 1 - deg[u]
                term *= comb(m, k) if m >= k else 0
            total += term
    return total / (denom**h_vertices)


def expected_labeled_copies(
    n: int, k: int, h_edges: tuple[tuple[int, int], ...], h_vertices: int
) -> float:
    """Exact E[# labeled copies of H] in G(n, k-out)."""
    return falling(n, h_vertices) * p_H_in_kout(n, k, h_edges, h_vertices)


def expected_unlabeled_copies(
    n: int, k: int, h_edges: tuple[tuple[int, int], ...], h_vertices: int, aut_size: int
) -> float:
    """E[# unlabeled copies of H] = E[# labeled copies] / |Aut(H)|."""
    return expected_labeled_copies(n, k, h_edges, h_vertices) / aut_size


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------


def p_triangle_closed(n: int, k: int) -> float:
    """Theorem 2 closed form: P(fixed triple is a triangle).

    With A = C(n-2,k)/C(n-1,k) and B = C(n-3,k)/C(n-1,k):
        P = 1 - 3 A^2 + 3 A^2 B - B^3.
    """
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    A = comb(n - 2, k) / comb(n - 1, k)
    B = comb(n - 3, k) / comb(n - 1, k) if n - 3 >= k else 0.0
    return 1.0 - 3.0 * A * A + 3.0 * A * A * B - B**3


def p_C4_closed(n: int, k: int) -> float:
    """Closed form for the 4-cycle (derived from Theorem 1):
    P = 1 - 4 A^2 + 4 A^2 B + 2 A^4 - 4 A^2 B^2 + B^4,
    A = C(n-2,k)/C(n-1,k), B = C(n-3,k)/C(n-1,k).
    """
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    A = comb(n - 2, k) / comb(n - 1, k)
    B = comb(n - 3, k) / comb(n - 1, k) if n - 3 >= k else 0.0
    return 1.0 - 4.0 * A * A + 4.0 * A * A * B + 2.0 * A**4 - 4.0 * A * A * B * B + B**4


def p_star_closed(n: int, k: int, m: int) -> float:
    """Closed form for the star K_{1,m} centered at the first vertex:
    P = sum_{j=0}^{m} (-1)^j C(m, j) C(n-1-j, k) [C(n-2,k)/C(n-1,k)]^j / C(n-1,k).
    """
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    A = comb(n - 2, k) / comb(n - 1, k)
    total = 0.0
    for j in range(m + 1):
        cval = comb(n - 1 - j, k) if n - 1 - j >= k else 0
        total += ((-1) ** j) * comb(m, j) * cval * (A**j)
    return total / comb(n - 1, k)


def expected_edges_closed(n: int, k: int) -> float:
    """Exact expected number of edges: C(n,2) * (2p - p^2), p = k/(n-1)."""
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    p = k / (n - 1)
    return comb(n, 2) * (2.0 * p - p * p)


def asymptotic_triangle_constant(k: int) -> float:
    """lim_{n->inf} E[#triangles] = (4k^3 - 3k^2)/3."""
    return (4.0 * k**3 - 3.0 * k**2) / 3.0


def er_comparison_limit(k: int) -> float:
    """lim of E[#triangles (k-out)] / E[#triangles (ER, matched density)]
    = 1 - 3/(4k)."""
    return 1.0 - 3.0 / (4.0 * k)


# ---------------------------------------------------------------------------
# Exact rational (Fraction) variants — independent arithmetic path
# ---------------------------------------------------------------------------


def p_H_exact(
    n: int, k: int, h_edges: tuple[tuple[int, int], ...], h_vertices: int
) -> Fraction:
    """Theorem 1 as an exact rational number (integer binomials only)."""
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    if h_vertices > n:
        return Fraction(0)
    denom = Fraction(comb(n - 1, k), 1) ** h_vertices
    total = Fraction(0)
    for r in range(len(h_edges) + 1):
        for F in itertools.combinations(h_edges, r):
            deg = [0] * h_vertices
            for u, w in F:
                deg[u] += 1
                deg[w] += 1
            term = Fraction((-1) ** r)
            for u in range(h_vertices):
                m = n - 1 - deg[u]
                term *= Fraction(comb(m, k), 1) if m >= k else Fraction(0)
            total += term
    return total / denom


def p_triangle_exact(n: int, k: int) -> Fraction:
    """Theorem 2 closed form as an exact rational number."""
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    A = Fraction(comb(n - 2, k), comb(n - 1, k))
    B = Fraction(comb(n - 3, k), comb(n - 1, k)) if n - 3 >= k else Fraction(0)
    return 1 - 3 * A * A + 3 * A * A * B - B**3


def p_C4_exact(n: int, k: int) -> Fraction:
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    A = Fraction(comb(n - 2, k), comb(n - 1, k))
    B = Fraction(comb(n - 3, k), comb(n - 1, k)) if n - 3 >= k else Fraction(0)
    return 1 - 4 * A * A + 4 * A * A * B + 2 * A**4 - 4 * A * A * B * B + B**4


def p_star_exact(n: int, k: int, m: int) -> Fraction:
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    A = Fraction(comb(n - 2, k), comb(n - 1, k))
    total = Fraction(0)
    for j in range(m + 1):
        cval = Fraction(comb(n - 1 - j, k), 1) if n - 1 - j >= k else Fraction(0)
        total += ((-1) ** j) * comb(m, j) * cval * (A**j)
    return total / Fraction(comb(n - 1, k), 1)


def expected_edges_exact(n: int, k: int) -> Fraction:
    if k > n - 1:
        raise ValueError("k must be at most n-1")
    p = Fraction(k, n - 1)
    return Fraction(comb(n, 2), 1) * (2 * p - p * p)


def triangle_moments_exact_frac(n: int, k: int) -> tuple[Fraction, Fraction, Fraction]:
    """Exact (E[T^2], E[T], Var(T)) as Fractions."""
    P3 = p_triangle_exact(n, k)
    N0 = Fraction(comb(n, 3), 1)
    PK4e = p_H_exact(n, k, K4_MINUS_EDGE[1], K4_MINUS_EDGE[0])
    Pbow = p_H_exact(n, k, BOWTIE[1], BOWTIE[0])
    P2t = p_H_exact(n, k, TWO_TRIANGLES[1], TWO_TRIANGLES[0])
    N2 = Fraction(12 * comb(n, 4), 1)
    N1 = Fraction(n * comb(n - 1, 2) * comb(n - 3, 2), 1)
    Nd = Fraction(comb(n, 3) * comb(n - 3, 3), 1)
    E2 = N0 * P3 + N2 * PK4e + N1 * Pbow + Nd * P2t
    E1 = N0 * P3
    return E2, E1, E2 - E1 * E1


# ---------------------------------------------------------------------------
# Exact first and second moments of the triangle count
# ---------------------------------------------------------------------------

# Union graphs of two triangles, classified by the overlap of the two
# (ordered) triples.  For s = t the union is a single triangle.
K4_MINUS_EDGE = (4, ((0, 1), (0, 2), (1, 2), (0, 3), (1, 3)))  # share an edge
BOWTIE = (5, ((0, 1), (1, 2), (0, 2), (0, 3), (3, 4), (0, 4)))  # share a vertex
TWO_TRIANGLES = (6, ((0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)))  # disjoint


def triangle_moments_exact(n: int, k: int) -> tuple[float, float, float]:
    """Exact (E[T^2], E[T], Var(T)) for the triangle count T of G(n,k-out).

    Decomposition over ordered pairs of triples (s, t) of vertices:
      E[T^2] = sum_s sum_t P(union edges of s and t all present),
    classified by |s & t| in {3, 2, 1, 0} via Theorem 1 applied to the
    union graph (K4 - e, bowtie, two disjoint triangles).
    """
    P3 = p_triangle_closed(n, k)
    N0 = comb(n, 3)
    PK4e = p_H_in_kout(n, k, K4_MINUS_EDGE[1], K4_MINUS_EDGE[0])
    Pbow = p_H_in_kout(n, k, BOWTIE[1], BOWTIE[0])
    P2t = p_H_in_kout(n, k, TWO_TRIANGLES[1], TWO_TRIANGLES[0])
    N2 = 12 * comb(n, 4)
    N1 = n * comb(n - 1, 2) * comb(n - 3, 2)
    Nd = comb(n, 3) * comb(n - 3, 3)
    E2 = N0 * P3 + N2 * PK4e + N1 * Pbow + Nd * P2t
    E1 = N0 * P3
    return E2, E1, E2 - E1 * E1


def triangle_moments_enum_full(n: int, k: int) -> tuple[float, float, float]:
    """(E[T^2], E[T], Var(T)) by exhaustive enumeration over ALL choice
    profiles (small n): ground truth for the moment decomposition."""
    mask_lists = [_ksubset_masks_excluding(n, k, v) for v in range(n)]
    s1 = 0.0
    s2 = 0.0
    count = 0
    for profile in itertools.product(*[range(len(ml)) for ml in mask_lists]):
        choices = np.array(
            [mask_lists[v][j] for v, j in enumerate(profile)], dtype=np.int64
        )
        adj = choices_to_adjacency(choices, n)
        T = count_triangles_adj(adj)
        s1 += T
        s2 += T * T
        count += 1
    return s2 / count, s1 / count, s2 / count - (s1 / count) ** 2


# ---------------------------------------------------------------------------
# Verification by exact enumeration
# ---------------------------------------------------------------------------


def p_triangle_enum_triple(n: int, k: int) -> float:
    """P(fixed triple forms a triangle) by exact enumeration of the three
    choosers' k-subset choices restricted to the relevant two targets each."""
    denom = comb(n - 1, k)
    # state encoding per chooser: which of its two specific targets chosen
    # 0: neither, 1: first, 2: second, 3: both
    pr = np.zeros(4)
    pr[0] = comb(n - 3, k) if n - 3 >= k else 0
    pr[1] = comb(n - 3, k - 1) if k >= 1 and n - 3 >= 0 else 0
    pr[2] = pr[1]
    pr[3] = comb(n - 3, k - 2) if k >= 2 and n - 3 >= 0 else 0
    pr = pr / denom
    total = 0.0
    for sa, sb, sc in itertools.product(range(4), repeat=3):
        # chooser a: targets (b, c): bit0 -> chose b, bit1 -> chose c
        xa, ya = bool(sa & 1), bool(sa & 2)
        # chooser b: targets (a, c)
        xb, yb = bool(sb & 1), bool(sb & 2)
        # chooser c: targets (a, b)
        xc, yc = bool(sc & 1), bool(sc & 2)
        if (xa or xb) and (yb or yc) and (xc or ya):
            total += pr[sa] * pr[sb] * pr[sc]
    return float(total)


def expected_triangles_enum_full(n: int, k: int) -> float:
    """E[#triangles] by direct enumeration over ALL choice profiles (small n)."""
    mask_lists = [_ksubset_masks_excluding(n, k, v) for v in range(n)]
    total = 0.0
    count = 0
    for profile in itertools.product(*[range(len(ml)) for ml in mask_lists]):
        choices = np.array(
            [mask_lists[v][j] for v, j in enumerate(profile)], dtype=np.int64
        )
        adj = choices_to_adjacency(choices, n)
        total += count_triangles_adj(adj)
        count += 1
    return total / count


def expected_triangles_enum_triple(n: int, k: int) -> float:
    """E[#triangles] = C(n,3) * P(triple), with P(triple) verified by direct
    enumeration over the triple choosers' full choice sets (subset-level, not
    the 4-state marginal)."""
    a0, b0, c0 = 0, 1, 2
    masks_a = _ksubset_masks_excluding(n, k, a0)
    masks_b = _ksubset_masks_excluding(n, k, b0)
    masks_c = _ksubset_masks_excluding(n, k, c0)

    def states(masks: np.ndarray, t1: int, t2: int) -> np.ndarray:
        return np.array(
            [((int(m) >> t1) & 1) | (((int(m) >> t2) & 1) << 1) for m in masks],
            dtype=np.int64,
        )

    sa_states = states(masks_a, b0, c0)
    sb_states = states(masks_b, a0, c0)
    sc_states = states(masks_c, a0, b0)
    # The three choosers are independent: the joint state histogram factorizes.
    na = np.bincount(sa_states, minlength=4).astype(np.int64)
    nb = np.bincount(sb_states, minlength=4).astype(np.int64)
    nc = np.bincount(sc_states, minlength=4).astype(np.int64)
    counts = na[:, None, None] * nb[None, :, None] * nc[None, None, :]
    total = 0
    for sa, sb, sc in itertools.product(range(4), repeat=3):
        xa, ya = bool(sa & 1), bool(sa & 2)
        xb, yb = bool(sb & 1), bool(sb & 2)
        xc, yc = bool(sc & 1), bool(sc & 2)
        if (xa or xb) and (yb or yc) and (xc or ya):
            total += counts[sa, sb, sc]
    return comb(n, 3) * total / (len(masks_a) * len(masks_b) * len(masks_c))


# ---------------------------------------------------------------------------
# Monte Carlo
# ---------------------------------------------------------------------------


def simulate_kout(n: int, k: int, reps: int, seed: int = 0) -> np.ndarray:
    """Monte Carlo sampling of triangle counts in G(n, k-out).

    Returns an array of triangle counts, one per replicate.
    """
    rng = np.random.default_rng(seed)
    counts = np.empty(reps, dtype=np.int64)
    for r in range(reps):
        perm = rng.permuted(np.tile(np.arange(n), (n, 1)), axis=1)
        picked = _drop_self_targets(perm, k)
        adj = np.zeros((n, n), dtype=np.int64)
        rows = np.repeat(np.arange(n), k)
        adj[rows, picked.ravel()] = 1
        adj = adj | adj.T
        ai = adj.astype(np.int64)
        counts[r] = int(np.trace(ai @ ai @ ai) // 6)
    return counts


def simulate_triangles_fast(
    n: int, k: int, reps: int, seed: int = 0, chunk: int = 32
) -> np.ndarray:
    """Vectorized Monte Carlo (batched, float32 BLAS) of triangle counts.

    For a batch of replicates the adjacency is built at once, symmetrized,
    and the triangle count is obtained as
        T = (1/3) * sum_{i<j} adj[i,j] * (A @ A)[i,j],
    which counts each triangle exactly once (via its three edges).  Integer
    values stay exact in float32 for n up to a few hundred; the final row
    sum is accumulated in float64.
    """
    rng = np.random.default_rng(seed)
    counts = np.empty(reps, dtype=np.int64)
    done = 0
    iu0, iu1 = np.triu_indices(n, 1)
    while done < reps:
        m = min(chunk, reps - done)
        perm = rng.permuted(np.tile(np.arange(n), (n * m, 1)), axis=1)
        picked = _drop_self_targets(perm, k)  # (n*m, k)
        rows = np.tile(np.repeat(np.arange(n), k), m)
        cols = picked.ravel()
        batch = np.zeros((m, n, n), dtype=np.float32)
        rep_idx = np.repeat(np.arange(m), n * k)
        batch[rep_idx, rows, cols] = 1.0
        batch = np.maximum(batch, np.transpose(batch, (0, 2, 1)))  # OR: 0/1
        common = np.matmul(batch, batch)  # (m, n, n)
        prod = batch[:, iu0, iu1] * common[:, iu0, iu1]
        counts[done : done + m] = np.rint(prod.sum(axis=1, dtype=np.float64) / 3.0)
        done += m
    return counts
