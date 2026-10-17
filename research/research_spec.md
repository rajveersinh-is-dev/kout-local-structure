# Research Specification — The local structure of the random k-out graph

Date: 2026-10-17. Author of record: Rajveersinh Vishal Pardeshi (Independent Researcher).

## 1. Exact problem statement

Let n ≥ 2 and 1 ≤ k ≤ n-1. In the classical **random k-out graph** G(n, k-out)
(the model in Bollobás's *Random Graphs*, §"regular graphs", and Frieze & Karoński,
*Introduction to Random Graphs*, Ch. 22): every vertex v independently and uniformly
selects a k-subset S_v of the other n-1 vertices; the undirected simple graph has
edge {u,v} iff u chose v or v chose u (the underlying simple graph of the directed
union; equivalent to the standard definition via the symmetric difference rule used
in the literature, which the formal model section will spell out).

**Research question.** For a fixed finite simple graph H, what is the *exact*
expected number of copies of H in G(n, k-out) for finite n and k, and what is its
sharp asymptotic behaviour as n → ∞ with k (and H) fixed? In particular, how does
the triangle density of the k-out model compare with that of the Erdős–Rényi model
at the same expected edge density?

## 2. Strongest relevant existing results

- Connectivity and Hamiltonicity of G(n, k-out): Fenner–Frieze (1982); Frieze (2025
  survey chapter "Fast Construction on a Restricted Budget") — the k-out graph is
  Hamiltonian/connected for k ≥ 2 whp.
- Random mappings (the k=1 directed functional graph): expected # of m-cycles = 1/m
  (classical, e.g., Harris 1960; Flajolet–Odlyzko).
- Subgraph densities in G(n,p): Ruciński (1988); in the configuration model:
  standard generating-function counts.
- Clustering of specific models: Newman (2009) clustered random graphs; Holme–Kim;
  Işlak–Yeşiloğlu (J. Applied Probability 2025) analyze a *different* "clustering
  index" for ER/regular/BA/WS models — k-out is not treated there (verified 2026-10-17).

## 3. The gap

For the *symmetrized* k-out model (the standard one used in the literature), the
expected number of copies of a fixed small graph H — the model's "local structure" —
has apparently not been computed in closed form for finite n. The obstruction that
makes the question non-trivial (and distinct from G(n,p)) is the *within-chooser
negative correlation*: for a fixed chooser u and distinct targets v ≠ w,
P(v ∈ S_u ∧ w ∈ S_u) = k(k-1)/((n-1)(n-2)) < (k/(n-1))².
Naive independence computations give the wrong answer (demonstrated in the paper).

## 4. Proposed contribution (updated 2026-10-18 after symbolic verification)

1. **Theorem 1 (general subgraph density).** For any fixed simple graph H with
   vertex set V(H), the expected number of labeled copies of H in G(n, k-out) is
      E[N_H] = (n)_{v} · Σ_{F ⊆ E(H)} (-1)^{|F|} ∏_{u ∈ V(H)} C(n-1-d_u(F), k) / C(n-1, k)^{v},
   where d_u(F) is the degree of u in the spanning subgraph ([v], F). Proof via
   ∏_e (1 - X̄_uv X̄_vu) expansion and the exact hypergeometric "avoidance"
   probability per chooser, with independence across choosers. Verified:
   symbolic (sympy) vs closed forms; exhaustive enumeration over ALL choice
   profiles (triangles, C4, stars, small n); Monte Carlo (36 pytest checks pass).
2. **Theorem 2 (triangles).** Exact closed form
   P(T) = 1 - 3A² + 3A²B - B³, A = C(n-2,k)/C(n-1,k), B = C(n-3,k)/C(n-1,k) = A·C
   with C = (n-2-k)/(n-2); equivalent identity
   P(T) = (1-A²)³ + A³(C-A)(3-A²-AC-C²) = naive-independent value + exact
   (negative) correction. Limit E[#triangles] = (4k³-3k²)/3, with corrections
   e1 = -k²(k-1)(4k-3)/2, e2 = k²(6k³-39k²+46k-15)/6 (machine-exact series).
3. **Theorem 3 (E-R comparison).** At matched expected edge density the ratio of
   expected triangle counts tends to 1 - 3/(4k), with 1/n-correction
   r1 = 3(4k-3)/(8k). For C4 the ratio limit is 1 - 1/k + 1/(8k²).
4. **Theorem 4 (stars, cancellation).** Exact closed form for stars K_{1,m}
   (alternating sum), and the first nonzero term of P(K_{1,m} at fixed labeled
   vertices) in 1/n is exactly of order n^{-m} (machine-verified m ≤ 6; proof
   via finite differences: the t^ℓ coefficient of the summand is a polynomial
   in j of degree ≤ ℓ, annihilated by the m-th difference for ℓ < m). Hence
   E[#K_{1,m}] = (c_m(k)/m!)·n(1+o(1)) with explicit c_m(k) for m ≤ 6.
5. **Theorem 5 (general asymptotics, orientation invariant).** For any fixed
   simple graph H with e edges, lim E[#labeled copies]/n^{v-e} = Ω_k(H) :=
   Σ_σ ∏_x (k)_{t_x(σ)} over the 2^e edge-owner assignments σ (each edge
   assigned to exactly one of its endpoints), t_x(σ) = number of edges owned
   by x. Ω_k(H) = 0 iff no orientation with max owned-degree ≤ k exists.
   Hand+machine verified: Ω_k(K₃) = 8k³-6k², Ω_k(C₄) = 2k²(8k²-8k+1),
   Ω_k(K_{1,m}) = Σ_j C(m,j)(k)_j k^{m-j}.
6. **Theorem 6 (exact second moment of the triangle count).** E[T²] via Theorem 1
   on the three union graphs (K4-e, bowtie, two disjoint triangles); exact
   Var(T) for all n; Var(T) → (4k³-3k²)/3 = E[T] (limit), with underdispersion
   v1 = -k²(16k²-21k+7)/2 < 0. Poisson-type limit conjecture stated as
   conjecture (two moments proven only; MC supports).
7. **k=1 special cases.** P(T) = 2/(n-1)³ and P(C4) = 2/(n-1)⁴ exactly for all n
   (symbolically verified); naive-independence gap has the exact compact form
   q³ - P(T) = -A³(C-A)(3-A²-AC-C²) with g3 = 6k², g4 = -3k²(7k-9).
8. **Edges.** E[#edges] = nk - nk²/(2(n-1)) exactly (also = C(n,2)(2p-p²)).

## 5. Assumptions required

k ≤ n-1; H finite simple; asymptotics with k and H fixed as n → ∞. All
combinatorial identities elementary; no hidden assumptions.

## 6. Evidence required

- Complete written proofs (inclusion–exclusion + hypergeometric identity).
- Exhaustive- enumeration verification: for a fixed triple, P(triangle) depends
  only on (n, k) and can be verified exactly by enumerating all (C(n-1,k))³
  choice profiles of the three triple vertices — feasible for n ≤ ~15, k ≤ ~4.
- Monte Carlo verification of E[#triangles] for larger n and of the asymptotics.
- Exact symbolic cross-check of all closed forms with sympy/mpmath.

## 7. Falsification tests

- Formula vs exhaustive enumeration for n ∈ {3,...,9}, several k: any mismatch
  falsifies Theorem 1/2.
- Formula vs Monte Carlo (10^6+ samples) for large n: statistical consistency.
- E[#triangles] limit vs (4k³-3k²)/3 for k = 1..5.

## 8. Known limitations (declared up front)

- Expectations only (not full distribution) for most H; the distribution of the
  triangle count is addressed only via Poisson-limit *heuristics/verification*,
  or flagged as open if not rigorously established in-session.
- Asymptotics for fixed k, fixed H only.

## 9. Minimum evidence for the abstract's claims

Theorems 1–3 proven and verified by (i) exhaustive enumeration cross-checks and
(ii) Monte Carlo; all numeric claims reproducible from the repository code.

## 10. Novelty classification (pre-deep-dive)

**Potentially novel** (promising differentiation; meaningful uncertainty remains —
to be confirmed by a deeper targeted search during execution; closest neighbours
checked 2026-10-17: Bollobás survey (BA-model triangles, not k-out); Işlak–Yeşiloğlu
(no k-out); Frieze book (connectivity/Hamiltonicity focus)).
