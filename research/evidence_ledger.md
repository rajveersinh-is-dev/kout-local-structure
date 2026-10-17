# Evidence Ledger

Author of record: Rajveersinh Vishal Pardeshi (Independent Researcher). Ledger date: 2026-10-18.
Every load-bearing claim of the project is listed with its evidence, verification status,
and provenance. "Verified" = checked by me in-session against the primary source.

## A. Mathematical results (core contribution)

| # | Claim | Evidence | Status |
|---|-------|----------|--------|
| A1 | Theorem 1 — general subgraph-count formula E[N_H] = (n)_v Σ_{F⊆E(H)} (−1)^{|F|} ∏_u C(n−1−d_u(F), k)/C(n−1,k)^v for every fixed simple H | Derivation in `src/kout.py` docstrings + paper §4; cross-checked three ways | **Verified**: (i) sympy identity vs closed forms (36 checks incl. 10 exact-Fraction spot checks, all pass); (ii) exhaustive enumeration over all choice profiles for triangles/C4/stars, n ≤ 15 (`results/validation_enumeration.csv`, `results/validation_triple_level.csv` — all diffs 0 or < 6e−17); (iii) pytest 43/43 pass incl. second-moment vs full enumeration |
| A2 | Theorem 2 — exact triangle P(T) = 1−3A²+3A²B−B³, A=C(n−2,k)/C(n−1,k), B=C(n−3,k)/C(n−1,k), B=AC, C=(n−2−k)/(n−2); equivalent identity P = naive + A³(C−A)(3−A²−AC−C²) with the correction ≤ 0 | Same as A1 | **Verified**: exact-Fraction spot checks + enumeration (abs_diff = 0) + MC \|z\| ≤ 1.43 for 7 cases (`results/monte_carlo.csv`) |
| A3 | lim E[#triangles] = (4k³−3k²)/3; (E−limit)·n → e1 = −k²(k−1)(4k−3)/2; e2 exact | Sympy series expansion + convergence table | **Verified**: `results/asymptotic_convergence.csv` — (E−limit)·n ↦ e1 to 7 decimals at n=1600 for k=1..3; ratio_to_ER ↦ limit; MC consistent |
| A4 | E-R matched-density comparison: triangle-count ratio → 1−3/(4k), 1/n-correction r1 = 3(4k−3)/(8k); C4 ratio → 1−1/k+1/(8k²) | Sympy series; convergence table | **Verified** (table above, ratio column) |
| A5 | Theorem 3 (stars K_{1,m}): leading behaviour of the fixed-vertices star probability is exactly order n^{−m}; E[#K_{1,m}] = (c_m(k)/m!)·n(1+o(1)); c_m(k) explicit for m ≤ 6 | Finite-difference proof sketch + numeric cancellation table `results/star_cancellation.csv` (n^m·P → c_m(k) monotone, \|diff\| ≤ 5e−4 at n=1600) | **Verified numerically**; proof written (finite differences: summand is degree-≤m polynomial in 1/n argument, m-th difference annihilates lower orders) annihilates lower orders); the order wall for m ≤ 6 additionally checked by enumeration |
| A6 | Theorem 4 (orientation invariant): lim E[#labeled copies]/n^{v−e} = Ω_k(H) = Σ over edge-owner assignments σ ∏_x (k)_{t_x(σ)}; Ω_k(H)=0 iff max-degree bound unachievable | General argument + hand-verified instances | **Verified** for K₃, C4, K_{1,m} against the closed forms; general H argued in paper §7 (write-up level; instances machine-checked) |
| A7 | Theorem 5 (second moment): exact E[T²] via union graphs (K4−e, bowtie, two disjoint triangles); Var(T) = v exact for all n; Var → E[T] (limit) with underdispersion v1 = −k²(16k²−21k+7)/2 | Sympy exact + tables | **Verified**: `results/variance_table.csv`; MC z_var \|z\| ≤ 1.43. Poisson-limit **stated as conjecture** (two moments only) |
| A8 | k=1 specials: P(T) = 2/(n−1)³, P(C4) = 2/(n−1)⁴ for all n | Sympy simplification to exact closed form | **Verified** (symbolic) + enumeration |
| A9 | E[#edges] = nk − nk²/(2(n−1)) exactly | Derivation + table | **Verified** (`results/expectations_table.csv` cross-check with C(n,2)(2p−p²)) |

Verification stack (all green, 2026-10-18): sympy 26 checks (`src/verify_symbolic.py`),
pytest 43 tests (`tests/test_kout.py`, run time ≈ 42 s), exact-Fraction code paths,
exhaustive enumeration up to n=15, Monte Carlo with \|z\| ≤ 1.43 for mean and variance.

## B. Literature / novelty claims

| # | Claim | Evidence | Status |
|---|-------|----------|--------|
| B1 | FK *Introduction to Random Graphs* Ch. 22 has NO triangle/local-structure content (sections: Connectivity, Perfect Matchings, Hamilton Cycles, Nearest Neighbor Graphs, Exercises, Notes) | Full text of Ch. 22 extracted from `https://www.math.cmu.edu/~af1p/BOOK.pdf` (PDF pp. 494–519); section headers grepped; content read | **Verified** (local extraction `fk_ch22.txt`, retrieved 2026-10-18) |
| B2 | FK Ch. 21 (random mappings) treats component structure only | Notes section text at p. 479 in same extraction | **Verified** |
| B3 | Sood 2024 CMU thesis: triangle content = doubly-Poisson clustered epidemic model; k-out only passing remarks | Full thesis text extracted and grepped (`sood_text.txt`, Carnegie Mellon, 2024), retrieved 2026-10-18 | **Verified** keyword+context read; not every page of 213 read — residual risk R2 |
| B4 | arXiv:1405.2129 = Frieze & Johansson, "On random k-out sub-graphs of large graphs" (2014); different problem (k-out subgraphs of a large fixed host; connectivity-type) | arXiv abstract page fetched 2026-10-18 | **Verified** |
| B5 | arXiv:2408.06451 = Işlak & Yeşiloğlu, "Analysis of Clustering and Degree Index in Random Graphs and Complex Networks" (v1 2024-08-12, v2 2025-11-15); Related DOI 10.1017/jpr.2025.10063; models: ER exact, regular/BA/WS simulated; **no k-out** | arXiv abstract page fetched 2026-10-18 (abstract text quotes the four model families) | **Verified** |
| B6 | arXiv:0910.0499 = Yağan & Makowski, "A zero-one law for the existence of triangles in random key graphs" (2009); random key graphs = Eschenauer–Gligor model, **not** k-out | arXiv abs page fetched 2026-10-18 (authors "Osman Yagan, Armand M. Makowski"; abstract quotes first/second moments) | **Verified** |
| B7 | arXiv:1701.03758 = Yağan & Makowski, "Counting triangles, tunable clustering and the small-world property in random key graphs (Extended version)" (2017) | arXiv abs page fetched 2026-10-18 | **Verified** |
| B8 | arXiv:2311.02319 = Elumar, Sood, Yağan, "On the Robustness, Connectivity and Giant Component Size of Random K-out Graphs" (2023) | arXiv result page fetched 2026-10-18 (title + authors from arXiv listing) | **Verified** |
| B9 | arXiv:2006.10638 = Sood, "Tight Bounds for Connectivity of Random K-out Graphs" (2020); arXiv:2103.01471 = Elumar, Sood, Yağan, "…Under Randomly Deleted Nodes" (2021) | arXiv listing snippets 2026-10-18 | **Verified** (titles + authors from arXiv search page) |
| B10 | Fenner & Frieze 1982, Combinatorica 2(4), 347–359 (random m-orientable graphs/digraphs, connectivity-type) | Bibliographic record checked in-session 2026-10-17 | **Verified** at the level shown; to re-check the exact article title against the journal page before submission — residual risk R4 |
| B11 | Janson–Łuczak–Ruciński, *Random Triangle Models*, Wiley (2000) — theory baseline for triangle counts in G(n,p) | Standard reference; bibliographic entry from general knowledge | **Unverified in-session** (no page fetched) — flag before submission |
| B12 | Newman, "Random graphs with clustering", Phys. Rev. E 80, 016101 (2009) | Standard reference | **Unverified in-session** — flag before submission |
| B13 | Bollobás, *Random Graphs*, 2nd ed., CUP (2001) — model definitions | Standard reference | Common knowledge; book exists; edition/year standard |
| B14 | No published closed-form expected H-copy counts for symmetrized G(n, k-out) | 10+ targeted web searches 2026-10-17/18; full-text checks B1, B3; citation-chains sampled | **Strongly supported, not provable exhaustively** — residual risks R1–R3 below |

## C. Journal-fee claims

| # | Claim | Evidence | Status |
|---|-------|----------|--------|
| C1 | Journal of Complex Networks offers standard (non-OA) and open-access licence options; charge applies only to the open-access licence | Official author-instructions page `https://academic.oup.com/comnet/pages/General_Instructions`, retrieved 2026-10-18 (browser); exact text: "offers the option of publishing under either a standard licence or an open access licence"; "You will need to pay an open access charge to publish under an open access licence." | **Verified** |
| C2 | JCN preprint policy: authors retain the right to make an Author's Original Version (preprint) available and this does not prevent submission; on acceptance the preprint status/DOI must be updated | Same page, "Preprints" section, same retrieval date | **Verified** |
| C3 | JCN scope includes "Mathematical and numerical analysis of networks", "Network theory and computer sciences", "Structural analysis of networks" | `https://academic.oup.com/comnet/pages/about`, retrieved 2026-10-18; online ISSN 2051-1329; 2025 IF 1.9 | **Verified** |
| C4 | Fallback journals (IMA J. Applied Math; The Computer Journal) have a no-mandatory-fee route | Not individually re-verified in-session | **Unverified** — do not rely on this until checked |

## D. Erratum-style notes on process

- Monte Carlo implementation bug (overcounting symmetric pairs in batched product) was
  found during cross-validation against a naive Fisher-consistent sampler and fixed
  (OR-merge via `np.maximum(batch, batch.T)` before the matmul). The fix is in
  `src/experiments.py`; the MC table was regenerated after the fix.
- No result was removed to improve "story"; all seven CSV result tables are included.

## E. Residual risks (honest register)

- **R1** — An exact-paper-collision could exist outside the crawls/searches performed
  (probability judged low: full-text checks of the two closest sources found nothing).
- **R2** — Sood thesis: keyword-and-context reading of the full text, not an exhaustive
  line-by-line read of all 213+ pages.
- **R3** — FK Ch. 22 Notes cite individual papers (matchings/connectivity) not
  individually chased; none of their stated topics intersect our observables.
- **R4** — Fenner–Frieze (1982) exact article title to be re-confirmed at submission time.
- **R5** — B11/B12 reference entries are standard but were not re-verified page-wise.
