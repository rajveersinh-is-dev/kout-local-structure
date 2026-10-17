# Literature Review — Local structure of the random k-out graph

Date of last search: 2026-10-18. Author of record: Rajveersinh Vishal Pardeshi (Independent Researcher).
Search log and per-item evidence notes are duplicated in `evidence_ledger.md`; this file
summarizes what the literature establishes and where the gap is.

## 1. Scope of the review

Question under study: for a fixed finite simple graph H, the *exact expected number of
labeled copies* of H in the symmetrized random k-out graph G(n, k-out) at finite (n, k),
and its sharp asymptotics for fixed k, H as n → ∞.

Reviewed bodies of work:
(a) the k-out graph literature (connectivity, Hamiltonicity, matchings, subgraphs of large hosts);
(b) the triangle/clustering literature for ER, key graphs, and clustered network models;
(c) random mappings (the k = 1 directed ancestor);
(d) subgraph-count asymptotics frameworks (Janson–Łuczak–Ruciński-type triangle models).

## 2. What the literature establishes

### 2.1 Global properties of G(n, k-out) — well studied

- **Fenner & Frieze (1982)** [Combinatorica 2(4), 1982, pp. 347–359] study connectivity-type
  questions for random m-orientable graphs/digraphs (an orientation model equivalent in
  spirit to k-out digraphs). Verified bibliographic record: Combinatorica, vol. 2, issue 4,
  1982, pp. 347–359.
- **Frieze & Karoński, *Introduction to Random Graphs* (Cambridge University Press), Ch. 22**
  is the standard textbook treatment of G(n, k-out). Evidence: I downloaded the
  authors' full text (`https://www.math.cmu.edu/~af1p/BOOK.pdf`, retrieved 2026-10-18) and
  extracted the complete Ch. 22 text (PDF pp. 494–519 of the 854-page file). Its sections
  are exactly: 22.1 Connectivity, 22.2 Perfect Matchings, 22.3 Hamilton Cycles,
  22.4 Nearest Neighbor Graphs, plus Exercises and Notes. **There is no section or result on
  triangles, counts of small fixed subgraphs, clustering, or any other local-structure
  observable.** Ch. 21 (random mappings) likewise contains only connectivity/component
  structure (Notes at p. 479 confirm: mappings are treated for component structure).
- **Frieze & Johansson, arXiv:1405.2129 (2014), "On random k-out sub-graphs of large graphs"**:
  considers a k-out *subgraph* of a fixed host graph G with large minimum degree and
  establishes connectivity-type results. The model class (k-out subgraphs of arbitrary
  large hosts) and the observable (connectivity) are both different from our question
  (complete host K_n; all fixed H; exact finite-n expectations).
- **Elumar–Sood–Yağan**: arXiv:2311.02319 "On the Robustness, Connectivity and Giant
  Component Size of Random K-out Graphs" (2023); see also arXiv:2103.01471 (2021) and,
  for tight connectivity bounds, Sood arXiv:2006.10638 (2020). All observable-focused on
  connectivity/robustness/giant-component scaling, not local structure.

### 2.2 Triangle/clustering content in *other* random graph models — well studied

- **Erdős–Rényi G(n, p):** expected subgraph counts are classical; the general theory of
  triangle counts (when the triangle count is asymptotically Poisson, normal, or in between)
  is the subject of **Janson–Łuczak–Ruciński, *Random Triangle Models* (Wiley, 2000)**.
  This is our comparison baseline: E[#T] = C(n,3)p³ exactly for G(n,p).
- **Random key graphs** (Eschenauer–Gligor scheme): **Yağan–Makowski, arXiv:0910.0499
  (2009), "A zero-one law for the existence of triangles in random key graphs"** prove a
  zero-one law for triangle appearance via first/second-moment methods; the extended
  treatment is **arXiv:1701.03758 (2017), "Counting triangles, tunable clustering and the
  small-world property in random key graphs (Extended version)"**, which also computes
  clustering coefficients and compares to ER. Key graphs are a *different* model
  (uniform key-pool intersection, zero-one threshold behaviour), but the paper's
  programme — triangles as a local-structure probe, contrasted against ER at matched
  density — is the nearest methodological sibling of ours. Their observables are
  zero-one laws and clustering coefficients of the key model, not exact finite-n expected
  copy counts in k-out.
- **Clustered/exponential models:** **Newman, "Random graphs with clustering"
  (Phys. Rev. E 80, 016101, 2009)** introduces creation of triangles via correlated
  stub-matching; Sood's **2024 CMU PhD thesis** (full text retrieved and searched,
  2026-10-18) treats triangles as *building blocks* ("single-edges and triangles" attached
  per node) in a doubly-Poisson multilayer epidemic model. Its k-out content is limited to
  passing remarks on connectivity; no intermediate-structure expectations for k-out.
- **Işlak–Yeşiloğlu, arXiv:2408.06451 "Analysis of Clustering and Degree Index in Random
  Graphs and Complex Networks"** (v1 2024, v2 Nov 2025; Related DOI 10.1017/jpr.2025.10063,
  confirmed on the arXiv abstract page): introduces a "clustering index" and computes
  upper bounds for ER, with simulations for random regular, Barabási–Albert, and
  Watts–Strogatz graphs. Verified from the abstract: **no k-out model appears anywhere**
  (the models named are ER, random regular, BA, WS).

### 2.3 Random mappings (k = 1 directed ancestor)

FK Ch. 21 covers uniform random mappings (connectivity and component structure, notes at
p. 479). The classical expected number of directed m-cycles in a random mapping is ∝ 1/m,
but this does not lift to our symmetrized model: in G(n, 1-out) the triangle requires an
undirected cycle on the triple, and we prove P(T) = 2/(n−1)³ exactly (both orientations),
which is different in structure from counting cycles of the functional digraph (e.g. the
directed 3-cycles are only one of the two orientations that our symmetrization makes
equivalent — and neither fact appears as an undirected-copy expectation in the references
above).

## 3. The gap this paper fills

No work found — after full-text checks of the two closest sources (FK Ch. 22; Sood thesis)
and 10+ targeted web searches (log in evidence_ledger.md) — computes:
1. any **exact finite-n** expected copy count of a fixed graph H in G(n, k-out);
2. the **within-chooser negative-correlation** obstruction P(v∈S_u ∧ w∈S_u) =
   k(k−1)/((n−1)(n−2)) < (k/(n−1))², which is what breaks naive G(n,p)-style shortcuts;
3. the **matched-density ER comparison** (engineering signal: k-out has fewer triangles
   than ER at equal expected edge density, ratio → 1 − 3/(4k));
4. a **second moment / dispersion** analysis of the triangle count in k-out.

## 4. Residual risks

- The FK book "Notes" for Ch. 22 may cite papers I have not individually opened
  (e.g. on matchings in k-out). I checked the notes' *topics* (they point to
  connectivity/nearest-neighbor references) but an exhaustive per-reference chase was not
  feasible in-session; recorded as residual risk R3 in evidence_ledger.md.
- Citation-search tools (Semantic Scholar "cited by" chains) were sampled, not exhausted.
- A paper very close in spirit could exist behind paywalls of obscure venues; mitigated
  by the fact that our contribution is exact finite-n formulas (a checkable artifact,
  unlikely to be replicated silently).

## 5. Search log summary

Queries run 2026-10-17/18 (details in evidence_ledger.md): "random k-out graph triangles",
"random k-out clustering coefficient", "expected number of triangles random k-out",
"subgraph density k-out graph", "random k-out graph local structure", "Frieze k-out
sub-graphs", "clustering index random graphs", "random key graphs triangles zero-one law",
"random neighbor graph clustering", "Newman random graphs with clustering k-out",
"Işlak Yeşiloğlu clustering", "Sood thesis k-out triangles", plus full-text greps of the
FK chapter and thesis extractions.
