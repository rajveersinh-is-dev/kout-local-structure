# Novelty Assessment

Author of record: Rajveersinh Vishal Pardeshi (Independent Researcher). Date: 2026-10-18.
Classification scheme: Strongly supported originality / Potentially novel / Incremental /
Not sufficiently novel.

## Verdict: **Strongly supported originality, with documented residual risk**

Supporting rationale (all items mirrored in `evidence_ledger.md` §B):

1. **The two closest sources were checked at full-text level, not by snippet.**
   - Frieze–Karoński, *Introduction to Random Graphs*, Ch. 22 (the canonical k-out
     chapter): complete text extracted from the authors' posted PDF
     (`https://www.math.cmu.edu/~af1p/BOOK.pdf`, retrieved 2026-10-18). Sections:
     Connectivity / Perfect Matchings / Hamilton Cycles / Nearest Neighbor Graphs /
     Exercises / Notes. No triangles, no counts of fixed subgraphs, no clustering.
   - Sood (2024 CMU thesis, the most recent large-scale treatment of k-out graphs):
     full text extracted and searched; its triangle content belongs to a doubly-Poisson
     clustered epidemic model in which triangles are inserted as *explicit building
     blocks*; k-out appears only in connectivity/robustness discussions.
2. **All nearest neighbours identified and differentiated.**
   - Frieze–Johansson arXiv:1405.2129: k-out subgraphs of a *large fixed host graph*,
     connectivity-type results — different model scope and different observable.
   - Işlak–Yeşiloğlu arXiv:2408.06451 (J. Applied Probability, Related DOI
     10.1017/jpr.2025.10063): clustering *index* with upper bounds for ER; simulations
     for regular/BA/WS; **k-out absent** (verified from abstract 2026-10-18).
   - Yağan–Makowski arXiv:0910.0499 and arXiv:1701.03758: triangle zero-one laws and
     clustering for *random key graphs* (Eschenauer–Gligor) — a different model with
     fundamentally different local dependence; methodologically the closest neighbour
     (triangles as local-structure probe, ER contrast at matched density), but its
     observables and proofs do not transfer to our questions.
   - Elumar–Sood–Yağan arXiv:2311.02319 / arXiv:2103.01471 / Sood arXiv:2006.10638:
     connectivity/robustness/giant-component — different observables.
3. **Targeted searches returned no collision**: 10+ distinct query families (see
   literature_review.md §5, log in evidence_ledger.md). Nothing surfaced computing
   exact finite-n expected copy counts for fixed H in symmetrized G(n, k-out).
4. **The contribution is a checkable artifact, which lowers "silent collision" risk**:
   exact formula (Theorem 1) + closed forms + exact 1/n-series (e1, e2, r1, c_m(k),
   v1) + orientation invariant. Any true prior duplicate would imply the same constants;
   none found.

## What the classification is NOT

- Not "first in the world" claims of the type banned by the honesty policy: we say
  "no collision found after the described search" rather than "first ever".
- Not claiming the *G(n,p)* part (expectations, asymptotics, Janson-type regimes) is
  new — that is standard (Janson–Łuczak–Ruciński-type triangle models; cited as
  background).

## Residual risk statement

Classification rests on (a) full-text checks of two closest sources, (b) titles +
abstracts of four nearest arXiv neighbours, (c) web-search sampling. It cannot rule out
an unindexed or paywalled paper computing the same quantities. Residual risks R1–R5 are
logged in `evidence_ledger.md` §E. If a collision is later found, the *method
(inclusion–exclusion + exact avoidance probabilities + orientation invariant) and the
verifiable artifact (code/tests/data) remain a contribution, at reduced aesthetic novelty.

## Novelty grade in rubric terms

Stage A rubric category "Novelty 25": initially scored at the "Potentially novel" level;
after the deep-dive, upgraded to the **strongly-supported-originality** level (score 23/25;
justification: full-text checks rather than snippets, no collision, differentiated from all
named neighbours). Significance scored 14/20: the triangle
comparison is a clean, quotable constant (1 − 3/(4k)) with applied meaning (fewer
triangles than ER at equal density), but it is one model family among many.
