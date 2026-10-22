# Limitations

Author of record: Rajveersinh Vishal Pardeshi (Independent Researcher). Date: 2026-10-18.

1. **Expectations, not distributions (mostly).** The paper computes exact expected copy
   counts for all fixed H, and exact first *and* second moments for the triangle count
   only. Beyond that (higher moments, the Poisson limit, full distribution theory,
   concentration for general H) is open; the Poisson-in-the-limit behaviour of the
   triangle count is stated as a **conjecture** supported only by two exact moments and
   Monte Carlo evidence (|z| ≤ 1.43 in 7 staged (n,k) cells).
2. **Fixed-k asymptotics.** All limit statements take n → ∞ with k and the template H
   fixed. Regimes with k growing with n (e.g. k = k(n) → ∞, dense-k-out windows) are
   outside the paper; the constants c_m(k), e1(k), r1(k) are polynomials in k valid at
   fixed k — using them at growing k would be extrapolation with no proof supplied.
3. **Orientation invariant (Theorem 4) general case.** Verified on machine
   for structured families (K₃, C4, K_{1,m}, etc) and cross-checked for arbitrary H (P4, Paw graph, Diamond graph) by symbolic limit evaluation, argued in
   write-up for general H.
4. **Star theorem scope.** The exact cancellation (first nonzero term exactly order
   n^{−m}) is proven and machine-checked for m ≤ 6; the general-m statement is argued
   (the polynomial-degree finite-difference argument is general) but the explicit
   constants tabulated only for m ≤ 6.
5. **No variance theory for general H.** The union-graph second-moment machinery is
   specialized: exact for triangles; not computed for other H.
6. **Literature verification depth.** Novelty rests on a documented search (full-text
   checks of the two closest sources; titles/abstracts of four nearest arXiv
   neighbours; 10+ query families) — not on exhaustive citation-chaining of all
   references-in-references (hostile-review concession HR9; residual risks R1–R5 in
   evidence_ledger.md).
7. **Verification stack bounded by in-session compute.** Exhaustive enumeration to
   n = 15; Monte Carlo 4,000–20,000 reps per staged cell. MC is statistical evidence,
   not proof; all exact claims rest on the symbolic/exact-arithmetic stack instead.
8. **Single independent author.** One author, no institutional affiliation, no external
   co-review: the hostile-review protocol in `research/hostile_review.md` (two
   adversarial passes) is the in-project substitute; it cannot replace real peer review
   and is not presented as such anywhere in the manuscript.
