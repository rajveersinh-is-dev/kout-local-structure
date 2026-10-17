# Hostile Review — two adversarial passes against the paper and project

Author of record: Rajveersinh Vishal Pardeshi (Independent Researcher). Dates: pass 1
2026-10-17, pass 2 2026-10-18 (after MC bug fix and doc write-up). The point of this file
is to try to *break* the contribution the way a skeptical referee would, and to record
what survived and what was repaired.

## Pass 1 — "Is the math even right?" (played by adversarial persona: an RSG referee)

**HR1. "Theorem 1 looks like boilerplate inclusion–exclusion sprinkled with
hypergeometric terms; maybe a triviality."**
*Attempted demolition.* If P(triple omitted) factorizes per chooser, why is any of this
new — isn't this exactly the classical decomposable-structure computation for G(n,p)?
*Result of the attack.* Naive independence is *false* here: for fixed chooser u and
distinct v ≠ w, P(v∈S_u ∧ w∈S_u) = k(k−1)/((n−1)(n−2)) which is strictly less than
(k/(n−1))². A referee computing P(T) "like for G(n,p)" gets the naive value
(1−A²)³ — which our table shows is up to 2.97× the truth at k=1, n=200. The paper's
closed form `1−3A²+3A²B−B³` differs by an exact negative correction, and the equivalence
`P = naive + A³(C−A)(3−A²−AC−C²)` with (C−A)<0 makes the sign of the correction manifest.
**Verdict: survived** — the obstruction is real, quantitative, and showcased (naive
relative error 0.59–2.97 in staged table).

**HR2. "You claim 'variance equals mean in the limit', so is the count Poisson? You
only have two moments."**
*Attempted demolition.* Referees hate moment-tiered conjecture stated as fact.
*Response.* We state the Poisson limit explicitly as a **conjecture** backed by two exact
moments + Monte Carlo (all |z| ≤ 1.43) and by the G(n,p) precedent where two moments
suffice only with additional work. The text flags it, the theorem numbering keeps it out
of the "Theorem" namespace. **Verdict: survived as labelled conjecture**; cannot be
upgraded in-session (that would need third-moment + method-of-moments machinery — noted
as future work).

**HR3. "Your enumeration 'verification' covers only n ≤ 15 — is that teeth enough?"**
*Response.* For a triple, P(T) depends only on (n, k): we enumerate ALL choice profiles
(C(n−1,k))³ — exact Fractions, abs_diff 0 across n = 3..15 and k up to the feasibility
wall; plus triple-level enumeration (n ≤ 15) and second-moment cross-check vs full joint
enumeration for small n. Deterministic exhaustive checks at exact-precision are stronger
than loose MC for small n; MC covers the large-n regime (7 (n,k) cells, |z| ≤ 1.43 for
both mean and var). **Verdict: survived.**

**HR4. "Stars: did you confuse the probability of a star at fixed vertices with
expectations over all vertices?"**
*Response.* Separated: fixed-labeled-vertices probability (order n^{−m} exactly, with
explicit c_m(k) constants and finite-difference proof) vs E[#K_{1,m} at ALL vertices] =
(c_m(k)/m!)·n(1+o(1)). The cancellation table shows n^m·P → c_m(k) within 5·10⁻⁴ at
n = 1600 for m ≤ 6 — no order confusion. **Verdict: survived.**

**HR5. "The 'orientation invariant' Ω_k(H): is it actually a *theorem* or a
plausibility sketch?"**
*Honest concession.* The general-H statement is argued in write-up and hand+machine
verified on structured instances (K₃, C4, K_{1,m}); the general proof is given but a
referee may want to scrutinize the slot-allocation argument for subtle
double-counting. Mitigation in the paper: stated as Theorem 4 with a full proof, its
special cases machine-checked; flagged in limitations. **Verdict: survived with
declared scrutiny-need** (this is the least bulletproof claim — we do not hide it).

**HR6. "Dependencies between choosers: u's choices affect v's, error might accumulate."**
*Attempted demolition.* Choosers ARE independent sampler-by-sampler, but the
*symmetrized edges* create coupling across pairs. Theorem 1 handles it: expectations are
linear; only the ∏_u per-chooser avoidance factors enter, coupling arises from the
indicator product over edges, not from marginal independence claims. The one place
coupling could have corrupted the second moment — matrix product shortcuts — was caught
and *the implementation bug fixed* (batched matmul overcounting; corrected via OR-merge
`np.maximum(batch, batch.T)`, verified Fisher-consistency against a naive sampler).
**Verdict: bug did not propagate into results** (MC table regenerated).

## Pass 2 — "Is this publishable as claimed?" (played by a second persona: an
applications-oriented reviewer at Journal of Complex Networks)

**HR7. "Why should a networks-applications audience care about one more random graph
model?"**
*Response.* G(n, k-out) is the standard *sparse O(nk)-edge* alternative to ER, used in
practice for random pairwise key predistribution in ad-hoc networks and related
degree-bounded designs (as the Elumar–Sood–Yağan abstract itself notes: "Random K-out
graphs are known to achieve connectivity easily, i.e., with far fewer edges (O(n)) as
compared to classical random graph models"). Its *local* structure differs sharply from
ER at matched density, and we quantify it exactly: matched-density triangle ratio →
1 − 3/(4k) at fixed k (for k=1, k-out has ~25% of ER triangles; at k=2, ~37.5%; and it
varies with k), plus exact 1/n corrections for finite n. Applications that
that *simulate* on k-out and assume ER local structure are quantifiably wrong; the
exact-Fraction case k=1, n=200 shows naive assumption overshoots triangles ~×4.
**Verdict: survived** (relevance is honest but moderate; the significance score stays
moderate too).

**HR8. "Independent researcher, single-author, no affiliation — how do I know the code
is what the paper says?"**
*Response.* Full stack in the public repo: 43 pytest tests (will re-run on CI-less
machines), exhaustive-enumeration CSVs with exact Fractions, MC z-scores table,
regenerable one-command experiments. Reproduction reproducibility claims are limited to
what is represented in the repo. **Verdict: survived at the level of verifiability;
independent-author status itself is not a defect the repo can fix.**

**HR9. "The novelty claim — is 'not found in a search' really novelty?"**
*Response.* Project standard: classification "strongly supported originality with
documented residual risk" — never "first ever". Paper phrasing mirrors this: "we are not
aware of any published closed-form… after the search described in §2". **Verdict:
survived as-modestly-claimed.**

**HR10. "Two moments, expected values only — no distribution theory for general H;Aren't
you selling asymptotic tables as 'local structure'?"**
*Concession.* Honest limitation: expectations + second moments only; higher distribution
theory for general H is open. Limitations file and paper discussion both carry this.
**Verdict: survived only because the claim is explicitly scoped in the abstract — the
paper must not oversell distributional statements.**

## Post-review defect log

- Defect D1 (found pass 2, fixed): MC batched-product overcount (HR6) — fixed and
  regenerated; z-scores re-checked.
- Defect D2 (found pass 1): abbreviation of the corrections e1/e2/r1 signs in an early
  draft table — signs unified with `research_spec.md` and cross-checked against the
  sympy outputs during pass 2 restatement (no code change implied).
- Defect D3 (open, declared): orientation invariant scrutiny-need (HR5).
