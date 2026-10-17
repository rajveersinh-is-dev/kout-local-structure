# Exact expected subgraph counts in the random *k*-out graph

**Repository:** https://github.com/rajveersinh-is-dev/kout-local-structure

Author of record: **Rajveersinh Vishal Pardeshi** (Independent Researcher).

This repository accompanies the manuscript
*"Exact expected subgraph counts in the random k-out graph: closed forms, a
matched-density Erdős–Rényi contrast, and an orientation invariant"*
(`paper/main.tex`, compiled PDF `paper/main.pdf`).

## Research question

In the symmetrized random *k*-out graph `G(n, k-out)` every vertex independently
selects a *k*-subset of the other `n−1` vertices and two vertices are adjacent if
**either** selected the other (OR-merge). What is the **exact expected number of
copies of a fixed graph H** in this model — exactly for finite `(n, k)`, and
asymptotically for fixed `H`, fixed `k`?

## Contribution (summary)

1. **Theorem 1 (general).** An exact closed formula for `E[#copies of H]` for every
   fixed finite simple graph `H` and all finite `(n, k)`, built from a per-chooser
   hypergeometric avoidance probability and chooser independence; within-chooser
   negative correlation `Pr(v∈S_u, w∈S_u) = k(k−1)/((n−1)(n−2))` is quantified.
2. **Triangles.** Exact closed form `P(△) = 1 − 3A² + 3A²B − B³`; limiting mean
   `E T → (4k³−3k²)/3`, limiting variance `Var(T) → (4k³−3k²)/3`, exact `1/n`
   terms, and `k=1` specials `P(△) = 2/(n−1)³`.
3. **ER contrast.** Matched-density ratio `E_tri^kout / E_tri^ER → 1 − 3/(4k)`
   (four-cycle: `1 − 1/k + 1/(8k²)`).
4. **Stars.** Exact cancellation theorem: the first non-zero term of the
   fixed-vertex-set star probability is exactly of order `n^(−m)`, and
   `E[#unlabelled K_{1,m}] = c_m(k)/m! · n(1+o(1))` with explicit
   `c_m(k) = Σ_j C(m,j)(k)_j k^(m−j)`; polynomials tabulated for `m ≤ 6`.
5. **Orientation invariant.** For any fixed `H`,
   `E[#labelled copies of H] ~ Ω_k(H) · n^(v−e)` where
   `Ω_k(H) = Σ_σ ∏_x (k)_{t_x(σ)}` sums over assignments σ of each edge of `H`
   to exactly one endpoint — an orientation invariant; `Ω_k(H) = 0` iff no
   assignment keeps every vertex's owned degree at most `k`. Complete
   two-sided Bonferroni bounds are proven for the truncation.
6. The Poisson-type limit of the triangle count is stated as a **conjecture**
   supported by two exact moments, and is not claimed as a theorem.

`Ω_k` examples: `Ω(K3) = 8k³−6k²`, `Ω(C4) = 8k²·d0` with
`d0 = 8k²−8k+1`, `Ω(K4−e) = 0` at `k=1`.

## Verification status (as generated 2026-10-18)

- **Unit tests:** `pytest` 43/43 pass (`tests/test_kout.py`).
- **Symbolic checks:** `src/verify_symbolic.py` — 36/36 exact identities pass
  (sympy, incl. orientation-invariant instances and exact-Fraction asymptotic
  spot checks with truncation bounds).
- **Exhaustive enumeration:** 25 staged small-`n` cases enumerate **all** choice
  profiles; max absolute difference `2.2e−15` (float resolution), exact Fractions
  on the closed-form side.
- **Monte Carlo:** 7 cells (4,000–20,000 replicates), worst `|z|` = 1.43 (mean),
  1.41 (variance).

## Reproduce

```bash
pip install -r requirements.txt   # numpy 2.4.6, sympy 1.14.0, pytest 9.1.1
python src/experiments.py         # regenerates every CSV in results/ (seeded MC)
python src/verify_symbolic.py     # 36 symbolic checks
python -m pytest -q               # 43 tests
```

`tests/` imports the module as `kout`; run pytest from the repository root
(`src/` is placed on the path as needed by the test file).

PDF: `cd paper && tectonic main.tex` (Tectonic 0.17.0 was used; any equivalent
LaTeX toolchain works).

## Repository layout

| Path | Contents |
|---|---|
| `paper/` | manuscript (`main.tex`, `references.bib`) and compiled `main.pdf` |
| `src/` | `kout.py` (formulae, closed forms, exact-Fraction paths, enumeration, MC), `experiments.py` (regenerates every CSV), `verify_symbolic.py` (36 checks) |
| `tests/` | pytest suite (43 tests) |
| `results/` | 8 CSV files + `summary.md`; **every number in the paper's tables comes from these files** |
| `research/` | research-process documents: literature review, evidence ledger (claim → evidence map, residual risks R1–R5), novelty assessment, two hostile-review passes, limitations, journal-selection record |

## Data provenance and licence

No external datasets are used. `results/*.csv` are generated entirely by
`src/experiments.py` from first principles over synthetic random graphs; the
enumeration and Monte Carlo code is deterministic up to the documented seed in
`src/experiments.py`. All committed CSVs were produced by this repository's own
code on 2026-10-18 (Python 3.14.6, numpy 2.4.6, sympy 1.14.0). No personal data,
no third-party datasets, no restricted-access sources.

Literature-verification documents quote publicly retrievable bibliographic
sources (arXiv abstract pages, publisher pages) with retrieval dates recorded.

## Limitations

See `research/limitations.md` and the Limitations section of the manuscript:
(a) distributional results are first/second-moment only (Poisson limit is a
conjecture); (b) fixed-`k`, fixed-`H` regimes only; (c) star constants tabulated
for `m ≤ 6` while the general-`m` formula is proven; (d) the orientation-invariant
implementation checks cover structured families (`K3`, `C4`, stars, `K4−e`), not
an arbitrary-`H` cross-check; (e) the novelty statement is a negative claim from a
targeted literature search, not a proof of no-overlap — residual risks R1–R5 are
recorded in `research/evidence_ledger.md`.

## Citation

`CITATION.cff` provides the citation metadata. Please cite:

> R. V. Pardeshi. *Exact expected subgraph counts in the random k-out graph:
> closed forms, a matched-density Erdős–Rényi contrast, and an orientation
> invariant.* Preprint, 2026.

## Licence

- **Code** (`src/`, `tests/`): MIT — see `LICENSE`.
- **Manuscript, research/reproducibility documents, and research data files**
  (`paper/`, `research/`, `results/`): Creative Commons Attribution 4.0
  International (CC BY 4.0) — see `LICENSE-DOCS`.

## Contact

Prepared as an independent-researcher project; no institutional affiliation and
no funding. Competing interests: none.
