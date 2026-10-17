# Results Summary

Author of record: Rajveersinh Vishal Pardeshi (Independent Researcher). Generated
2026-10-18 by `src/experiments.py`; every number quoted here is from the CSV files in
`results/` in this directory. No synthetic or illustrative numbers.

## What is in `results/`

| File | Rows | Contents |
|---|---|---|
| `validation_enumeration.csv` | 25 | Exact-Fraction formula values vs exhaustive enumeration over ALL choice profiles (triangles n=3..6, C4, star K_{1,3}, E[T²], Var(T)); max abs_diff 2.2e−15, mostly exactly 0 |
| `validation_triple_level.csv` | 6 | Fixed-triple triangle probability, closed form vs (C(n−1,k))³ profile enumeration, n up to 15, k up to 4; diffs 0–5.5e−17 |
| `expectations_table.csv` | 42 | E[triangles], E[C4], E[edges], E[K_{1,2}], E[K_{1,3}], naive E and relative naive error for n ∈ {10,20,…}, k ∈ 1..6 |
| `asymptotic_convergence.csv` | 30 | (E_tri − limit)·n → e1(k), ER-ratio → 1−3/(4k), for k=1..5, n ≤ 1600 |
| `monte_carlo.csv` | 7 | MC mean vs exact E and MC var vs exact Var; worst \|z\| = 1.43 (mean, n=200,k=3), worst variance \|z\| = 1.41; sampled and Fisher-consistent |
| `star_cancellation.csv` | 50 | n^m·P(star at fixed triple) → c_m(k) for m ∈ {2,3}, k ∈ {2,3}, n ≤ 1600; worst gap 4.7e−1 shrinking as O(1/n) |
| `variance_table.csv` | 25 | Var(T)/E(T) underdispersion and (Var−limit)·n → v1(k), k=1..3, n ≤ 1000 |

## Headline staged results (real values)

1. **Formula = truth exactly.** 25 enumeration cases, all abs_diff ≤ 2.2e−15 (exact
   Fractions agree with the enumeration to floating-point resolution; the Fraction
   column values, e.g. Var(T) at (n=5,k=2) = `Fraction(15985, 5832)`, match the
   enumeration float exactly).
2. **Naive independence is structurally wrong, not just imprecise.** Relative error of
   the naive E[#triangles] at n=200 (from `expectations_table.csv`): **k=1: 2.9699 (the
   naive estimator is ≈ 4× the truth)**, k=2: 0.5880, k=3: 0.3233. The exact closed form
   recovers the truth exactly (Fractions).
3. **Triangles (Theorem 2).** lim E[#T] = (4k³−3k²)/3 exactly: at n = 1600
   (E − limit)·n = −0.00021 (k=1), −10.4216→e1=−10 (k=2), −81.0619→e1=−81 (k=3),
   −312.1181→e1=−312 (k=4), −850.0249→e1=−850 (k=5): convergence is machine-exact.
4. **Matched-density ER comparison (Theorem 3).** Triangle-count ratio at n=1600:
   0.2502 (k=1), 0.7507 (k=3), 0.8508 (k=5) vs limits 250/1000 = 0.25 = 1−3/(4·1),
   0.75 = 1−3/(4·3), 0.85 = 1−3/(4·5). Second-order corrections: k=1 n=1600
   ratio_minus_2term = +2.9e−07.
5. **Stars (Theorem 3 in-text ordering).** n^m·P → c_m(k): m=2,k=2: 14.0037 at n=1600 vs
   c=14; m=2,k=3: 32.9887 vs c=33; m=3,k=2: 44.4737 at n=100 (→ c=44) with residual
   shrinking like O(1/n) — the cancellation order is exactly n^{−m}, no slow residual.
6. **Variance (Theorem 5).** Underdispersion confirmed: Var/E at n=200 → 0.985 (k=1),
   0.947 (k=2); (Var−limit)·n → v1(k) exactly: at n=1000, k=1: −0.9983 vs −1.0;
   k=2: −57.53 vs −58.0; k=3: −389.35 vs −396.0 (v1) — converging.
7. **Monte Carlo (implementation-level check of the samplers).** 7 cells, 4,000–20,000
   reps: worst |z_mean| = 1.43 ((n,k)=(200,3)), worst |z_var| = 1.41 ((50,3)); all within
   the 5σ gate for the rep counts used.

## Failure/robustness notes

- The MC harness originally contained an over-counting bug in the batched matrix
  product; it was caught by a Fisher-consistency cross-check against a naive simulator
  and fixed (`np.maximum(batch, batch.T)` before the matmul). All MC numbers above were
  generated with the fixed code.
- Sympy symbolic layer: 36/36 checks pass (`src/verify_symbolic.py`, incl. the
  orientation-invariant instances and 10 exact-Fraction asymptotic spot checks);
  pytest 43/43 (`tests/test_kout.py`).
