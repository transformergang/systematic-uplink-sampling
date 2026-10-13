# Exact checks of the statements

Requirements: Python 3, NumPy, SciPy (SciPy is used for the linear programmes in `enum_s.py` only). `sh run_all.sh` runs everything; it takes about twenty seconds.

| file | what it verifies |
|---|---|
| `enum_s.py` | Proposition 1 and Corollary 1 by exact enumeration of all `2^N` outcome vectors (small instances of up to a handful of units) and linear programming over all joint laws of the uploads, on 200 random instances per check (random probabilities, block partitions, events, closed-loop probability tables). |
| `madow_exact.py` | Exact variance of the single-start design (scan of the cumulative-probability axis, `O(n log n)`), the covariogram and the identity of Proposition 2(ii), and the exact variance of EXCL. Library used by the test below. |
| `test_madow_exact.py` | Proposition 2 and the variance comparison between the single start and EXCL: enumeration on rational probabilities, floor/ceiling counts, the identity `N^2 Var = sum_j gamma(j) - int gamma`, the criterion `Delta >= 0`, the sufficient condition, the periodic counterexample (ratio `L`), the counterexample with alternating probabilities 0.9 and 0.1 (12.52 and 67.8 against 11.11), the closed form of `N^2 (Var_Madow - Var_EXCL)`, and a scale test with 20000 units. |
| `../code/python/test_reference.py` | Properties of the reference implementation (one upload per block, floor/ceiling counts, marginals, unbiasedness, `V1 = V-hat`). |

## The checks of `enum_s.py`

| check | statement |
|---|---|
| E1 | The largest probability of hitting an event over all joint laws that have exact conditional marginals in every block (linear programme) equals `1 - prod_k (1 - m_k)`; EXCL attains it and has no repeated hits within a block (`W_e = 0`). |
| E2 | `E sum_k M_k - P(e hit) = R_e + W_e` with `R_e, W_e >= 0`, on vertices of the linear programme, open and closed loop. |
| E3 | With only the unconditional marginals the linear programme gives `min(1, sum p)` (the single start attains it); with delayed information the bound of E1 can be exceeded. |
| E4 | A design with correct unconditional but wrong conditional marginals is biased for some predictor with values in `{0,1}` that depends on earlier blocks. |
| E5 | The numbers of the closed-loop counterexample (EXCL is not uniformly optimal when the masses depend on earlier outcomes). |
| E6 | Corollary 1: marginals and variance of the budget, event hit probability `1 - prod_k (1 - m_k) >= 1 - prod_t (1 - p_t)`, the bracket for the gain, the variance change `-sum_k sum_{s != t in B_k} r_s r_t`, and the expected variance estimate. |
| E7 | Monte-Carlo agreement of `code/python/samplers.simulate_excl_openloop` (unit-level and event-level hit frequencies) with the exact law. |

## Expected output

`python3 enum_s.py` ends with a table; every row must read `PASS`:

```
| check | n inst | max abs err | strict failures | pass | s | what |
| E1 | 200 | 1.78e-15 | 0 | PASS | ... |
| E2 | 200 | 9.09e-16 | 0 | PASS | ... |
| E3 | 200 | 1.11e-15 | 0 | PASS | ... |
| E4 | 200 | 1.11e-16 | 0 | PASS | ... |
| E5 | 201 | 4.28e-16 | 0 | PASS | ... |
| E6 | 200 | 3.55e-15 | 0 | PASS | ... |
| E7 | 2 | 2.15 | 0 | PASS | ... |
```

(For E7 the printed value is the largest z-score, accepted below 4.5.) `python3 test_madow_exact.py` ends with

```
=== summary
test_M1 OK
test_M2 OK
test_M3 OK
test_scale OK
```

(`M1`, `M2`, `M3` are the three groups of tests: the basic properties of the single start, the variance identity and criterion, and the comparison with EXCL.) The result file of `enum_s.py` is written to `results/enum_results.json` (`--smoke` for a quick run with 5 instances).
