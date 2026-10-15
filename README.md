# Systematic unequal-probability sampling under an uplink budget

Companion material for the paper

> *Systematic Unequal-Probability Sampling for Active Statistical Inference and Event Recall under an Uplink Budget on Passive Acoustic Nodes*
> (submitted to IEEE ISCAS 2027).

The paper has four pages. This repository holds what does not fit there: the proofs, the exact-enumeration checks, the full result tables, a reference implementation of the draws and the estimator, and the integer C code that ran on the microcontroller boards.

## Contents

| path | what it is |
|---|---|
| `proofs/proofs.pdf` (`proofs.tex`) | Assumptions and proofs of Lemma 1, Propositions 1 and 2 and Corollary 1, the asymptotic statement behind Lemma 1, the remarks made without proof, and two supplementary statements (exact variance difference between the single start and EXCL; sufficiency of the uplink sketch). Defines the successive-difference variance estimate. |
| `enumeration/` | Exact enumeration and linear-programming checks of the statements (`enum_s.py`, `test_madow_exact.py`, `run_all.sh`) with a README that lists what each check verifies and the expected output. |
| `tables/main_results.md` | Every dataset x budget cell of the main comparison: coverage, matched width, event recall and upload rate for the closed-loop rule (Trigger, Uniform+PPI, Duty cycle, Independent, EXCL, EXCL without onset alignment, Madow with open-loop and closed-loop price, union bound) and for Robust Sampling and Active Measurement with independent, EXCL and Madow draws; bias of the point estimates (including Trigger and the ordered pivotal design) and the decision delay of the ordered pivotal design; summary of the test deployments. |
| `tables/score_quality.md` | Seven score-quality levels (random score to near-oracle score): coverage and event recall of EXCL and of independent draws, recall ratio, matched-width ratio. |
| `tables/onboard_score.md` | The vessel result with the 14-band score and with the integer on-board score chain. |
| `tables/robustness.md` | Event definition (gaps merged: 0, 1, 2); block statistics; seven controller configurations and the block-mass cap; variance estimates of EXCL (the estimate of the paper and the successive-difference estimate); the check of the single-start condition on development and test data and the periodic stress test. |
| `tables/balanced_active_inference.md` | The comparison with balanced active inference (cube method balanced on model uncertainty, run offline on the whole test deployment) at the probabilities of Robust Sampling and Active Measurement: coverage, matched width and event recall next to independent, EXCL and single-start draws, two runs with different seeds. |
| `tables/node_cost.md` | Cycles per decision on a Cortex-M3 and a Cortex-M0+, energy, band-level score cost, object-code size, and the uplink share of the sketch. States what was and what was not measured. |
| `code/python/` | NumPy reference implementation: probability rules (closed-loop rule with price and token bucket; causal Robust Sampling and Active Measurement), the three draws (independent, EXCL, Madow), the estimator with `V-hat`, `V_plug` and the successive-difference estimate, events and event recall, coverage and matched width. `demo_synthetic.py` runs the three draws at the same probabilities on a synthetic clustered population. |
| `code/node/` | Portable integer C code of the decision rules, the test vectors that fix their output bit by bit, the native test, the band-level score, and the RIOT benchmark application. |

## The method in brief

A node sees units `t = 1..N` (segments of audio) and gives each an upload probability `p_t` (a function of an on-board score; the rule is an input). The label `y_t` becomes known only for uploaded units. The paper changes the draw, not the rule: systematic unequal-probability sampling at the same probabilities.

* **Single start (Madow), for score-only (open-loop) rules.** One uniform start `u` per recording; unit `t` is uploaded iff `u + integer` falls in `(C_{t-1}, C_t]`, where `C` is the cumulative probability within the recording.
* **Blocked restart (EXCL), for any rule, in particular a closed-loop rule that adapts a price to its own uploads.** Time is cut into blocks of mass `sum p_t <= 1`. A new block starts at a recording start, when the mass would exceed 1, or at an onset alignment (the bin prevalence reaches the development prevalence after at least two consecutive units below it); price and token bucket are updated only when a block closes. Each block draws one uniform `u_k`; unit `t` of the block is uploaded iff `A_{t-1} < u_k <= A_t`, `A` the within-block cumulative probability. At most one unit per block is uploaded, the conditional marginals are exact, and the restart costs some event recall (the price of the conditional marginals).
* **Estimator and interval.** With a predictor `f_t` (here a prequential predictor) and `r_t = y_t - f_t`:
  `theta-hat = (1/N) sum_t [f_t + xi_t r_t / p_t]`, `V-hat = (1/N^2) sum_t xi_t r_t^2 / p_t^2`, interval `theta-hat +- 1.96 sqrt(max(V-hat, V_plug))` with `V_plug = N^-2 sum_t v~_t / p_t`. `V-hat` is conservative in expectation for EXCL (Lemma 1); the single start has no such guarantee, but the condition under which it holds can be checked on development data.

## How to run

Requirements: Python 3 with NumPy (SciPy for the linear programmes in `enumeration/enum_s.py`), and a C compiler for the node code.

```
sh enumeration/run_all.sh                  # exact checks of the statements, about twenty seconds
cd code/python && python3 demo_synthetic.py   # three draws at the same probabilities, a few seconds
cd code/python && python3 test_reference.py   # self-checks of the reference implementation
cd code/node && make test                  # or: CC=gcc sh test.sh ; expects the line PASS: 66/66 ...
```

The demonstration prints, for the independent draw, EXCL and the single start at the same Active Measurement probabilities on a synthetic population (400 replications, `b = 0.05`), the coverage of the nominal 95% interval, event recall, nominal and matched width relative to independent draws, and the upload rate. A typical output:

```
design                  coverage  event recall  nominal width  matched width    rate/b
Independent                0.940         0.480          1.000          1.000     0.992
EXCL                       0.980         0.565          1.030          0.800     0.996
Madow (single start)       0.978         0.583          1.028          0.890     0.994
```

The native node test ends with `PASS: 66/66 vectors bit-exact (34 EXCL-family + 32 Madow), selfchecks 4/4`.

## Data

No data are redistributed here; the tables contain derived summary numbers only. The experiments of the paper use

* the BioDCASE 2026 Task 2 development set, *2026 BioDCASE Task 2: Development set*, Zenodo, doi:10.5281/zenodo.18832958, built from the open-access annotated Antarctic recordings of B. S. Miller et al., *Scientific Reports* 11:806, 2021, doi:10.1038/s41598-020-78995-8;
* the Belgian North Sea vessel-noise dataset of W. Decrop et al., *IEEE J. Sel. Topics Appl. Earth Observ. Remote Sens.* 18:19710--19726, 2025, doi:10.1109/JSTARS.2025.3593779 (clips paired with the Automatic Identification System);
* NOAA SanctSound raw passive acoustic data (Gray's Reef), NOAA National Centers for Environmental Information, doi:10.25921/saca-sp25.

Obtain them from their providers under their terms. The test vectors in `code/node/vectors/` contain score-bin indices only (no audio).

## Note for implementers

The decision code draws from a xorshift32 generator. Seed each node with a well-mixed non-zero 32-bit value, for example a hashed device identifier. Small consecutive integer seeds (1, 2, 3, ...) give strongly correlated first outputs, and the first output is the one that fixes the start of the first block or recording.

## What is not in this repository

* The paper itself, the audio, annotations or scores of the datasets above, and any third-party code.
* The scripts that run the full Monte-Carlo experiments on the datasets (they read the datasets); the tables contain their outputs.
