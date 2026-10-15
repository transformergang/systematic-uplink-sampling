# Balanced active inference at the same probabilities

The paper (Results, first paragraph) compares the draws with balanced active inference (B. Chen, Z. Zhou, L. Peng and Z. Wang, NeurIPS 2025), which replaces independent draws by the cube method and balances the sample on model uncertainty. This file gives the full comparison.

**Setting.** Test deployments of the three datasets, the two open-loop rules of the paper (Robust Sampling, Active Measurement) with exactly the probabilities `p_t` of the main runs, `b = 0.05` and `0.10`. The balanced draw is the cube method (flight phase, then landing phase) at inclusion probabilities `p_t` with one balancing variable `u_t = sqrt(v_j(t))`, the bin-level residual standard deviation (the analogue of the model uncertainty used by balanced active inference). It selects from the whole test deployment at once, which a node cannot do; it is a reference, not a competitor on equal terms. Estimator, predictor, interval, events and the Uniform+PPI reference are those of the main runs. The cube method is the C++ code of the R package BalancedSampling (version 2.1.1), third-party code not included here. Independent, EXCL and single-start (Madow) draws are re-run in the same job. 500 replications per cell, run twice with different random seeds; the main table of the paper uses 2000 replications.

Cov: coverage of the nominal 95% interval; W: matched width relative to Uniform+PPI; Rec: event recall.

**Result.** At the same probabilities the balanced draw returns the same share of events as independent draws in every cell (paired difference below 0.002), and fewer than EXCL and the single start. Its matched width is close to that of independent draws (the sign of the difference changes between the two runs); against EXCL and the single start it is wider under Robust Sampling in both runs, and mixed under Active Measurement on whale and vessel data. Balancing on uncertainty does not spread the uploads in time, so it does not raise event recall.

## Run 1

`b = 0.05`

| Rule | Draw | Whale Cov | Whale W | Whale Rec | Vessel Cov | Vessel W | Vessel Rec | Gray's Reef Cov | Gray's Reef W | Gray's Reef Rec |
|---|---|---|---|---|---|---|---|---|---|---|
| Robust Sampling | independent | 0.95 | 0.93 | 0.185 | 0.95 | 0.94 | 0.197 | 0.96 | 0.93 | 0.083 |
| Robust Sampling | EXCL | 0.97 | 0.86 | 0.204 | 0.97 | 0.88 | 0.218 | 0.97 | 0.88 | 0.086 |
| Robust Sampling | Madow | 0.98 | 0.84 | 0.206 | 0.99 | 0.82 | 0.225 | 0.97 | 0.89 | 0.086 |
| Robust Sampling | balanced (cube) | 0.94 | 0.97 | 0.185 | 0.95 | 0.98 | 0.197 | 0.94 | 0.98 | 0.082 |
| Active Measurement | independent | 0.96 | 1.27 | 0.228 | 0.96 | 1.15 | 0.285 | 0.96 | 0.88 | 0.187 |
| Active Measurement | EXCL | 0.96 | 1.26 | 0.260 | 0.96 | 1.09 | 0.331 | 0.97 | 0.82 | 0.198 |
| Active Measurement | Madow | 0.97 | 1.21 | 0.262 | 0.97 | 1.04 | 0.342 | 0.99 | 0.75 | 0.203 |
| Active Measurement | balanced (cube) | 0.95 | 1.31 | 0.229 | 0.94 | 1.18 | 0.286 | 0.95 | 0.90 | 0.187 |

`b = 0.10`

| Rule | Draw | Whale Cov | Whale W | Whale Rec | Vessel Cov | Vessel W | Vessel Rec | Gray's Reef Cov | Gray's Reef W | Gray's Reef Rec |
|---|---|---|---|---|---|---|---|---|---|---|
| Robust Sampling | independent | 0.96 | 0.92 | 0.327 | 0.95 | 0.99 | 0.343 | 0.94 | 0.98 | 0.158 |
| Robust Sampling | EXCL | 0.98 | 0.84 | 0.382 | 0.99 | 0.84 | 0.394 | 0.98 | 0.82 | 0.169 |
| Robust Sampling | Madow | 0.98 | 0.83 | 0.386 | 0.99 | 0.82 | 0.411 | 0.99 | 0.69 | 0.171 |
| Robust Sampling | balanced (cube) | 0.96 | 0.91 | 0.327 | 0.96 | 0.94 | 0.344 | 0.95 | 0.99 | 0.159 |
| Active Measurement | independent | 0.96 | 1.25 | 0.388 | 0.95 | 1.22 | 0.469 | 0.95 | 0.93 | 0.336 |
| Active Measurement | EXCL | 0.97 | 1.23 | 0.472 | 0.98 | 1.05 | 0.535 | 0.98 | 0.82 | 0.351 |
| Active Measurement | Madow | 0.97 | 1.21 | 0.481 | 0.98 | 1.06 | 0.564 | 0.98 | 0.81 | 0.365 |
| Active Measurement | balanced (cube) | 0.96 | 1.24 | 0.388 | 0.94 | 1.24 | 0.469 | 0.96 | 0.88 | 0.337 |

## Run 2 (different seeds)

`b = 0.05`

| Rule | Draw | Whale Cov | Whale W | Whale Rec | Vessel Cov | Vessel W | Vessel Rec | Gray's Reef Cov | Gray's Reef W | Gray's Reef Rec |
|---|---|---|---|---|---|---|---|---|---|---|
| Robust Sampling | independent | 0.95 | 0.94 | 0.185 | 0.94 | 1.01 | 0.198 | 0.96 | 0.89 | 0.083 |
| Robust Sampling | EXCL | 0.97 | 0.89 | 0.204 | 0.97 | 0.90 | 0.218 | 0.96 | 0.90 | 0.086 |
| Robust Sampling | Madow | 0.98 | 0.87 | 0.206 | 0.98 | 0.80 | 0.225 | 0.97 | 0.86 | 0.087 |
| Robust Sampling | balanced (cube) | 0.96 | 0.90 | 0.185 | 0.95 | 0.95 | 0.197 | 0.94 | 0.97 | 0.082 |
| Active Measurement | independent | 0.94 | 1.37 | 0.228 | 0.95 | 1.22 | 0.286 | 0.96 | 0.85 | 0.187 |
| Active Measurement | EXCL | 0.96 | 1.23 | 0.260 | 0.96 | 1.17 | 0.330 | 0.97 | 0.84 | 0.200 |
| Active Measurement | Madow | 0.95 | 1.28 | 0.262 | 0.97 | 1.09 | 0.342 | 0.99 | 0.76 | 0.202 |
| Active Measurement | balanced (cube) | 0.95 | 1.27 | 0.229 | 0.97 | 1.11 | 0.286 | 0.95 | 0.91 | 0.187 |

`b = 0.10`

| Rule | Draw | Whale Cov | Whale W | Whale Rec | Vessel Cov | Vessel W | Vessel Rec | Gray's Reef Cov | Gray's Reef W | Gray's Reef Rec |
|---|---|---|---|---|---|---|---|---|---|---|
| Robust Sampling | independent | 0.95 | 0.93 | 0.327 | 0.95 | 0.98 | 0.345 | 0.94 | 0.98 | 0.159 |
| Robust Sampling | EXCL | 0.98 | 0.84 | 0.382 | 0.98 | 0.88 | 0.395 | 0.98 | 0.86 | 0.169 |
| Robust Sampling | Madow | 0.98 | 0.81 | 0.387 | 0.99 | 0.82 | 0.412 | 0.99 | 0.70 | 0.170 |
| Robust Sampling | balanced (cube) | 0.96 | 0.88 | 0.327 | 0.94 | 1.00 | 0.344 | 0.95 | 0.96 | 0.158 |
| Active Measurement | independent | 0.95 | 1.30 | 0.388 | 0.96 | 1.17 | 0.470 | 0.96 | 0.90 | 0.338 |
| Active Measurement | EXCL | 0.97 | 1.24 | 0.472 | 0.97 | 1.13 | 0.536 | 0.98 | 0.77 | 0.352 |
| Active Measurement | Madow | 0.97 | 1.22 | 0.481 | 0.98 | 1.04 | 0.565 | 0.99 | 0.82 | 0.365 |
| Active Measurement | balanced (cube) | 0.96 | 1.26 | 0.388 | 0.95 | 1.20 | 0.470 | 0.97 | 0.86 | 0.338 |
